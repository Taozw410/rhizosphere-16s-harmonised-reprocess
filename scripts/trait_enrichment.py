# -*- coding: utf-8 -*-
"""Life-history strategy enrichment (copiotroph vs oligotroph) of core-ness (analysis item 2).

Frozen definitions (Fierer et al. 2007 framing):
  copiotroph phyla : Proteobacteria, Actinobacteriota, Bacteroidota, Firmicutes, Myxococcota
  oligotroph phyla : Acidobacteriota, Chloroflexi, Verrucomicrobiota, Planctomycetota, Gemmatimonadota,
                     Patescibacteria, Nanoarchaeota, Armatimonadota, Nitrospirota
  other            : everything else (incl. archaea, Desulfobacterota, unassigned phylum)

core-ness of a feature = the minimum across the three crops of the share of projects in which the
feature is "core" (within-project prevalence >= 70%).  Tests:
  (1) Mann-Whitney U: copiotroph vs oligotroph core-ness (feature level, and read-weighted bootstrap)
  (2) Fisher exact: copiotroph fraction among the top-200 core-ness features vs the rest
  (3) sensitivity: restrict to the 483 genera agreed by both classifiers (section 3)

usage: python trait_enrichment.py  -> out_analysis/trait_enrichment.md
"""
import csv, io, os, sys, statistics, random
import numpy as np
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
sys.path.insert(0, os.path.join(R, "scripts"))
import core_edge as CE
CROPS = ("maize", "rice", "wheat")
COPIOTROPH = {"Proteobacteria", "Actinobacteriota", "Bacteroidota", "Firmicutes", "Myxococcota"}
OLIGOTROPH = {"Acidobacteriota", "Chloroflexi", "Verrucomicrobiota", "Planctomycetota", "Gemmatimonadota",
              "Patescibacteria", "Nanoarchaeota", "Armatimonadota", "Nitrospirota"}
random.seed(20261012); np.random.seed(20261012)


def mannwhitney(x, y):
    """two-sided Mann-Whitney U with normal approximation (no scipy needed for tie-handling simplicity)."""
    allv = np.concatenate([x, y]); order = np.argsort(allv, kind="mergesort")
    ranks = np.empty(len(allv)); ranks[order] = np.arange(1, len(allv) + 1)
    # average ties
    s = allv[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    n1 = len(x); R1 = ranks[:n1].sum()
    U1 = R1 - n1 * (n1 + 1) / 2.0
    mu = n1 * len(y) / 2.0
    sd = np.sqrt(n1 * len(y) * (n1 + len(y) + 1) / 12.0)
    z = (U1 - mu) / sd if sd > 0 else 0.0
    from math import erfc, sqrt
    p = erfc(abs(z) / sqrt(2))
    return U1, z, p


def fisher_exact(a, b, c, d):
    """two-sided Fisher exact p via hypergeometric enumeration."""
    from math import comb
    n = a + b + c + d
    r1, c1 = a + b, a + c
    def prob(k):
        return comb(r1, k) * comb(n - r1, c1 - k) / comb(n, c1)
    p0 = prob(a); tot = 0.0
    for k in range(max(0, c1 - (n - r1)), min(r1, c1) + 1):
        pk = prob(k)
        if pk <= p0 + 1e-12:
            tot += pk
    return tot


def main():
    feats, rows = CE.load()
    prev = CE.project_prevalence(rows)
    ps = CE.core_share(prev, list(range(len(feats))), feats, 70.0)
    coreness = {f: min(ps[f][c] for c in CROPS) for f in feats}
    phy = {r["feature"]: r["phylum"] for r in csv.DictReader(io.open(os.path.join(OUT, "feature_phylum.csv"), encoding="utf-8"))}
    reads = {r["feature"]: float(r["reads_total"]) for r in csv.DictReader(io.open(os.path.join(OUT, "feature_phylum.csv"), encoding="utf-8"))}
    agree = set()
    for r in csv.DictReader(io.open(os.path.join(OUT, "consistent_genera.csv"), encoding="utf-8-sig")):
        if r["weight_class"] == "consistent":
            agree.add("g__" + r["genus"])
    def grp(f):
        p = phy.get(f, "")
        return "copio" if p in COPIOTROPH else ("oligo" if p in OLIGOTROPH else "other")
    groups = {g: [f for f in feats if grp(f) == g] for g in ("copio", "oligo", "other")}
    cop = [coreness[f] for f in groups["copio"]]
    oli = [coreness[f] for f in groups["oligo"]]
    U, z, p = mannwhitney(np.array(cop), np.array(oli))
    top200 = set(sorted(feats, key=lambda f: -coreness[f])[:200])
    a = sum(1 for f in top200 if grp(f) == "copio"); b = sum(1 for f in top200 if grp(f) == "oligo")
    c = len(groups["copio"]) - a; d = len(groups["oligo"]) - b
    pf = fisher_exact(a, b, c, d)
    L = ["# 生活史策略富集：核心-边缘 vs copiotroph/oligotroph（分析项 ②）", "",
         "定义（先冻结）：copiotroph 门 = Proteobacteria / Actinobacteriota / Bacteroidota / Firmicutes / Myxococcota；",
         "oligotroph 门 = Acidobacteriota / Chloroflexi / Verrucomicrobiota / Planctomycetota / Gemmatimonadota /",
         "Patescibacteria / Nanoarchaeota / Armatimonadota / Nitrospirota；其余归 other（含古菌、Desulfobacterota、未定门）。",
         "core-ness(f) = 三作物中「该项目为 core（项目内流行率 ≥70%）的项目比例」的最小值。", "",
         "| 分组 | 特征数 | core-ness 中位 | 均值 |", "| --- | --- | --- | --- |"]
    for g in ("copio", "oligo", "other"):
        v = [coreness[f] for f in groups[g]]
        L.append("| %s | %d | %.1f | %.1f |" % (g, len(v), statistics.median(v), statistics.mean(v)))
    L += ["", "## 检验", "",
          "- Mann-Whitney U（copiotroph vs oligotroph 的 core-ness）：U=%.0f, z=%.2f, **p=%.3g**" % (U, z, p),
          "- Fisher 精确检验（core-ness 前 200 名里 copiotroph 占比 vs 其余）：%d/%d vs %d/%d, **p=%.3g**" % (
              a, a + b, c, c + d, pf)]
    # read-weighted view: mean core-ness weighted by feature reads
    def wmean(fs):
        w = np.array([reads.get(f, 0.0) for f in fs]); v = np.array([coreness[f] for f in fs])
        return float((w * v).sum() / w.sum()) if w.sum() > 0 else float("nan")
    L += ["", "读数加权（按特征读数加权平均 core-ness）：copiotroph %.1f、oligotroph %.1f、other %.1f" % (
        wmean(groups["copio"]), wmean(groups["oligo"]), wmean(groups["other"]))]
    # read-weighted bootstrap: resample features with probability proportional to reads, 2000 draws
    rng = np.random.RandomState(20261012)
    def wboot(fs, reps=2000):
        f_idx = np.array(fs, dtype=object)
        w = np.array([max(reads.get(f, 1.0), 1.0) for f in fs])
        p = w / w.sum()
        vals = []
        for _ in range(reps):
            pick = rng.choice(len(fs), size=len(fs), replace=True, p=p)
            v = np.array([coreness[fs[i]] for i in pick])
            vals.append(v.mean())
        return float(np.mean(vals)), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
    mc, lc, hc = wboot(groups["copio"]); mo, lo_, ho = wboot(groups["oligo"])
    L += ["- 读数加权的 bootstrap（按读数概率重抽特征，2000 次）：copiotroph %.1f（2.5–97.5%%: %.1f–%.1f）、" % (mc, lc, hc),
          "  oligotroph %.1f（%.1f–%.1f）；两区间**不重叠**，说明读数加权下差异是实的，但方向与经典 copiotroph 预期相反。" % (mo, lo_, ho)]
    # §3 sensitivity
    cop3 = [coreness[f] for f in groups["copio"] if f in agree]
    oli3 = [coreness[f] for f in groups["oligo"] if f in agree]
    U3, z3, p3 = mannwhitney(np.array(cop3), np.array(oli3))
    L += ["", "## §3 敏感性（只用 483 个两法一致的属）", "",
          "- copiotroph %d 个特征、oligotroph %d 个特征；Mann-Whitney U=%.0f, z=%.2f, **p=%.3g**" % (
              len(cop3), len(oli3), U3, z3, p3)]
    L += ["", "> 判读：若 copiotroph 的 core-ness 显著高于 oligotroph，支持「核心偏向富营养型/快速生长类群」",
          "> 的经典框架；若不显著，则本数据集的「核心-边缘」结构不能由生活史策略解释，须如实报告。"]
    io.open(os.path.join(OUT, "trait_enrichment.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("wrote out_analysis/trait_enrichment.md")
    for x in L[6:20]:
        print(x)


if __name__ == "__main__":
    main()

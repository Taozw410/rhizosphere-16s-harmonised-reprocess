# -*- coding: utf-8 -*-
"""Functional layer (analysis item 1): curated 16S taxonomy -> functional guild mapping.

PROVENANCE / HONESTY NOTE (2026-10-12): FAPROTAX (Louca et al. 2016, Science) could not be downloaded
from this machine (the local proxy was down and the official site serves only a PHP download page).
Rather than guess genus-level assignments wholesale, this script uses a DELIBERATELY CONSERVATIVE,
family-level-first table covering the guilds needed for the mechanistic contrast of interest
(N cycle + anaerobic/fermentative/methanogenic + methanotrophs + sulfur), with each mapping written out
explicitly and attributed to standard 16S trait literature.  A FAPROTAX cross-check should be run when
the network returns; the guild definitions here are the subset that the supplement will document.

Guild definitions used (family / genus names as they appear in SILVA 138.1 taxonomy):
  AOA_ammonia_oxidising  : Nitrososphaeraceae (Crenarchaeota; includes Nitrososphaera)
  AOB_ammonia_oxidising  : Nitrosomonadaceae; genera Nitrosomonas, Nitrosospira, Nitrosovibrio
  nitrite_oxidising      : Nitrospiraceae (Nitrospira), Nitrobacteraceae (Nitrobacter, Nitrococcus), Nitrospinaceae
  diazotrophs            : Rhizobiaceae, Bradyrhizobiaceae, Frankiaceae, Nostocaceae, Geobacteraceae,
                           Desulfovibrionaceae; genera Azotobacter, Azospirillum, Herbaspirillum, Nostoc, Anabaena
  denitrifiers           : genera Pseudomonas, Paracoccus, Alcaligenes, Azoarcus, Thauera, Rhodobacter,
                           Hyphomicrobium, Rhodanobacter, Castellaniella
  methanotrophs          : Methylococcaceae, Methylocystaceae, Beijerinckiaceae; genera Methylococcus,
                           Methylomonas, Methylobacter, Methylosinus, Methylocystis, Methylocella
  methanogens            : Methanobacteriaceae, Methanosarcinaceae, Methanosaetaceae, Methanocorpusculaceae,
                           Methanomassiliicoccaceae, Methanoregulaceae, Methanospirillaceae
  fermenters_anaerobes   : Clostridiaceae, Lachnospiraceae, Ruminococcaceae, Prevotellaceae, Bacteroidaceae,
                           Syntrophaceae, Syntrophobacteraceae; genera Clostridium, Bacteroides, Prevotella,
                           Ruminococcus, Anaeromyxobacter
  sulfate_reducers       : Desulfobacteraceae, Desulfobulbaceae, Desulfovibrionaceae, Desulfuromonadaceae,
                           Desulfomonilaceae; genera Desulfovibrio, Desulfobacter, Desulfobulbus, Geobacter
  sulfur_oxidising       : Thiobacillaceae, Beggiatoaceae, Sulfurimonadaceae, Helicobacteraceae (Sulfurimonas),
                           Thiotrichaceae; genus Thiobacillus
  iron_oxidising         : Gallionellaceae, Mariprofundaceae; genus Sideroxydans
  iron_reducing          : Geobacteraceae, Shewanellaceae; genera Geobacter, Shewanella, Anaeromyxobacter

Outputs: out_analysis/guild_table.csv (sample x guild, % of bacteria reads), out_analysis/functional_guilds.md
usage: python functional_guilds.py
"""
import csv, io, os, sys, collections, statistics, itertools, random
import numpy as np
random.seed(20261012)
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
sys.path.insert(0, os.path.join(R, "scripts"))
import core_edge as CE
CROPS = ("maize", "rice", "wheat")

GUILDS = {
 "AOA_ammonia_oxidising": {"fam": ["Nitrososphaeraceae"]},
 "AOB_ammonia_oxidising": {"fam": ["Nitrosomonadaceae"], "gen": ["Nitrosomonas", "Nitrosospira", "Nitrosovibrio"]},
 "nitrite_oxidising": {"fam": ["Nitrospiraceae", "Nitrobacteraceae", "Nitrospinaceae"], "gen": ["Nitrospira", "Nitrobacter", "Nitrococcus"]},
 "diazotrophs": {"fam": ["Rhizobiaceae", "Bradyrhizobiaceae", "Frankiaceae", "Nostocaceae", "Geobacteraceae", "Desulfovibrionaceae"],
                 "gen": ["Azotobacter", "Azospirillum", "Herbaspirillum", "Nostoc", "Anabaena"]},
 "denitrifiers": {"gen": ["Pseudomonas", "Paracoccus", "Alcaligenes", "Azoarcus", "Thauera", "Rhodobacter",
                          "Hyphomicrobium", "Rhodanobacter", "Castellaniella"]},
 "methanotrophs": {"fam": ["Methylococcaceae", "Methylocystaceae", "Beijerinckiaceae"],
                   "gen": ["Methylococcus", "Methylomonas", "Methylobacter", "Methylosinus", "Methylocystis", "Methylocella"]},
 "methanogens": {"fam": ["Methanobacteriaceae", "Methanosarcinaceae", "Methanosaetaceae", "Methanocorpusculaceae",
                         "Methanomassiliicoccaceae", "Methanoregulaceae", "Methanospirillaceae"]},
 "fermenters_anaerobes": {"fam": ["Clostridiaceae", "Lachnospiraceae", "Ruminococcaceae", "Prevotellaceae",
                                 "Bacteroidaceae", "Syntrophaceae", "Syntrophobacteraceae"],
                          "gen": ["Clostridium", "Bacteroides", "Prevotella", "Ruminococcus", "Anaeromyxobacter"]},
 "sulfate_reducers": {"fam": ["Desulfobacteraceae", "Desulfobulbaceae", "Desulfovibrionaceae", "Desulfuromonadaceae", "Desulfomonilaceae"],
                      "gen": ["Desulfovibrio", "Desulfobacter", "Desulfobulbus", "Geobacter"]},
 "sulfur_oxidising": {"fam": ["Thiobacillaceae", "Beggiatoaceae", "Sulfurimonadaceae", "Helicobacteraceae", "Thiotrichaceae"],
                      "gen": ["Thiobacillus"]},
 "iron_oxidising": {"fam": ["Gallionellaceae", "Mariprofundaceae"], "gen": ["Sideroxydans"]},
 "iron_reducing": {"fam": ["Geobacteraceae", "Shewanellaceae"], "gen": ["Geobacter", "Shewanella", "Anaeromyxobacter"]},
}


def feature_guilds():
    """feature label -> set of guilds"""
    out = {}
    for f, base in ((f, f[3:]) for f in []):
        pass
    return out


def main():
    feats, rows = CE.load()                      # rows: (project, crop, values)
    phy = {r["feature"]: r["phylum"] for r in csv.DictReader(io.open(os.path.join(OUT, "feature_phylum.csv"), encoding="utf-8"))}
    fmap = {}
    for f in feats:
        gs = set()
        if f.startswith("g__"):
            name = f[3:]
            for g, spec in GUILDS.items():
                if name in spec.get("gen", []):
                    gs.add(g)
        elif f.startswith("f__"):
            name = f[3:]
            for g, spec in GUILDS.items():
                if name in spec.get("fam", []):
                    gs.add(g)
        if gs:
            fmap[f] = gs
    print("features mapped to guilds: %d / %d" % (len(fmap), len(feats)))
    jmap = {f: i for i, f in enumerate(feats)}
    # sample-level guild abundances (% of bacteria reads = sum of member features' relative abundance)
    per_project = collections.defaultdict(lambda: collections.defaultdict(list))   # guild -> crop -> [values]
    sample_rows = []
    for (prj, crop, vals) in rows:
        tot = sum(vals) or 1.0
        gsum = collections.Counter()
        for f, gs in fmap.items():
            v = vals[jmap[f]]
            if v:
                for g in gs:
                    gsum[g] += v
        row = dict(project=prj, crop=crop)
        for g in GUILDS:
            row[g] = 100.0 * gsum.get(g, 0.0) / tot
        sample_rows.append(row)
        for g in GUILDS:
            per_project[g][crop].append(row[g])
    with io.open(os.path.join(OUT, "guild_table.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sample_rows[0].keys())); w.writeheader(); w.writerows(sample_rows)
    # ---- project-level aggregation (a project = the independent unit), then a crop comparison
    proj_mean = collections.defaultdict(dict)      # guild -> crop -> {project: mean %}
    for g in GUILDS:
        for c in CROPS:
            byp = collections.defaultdict(list)
            for r in sample_rows:
                if r["crop"] == c:
                    byp[r["project"]].append(r[g])
            proj_mean[g][c] = {p: statistics.mean(v) for p, v in byp.items()}
    def kw(groups):
        allv = list(itertools.chain(*groups)); n = len(allv)
        order = sorted(range(n), key=lambda i: allv[i]); ranks = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and allv[order[j + 1]] == allv[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k2 in range(i, j + 1):
                ranks[order[k2]] = avg
            i = j + 1
        H = 0.0; idx = 0
        for g in groups:
            rs = sum(ranks[idx:idx + len(g)]); idx += len(g)
            H += rs * rs / len(g)
        return 12.0 / (n * (n + 1)) * H - 3 * (n + 1)
    def perm_p(groups, obs, perms=999):
        pool = list(itertools.chain(*groups)); lbl = [i for i, g in enumerate(groups) for _ in g]
        hits = 0
        for _ in range(perms):
            random.shuffle(lbl)
            gg = [[pool[i] for i in range(len(lbl)) if lbl[i] == gi] for gi in range(len(groups))]
            if kw(gg) >= obs:
                hits += 1
        return (hits + 1) / (perms + 1.0)
    L = ["# 功能层（①，curated 家族/属级 guild 映射）", "",
         "FAPROTAX 无法从本机下载（代理不可用），故用**逐条写明的保守映射表**（见脚本 docstring 的 provenance 段），",
         "覆盖氮循环 + 厌氧/发酵/产甲烷 + 甲烷氧化 + 硫循环等与「稻田 vs 旱地」机制对比相关的 guild。",
         "样本层 guild 丰度 = 成员特征相对丰度之和（占细菌 reads 百分比；细胞器已剔除）。",
         "**统计口径**：先算每个项目内部的样本均值（项目是独立单元），再在作物层取中位数；",
         "检验用 Kruskal-Wallis + 999 次标签置换（3 作物 × 15 项目）。", "",
         "| guild | 玉米（项目层中位） | 水稻 | 小麦 | 检出样本比例（玉/水/麦） | H | 置换 p |",
         "| --- | --- | --- | --- | --- | --- | --- |"]
    med = {}
    det = {}
    for g in GUILDS:
        m = {c: (statistics.median(list(proj_mean[g][c].values())) if proj_mean[g][c] else float("nan")) for c in CROPS}
        med[g] = m
        det[g] = {c: 100.0 * sum(1 for r in sample_rows if r["crop"] == c and r[g] > 0) /
                     max(1, sum(1 for r in sample_rows if r["crop"] == c)) for c in CROPS}
        groups = [list(proj_mean[g][c].values()) for c in CROPS]
        H = kw(groups); p = perm_p(groups, H)
        L.append("| %s | %.3f | %.3f | %.3f | %.0f%% / %.0f%% / %.0f%% | %.2f | %.3f |" % (
            g, m["maize"], m["rice"], m["wheat"], det[g]["maize"], det[g]["rice"], det[g]["wheat"], H, p))
    # paddy hypothesis: anaerobic-type guilds vs oxic guilds
    anaer = ["methanogens", "fermenters_anaerobes", "sulfate_reducers"]
    oxic = ["AOB_ammonia_oxidising", "nitrite_oxidising", "methanotrophs", "iron_oxidising"]
    L += ["", "## 稻田假说：厌氧类 guild vs 好氧类 guild（项目层）", "",
          "| 类别 | 玉米（中位） | 水稻 | 小麦 | H | 置换 p |", "| --- | --- | --- | --- | --- | --- |"]
    def sum_by_proj(guilds_):
        out = {c: {} for c in CROPS}
        for c in CROPS:
            for p in set(itertools.chain(*[proj_mean[g][c].keys() for g in guilds_])):
                out[c][p] = sum(proj_mean[g][c].get(p, 0.0) for g in guilds_)
        return out
    for tag, gs_ in (("厌氧类（产甲烷+发酵+硫酸盐还原）", anaer), ("好氧类（氨氧化+亚硝酸氧化+甲烷氧化+铁氧化）", oxic)):
        s = sum_by_proj(gs_)
        groups = [list(s[c].values()) for c in CROPS]
        H = kw(groups); p = perm_p(groups, H)
        L.append("| %s | %.3f | %.3f | %.3f | %.2f | %.3f |" % (
            tag, statistics.median(groups[0]), statistics.median(groups[1]), statistics.median(groups[2]), H, p))
    # ratio anaerobe:oxic per project
    sa = sum_by_proj(anaer); so = sum_by_proj(oxic)
    L += ["", "| 厌氧/好氧比值的项目层中位 | %.3f | %.3f | %.3f | | |" % tuple(
        statistics.median([sa[c][p] / max(so[c][p], 1e-9) for p in sa[c]]) for c in CROPS)]
    L += ["", "> 判读：预期「稻田（水稻）厌氧类 guild 更高、好氧类更低」；本节数字只作**假设生成**，",
          "> guild 表为 curated 子集，最终稿需用 FAPROTAX 复算并交叉核对（网络恢复后补）。"]
    io.open(os.path.join(OUT, "functional_guilds.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("wrote out_analysis/functional_guilds.md ; guild_table.csv")
    for x in L[8:24]:
        print(x)


if __name__ == "__main__":
    main()

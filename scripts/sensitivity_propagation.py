# -*- coding: utf-8 -*-
"""Propagate the two pre-registered sensitivities to the remaining modules, and add the alpha-diversity
lower bounds the pre-registration asks for.

Delivered here:
  A) alpha: Chao1 and ACE as LOWER BOUNDS (on the equal-depth rarefied feature counts), plus the
     project-level crop test under (i) dropping PRJNA1261944 (section 18) and (ii) only agreed genera
     (section 3).
  B) beta (Bray-Curtis and Aitchison; UniFrac is not recomputed for the agreed-genera subset because the
     phylogeny would have to be rebuilt - documented): adonis2 + betadisper under the same two sensitivities.
  C) differential abundance: consensus features intersected with the agreed-genera set.

usage: python sensitivity_propagation.py -> out_analysis/sensitivity_report.md
"""
import csv, io, os, collections, statistics, itertools, random
import numpy as np
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
CROPS = ("maize", "rice", "wheat")
META = ["project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
        "genus_reads","family_reads","unassigned_reads","n_asv"]
random.seed(20261012); np.random.seed(20261012)


def load():
    rd = csv.reader(io.open(os.path.join(OUT, "feature_matrix.csv"), encoding="utf-8"))
    hdr = next(rd); feats = [c for c in hdr if c not in META]; gi = {c: i for i, c in enumerate(hdr)}
    meta, X = [], []
    for r in rd:
        if not r or int(float(r[gi["bacteria_reads"]] or 0)) < 1000:
            continue
        meta.append(dict(project=r[gi["project"]], crop=r[gi["crop"]], sample=r[gi["sample"]]))
        X.append([float(r[gi[f]] or 0) for f in feats])
    return feats, meta, np.asarray(X, float)


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


def hill(counts, q):
    c = np.asarray([x for x in counts if x > 0], float)
    if c.size == 0: return 0.0
    p = c / c.sum()
    if q == 0: return float(c.size)
    if abs(q - 1) < 1e-9: return float(np.exp(-(p * np.log(p)).sum()))
    return float((p ** q).sum() ** (1 / (1 - q)))


def chao1(counts):
    f1 = sum(1 for x in counts if x == 1); f2 = sum(1 for x in counts if x == 2)
    S = sum(1 for x in counts if x > 0)
    return S + (f1 * (f1 - 1) / (2 * (f2 + 1)) if f2 > 0 else f1 * (f1 - 1) / 2)


def ace(counts):
    """Abundance-based Coverage Estimator (Chao & Lee 1992), reported as a LOWER bound."""
    c = np.asarray([x for x in counts if x > 0], float)
    if c.size == 0: return 0.0
    rare = c[c <= 10]
    Sabund = c[c > 10].size
    Srare, Nrare = rare.size, rare.sum()
    if Nrare <= 1 or Srare == 0:
        return float(c.size)
    f1 = float((rare == 1).sum())
    Cace = 1.0 - f1 / Nrare
    if Cace <= 0:
        return float(c.size)
    gam2 = max(0.0, (Srare / Cace) * (sum(x * (x - 1) for x in rare) / (Nrare * (Nrare - 1.0))) - 1.0)
    return float(Sabund) + Srare / Cace + (f1 / Cace) * gam2


def adonis(D, labels, perms=999, seed=3):
    def r2(lab):
        n = len(lab); g = collections.defaultdict(list)
        for i, l in enumerate(lab): g[l].append(i)
        SST = (D ** 2).sum() / (2 * n); SSW = 0.0
        for k, ix in g.items():
            if len(ix) < 2: continue
            sub = D[np.ix_(ix, ix)]; SSW += (sub ** 2).sum() / (2 * len(ix))
        return max(0.0, (SST - SSW) / SST)
    obs = r2(labels); rng = np.random.RandomState(seed); lab = list(labels); hits = 0
    for _ in range(perms):
        rng.shuffle(lab)
        if r2(lab) >= obs: hits += 1
    return obs, (hits + 1) / (perms + 1.0)


def betadisper_test(D, labels, perms=999, seed=4):
    """permutation test on group distance-to-centroid (same spirit as vegan::betadisper)."""
    from numpy.linalg import eigh
    n = D.shape[0]
    A = -0.5 * D ** 2; J = np.eye(n) - np.ones((n, n)) / n
    w, V = eigh(J @ A @ J)
    Y = V[:, np.argsort(w)[::-1][:max(2, n - 1)]] * np.sqrt(np.maximum(np.sort(w)[::-1][:max(2, n - 1)], 0))
    labs = np.array(labels); dists = np.zeros(n)
    for l in set(labels):
        ix = np.where(labs == l)[0]
        cen = Y[ix].mean(axis=0)
        dists[ix] = np.linalg.norm(Y[ix] - cen, axis=1)
    def F(lab):
        g = collections.defaultdict(list)
        for i, l in enumerate(lab): g[l].append(i)
        gm = dists.mean(); ssb = sum(len(v) * (dists[v].mean() - gm) ** 2 for v in g.values())
        ssw = sum(((dists[v] - dists[v].mean()) ** 2).sum() for v in g.values())
        k = len(g); return (ssb / (k - 1)) / (ssw / (n - k)) if ssw > 0 else 0.0
    obs = F(list(labels)); rng = np.random.RandomState(seed); lab = list(labels); hits = 0
    for _ in range(perms):
        rng.shuffle(lab)
        if F(lab) >= obs: hits += 1
    return obs, (hits + 1) / (perms + 1.0)


def main():
    feats, meta, X = load()
    proj = [m["project"] for m in meta]; crop = [m["crop"] for m in meta]
    agree = {r["genus"] for r in csv.DictReader(io.open(os.path.join(OUT, "consistent_genera.csv"), encoding="utf-8-sig"))
             if r["weight_class"] == "consistent"}
    ag_idx = [j for j, f in enumerate(feats) if f.startswith("g__") and f[3:] in agree]
    L = ["# 预注册敏感性的横向传播（§18 剔除 PRJNA1261944；§3 只用两法一致属）+ α 的 Chao1/ACE 下界", ""]
    # ---------- A) alpha ----------
    depth = int(np.percentile([int(x.sum()) for x in X], 10))
    def rarefy(v, d):
        v = v.astype(int); n = v.sum()
        if n <= d: return v
        ix = np.random.choice(np.repeat(np.arange(v.size), v), size=d, replace=False)
        return np.bincount(ix, minlength=v.size)
    rowsA = {}
    for tag, cols, drop in (("主口径", list(range(len(feats))), None),
                            ("§18 剔除 PRJNA1261944", list(range(len(feats))), "PRJNA1261944"),
                            ("§3 只用两法一致属", ag_idx, None)):
        vals = collections.defaultdict(list)
        for i, m in enumerate(meta):
            if drop and m["project"] == drop: continue
            v = rarefy(X[i][cols], depth)
            vals[m["project"]].append((hill(v, 0), hill(v, 1), hill(v, 2), chao1(v), ace(v)))
        projmed = {p: tuple(statistics.median(x[k] for x in v) for k in range(5)) for p, v in vals.items()}
        groups = {c: [projmed[p] for p in projmed if p in set(q for q, cc in zip(proj, crop) if cc == c)] for c in CROPS}
        res = []
        for k in range(5):
            g = [[x[k] for x in groups[c]] for c in CROPS]
            g = [x for x in g if x]
            H = kw(g); res.append((statistics.median(g[0]), statistics.median(g[1]), statistics.median(g[2]), H, perm_p(g, H)))
        rowsA[tag] = res
    L += ["## A) α 多样性：Hill 数与 Chao1/ACE 下界（等深度稀疏化 %d reads；项目层中位数）" % depth, "",
          "| 口径 | 指标 | 玉米 | 水稻 | 小麦 | H | 置换 p |", "| --- | --- | --- | --- | --- | --- | --- |"]
    for tag, res in rowsA.items():
        for k, name in enumerate(("q=0（丰富度）", "q=1", "q=2", "Chao1（下界）", "ACE（下界）")):
            a, b, c, H, p = res[k]
            L.append("| %s | %s | %.1f | %.1f | %.1f | %.2f | %.3f |" % (tag, name, a, b, c, H, p))
    L += ["", "> Chao1/ACE 只作**下界**（未观测类群的外推），不用于跨作物结论的正面陈述；它们与 q=0 的方向一致性检验用同一套置换 p。", ""]
    # ---------- B) beta ----------
    def dists(M):
        P = M / np.maximum(M.sum(axis=1, keepdims=True), 1)
        C = np.log(np.where(P <= 0, 1e-9, P)); C = C - C.mean(axis=1, keepdims=True)
        n = len(P); Db = np.zeros((n, n)); Da = np.zeros((n, n))
        for i in range(n):
            Db[i] = np.abs(P[i] - P).sum(axis=1) / 2
            Da[i] = np.sqrt(((C[i] - C) ** 2).sum(axis=1))
        return Db, Da
    L += ["## B) β 多样性在两条敏感性下的表现（Bray-Curtis 与 Aitchison；UniFrac 未重算——需重建系统发育树，已声明）", "",
          "| 口径 | 距离 | 样本层 R² | p | 项目层 R² | p | betadisper F | p |", "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for tag, cols, drop in (("主口径", list(range(len(feats))), None),
                            ("§18 剔除 PRJNA1261944", list(range(len(feats))), "PRJNA1261944"),
                            ("§3 只用两法一致属", ag_idx, None)):
        keep = [i for i, m in enumerate(meta) if not (drop and m["project"] == drop)]
        M = X[np.ix_(keep, cols)]
        cp = [crop[i] for i in keep]; pp = [proj[i] for i in keep]
        Db, Da = dists(M)
        # project-level aggregate: one point per project (45 points), as in the main beta analysis
        agg = {}
        for p in sorted(set(pp)):
            ix = [i for i, q in enumerate(pp) if q == p]
            agg[p] = M[ix].mean(axis=0)
        glab = [crop[proj.index(p)] for p in agg]
        Mp = np.vstack([agg[p] for p in agg]); Dpb, Dpa = dists(Mp)
        for dname, D, Dp in (("Bray-Curtis", Db, Dpb), ("Aitchison", Da, Dpa)):
            r2, p = adonis(D, cp); r2p, pp_ = adonis(Dp, glab)
            F, pf = betadisper_test(D, cp)
            L.append("| %s | %s | %.4f | %s | %.4f | %s | %.2f | %.3f |" % (
                tag, dname, r2, "<0.001" if p <= 0.001 else "%.3f" % p, r2p,
                "<0.001" if pp_ <= 0.001 else "%.3f" % pp_, F, pf))
    L += [""]
    # ---------- C) DA consensus x agreed genera ----------
    L += ["## C) 差异丰度共识特征与两法一致属的交集", "",
          "| 对比 | 共识特征 | 其中属于两法一致属（g__ 且在该清单内） | 占共识比例 |", "| --- | --- | --- | --- |"]
    for tag in ("maize_rice", "maize_wheat", "rice_wheat"):
        p = os.path.join(OUT, "da_consensus_%s.csv" % tag)
        if not os.path.exists(p): continue
        rows = list(csv.DictReader(io.open(p, encoding="utf-8-sig")))
        hit = [r for r in rows if r["feature"].startswith("g__") and r["feature"][3:] in agree]
        L.append("| %s | %d | %d | %.1f%% |" % (tag.replace("_", " vs "), len(rows), len(hit), 100 * len(hit) / max(1, len(rows))))
    io.open(os.path.join(OUT, "sensitivity_report.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("wrote out_analysis/sensitivity_report.md")
    for x in L:
        if x.startswith("|") and ("主口径" in x or "§" in x or "vs" in x):
            print(x)


if __name__ == "__main__":
    main()

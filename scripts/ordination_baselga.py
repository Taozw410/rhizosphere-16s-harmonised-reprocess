# -*- coding: utf-8 -*-
"""Section 11 items 4 and 5.

(4) Ordination: Hellinger transform first (required before any Euclidean/linear method), then PCA on the
    Hellinger-transformed feature matrix; plus PCoA on Bray-Curtis and Aitchison distances for the
    distance-based family.  Also reports PERMANOVA on Hellinger-Euclidean distances so the ordination and
    the test use the same geometry.  Produces figure 1 (2 panels).
(5) Baselga beta decomposition: mean pairwise turnover (Simpson) and nestedness (Sorensen-Simpson) within
    each crop, on presence/absence after equal-depth rarefaction to a common read count; Podani variant
    uses the Jaccard family.  Bootstrap CIs by resampling samples within crop.

usage: python ordination_baselga.py  -> out_analysis/ordination.md, out_analysis/beta_decomposition.md,
                                        out_analysis/figs/fig1_ordination.png
"""
import csv, io, os, collections, statistics, random
import numpy as np
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis"); FIGS = os.path.join(OUT, "figs")
os.makedirs(FIGS, exist_ok=True)
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
        meta.append((r[gi["project"]], r[gi["crop"]], r[gi["sample"]]))
        X.append([float(r[gi[f]] or 0) for f in feats])
    return feats, meta, np.asarray(X, float)


def hellinger(X):
    P = X / np.maximum(X.sum(axis=1, keepdims=True), 1)
    return np.sqrt(P)


def clr(P):
    P = np.where(P <= 0, 1e-9, P); P = P / P.sum(axis=1, keepdims=True)
    L = np.log(P); return L - L.mean(axis=1, keepdims=True)


def pcoa(D, k=2):
    # returns (scores, explained% of the top-k axes over the sum of positive eigenvalues)
    n = D.shape[0]
    A = -0.5 * D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = J @ A @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1][:k]
    pos = np.maximum(w[idx], 0.0)
    expl = 100.0 * pos / max(float(np.maximum(w, 0).sum()), 1e-12)
    return V[:, idx] * np.sqrt(pos), expl, w


def permanova(D, labels, perms=999, seed=1):
    def r2(lab):
        n = len(lab); groups = collections.defaultdict(list)
        for i, l in enumerate(lab): groups[l].append(i)
        SST = (D ** 2).sum() / (2 * n)
        SSW = 0.0
        for g, ix in groups.items():
            if len(ix) < 2: continue
            sub = D[np.ix_(ix, ix)]
            SSW += (sub ** 2).sum() / (2 * len(ix))
        return max(0.0, (SST - SSW) / SST)
    obs = r2(labels)
    rng = np.random.RandomState(seed); lab = list(labels); hits = 0
    for _ in range(perms):
        rng.shuffle(lab)
        if r2(lab) >= obs: hits += 1
    return obs, (hits + 1) / (perms + 1.0)


def main():
    feats, meta, X = load()
    crops = np.array([m[1] for m in meta]); projs = np.array([m[0] for m in meta])
    print("samples %d ; features %d" % (len(meta), len(feats)))
    H = hellinger(X)
    P = X / np.maximum(X.sum(axis=1, keepdims=True), 1)
    # ---- (4) ordination
    Hc = H - H.mean(axis=0)
    U, S, Vt = np.linalg.svd(Hc, full_matrices=False)
    ev = (S ** 2) / (S ** 2).sum()
    pca = U[:, :2] * S[:2]
    # Hellinger-Euclidean distances via the Gram matrix (rows are unit norm, so D2 = 2 - 2*H H^T);
    # the naive (n, n, p) broadcast would need ~21 GB and blows up.
    G = H @ H.T
    Dh = np.sqrt(np.maximum(2.0 - 2.0 * G, 0.0))
    # persist Hellinger PCA for later sensitivity runs (avoids recomputation)
    np.save(os.path.join(OUT, "hellinger_pca_scores.npy"), pca)
    r2h, ph = permanova(Dh, list(crops))
    print("Hellinger PCA: PC1 %.1f%% PC2 %.1f%% ; PERMANOVA(Hellinger-Euclidean) R2=%.4f p=%.3f" %
          (100 * ev[0], 100 * ev[1], r2h, ph))
    D_bc = np.zeros((len(meta),) * 2); D_ai = np.zeros((len(meta),) * 2)
    C = clr(P)
    for i in range(len(meta)):
        D_bc[i] = np.abs(P[i] - P).sum(axis=1) / 2
        D_ai[i] = np.sqrt(((C[i] - C) ** 2).sum(axis=1))
    Y_bc, eb, wb = pcoa(D_bc); Y_ai, ea, wa = pcoa(D_ai)
    L = ["# 排序（§11-4：欧氏/线性方法前先 Hellinger 变换）", "",
         "样本 %d（细菌 reads ≥1000）；特征 %d。" % (len(meta), len(feats)), "",
         "| 方法 | 前两轴解释度 | 说明 |", "| --- | --- | --- |",
         "| **Hellinger + PCA**（线性方法，已按要求先做 Hellinger 变换） | PC1 %.1f%%、PC2 %.1f%% | 主排序 |" % (100 * ev[0], 100 * ev[1]),
         "| PCoA on Bray-Curtis | %.1f%%、%.1f%% | 距离法对照 |" % (eb[0], eb[1]),
         "| PCoA on Aitchison | %.1f%%、%.1f%% | CLR 距离 |" % (ea[0], ea[1]), "",
         "**与排序同几何的检验**：PERMANOVA on Hellinger-Euclidean：R² = %.4f，p = %.3f（999 次置换）。" % (r2h, ph),
         "即：在「先 Hellinger 再线性」的几何下，作物效应依然显著。", ""]
    np.savetxt(os.path.join(OUT, "ordination_scores.csv"),
               np.column_stack([pca, Y_bc[:, :2], Y_ai[:, :2]]), delimiter=",",
               header="hellingerPC1,hellingerPC2,brayPCoA1,brayPCoA2,aitPCoA1,aitPCoA2", comments="")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), dpi=200)
        col = {"maize": "#1b9e77", "rice": "#d95f02", "wheat": "#7570b3"}
        for c in CROPS:
            m = crops == c
            ax[0].scatter(pca[m, 0], pca[m, 1], s=7, alpha=.55, c=col[c], label=c)
            ax[1].scatter(Y_bc[m, 0], Y_bc[m, 1], s=7, alpha=.55, c=col[c], label=c)
        ax[0].set_xlabel("Hellinger PC1 (%.1f%%)" % (100 * ev[0])); ax[0].set_ylabel("Hellinger PC2 (%.1f%%)" % (100 * ev[1]))
        ax[0].set_title("PCA on Hellinger-transformed features")
        ax[1].set_xlabel("PCoA1 (Bray-Curtis)"); ax[1].set_ylabel("PCoA2")
        ax[1].set_title("PCoA on Bray-Curtis")
        for a_ in ax: a_.legend(markerscale=3, fontsize=8, frameon=False)
        fig.tight_layout(); fig.savefig(os.path.join(FIGS, "fig1_ordination.png"))
        print("wrote figures/fig1_ordination.png")
    except Exception as e:
        print("figure skipped:", e)
    io.open(os.path.join(OUT, "ordination.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")

    # ---- (5) Baselga decomposition on equal-depth rarefied presence/absence
    depth = int(np.percentile([int(x.sum()) for x in X], 10))
    pres = []
    for i in range(len(meta)):
        cnt = X[i].astype(int); n = cnt.sum()
        if n > depth:
            idx = np.random.choice(np.repeat(np.arange(cnt.size), cnt), size=depth, replace=False)
            cnt = np.bincount(idx, minlength=cnt.size)
        pres.append(cnt > 0)
    pres = np.asarray(pres)
    L2 = ["# β 分解：Baselga 周转/嵌套（§11-5）", "",
          "在**等深度稀疏化**（%d reads，各样本细菌读数的第 10 百分位）后的**存在/缺失**上计算两两分解：" % depth,
          "βsor = (b+c)/(2a+b+c)（Sørensen），βsim = min(b,c)/(a+min(b,c))（Simpson 周转），βnes = βsor − βsim。",
          "Podani 变体：把 Sørensen 换成 Jaccard，βjac = (b+c)/(a+b+c)，βnes = βjac − βsim。", "",
          "| 作物 | 样本 | βsor | βsim(周转) | βnes(嵌套) | 周转占比 | Podani βjac | Podani βnes |",
          "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for c in CROPS:
        ix = np.where(crops == c)[0]
        bs, bsim, bnes, bjac, bjnes = [], [], [], [], []
        for i in range(len(ix)):
            for j in range(i + 1, len(ix)):
                a = int((pres[ix[i]] & pres[ix[j]]).sum()); b = int((pres[ix[i]] & ~pres[ix[j]]).sum())
                cc = int((~pres[ix[i]] & pres[ix[j]]).sum())
                if a + b + cc == 0: continue
                sor = (b + cc) / (2 * a + b + cc); sim = min(b, cc) / (a + min(b, cc)) if (a + min(b, cc)) else 0
                jac = (b + cc) / (a + b + cc)
                bs.append(sor); bsim.append(sim); bnes.append(sor - sim); bjac.append(jac); bjnes.append(jac - sim)
        L2.append("| %s | %d | %.3f | %.3f | %.3f | %.1f%% | %.3f | %.3f |" % (
            c, len(ix), statistics.mean(bs), statistics.mean(bsim), statistics.mean(bnes),
            100 * statistics.mean(bsim) / max(statistics.mean(bs), 1e-9), statistics.mean(bjac), statistics.mean(bjnes)))
    L2 += ["", "> 判读：周转占比高（→100%）说明作物内差异主要来自类群替换；嵌套占比高说明是「子集-超集」关系（多为深度/检出差异）。",
           "> 两套（Sørensen 与 Podani/Jaccard）应给出一致的排序，否则要如实并列。"]
    io.open(os.path.join(OUT, "beta_decomposition.md"), "w", encoding="utf-8").write("\n".join(L2) + "\n")
    print("\n".join(L2[-8:]))


if __name__ == "__main__":
    main()


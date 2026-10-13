# -*- coding: utf-8 -*-
"""Remaining figures (fig1 ordination and fig2 calibration are produced by their own scripts).

fig3 core-edge: min crop core-share vs mean relative abundance, near-core highlighted, threshold grid inset
fig4 functional layer: project-level boxplots of the focus FAPROTAX guilds by crop
fig5 variance decomposition: study-level variance fractions (crop / shared / technical / residual) plus the
     sample-level project-vs-crop R2 for context.

usage: python figures_extra.py -> out_analysis/figs/fig3_core_edge.png, fig4_functional.png, fig5_variance.png
"""
import csv, io, os, statistics, collections, itertools, random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis"); FIGS = os.path.join(OUT, "figs")
os.makedirs(FIGS, exist_ok=True)
CROPS = ("maize", "rice", "wheat")
COL = {"maize": "#1b9e77", "rice": "#d95f02", "wheat": "#7570b3"}
META = ["project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
        "genus_reads","family_reads","unassigned_reads","n_asv"]
random.seed(20261012)


def load_fm():
    rd = csv.reader(io.open(os.path.join(OUT, "feature_matrix.csv"), encoding="utf-8"))
    hdr = next(rd); feats = [c for c in hdr if c not in META]; gi = {c: i for i, c in enumerate(hdr)}
    meta, X = [], []
    for r in rd:
        if not r or int(float(r[gi["bacteria_reads"]] or 0)) < 1000:
            continue
        meta.append((r[gi["project"]], r[gi["crop"]], r[gi["sample"]]))
        X.append([float(r[gi[f]] or 0) for f in feats])
    return feats, meta, np.asarray(X, float)


def main():
    feats, meta, X = load_fm()
    proj = [m[0] for m in meta]; crop = [m[1] for m in meta]
    # ---------- fig3: core-edge ----------
    ce = list(csv.DictReader(io.open(os.path.join(OUT, "core_edge.csv"), encoding="utf-8")))
    # core-ness from the project-level shares: recompute quickly (min across crops of the share of projects
    # in which the feature reaches >=70% within-project prevalence)
    byp = collections.defaultdict(lambda: collections.defaultdict(list))
    for i, (p, c, s) in enumerate(meta):
        byp[c][p].append(X[i])
    feats_j = {f: j for j, f in enumerate(feats)}
    coreness, abund = {}, {}
    for f, j in feats_j.items():
        shares = []
        for c in CROPS:
            n = 0
            for p, rows in byp[c].items():
                M = np.asarray(rows)
                if (M[:, j] > 0).mean() * 100 >= 70:
                    n += 1
            shares.append(100.0 * n / max(1, len(byp[c])))
        coreness[f] = min(shares)
        P = X[:, j] / np.maximum(X.sum(axis=1), 1)
        abund[f] = 100.0 * P.mean()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), dpi=200)
    fs = list(feats_j)
    xs = [abund[f] for f in fs]; ys = [coreness[f] for f in fs]
    kind = ["genus" if f.startswith("g__") else ("family" if f.startswith("f__") else "unassigned") for f in fs]
    for k, col, lab in (("genus", "#377eb8", "genus"), ("family", "#e41a1c", "family fallback"), ("unassigned", "#4daf4a", "unassigned")):
        m = [i for i, q in enumerate(kind) if q == k]
        ax[0].scatter([xs[i] for i in m], [ys[i] for i in m], s=8, alpha=.5, c=col, label=lab)
    ax[0].axhline(70, ls="--", c="grey", lw=1); ax[0].set_xscale("symlog", linthresh=0.01)
    ax[0].set_xlabel("mean relative abundance (%)"); ax[0].set_ylabel("min crop core-share (%)")
    ax[0].set_title("Core-edge structure (project-level core = prevalence ≥70%)")
    ax[0].legend(fontsize=8, frameon=False, markerscale=2)
    grid = [[27, 3, 1], [6, 1, 1], [1, 1, 0]]
    im = ax[1].imshow(grid, cmap="viridis", origin="lower")
    for i in range(3):
        for j in range(3):
            ax[1].text(j, i, grid[i][j], ha="center", va="center", color="w", fontsize=9)
    ax[1].set_xticks(range(3)); ax[1].set_xticklabels(["≥50%", "≥70%", "≥90%"])
    ax[1].set_yticks(range(3)); ax[1].set_yticklabels(["50%", "70%", "90%"])
    ax[1].set_xlabel("share of projects that must be core"); ax[1].set_ylabel("within-project prevalence threshold")
    ax[1].set_title("Universal-core feature count")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "fig3_core_edge.png")); plt.close(fig)
    # ---------- fig4: functional ----------
    rd = csv.reader(io.open(os.path.join(OUT, "faprotax_functions.tsv"), encoding="utf-8"), delimiter="\t")
    h = next(rd); cols = h[2:]
    idx = {s: i for i, s in enumerate(cols)}
    FOCUS = ["methanotrophy", "methylotrophy", "aerobic_chemoheterotrophy", "nitrate_respiration", "fermentation"]
    vals = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rd:
        g = r[0].strip()
        if g not in FOCUS: continue
        for i, (p, c, s) in enumerate(meta):
            j = idx.get("%s|%s" % (p, s))
            if j is None: continue
            vals[g][(c, p)].append(100.0 * float(r[2 + j] or 0) / max(1.0, X[i].sum()))
    fig, axes = plt.subplots(1, len(FOCUS), figsize=(13, 3.8), dpi=200)
    for ax_, g in zip(axes, FOCUS):
        data = [[statistics.mean(v) for (c, p), v in vals[g].items() if c == cc] for cc in CROPS]
        bp = ax_.boxplot(data, tick_labels=list(CROPS), patch_artist=True, showfliers=True,
                         flierprops=dict(marker="o", ms=2, alpha=.4))
        ax_.set_xticks(range(1, 1 + len(CROPS)))
        ax_.set_xticklabels(list(CROPS))
        for patch, cc in zip(bp["boxes"], CROPS): patch.set_facecolor(COL[cc]); patch.set_alpha(.6)
        ax_.set_title(g, fontsize=8); ax_.tick_params(labelsize=7)
        ax_.set_ylabel("% of bacterial reads" if g == FOCUS[0] else "")
    fig.suptitle("FAPROTAX functional groups by crop (project-level means; n=15 per crop)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "fig4_functional.png")); plt.close(fig)
    # ---------- fig5: variance decomposition ----------
    vd = {"crop only": [0.058, 0.080], "shared": [0.157, 0.080], "technical only": [0.033, 0.016], "residual": [0.752, 0.824]}
    labels = ["Bray-Curtis", "Aitchison"]
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 4.0), dpi=200)
    bottom = np.zeros(2)
    for k, col in (("crop only", "#d95f02"), ("shared", "#7570b3"), ("technical only", "#1b9e77"), ("residual", "#999999")):
        ax[0].bar(labels, vd[k], bottom=bottom, label=k, color=col)
        bottom += np.array(vd[k])
    ax[0].set_ylabel("adjusted R² fraction"); ax[0].set_title("Study-level variance partitioning (45 projects)")
    ax[0].legend(fontsize=7, frameon=False)
    ax[1].bar(["crop", "project"], [0.0777, 0.6246], color=["#d95f02", "#333333"])
    for i, v in enumerate([0.0777, 0.6246]): ax[1].text(i, v + .01, "%.3f" % v, ha="center", fontsize=8)
    ax[1].set_ylabel("R²"); ax[1].set_title("Sample level: crop vs study (Bray-Curtis)")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "fig5_variance.png")); plt.close(fig)
    print("wrote figs: fig3_core_edge.png, fig4_functional.png, fig5_variance.png")
    for f in sorted(os.listdir(FIGS)):
        print("   ", f, os.path.getsize(os.path.join(FIGS, f)) // 1024, "KB")


if __name__ == "__main__":
    main()

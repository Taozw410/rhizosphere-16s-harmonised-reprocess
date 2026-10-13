# -*- coding: utf-8 -*-
"""Section 11 item 7, completed: the ML block with everything the pre-registration asks for.

Adds to the earlier run:
  (a) apparent (in-sample) AUC vs out-of-fold AUC -> the optimism gap;
  (b) probability calibration table + Brier score + calibration figure;
  (c) shortcut check WITH the amplicon-region covariate (not just depth/ASV/sample counts);
  (d) REGION-STRATIFIED leave-one-study-out (train and test inside one region class only) - the decisive
      check that the classifier is not reading the amplicon region;
  (e) leave-one-CLUSTER-out using the section-2 platform clusters (grouping the paired studies);
  plus the two pre-registered sensitivities (agreed genera; drop PRJNA1261944).
usage: python ml_leaveout2.py -> out_analysis/ml_report2.md, out_analysis/figs/fig2_calibration.png
"""
import csv, io, os, collections
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis"); FIGS = os.path.join(OUT, "figs")
os.makedirs(FIGS, exist_ok=True)
META = ["project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
        "genus_reads","family_reads","unassigned_reads","n_asv"]
CROPS = ("maize", "rice", "wheat")
rng = np.random.RandomState(20261012)


def load():
    rd = csv.reader(io.open(os.path.join(OUT, "feature_matrix.csv"), encoding="utf-8"))
    hdr = next(rd); feats = [c for c in hdr if c not in META]; gi = {c: i for i, c in enumerate(hdr)}
    meta, X = [], []
    for r in rd:
        if not r or int(float(r[gi["bacteria_reads"]] or 0)) < 1000:
            continue
        meta.append(dict(project=r[gi["project"]], crop=r[gi["crop"]], sample=r[gi["sample"]],
                         depth=float(r[gi["bacteria_reads"]] or 0), n_asv=float(r[gi["n_asv"]] or 0),
                         unassigned=float(r[gi["unassigned_reads"]] or 0)))
        X.append([float(r[gi[f]] or 0) for f in feats])
    return feats, meta, np.asarray(X, float)


def clr(P):
    P = np.where(P <= 0, 1e-9, P); P = P / P.sum(axis=1, keepdims=True)
    L = np.log(P); return L - L.mean(axis=1, keepdims=True)


def fit_predict(Xtr, ytr, Xte):
    m = LogisticRegression(max_iter=2000, C=1.0)
    m.fit(Xtr, ytr)
    return m.predict_proba(Xte), m.classes_


def loo(X, y, groups, subset=None):
    """leave-one-group-out; returns pooled probabilities (n x 3) in the fixed crop order, plus apparent probs."""
    classes = list(CROPS)
    yi = np.array([classes.index(v) for v in y])
    oof = np.zeros((len(y), 3)); app = np.zeros((len(y), 3))
    idx_all = np.arange(len(y)) if subset is None else np.asarray(subset)
    for g in sorted(set(np.asarray(groups)[idx_all])):
        te = np.array([(q == g) and (subset is None or i in set(idx_all.tolist()))
                       for i, q in enumerate(groups)])
        tr = np.array([(q != g) and (subset is None or i in set(idx_all.tolist()))
                       for i, q in enumerate(groups)])
        if te.sum() == 0 or len(set(yi[tr])) < 2:
            continue
        pr, cls = fit_predict(X[tr], yi[tr], X[te])
        for j, c in enumerate(cls):
            oof[te, c] = pr[:, j]
        pra, _ = fit_predict(X[tr], yi[tr], X[tr])
        for j, c in enumerate(cls):
            app[tr, c] = pra[:, j]
    keep = np.ones(len(y), bool) if subset is None else np.zeros(len(y), bool); keep[idx_all] = True
    auc = roc_auc_score(yi[keep], oof[keep], multi_class="ovr", average="macro")
    ap = roc_auc_score(yi[app.sum(1) > 0], app[app.sum(1) > 0], multi_class="ovr", average="macro")
    return auc, oof, yi, keep, ap


def boot(yi, groups, oof, keep, reps=30):
    ps = sorted(set(np.asarray(groups)[keep])); vals = []
    for _ in range(reps):
        pick = rng.choice(ps, size=len(ps), replace=True)
        ix = np.concatenate([np.where((np.asarray(groups) == p) & keep)[0] for p in pick])
        try:
            vals.append(roc_auc_score(yi[ix], oof[ix], multi_class="ovr", average="macro"))
        except Exception:
            pass
    return (float(np.percentile(vals, 5)), float(np.percentile(vals, 95)), len(vals)) if vals else (None, None, 0)


def main():
    feats, meta, X = load()
    proj = [m["project"] for m in meta]; y = [m["crop"] for m in meta]
    cl = {r["project"]: r["platform_cluster"] for r in csv.DictReader(io.open(os.path.join(OUT, "platform_cluster.csv"), encoding="utf-8"))}
    cov = {r["project"]: r for r in csv.DictReader(io.open(os.path.join(OUT, "study_covariates.csv"), encoding="utf-8"))}
    region = {p: ("V4-like" if float(c["asv_len_mode"]) <= 320 else "V3V4-like") for p, c in cov.items()}
    agree = {r["genus"] for r in csv.DictReader(io.open(os.path.join(OUT, "consistent_genera.csv"), encoding="utf-8-sig"))
             if r["weight_class"] == "consistent"}
    keep_j = [j for j in range(len(feats)) if X[:, j].sum() > 0]
    R_ = X / np.maximum(X.sum(axis=1, keepdims=True), 1)
    Xc = clr(R_[:, keep_j]); feats_k = [feats[j] for j in keep_j]
    L = ["# 留一研究 ML 判别（§11-7 完整版）", "",
         "样本 %d；折 = 留出一个完整项目（或一个 platform_cluster）；特征 = §17 特征空间 CLR；分类器 = L2 多项逻辑回归。" % len(meta), ""]
    # 1 main
    auc, oof, yi, keep, app_auc = loo(Xc, y, proj)
    lo, hi, n = boot(yi, proj, oof, keep)
    per = {c: roc_auc_score((yi == i).astype(int), oof[:, i]) for i, c in enumerate(CROPS)}
    L += ["## 1 主口径（留一项目）", "", "| 指标 | 值 |", "| --- | --- |",
          "| 折外合并 AUC（ovr macro） | **%.3f** |" % auc,
          "| 项目 bootstrap 5–95%% | %.3f–%.3f（%d 次） |" % (lo, hi, n),
          "| 每类 AUC（玉米/水稻/小麦） | %.3f / %.3f / %.3f |" % (per["maize"], per["rice"], per["wheat"]),
          "| **表观（训练内）AUC** | %.3f |" % app_auc,
          "| **乐观度 = 表观 − 折外** | **%.3f** |" % (app_auc - auc), ""]
    # 2 calibration + Brier
    bins = np.linspace(0, 1, 6)
    rows = []
    for i, c in enumerate(CROPS):
        ybin = (yi == i).astype(int); p = oof[keep, i]
        for b in range(len(bins) - 1):
            m = (p >= bins[b]) & (p < bins[b + 1])
            if m.sum() >= 10:
                rows.append((c, bins[b], bins[b + 1], m.sum(), p[m].mean(), ybin[keep][m].mean()))
    brier = np.mean([brier_score_loss((yi[keep] == i).astype(int), oof[keep, i]) for i in range(3)])
    L += ["## 2 概率校准", "",
          "| 类别 | 预测区间 | n | 平均预测概率 | 实际频率 |", "| --- | --- | --- | --- | --- |"]
    for c, a, b, nn, mp, mf in rows:
        L.append("| %s | %.1f–%.1f | %d | %.3f | %.3f |" % (c, a, b, nn, mp, mf))
    L += ["", "多分类 Brier 分数（三类平均，越小越准）= **%.3f**；与校准表一并看：预测概率与实际频率越贴近越好。" % brier, ""]
    # 3 shortcut with region
    T = np.column_stack([[m["depth"] for m in meta], [m["n_asv"] for m in meta],
                         [m["unassigned"] / max(1.0, m["depth"]) for m in meta],
                         [sum(1 for q in proj if q == m["project"]) for m in meta],
                         [1.0 if region.get(m["project"]) == "V3V4-like" else 0.0 for m in meta]])
    T = np.log1p(T); T = (T - T.mean(0)) / (T.std(0) + 1e-9)
    auc_s, oof_s, yi_s, keep_s, _ = loo(T, y, proj)
    L += ["## 3 捷径检验（技术协变量，**含扩增子区域**）", "",
          "协变量 = log 深度、log ASV 数、未定 reads 比例、项目样本数、扩增子区域类（V4 类 vs V3–V4 类）。",
          "**只用这些技术协变量**判别作物的折外 AUC = **%.3f**（≈0.5 表示无捷径可利用）。" % auc_s, ""]
    # 4 region-stratified LOSO
    L += ["## 4 区域分层留一（在同一区域类内部判别）", "",
          "| 区域类 | 项目数（玉米/水稻/小麦） | 折外 AUC |", "| --- | --- | --- |"]
    for reg in ("V3V4-like", "V4-like"):
        sub = [i for i, m in enumerate(meta) if region.get(m["project"]) == reg]
        cnt = collections.Counter(y[i] for i in sub)
        if len(set(y[i] for i in sub)) < 3 or min(cnt.values()) < 3:
            L.append("| %s | %d/%d/%d | 不足 3 类或某类 <3（跳过） |" % (reg, cnt["maize"], cnt["rice"], cnt["wheat"]))
            continue
        a2, o2, y2, k2, _ = loo(Xc, y, proj, subset=sub)
        L.append("| %s | %d/%d/%d | **%.3f** |" % (reg, cnt["maize"], cnt["rice"], cnt["wheat"], a2))
    # 5 leave-one-cluster-out
    cls_groups = [cl.get(p, p) for p in proj]
    a3, o3, y3, k3, _ = loo(Xc, y, cls_groups)
    L += ["", "## 5 按 platform_cluster 分组的留一（§2 集群）", "",
          "折 = 留出一个集群（6 个 2 成员集群 + 其余单例），共 %d 折；折外 AUC = **%.3f**。" % (len(set(cls_groups)), a3), ""]
    # 6 sensitivities
    idx_ag = [j for j, f in enumerate(feats_k) if f.startswith("g__") and f[3:] in agree]
    a4, *_ = loo(Xc[:, idx_ag], y, proj)
    sub = [i for i, m in enumerate(meta) if m["project"] != "PRJNA1261944"]
    a5, *_ = loo(Xc, y, proj, subset=sub)
    L += ["## 6 预注册敏感性", "", "| 口径 | 折外 AUC |", "| --- | --- |",
          "| 只用两法一致的属（§3） | %.3f |" % a4,
          "| 剔除 PRJNA1261944（§18） | %.3f |" % a5]
    io.open(os.path.join(OUT, "ml_report2.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.6, 4.4), dpi=200)
        for c in CROPS:
            xs = [r[4] for r in rows if r[0] == c]; ys = [r[5] for r in rows if r[0] == c]
            ax.plot(xs, ys, "o-", label=c, ms=4)
        ax.plot([0, 1], [0, 1], "--", color="grey", lw=1)
        ax.set_xlabel("predicted probability"); ax.set_ylabel("observed frequency")
        ax.set_title("Calibration (leave-one-study-out)"); ax.legend(fontsize=8, frameon=False)
        fig.tight_layout(); fig.savefig(os.path.join(FIGS, "fig2_calibration.png"))
    except Exception as e:
        print("figure skipped:", e)
    print("wrote ml_report2.md ; main AUC=%.3f (apparent %.3f) ; shortcut(with region)=%.3f ; cluster-LOO=%.3f" %
          (auc, app_auc, auc_s, a3))
    for x in L:
        if x.startswith("|") and ("AUC" in x or "乐观" in x or "Brier" in x or "剔除" in x or "一致" in x):
            print(x)


if __name__ == "__main__":
    main()

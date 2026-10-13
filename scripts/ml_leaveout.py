# -*- coding: utf-8 -*-
"""Leave-one-STUDY-out ML discrimination (section 11 binding requirement) + shortcut check + sensitivities.

- Features: CLR-transformed section-17 features (zero replacement inside run, documented).
- Split: leave one complete PROJECT out (45 folds); no sample from a held-out study is ever in training.
- Readout: pooled out-of-fold probabilities -> multiclass AUC (ovr macro) + per-crop AUC,
  uncertainty by bootstrapping PROJECTS (30 resamples).
- Sensitivities: (a) drop PRJNA1261944 (section 18); (b) use only the 483 agreed genera (section 3);
  (c) SHORTCUT CHECK: predict crop from technical covariates only (read depth, ASV count, sample count of
  the project, quality chars) - if that alone discriminates, crop signal may be confounded by study-level
  technical differences.

usage: python ml_leaveout.py
"""
import csv, io, os, collections, statistics
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
META = ["project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
        "genus_reads","family_reads","unassigned_reads","n_asv"]
CROPS = ("maize", "rice", "wheat")


def load():
    rd = csv.reader(io.open(os.path.join(OUT, "feature_matrix.csv"), encoding="utf-8"))
    hdr = next(rd); feats = [c for c in hdr if c not in META]; gi = {c: i for i, c in enumerate(hdr)}
    meta, X = [], []
    for r in rd:
        if not r or int(float(r[gi["bacteria_reads"]] or 0)) < 1000:
            continue
        meta.append(dict(project=r[gi["project"]], crop=r[gi["crop"]], sample=r[gi["sample"]],
                         depth=float(r[gi["bacteria_reads"]]), n_asv=float(r[gi["n_asv"]]),
                         unassigned=float(r[gi["unassigned_reads"]])))
        X.append([float(r[gi[f]] or 0) for f in feats])
    return feats, meta, np.asarray(X, float)


def clr(P):
    P = np.asarray(P, float)
    P = np.where(P <= 0, 1e-9, P)
    P = P / P.sum(axis=1, keepdims=True)
    L = np.log(P)
    return L - L.mean(axis=1, keepdims=True)


def loo_auc(Xf, y, proj, seed=0):
    """leave-one-project-out; returns pooled AUC + per-class AUC + out-of-fold probabilities."""
    y = np.asarray(y); proj = np.asarray(proj)
    classes = sorted(set(y))
    oof = np.zeros((Xf.shape[0], len(classes)))
    for p in sorted(set(proj)):
        te = np.array([q == p for q in proj]); tr = ~te
        if len(set(y[tr])) < 2 or te.sum() == 0:
            continue
        # sklearn >= 1.5 dropped the multi_class argument (multinomial is the default for multi-class)
        m = LogisticRegression(max_iter=2000, C=1.0)
        m.fit(Xf[tr], y[tr])
        pr = m.predict_proba(Xf[te])
        # align class columns
        for j, c in enumerate(m.classes_):
            oof[te, classes.index(c)] = pr[:, j]
    yi = np.array([classes.index(v) for v in y])
    auc = roc_auc_score(yi, oof, multi_class="ovr", average="macro")
    per = {}
    for c in classes:
        yy = (yi == classes.index(c)).astype(int)
        if 0 < yy.sum() < len(yy):
            per[c] = roc_auc_score(yy, oof[:, classes.index(c)])
    return auc, per, oof, classes


def boot_projects(y, proj, oof, classes, reps=30, seed=1):
    rnd = np.random.RandomState(seed)
    ps = sorted(set(proj)); yi = np.array([classes.index(v) for v in y]); vals = []
    for _ in range(reps):
        pick = rnd.choice(ps, size=len(ps), replace=True)
        idx = np.concatenate([np.where(np.array(proj) == p)[0] for p in pick])
        try:
            vals.append(roc_auc_score(yi[idx], oof[idx], multi_class="ovr", average="macro"))
        except Exception:
            pass
    return (np.percentile(vals, 5), np.percentile(vals, 95), len(vals)) if vals else (None, None, 0)


def main():
    feats, meta, X = load()
    proj = [m["project"] for m in meta]; y = [m["crop"] for m in meta]
    print("samples %d ; projects %d ; features %d" % (len(meta), len(set(proj)), len(feats)))
    agree = set()
    for r in csv.DictReader(io.open(os.path.join(OUT, "consistent_genera.csv"), encoding="utf-8-sig")):
        if r["weight_class"] == "consistent":
            agree.add(r["genus"])
    # relative abundance -> CLR, dropping all-zero features
    keep_j = [j for j in range(len(feats)) if X[:, j].sum() > 0]
    R_ = X / np.maximum(X.sum(axis=1, keepdims=True), 1)
    Xc = clr(R_[:, keep_j])
    feats_k = [feats[j] for j in keep_j]
    lines = []
    def run(tag, Xf, yy=y, pp=proj):
        auc, per, oof, classes = loo_auc(Xf, yy, pp)
        lo, hi, n = boot_projects(yy, pp, oof, classes)
        lines.append((tag, Xf.shape[1], auc, lo, hi, per))
        print("%-42s features=%-5d AUC=%.3f (%s-%s) per-class=%s" %
              (tag, Xf.shape[1], auc, ("%.3f" % lo) if lo else "-", ("%.3f" % hi) if hi else "-",
               {k: round(v, 3) for k, v in per.items()}))
    run("main (all features)", Xc)
    idx = [j for j, f in enumerate(feats_k) if f.startswith("g__") and f[3:] in agree]
    run("only agreed genera (section 3)", Xc[:, idx])
    keep_s = [i for i, p in enumerate(proj) if p != "PRJNA1261944"]
    run("without PRJNA1261944 (section 18)", Xc[keep_s], [y[i] for i in keep_s], [proj[i] for i in keep_s])
    # shortcut check: technical covariates only
    T = np.column_stack([[m["depth"] for m in meta], [m["n_asv"] for m in meta],
                         [m["unassigned"] / max(1.0, m["depth"]) for m in meta]])
    proj_n = collections.Counter(proj)
    T = np.column_stack([T, [proj_n[m["project"]] for m in meta]])
    T = np.log1p(T)
    run("shortcut: technical covariates only", (T - T.mean(0)) / (T.std(0) + 1e-9))
    rep = ["# 留一研究 ML 判别（§11 绑定要求）", "",
           "样本 %d（≥1000 细菌 reads）；折 = 留出一个完整项目（45 折，测试集从不进训练）。" % len(meta),
           "特征 = §17 特征空间的 CLR；分类器 = 多项逻辑回归（L2）。读数 = 折外概率合并后的多分类 AUC（ovr macro），",
           "不确定性 = 按**项目** bootstrap 30 次。", "",
           "| 口径 | 特征数 | 合并 AUC | AUC 的 5–95%（项目 bootstrap） | 每类 AUC（玉米/水稻/小麦） |",
           "| --- | --- | --- | --- | --- |"]
    for tag, nf, auc, lo, hi, per in lines:
        rep.append("| %s | %d | %.3f | %s | %s |" % (
            tag, nf, auc, ("%.3f–%.3f" % (lo, hi)) if lo else "—",
            " / ".join("%.3f" % per.get(c, float("nan")) for c in CROPS)))
    rep += ["", "> 判读：主口径 AUC 与其 5–95% 下界一起读；捷径检验（只用技术协变量）若也明显 >0.5，",
            "> 说明作物信号里混着研究层技术差异，必须如实写入限制；§3 与 §18 两条敏感性给出结论的边界。"]
    io.open(os.path.join(OUT, "ml_report.md"), "w", encoding="utf-8").write("\n".join(rep) + "\n")
    print("wrote out_analysis/ml_report.md")


if __name__ == "__main__":
    main()

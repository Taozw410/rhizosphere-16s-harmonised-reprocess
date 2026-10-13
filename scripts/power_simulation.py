# -*- coding: utf-8 -*-
"""Project-specific power / uncertainty simulation for the final paper's claims.

Builds a genus-level table from every finished project we already have, then:
  A) PERMANOVA power: how large a crop-driven compositional shift must be to be detected with
     k = 10 vs 15 projects per crop (Bray-Curtis, 199 permutations, real within-group dispersion).
  B) Leave-one-STUDY-out AUC and its uncertainty: what discrimination can be claimed when only
     ~10 projects per crop exist (logistic regression on CLR-transformed top genera).

Outputs: out50k/genus_table_combined.csv, out50k/power_simulation_report.md
"""
import csv, io, os, re, glob, math, random, statistics, collections, datetime, sys
import numpy as np

R = r"D:\CodexWork\mic-re"
RE_PRJ = re.compile(r"^PRJ[A-Z]*\d+$")
random.seed(20261009); np.random.seed(20261009)


def rd(p):
    return list(csv.DictReader(io.open(p, encoding="utf-8-sig"))) if os.path.exists(p) else []


def crop_map():
    m = {}
    for f in glob.glob(os.path.join(R, "meta", "*.csv")):
        for r in rd(f):
            k = (r.get("study") or r.get("project") or "").strip()
            c = (r.get("crop") or "").strip().lower()
            if k and c in ("maize", "rice", "wheat"):
                m.setdefault(k, c)
    return m


def genus_table_for(pdir):
    """-> (samples, {genus: [counts...]}) or None"""
    at, tt = os.path.join(pdir, "asv_table.csv"), os.path.join(pdir, "asv_taxonomy.csv")
    if not (os.path.exists(at) and os.path.exists(tt)):
        return None
    tax = {}
    with io.open(tt, encoding="utf-8-sig") as fh:
        r = csv.reader(fh); next(r, None)
        for row in r:
            if len(row) < 7: continue
            seq = row[0].strip().strip('"')
            g = (row[6] or "").strip()
            tax[seq] = g if g and g.upper() != "NA" else "unclassified"
    if not tax:
        return None
    samples = None
    agg = collections.defaultdict(list)
    with io.open(at, encoding="utf-8-sig") as fh:
        r = csv.reader(fh)
        hdr = next(r)
        samples = [x.strip().strip('"') for x in hdr[1:]]
        for row in r:
            if len(row) < 2: continue
            seq = row[0].strip().strip('"')
            g = tax.get(seq, "unclassified")
            vals = agg[g]
            if not vals:
                vals.extend([0.0] * len(samples))
            for i, x in enumerate(row[1:len(samples) + 1]):
                try: vals[i] += float(x or 0)
                except ValueError: pass
    return samples, agg


def build_table(only=None):
    cm = crop_map()
    meta_rows = []
    genera = collections.Counter()
    per_project = {}
    for base in ("out_stage2", "out50k"):
        d0 = os.path.join(R, base)
        if not os.path.isdir(d0): continue
        for p in sorted(os.listdir(d0)):
            if not RE_PRJ.match(p) or p in per_project: continue
            if only is not None and p not in only: continue   # 2026-10-12: 支持只用最终入选项目重跑（§8.3）
            crop = cm.get(p)
            if not crop: continue
            gt = genus_table_for(os.path.join(d0, p))
            if not gt: continue
            samples, agg = gt
            tot = np.zeros(len(samples))
            for g, v in agg.items():
                tot += np.array(v)
            keep = tot > 0
            # 2026-10-12: 原来是 <5 就剔除，结果把只有 4 个样本的 PRJNA753114 挡在门外，
            # 使"最终 45 个入选项目"的模拟只剩 44 个（玉米 14）。重采样本身用 min(12, n)，
            # 不要求 ≥5，所以门槛降到 4，让模拟与最终选集完全一致（15/15/15）。
            if keep.sum() < 4: continue
            per_project[p] = dict(crop=crop, samples=[s for s, k in zip(samples, keep) if k],
                                  mat={g: (np.array(v)[keep] / tot[keep] * 100.0) for g, v in agg.items()})
            for g in per_project[p]["mat"]: genera[g] += 1
    for p, d in per_project.items():
        for s in d["samples"]:
            meta_rows.append((s, p, d["crop"]))
    return per_project, meta_rows


# ---------- PERMANOVA (Bray-Curtis, pseudo-F, permutation p) ----------
def bray_curtis(X):
    n = X.shape[0]
    s = X.sum(axis=1, keepdims=True)
    Xn = X / np.where(s == 0, 1, s)
    ab = np.abs(Xn[:, None, :] - Xn[None, :, :]).sum(axis=2)
    tot = Xn.sum(axis=1)
    denom = tot[:, None] + tot[None, :]
    D = np.where(denom == 0, 0, ab / denom)
    np.fill_diagonal(D, 0.0)
    return D


def permanova_F(D2, labels, groups):
    n = D2.shape[0]
    SST = D2[np.triu_indices(n, 1)].sum() / n
    SSW = 0.0
    for g in groups:
        idx = np.where(labels == g)[0]
        ng = len(idx)
        if ng < 2: return None
        sub = D2[np.ix_(idx, idx)]
        SSW += sub[np.triu_indices(ng, 1)].sum() / ng
    a, b = len(groups) - 1, n - len(groups)
    if b <= 0: return None
    SSA = SST - SSW
    if SSW <= 1e-12: return None
    return (SSA / a) / (SSW / b), (SSA / SST if SST > 0 else 0.0)


def permanova(proj_stack, k, delta, reps, perms, signatures, groups):
    """project_stack: {crop: [(X_samples, project)]}; returns (power, mean R2 at delta=0)"""
    hits = 0; r2s = []
    pool = {c: [i for i in range(len(proj_stack[c]))] for c in groups}
    for rep in range(reps):
        Xs = []; labels = []
        for gi, c in enumerate(groups):
            idx = np.random.choice(pool[c], size=k, replace=True)
            for i in idx:
                X, _ = proj_stack[c][i]
                m = np.random.choice(len(X), size=min(12, len(X)), replace=False)
                Xs.append(X[m])
                labels += [gi] * len(m)
        X = np.vstack(Xs).astype(float)
        if delta > 0:
            Xn = X / np.where(X.sum(axis=1, keepdims=True) == 0, 1, X.sum(axis=1, keepdims=True))
            for gi, c in enumerate(groups):
                sig = signatures[c]
                Xn[labels == gi] = (1 - delta) * Xn[labels == gi] + delta * sig[:Xn.shape[1]]
            X = Xn
        D = bray_curtis(X)
        D2 = D ** 2
        lab = np.array(labels)
        obs = permanova_F(D2, lab, list(range(len(groups))))
        if obs is None: continue
        F, r2 = obs
        r2s.append(r2)
        cnt = 1
        for _ in range(perms):
            perm = np.random.permutation(lab)
            r = permanova_F(D2, perm, list(range(len(groups))))
            if r and r[0] >= F: cnt += 1
        if cnt / (perms + 1) <= 0.05: hits += 1
    return hits / reps, (statistics.mean(r2s) if r2s else float("nan"))


def main():
    only, out_path, label = None, os.path.join(R, "out50k", "power_simulation_report.md"), "用本项目已跑完数据"
    if "--selection" in sys.argv:
        sel = sys.argv[sys.argv.index("--selection") + 1]
        only = set()
        for r in csv.DictReader(io.open(sel, encoding="utf-8-sig")):
            if (r.get("rank") or "").strip():
                only.add(r["project"])
        label = "只用最终入选项目（%s）" % os.path.basename(sel)
        out_path = os.path.join(R, "out50k", "power_simulation_report_final45.md")
    if "--out" in sys.argv:
        out_path = sys.argv[sys.argv.index("--out") + 1]
    pp, meta_rows = build_table(only=only)
    crops = collections.Counter(d["crop"] for d in pp.values())
    log = ["# 功效与不确定性模拟（%s）" % label, "", "生成时间：%s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ""]
    log.append("参与建模的项目：%d（玉米 %d / 水稻 %d / 小麦 %d），样本 %d" % (
        len(pp), crops.get("maize", 0), crops.get("rice", 0), crops.get("wheat", 0), len(meta_rows)))
    # common genus space, top 120 by mean abundance for speed
    allg = collections.Counter()
    for d in pp.values():
        for g, v in d["mat"].items():
            allg[g] += float(np.mean(v))
    keepg = [g for g, _ in allg.most_common(120)]
    gidx = {g: i for i, g in enumerate(keepg)}
    groups = ["maize", "rice", "wheat"]
    stack = {c: [] for c in groups}
    for p, d in pp.items():
        n = len(d["samples"])
        X = np.zeros((n, len(keepg)))
        for g, v in d["mat"].items():
            if g in gidx: X[:, gidx[g]] = v
        stack[d["crop"]].append((X, p))
    missing = [c for c in groups if len(stack[c]) < 5]
    if missing:
        log.append("")
        log.append("> 注意：%s 可用项目不足 5 个，模拟只能作为示意。" % ", ".join(missing))
    sig = {}
    for c in groups:
        prof = np.mean([X.mean(axis=0) for X, _ in stack[c]], axis=0) if stack[c] else np.zeros(len(keepg))
        sig[c] = prof / prof.sum() if prof.sum() > 0 else prof
    # A) observed effect size + project-level bootstrap uncertainty (the honest question here is
    # "how precise is the crop effect when only k projects per crop exist", not an artificial power curve)
    log += ["", "## A. 作物效应量与项目级不确定性（Bray-Curtis PERMANOVA，99 次置换）", "",
            "做法：从三个作物的真实项目里按项目有放回重采样，每项目抽 12 个样本；对每个重采样计算「作物」因子的 R² 与置换 p 值。",
            "样本层 n≈360（k=10）或 540（k=15）；另给「项目层」（每项目取组成均值，k 个点/作物）的结果。", "",
            "| 情形 | 平均 R² | R² 的 5–95% | 检出率（p<0.05） |", "| --- | --- | --- | --- |"]

    def boot_once(k, agg_project, reps=40):
        r2s = []; hits = 0
        for rep in range(reps):
            Xs = []; labels = []
            for gi, c in enumerate(groups):
                pl = range(len(stack[c]))
                idx = np.random.choice(list(pl), size=k, replace=True)
                for i in idx:
                    X, _ = stack[c][i]
                    if agg_project:
                        Xs.append(X.mean(axis=0, keepdims=True)); labels += [gi]
                    else:
                        m = np.random.choice(len(X), size=min(12, len(X)), replace=False)
                        Xs.append(X[m]); labels += [gi] * len(m)
            X = np.vstack(Xs).astype(float)
            X = X / np.where(X.sum(axis=1, keepdims=True) == 0, 1, X.sum(axis=1, keepdims=True))
            D2 = bray_curtis(X) ** 2
            lab = np.array(labels)
            obs = permanova_F(D2, lab, list(range(len(groups))))
            if obs is None: continue
            F, r2 = obs
            r2s.append(r2)
            cnt = 1
            for _ in range(99):
                r = permanova_F(D2, np.random.permutation(lab), list(range(len(groups))))
                if r and r[0] >= F: cnt += 1
            if cnt / 100 <= 0.05: hits += 1
        return r2s, hits / max(1, len(r2s))

    obs_all = None
    Xs = []; labels = []
    for gi, c in enumerate(groups):
        for X, _ in stack[c]:
            Xs.append(X); labels += [gi] * len(X)
    X = np.vstack(Xs).astype(float)
    X = X / np.where(X.sum(axis=1, keepdims=True) == 0, 1, X.sum(axis=1, keepdims=True))
    o = permanova_F(bray_curtis(X) ** 2, np.array(labels), list(range(3)))
    if o:
        log.append("| 全部已跑完项目（样本层，参考） | %.3f | — | — |" % o[1])

    for k in (10, 15):
        r2s, hit = boot_once(k, False)
        if r2s:
            q = np.percentile(r2s, [5, 95])
            log.append("| k=%d / 作物（样本层） | %.3f | %.3f–%.3f | %.0f%% |" % (k, statistics.mean(r2s), q[0], q[1], 100 * hit))
        r2p, hitp = boot_once(k, True)
        if r2p:
            q = np.percentile(r2p, [5, 95])
            log.append("| k=%d / 作物（项目层：每项目 1 个点） | %.3f | %.3f–%.3f | %.0f%% |" % (k, statistics.mean(r2p), q[0], q[1], 100 * hitp))
    log.append("")
    log.append("> 判读：R² 的 5–95% 区间给出「换一批项目」的波动范围；项目层的检出率回答「若只看 k 个独立研究，能否稳住这个差异」。")

    io.open(out_path, "w", encoding="utf-8").write("\n".join(log) + "\n")
    print("PARTIAL (PERMANOVA) written", flush=True)
    # B) leave-one-study-out AUC
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    Xs = []; ys = []; ps = []
    for c in groups:
        for X, p in stack[c]:
            for i in range(len(X)):
                Xs.append(X[i]); ys.append(groups.index(c)); ps.append(p)
    X = np.array(Xs, dtype=float); y = np.array(ys); proj = np.array(ps)
    X = X + 1e-4
    X = np.log(X) - np.log(X).mean(axis=1, keepdims=True)
    X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-9)

    def pooled_lopo(X, y, proj, seed=0):
        """pool all out-of-fold probabilities (leave one project out), then one multiclass AUC"""
        rng = np.random.RandomState(seed)
        yt = []; yp = []
        for p in sorted(set(proj)):
            te = proj == p; tr = ~te
            if len(set(y[tr])) < 3:
                continue
            clf = LogisticRegression(max_iter=2000)
            clf.fit(X[tr], y[tr])
            yt.append(y[te]); yp.append(clf.predict_proba(X[te]))
        if not yt:
            return None, None
        yt = np.concatenate(yt); yp = np.vstack(yp)
        per = {}
        for gi, c in enumerate(groups):
            if (yt == gi).sum() and (yt != gi).sum():
                per[c] = roc_auc_score((yt == gi).astype(int), yp[:, gi])
        return roc_auc_score(yt, yp, multi_class="ovr", average="macro"), per

    auc, per = pooled_lopo(X, y, proj)
    log += ["", "## B. 留一研究（leave-one-study-out）判别 AUC", "",
            "做法：每次留出一个完整项目做测试，把**全部折外预测概率合并**后计算多分类 AUC（ovr macro）；"
            "不确定性用「按项目 bootstrap」估计（重抽样项目，重算合并 AUC）。", "",
            "| 情形 | 折数 | 合并 AUC | AUC 的 5–95% | 每类 AUC（玉米/水稻/小麦） |", "| --- | --- | --- | --- | --- |"]
    if auc is not None:
        log.append("| 全部可用项目 | %d | %.3f | — | %s |" % (
            len(set(proj)), auc,
            " / ".join("%s %.3f" % (c, per.get(c, float("nan"))) for c in groups)))
    boots = []
    projects = sorted(set(proj))
    for rep in range(30):
        rng = np.random.RandomState(100 + rep)
        pick = rng.choice(projects, size=len(projects), replace=True)
        keep = np.isin(proj, pick)
        if len(set(proj[keep])) < 5: continue
        a, _ = pooled_lopo(X[keep], y[keep], proj[keep], seed=rep)
        if a is not None: boots.append(a)
    if boots:
        q = np.percentile(boots, [5, 95])
        log.append("| 项目 bootstrap（重抽样研究） | %d 次 | %.3f | %.3f–%.3f | — |" % (
            len(boots), statistics.mean(boots), q[0], q[1]))
    # k=10 per crop variant
    means = []
    for rep in range(10):
        sel = []
        for c in groups:
            pl = [p for X_, p in stack[c]]
            if pl: sel += list(np.random.choice(pl, size=min(10, len(pl)), replace=False))
        keep = np.isin(proj, sel)
        a, _ = pooled_lopo(X[keep], y[keep], proj[keep], seed=rep)
        if a is not None: means.append(a)
    if means:
        log.append("| 每作物仅 10 个项目（10 次抽样） | — | %.3f | %.3f–%.3f | — |" % (
            statistics.mean(means), np.percentile(means, 5), np.percentile(means, 95)))
    log.append("")
    log.append("> 判读：折间 5–95% 区间宽度就是「换一批项目」会带来多大波动；若区间跨到 0.5 附近，说明判别结论不可靠。"
               "所有折均为「留出一个完整项目」，即不允许同项目样本同时出现在训练与测试中。")
    io.open(out_path, "w", encoding="utf-8").write("\n".join(log) + "\n")
    print("report written; see", out_path)


if __name__ == "__main__":
    main()

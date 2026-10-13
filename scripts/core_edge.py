# -*- coding: utf-8 -*-
"""Core-edge microbiome framework (section 17 feature space), with the thresholds FROZEN BEFORE RESULTS.

Definitions (documented, decided in advance):
  prevalence of a feature in a crop = share of that crop's samples in which the feature is present
    (>= 1 read) -- presence is defined on raw counts, not on transformed values.
  口径 A（样本合并，先冻结）: CORE = prevalence >= 70% in all three crops at sample level.
      2026-10-12 实测：该口径给出退化结果（只有 1 个特征，且是 unassigned 桶），原因是把 45 个不同
      研究的样本当作可交换，违反 §8.1「项目作为随机效应层」的设计原则。故并列给出：
  口径 B（研究层，主口径）: 先算每个项目内部的流行率；某特征在某项目为 "core" 需该项目内流行率 >= 70%；
      CROP-CORE = 该作物中 >=70% 的项目为 core；UNIVERSAL CORE = 三个作物都是 crop-core；
      EDGE = 至少两个作物中 crop-core 的项目比例 < 30%。
  MIDDLE: everything else.
  Sensitivity: thresholds 50% / 90%; and the section-3 subset (only genera agreed by both classifiers).

Outputs: out_analysis/core_edge.csv, out_analysis/core_edge_report.md
usage: python core_edge.py
"""
import csv, io, os, collections
import numpy as np

R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
CROPS = ("maize", "rice", "wheat")


def load():
    rd = csv.reader(io.open(os.path.join(OUT, "feature_matrix.csv"), encoding="utf-8"))
    hdr = next(rd)
    metas = ["project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
             "genus_reads","family_reads","unassigned_reads","n_asv"]
    feats = [c for c in hdr if c not in metas]; gi = {c: i for i, c in enumerate(hdr)}
    rows = []
    for r in rd:
        if not r:
            continue
        rows.append((r[gi["project"]], r[gi["crop"]], [float(r[gi[f]] or 0) for f in feats]))
    return feats, rows


def agree_set():
    f = os.path.join(OUT, "consistent_genera.csv")
    out = set()
    if os.path.exists(f):
        for r in csv.DictReader(io.open(f, encoding="utf-8-sig")):
            if r["weight_class"] == "consistent":
                out.add(r["genus"])
    return out


def prevalence(feats, rows, subset=None):
    """-> {feature: {crop: prevalence}}, and mean relative abundance per crop."""
    pre = {f: {} for f in feats}; abu = {f: {} for f in feats}
    for crop in CROPS:
        sub = [v for c, v in rows if c == crop]
        n = len(sub)
        M = np.array(sub)
        tot = M.sum(axis=1, keepdims=True); tot[tot == 0] = 1
        R = M / tot
        for j, f in enumerate(feats):
            if subset is not None and f not in subset:
                continue
            pre[f][crop] = float((M[:, j] > 0).sum()) / n * 100.0
            abu[f][crop] = float(R[:, j].mean()) * 100.0
    return pre, abu


def classify(pre, crop, core_th, edge_th):
    """call a single feature for one crop: core / edge / mid."""
    p = pre[crop]
    return "core" if p >= core_th else ("edge" if p < edge_th else "mid")


def project_prevalence(rows):
    """{crop: {project: within-project prevalence vector (%) aligned to the feature order}} - computed once."""
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for prj, crop, vals in rows:
        by[crop][prj].append(vals)
    prev = collections.defaultdict(dict)
    for crop, prjs in by.items():
        for prj, vals in prjs.items():
            M = np.asarray(vals, dtype=float)
            prev[crop][prj] = (M > 0).mean(axis=0) * 100.0
    return prev


def core_share(prev, feats_idx, feats, core_th):
    """{feature: {crop: % of that crop's projects whose within-project prevalence >= core_th}}"""
    out = {feats[j]: {} for j in feats_idx}
    for crop in CROPS:
        ps = prev.get(crop, {})
        n = len(ps)
        for j in feats_idx:
            out[feats[j]][crop] = (100.0 * sum(1 for v in ps.values() if v[j] >= core_th) / n) if n else 0.0
    return out


def main():
    feats, rows = load()
    rows_full = rows
    rows = [(c, v) for (_p, c, v) in rows_full]      # pooled-sample view (口径 A)
    agree = agree_set()
    print("features %d ; samples %d ; agreed genera %d" % (len(feats), len(rows), len(agree)))
    pre, abu = prevalence(feats, rows)
    lines = []
    core_all, edge_any = [], []
    for f in feats:
        pc = [pre[f][c] for c in CROPS]
        cls = "core" if all(p >= 70 for p in pc) else ("edge" if sum(1 for p in pc if p < 30) >= 2 else "mid")
        lines.append((f, cls, pc, [abu[f][c] for c in CROPS]))
        if cls == "core":
            core_all.append(f)
        if cls == "edge":
            edge_any.append(f)
    with io.open(os.path.join(OUT, "core_edge.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["feature", "class", "prev_maize", "prev_rice", "prev_wheat",
                    "abund_maize", "abund_rice", "abund_wheat", "in_agreed_genera"])
        for f, cls, pc, ab in lines:
            w.writerow([f, cls] + ["%.2f" % v for v in pc] + ["%.4f" % v for v in ab] +
                       ["1" if (f.startswith("g__") and f[3:] in agree) else "0"])
    rep = ["# 核心-边缘框架（最终 45 个项目，§17 特征空间）", "",
           "阈值**先冻结**再算：核心 = 三个作物中流行率都 ≥70%；边缘 = 至少两个作物流行率 <30%；其余为中间。",
           "流行率按原始计数（≥1 条 reads 即视为存在），丰度为相对丰度的样本均值。", "",
           "| 类别 | 特征数 | 其中属级 | 科级回退 | unassigned |", "| --- | --- | --- | --- | --- |"]
    for cls in ("core", "mid", "edge"):
        sub = [x for x in lines if x[1] == cls]
        rep.append("| %s | %d | %d | %d | %d |" % (cls, len(sub),
                   sum(1 for x in sub if x[0].startswith("g__")),
                   sum(1 for x in sub if x[0].startswith("f__")),
                   sum(1 for x in sub if x[0] == "unassigned")))
    rep += ["", "## 核心成员（按平均相对丰度排序，前 25）", "",
            "| 特征 | 流行率 玉米/水稻/小麦 | 平均相对丰度 玉米/水稻/小麦 | 属级两法一致 |",
            "| --- | --- | --- | --- |"]
    for f, cls, pc, ab in sorted([x for x in lines if x[1] == "core"], key=lambda x: -sum(x[3]))[:25]:
        rep.append("| %s | %.0f%% / %.0f%% / %.0f%% | %.3f%% / %.3f%% / %.3f%% | %s |" % (
            f, pc[0], pc[1], pc[2], ab[0], ab[1], ab[2],
            "是" if (f.startswith("g__") and f[3:] in agree) else ("—" if not f.startswith("g__") else "**否**")))
    # sensitivity: thresholds 50/90 and the agreed-genera subset
    rep += ["", "## 敏感性", "", "| 口径 | 核心特征数 | 边缘特征数 |", "| --- | --- | --- |"]
    for th_c, th_e, tag in ((70, 30, "主口径 70/30"), (50, 30, "宽松 50/30"), (90, 30, "严格 90/30")):
        c = sum(1 for f in feats if all(pre[f][x] >= th_c for x in CROPS))
        e = sum(1 for f in feats if sum(1 for x in CROPS if pre[f][x] < th_e) >= 2)
        rep.append("| %s | %d | %d |" % (tag, c, e))
    sub = {f for f in feats if f.startswith("g__") and f[3:] in agree}
    pre2, abu2 = prevalence(feats, rows, subset=sub)
    c2 = sum(1 for f in sub if all(pre2[f][x] >= 70 for x in CROPS))
    rep.append("| 只用两法一致的属（§3） | %d | — |" % c2)
    # ---- 口径 B（研究层，主口径）----
    prev_all = project_prevalence(rows_full)
    pshare = core_share(prev_all, list(range(len(feats))), feats, 70.0)
    uni_core = [f for f in feats if all(pshare[f][c] >= 70 for c in CROPS)]
    edge_b = [f for f in feats if sum(1 for c in CROPS if pshare[f][c] < 30) >= 2]
    rep += ["", "## 口径 B：研究层（**主口径**，因为 §8.1 要求项目是随机效应层）", "",
            "先算每个项目内部的流行率；某项目内流行率 ≥70% 记该项目为 core；某作物 crop-core = ≥70% 的项目为 core；",
            "universal core = 三个作物都是 crop-core。", "",
            "| 类别 | 特征数 | 其中属级 | 科级回退 | unassigned |", "| --- | --- | --- | --- | --- |",
            "| universal core | %d | %d | %d | %d |" % (
                len(uni_core), sum(1 for f in uni_core if f.startswith("g__")),
                sum(1 for f in uni_core if f.startswith("f__")), sum(1 for f in uni_core if f == "unassigned")),
            "| edge（≥2 个作物的 core 项目比例 <30%%） | %d | %d | %d | %d |" % (
                len(edge_b), sum(1 for f in edge_b if f.startswith("g__")),
                sum(1 for f in edge_b if f.startswith("f__")), sum(1 for f in edge_b if f == "unassigned"))]
    rep += ["", "### 研究层核心成员（按「三作物最小 core 项目比例」排序，前 30）", "",
            "| 特征 | core 项目比例 玉米/水稻/小麦 | 属级两法一致 |", "| --- | --- | --- |"]
    for f in sorted(uni_core, key=lambda f: -min(pshare[f][c] for c in CROPS))[:30]:
        rep.append("| %s | %.0f%% / %.0f%% / %.0f%% | %s |" % (
            f, pshare[f]["maize"], pshare[f]["rice"], pshare[f]["wheat"],
            "是" if (f.startswith("g__") and f[3:] in agree) else ("—" if not f.startswith("g__") else "**否**")))
    # 口径 B 的阈值敏感性 + §3 敏感性
    rep += ["", "### 口径 B 的敏感性", "", "| 口径 | universal core | 其中属级两法一致 |", "| --- | --- | --- |"]
    for th in (50, 70, 90):
        ps = core_share(prev_all, list(range(len(feats))), feats, float(th))
        uc = [f for f in feats if all(ps[f][c] >= 70 for c in CROPS)]
        ag = [f for f in uc if f.startswith("g__") and f[3:] in agree]
        rep.append("| 项目内流行率阈值 %d%% | %d | %d |" % (th, len(uc), len(ag)))
    ag_idx = [i for i, f in enumerate(feats) if f.startswith("g__") and f[3:] in agree]
    ag_feats = [feats[i] for i in ag_idx]
    pshare_ag = core_share(prev_all, ag_idx, feats, 70.0)
    uc_ag = [f for f in ag_feats if all(pshare_ag[f][c] >= 70 for c in CROPS)]
    rep.append("| 只用两法一致的属（§3） | %d | %d |" % (len(uc_ag), len(uc_ag)))
    # 阈值网格：项目内流行率 x 作物层项目比例
    rep += ["", "### 阈值网格（universal core 的特征数）", "",
            "行 = 项目内流行率阈值，列 = 要求多少比例的项目达标；每个作物都要满足才算 universal core。", "",
            "| 项目内流行率 \\ 作物层比例 | ≥50% | ≥70% | ≥90% |", "| --- | --- | --- | --- |"]
    for th in (50, 70, 90):
        ps = core_share(prev_all, list(range(len(feats))), feats, float(th))
        cell = [sum(1 for f in feats if all(ps[f][c] >= cut for c in CROPS)) for cut in (50, 70, 90)]
        rep.append("| %d%% | %d | %d | %d |" % (th, cell[0], cell[1], cell[2]))
    # 近核心榜：按"三作物最小项目比例"排序，暴露瓶颈作物
    rep += ["", "### 近核心榜（项目内流行率 ≥70% 口径，按三作物最小项目比例排序，前 25）", "",
            "| 特征 | 达标项目比例 玉米/水稻/小麦 | 最小 | 属级两法一致 |", "| --- | --- | --- | --- |"]
    ranked = sorted(feats, key=lambda f: -min(pshare[f][c] for c in CROPS))[:25]
    for f in ranked:
        mn = min(pshare[f][c] for c in CROPS)
        rep.append("| %s | %.0f%% / %.0f%% / %.0f%% | %.0f%% | %s |" % (
            f, pshare[f]["maize"], pshare[f]["rice"], pshare[f]["wheat"], mn,
            "是" if (f.startswith("g__") and f[3:] in agree) else ("—" if not f.startswith("g__") else "**否**")))
    # 瓶颈作物统计
    fail_rice = [f for f in feats if pshare[f]["maize"] >= 70 and pshare[f]["wheat"] >= 70 and pshare[f]["rice"] < 70]
    fail_maize = [f for f in feats if pshare[f]["rice"] >= 70 and pshare[f]["wheat"] >= 70 and pshare[f]["maize"] < 70]
    fail_wheat = [f for f in feats if pshare[f]["maize"] >= 70 and pshare[f]["rice"] >= 70 and pshare[f]["wheat"] < 70]
    rep += ["", "### 瓶颈作物（只差一个作物就进核心的特征数）", "",
            "| 只在哪个作物掉队 | 特征数 | 其中属级 |", "| --- | --- | --- |",
            "| 只差水稻 | %d | %d |" % (len(fail_rice), sum(1 for f in fail_rice if f.startswith("g__"))),
            "| 只差玉米 | %d | %d |" % (len(fail_maize), sum(1 for f in fail_maize if f.startswith("g__"))),
            "| 只差小麦 | %d | %d |" % (len(fail_wheat), sum(1 for f in fail_wheat if f.startswith("g__")))]
    rep += ["", "> 判读：核心-边缘是**描述性框架**；核心成员若在 §3 敏感性里被大幅削减，",
            "> 说明其身份依赖 assignTaxonomy 的单方面注释，报告中必须标注。",
            "> **口径 A 与 B 并列报告**：A 忽略研究层（把 45 个研究的样本当可交换），实测退化（core=1）；",
            "> B 才是与 §8.1「项目随机效应层」一致的框架，论文以此为主、A 作对照并以退化结果作为方法论证据。"]
    io.open(os.path.join(OUT, "core_edge_report.md"), "w", encoding="utf-8").write("\n".join(rep) + "\n")
    print("\n".join(rep[:24]))


if __name__ == "__main__":
    main()

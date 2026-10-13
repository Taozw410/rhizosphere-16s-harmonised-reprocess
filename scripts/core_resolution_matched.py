# -*- coding: utf-8 -*-
"""Resolution-matched core analysis: is the "rice bottleneck" biology or a taxonomic-resolution artifact?

Logic: rice projects have the highest unassigned (resolution-limited) share (22.7% vs maize 16.5% /
wheat 12.3%).  If within-project prevalence depends on a project's resolution, then comparing crops at
unequal resolution is confounded.  This script therefore:
  1) re-runs the project-level core definition on the full set (baseline);
  2) re-runs it on RESOLUTION-MATCHED subsets (top-k projects per crop by resolution; and a common
     resolution band), so all three crops are compared at equal resolution;
  3) quantifies, per feature, how strongly its within-project prevalence depends on the project's
     resolution (linear slope across the 45 projects).

usage: python core_resolution_matched.py   -> out_analysis/core_resolution_matched.md
"""
import csv, io, os, sys, collections, statistics
import numpy as np
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
sys.path.insert(0, os.path.join(R, "scripts"))
import core_edge as CE        # reuse load/project_prevalence/core_share (the frozen definitions)
CROPS = ("maize", "rice", "wheat")


def main():
    feats, rows = CE.load()
    prev = CE.project_prevalence(rows)                       # crop -> project -> prevalence vector
    cov = {r["project"]: r for r in csv.DictReader(io.open(os.path.join(OUT, "study_covariates.csv"), encoding="utf-8"))}
    crop_of = {p: c for (p, c, _v) in rows}
    res = {p: 100.0 - float(cov[p]["unassigned_pct"]) for p in crop_of if p in cov}
    L = ["# 分辨率匹配的核心分析（回答「水稻瓶颈是生物学还是分辨率伪影」）", "",
         "分辨率定义：resolution = 100 − 未定（unassigned）reads 占比。三作物的未定占比中位：",
         "水稻 22.7%、玉米 16.5%、小麦 12.3%（水稻最低，且这是本数据集里最直接的「分类学分辨率」指标）。", ""]

    def core_counts(prev_sub, feats_sub, th_in=70.0, th_crop=70.0):
        ps = CE.core_share(prev_sub, list(range(len(feats_sub))), feats_sub, th_in)
        uni = [f for f in feats_sub if all(ps[f][c] >= th_crop for c in CROPS)]
        return uni, ps

    # 1) baseline
    uni, ps = core_counts(prev, feats)
    L += ["## 1. 基线（全部 45 个项目，每作物 15 个）", "",
          "- universal core（项目内流行率 ≥70%%、且 ≥70%% 项目达标）：**%d** 个%s" % (
              len(uni), "（" + ", ".join(uni[:5]) + "）" if uni else ""),
          "- 只差一个作物就进核心的特征数：水稻 %d / 玉米 %d / 小麦 %d" % (
              sum(1 for f in feats if ps[f]["maize"] >= 70 and ps[f]["wheat"] >= 70 and ps[f]["rice"] < 70),
              sum(1 for f in feats if ps[f]["rice"] >= 70 and ps[f]["wheat"] >= 70 and ps[f]["maize"] < 70),
              sum(1 for f in feats if ps[f]["maize"] >= 70 and ps[f]["rice"] >= 70 and ps[f]["wheat"] < 70)), ""]

    # 2) resolution-matched subsets: top-k projects per crop by resolution
    L += ["## 2. 分辨率匹配后", "",
          "每个作物取**分辨率最高的 k 个项目**（三作物同一 k），再重算同一套核心定义。",
          "若水稻仍是瓶颈，说明瓶颈不是分辨率造成的；若瓶颈消失，则原结论是分辨率伪影。", "",
          "| k（每作物项目数） | 项目内流行率阈值 | universal core | 只差水稻 | 只差玉米 | 只差小麦 | 三作物分辨率中位 |",
          "| --- | --- | --- | --- | --- | --- | --- |"]
    for k in (15, 12, 10, 8, 6):
        keep = {}
        for c in CROPS:
            ps_ = [p for p in prev.get(c, {}) if p in res]
            keep[c] = set(sorted(ps_, key=lambda p: -res[p])[:k])
        sub_rows = [(p, c, v) for (p, c, v) in rows if p in keep[c]]
        sub_prev = CE.project_prevalence(sub_rows)
        for th_in in (70.0, 50.0):
            uni2, ps2 = core_counts(sub_prev, feats, th_in=th_in)
            fr = sum(1 for f in feats if ps2[f]["maize"] >= 70 and ps2[f]["wheat"] >= 70 and ps2[f]["rice"] < 70)
            fm_ = sum(1 for f in feats if ps2[f]["rice"] >= 70 and ps2[f]["wheat"] >= 70 and ps2[f]["maize"] < 70)
            fw = sum(1 for f in feats if ps2[f]["maize"] >= 70 and ps2[f]["rice"] >= 70 and ps2[f]["wheat"] < 70)
            med = " / ".join("%.1f" % statistics.median([res[p] for p in keep[c]]) for c in CROPS)
            L.append("| %d | %.0f%% | %d | %d | %d | %d | %s |" % (k, th_in, len(uni2), fr, fm_, fw, med))

    # 3) per-feature dependence of prevalence on resolution
    L += ["", "## 3. 逐特征：流行率对项目分辨率的依赖", "",
          "对近核心榜前 20 个特征，用 45 个项目的（项目分辨率, 该项目内流行率）做线性回归；",
          "斜率显著为正说明「分辨率越高的项目越容易检出该特征」，即核心与否部分取决于分辨率。", "",
          "| 特征 | 斜率（每 +1pp 分辨率带来的流行率变化，pp） | R² | 显著性（|t|>2） |", "| --- | --- | --- | --- |"]
    ranked = sorted(feats, key=lambda f: -min(ps[f][c] for c in CROPS))[:20]
    jmap = {f: i for i, f in enumerate(feats)}
    for f in ranked:
        j = jmap[f]
        xs, ys = [], []
        for c in CROPS:
            for p, vec in prev[c].items():
                if p in res:
                    xs.append(res[p]); ys.append(vec[j])
        xs = np.asarray(xs); ys = np.asarray(ys)
        if xs.std() == 0:
            continue
        b, a = np.polyfit(xs, ys, 1)
        yhat = a + b * xs
        ss_res = ((ys - yhat) ** 2).sum(); ss_tot = ((ys - ys.mean()) ** 2).sum()
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        t = b / (np.sqrt(ss_res / max(1, len(xs) - 2)) / (xs.std() * np.sqrt(len(xs)))) if ss_res > 0 else 0
        L.append("| %s | %+.2f | %.2f | %s |" % (f, b, r2, "是" if abs(t) > 2 else "否"))
    L += ["", "> 判读规则（先写死）：若 k=10 的分辨率匹配子集里「只差水稻」仍明显多于只差玉米/小麦，",
          "> 则水稻瓶颈**不是**分辨率伪影；反之则在论文里必须写明它是分辨率差异造成的。"]
    io.open(os.path.join(OUT, "core_resolution_matched.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L[:60]))


if __name__ == "__main__":
    main()

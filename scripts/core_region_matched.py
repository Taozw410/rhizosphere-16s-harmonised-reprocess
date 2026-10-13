# -*- coding: utf-8 -*-
"""Region-matched core analysis: the second structural confound is amplicon region.

Maize projects are mostly ~250-301 bp amplicons (V4-like) while rice/wheat are often 415-427 bp
(V3-V4-like).  This script splits projects by amplicon length mode and re-runs the project-level core
definition WITHIN the V4-like subset only, so all three crops are compared at the same region.

usage: python core_region_matched.py  -> out_analysis/core_region_matched.md
"""
import csv, io, os, sys, statistics
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
sys.path.insert(0, os.path.join(R, "scripts"))
import core_edge as CE
CROPS = ("maize", "rice", "wheat")


def main():
    feats, rows = CE.load()
    cov = {r["project"]: r for r in csv.DictReader(io.open(os.path.join(OUT, "study_covariates.csv"), encoding="utf-8"))}
    L = ["# 区域匹配的核心分析（处理扩增子区域混杂）", "",
         "玉米项目的扩增子长度主峰多在 250–301 bp（V4 类），水稻/小麦多在 415–427 bp（V3–V4 类）。",
         "为避免「区域」与「作物」共线造成的伪结论，下面只在**同一区域类**内比较三个作物。", "",
         "| 区域类 | 玉米项目数 | 水稻 | 小麦 | universal core | 只差水稻 | 只差玉米 | 只差小麦 |",
         "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for tag, lo, hi in (("V4 类（≤320 bp）", 0, 320), ("V3–V4 类（>320 bp）", 321, 10000)):
        keep = {c: set() for c in CROPS}
        for p, c, _v in rows:
            m = float(cov[p]["asv_len_mode"]) if p in cov else 0
            if lo <= m <= hi:
                keep[c].add(p)
        if min(len(keep[c]) for c in CROPS) < 4:
            L.append("| %s | %d | %d | %d | —（某作物项目 <4，无法比较） | — | — | — |" % (
                tag, len(keep["maize"]), len(keep["rice"]), len(keep["wheat"])))
            continue
        sub_rows = [(p, c, v) for (p, c, v) in rows if p in keep[c]]
        prev = CE.project_prevalence(sub_rows)
        ps = CE.core_share(prev, list(range(len(feats))), feats, 70.0)
        uni = [f for f in feats if all(ps[f][c] >= 70 for c in CROPS)]
        fr = sum(1 for f in feats if ps[f]["maize"] >= 70 and ps[f]["wheat"] >= 70 and ps[f]["rice"] < 70)
        fm = sum(1 for f in feats if ps[f]["rice"] >= 70 and ps[f]["wheat"] >= 70 and ps[f]["maize"] < 70)
        fw = sum(1 for f in feats if ps[f]["maize"] >= 70 and ps[f]["rice"] >= 70 and ps[f]["wheat"] < 70)
        L.append("| %s | %d | %d | %d | %d | %d | %d | %d |" % (
            tag, len(keep["maize"]), len(keep["rice"]), len(keep["wheat"]), len(uni), fr, fm, fw))
    # which crops have projects in each region class
    L += ["", "## 各作物落入两个区域类的项目", ""]
    for c in CROPS:
        v4 = [p for p, cc, _ in rows if cc == c and float(cov[p]["asv_len_mode"]) <= 320]
        v34 = [p for p, cc, _ in rows if cc == c and float(cov[p]["asv_len_mode"]) > 320]
        L.append("- **%s**：V4 类 %d 个（%s）；V3–V4 类 %d 个（%s）" % (
            c, len(set(v4)), ", ".join(sorted(set(v4))[:4]), len(set(v34)), ", ".join(sorted(set(v34))[:4])))
    L += ["", "> 判读：若某个区域类里三作物都有 ≥4 个项目，则该类内部的比较不受区域混杂影响；",
          "> 若区域内水稻仍是瓶颈，则「水稻瓶颈」在区域混杂之外依然成立。"]
    io.open(os.path.join(OUT, "core_region_matched.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("wrote out_analysis/core_region_matched.md")


if __name__ == "__main__":
    main()

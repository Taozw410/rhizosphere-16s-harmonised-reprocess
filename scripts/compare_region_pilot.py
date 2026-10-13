# -*- coding: utf-8 -*-
"""比较"现行全长参考"与"区域匹配参考"的属级注释率。"""
import csv, io, os, glob
R = r"D:\CodexWork\mic-re"
def genus_rate(path):
    if not os.path.exists(path): return None
    n = g = 0
    for row in csv.DictReader(io.open(path, encoding="utf-8-sig")):
        n += 1
        v = (row.get("Genus") or "").strip()
        if v and v.upper() != "NA": g += 1
    return (100.0 * g / n) if n else None
out = io.open(os.path.join(R, "out50k", "region_pilot", "comparison.md"), "w", encoding="utf-8")
out.write("# 区域匹配参考 vs 现行全长参考：属级注释率对比（试点）\n\n")
out.write("| 项目 | ASV 数 | 现行（全长 SILVA） | 区域参考 V4 | 区域参考 V3–V4 | 最好 | 结论 |\n")
out.write("| --- | --- | --- | --- | --- | --- | --- |\n")
for p in ("PRJDB19884", "PRJEB65459", "PRJEB40690", "PRJNA394069"):
    cur = genus_rate(os.path.join(R, "out_stage2", p, "asv_taxonomy.csv"))
    vals = {}
    for k, suffix in (("V4", "_V4_tax.csv"), ("V3V4", "_V3V4_tax.csv")):
        f = os.path.join(R, "out50k", "region_pilot", p + suffix)
        vals[k] = genus_rate(f)
    allv = [v for v in [vals["V4"], vals["V3V4"]] if v is not None]
    best = max(allv) if allv else None
    gain = (best - cur) if (best is not None and cur is not None) else None
    verdict = "-"
    if gain is not None:
        verdict = "提升 %.1f 个百分点 → %s" % (gain, "达到 ≥60%" if best >= 60 else ("接近门槛" if best >= 50 else "仍偏低"))
    out.write("| %s | - | %s | %s | %s | %s | %s |\n" % (
        p, ("%.1f%%" % cur) if cur is not None else "-",
        ("%.1f%%" % vals["V4"]) if vals["V4"] is not None else "-",
        ("%.1f%%" % vals["V3V4"]) if vals["V3V4"] is not None else "-",
        ("%.1f%%" % best) if best is not None else "-", verdict))
out.close()
print("wrote out50k/region_pilot/comparison.md")

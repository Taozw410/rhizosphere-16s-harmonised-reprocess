# -*- coding: utf-8 -*-
import csv, io, os
R = r"D:\CodexWork\mic-re"
def rd(p): return list(csv.DictReader(io.open(p, encoding="utf-8-sig"))) if os.path.exists(p) else []
rows = []
for p in sorted(os.listdir(os.path.join(R, "out_stage2"))):
    d = os.path.join(R, "out_stage2", p)
    if not (p.startswith("PRJ") and os.path.isdir(d)): continue
    tp = os.path.join(d, "asv_taxonomy.csv"); ap = os.path.join(d, "asv_table.csv"); rt = os.path.join(d, "read_tracking.csv")
    if not (os.path.exists(tp) and os.path.exists(ap)): continue
    ab = {}
    with io.open(ap, encoding="utf-8-sig") as fh:
        rr = csv.reader(fh); next(rr, None)
        for row in rr:
            if row: ab[row[0].strip().strip('"')] = sum(float(x) for x in row[1:] if x not in ("", "NA"))
    tot = g_reads = ph_reads = fam = 0.0; n = nass = 0
    for row in csv.DictReader(io.open(tp, encoding="utf-8-sig")):
        seq = (row.get("") or "").strip().strip('"'); a = ab.get(seq, 0.0)
        tot += a; n += 1
        g = (row.get("Genus") or "").strip(); f = (row.get("Family") or "").strip(); ph = (row.get("Phylum") or "").strip()
        if g and g.upper() != "NA": nass += 1; g_reads += a
        if f and f.upper() != "NA": fam += a
        if ph and ph.upper() != "NA": ph_reads += a
    m = rd(rt); merge = None
    if m:
        tf = sum(float(x.get("filtered") or 0) for x in m); tm = sum(float(x.get("merged") or 0) for x in m)
        if tf: merge = 100.0*tm/tf
    rows.append(dict(p=p, n=n, genus=100.0*nass/n if n else 0, g_reads=100.0*g_reads/tot if tot else 0,
                     ph=100.0*ph_reads/tot if tot else 0, fam=100.0*fam/tot if tot else 0, merge=merge))
rows.sort(key=lambda r: r["g_reads"])
print("%-14s %-7s %-10s %-12s %-12s %-8s %s" % ("项目","ASV数","属级(ASV%)","属级(reads%)","到门(reads%)","合并率","标记"))
cur = new = 0
for r in rows:
    ok_cur = (r["g_reads"] >= 60) and (r["merge"] is None or r["merge"] >= 40)
    ok_new = (r["ph"] >= 80) and (r["fam"] >= 50) and (r["merge"] is None or r["merge"] >= 40)
    cur += ok_cur; new += ok_new
    mark = ("[现行OK]" if ok_cur else "") + ("[新指标可救]" if (ok_new and not ok_cur) else "")
    print("%-14s %-7d %-10s %-12s %-12s %-8s %s" % (r["p"], r["n"], "%.1f%%" % r["genus"], "%.1f%%" % r["g_reads"],
          "%.1f%%" % r["ph"], "%.1f%%" % r["merge"] if r["merge"] is not None else "-", mark))
print("")
print("现行 §4（属级 reads>=60%% + 合并率>=40%%）：%d 个" % cur)
print("新指标（到门>=80%% + 到科>=50%% + 合并率>=40%%）：%d 个" % new)

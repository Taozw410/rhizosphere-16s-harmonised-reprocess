# -*- coding: utf-8 -*-
"""Apply a truncLen plan to a DADA2 primer table: empty primers (= no cutadapt trimming, so the
planned truncLen stays valid) + the planned truncLen. Rows whose note is not "ok" get 0,0.
Usage: python apply_trunclen_plan.py <plan.csv> <primer_table.csv>
"""
import io, os, sys, csv
plan, table = sys.argv[1], sys.argv[2]
primer_src = sys.argv[3] if len(sys.argv) > 3 else ""
psrc = {}
if primer_src and os.path.exists(primer_src):
    psrc = {r["project"]: r for r in csv.DictReader(io.open(primer_src, encoding="utf-8-sig"))}
rows = list(csv.DictReader(io.open(plan, encoding="utf-8-sig")))
planmap = {r["project"]: r for r in rows}
raw = io.open(table, encoding="utf-8-sig").read().replace("\r\n", "\n").split("\n")
header = raw[0].split(",")
cols = {n: i for i, n in enumerate(header)}
out = [header]
changed = []
for line in raw[1:]:
    if not line.strip():
        continue
    f = line.split(",")
    while len(f) < len(header):
        f.append("")
    p = f[cols["project"]]
    r = planmap.get(p)
    if r:
        ok = r["note"].strip() == "ok"
        tf = r["tf"] if ok else "0"
        tr = r["tr"] if ok else "0"
        for c, v in (("primer_f", (psrc.get(p, {}).get("primer_f") or "") if psrc else ""),
                     ("primer_r", (psrc.get(p, {}).get("primer_r") or "") if psrc else ""),
                     ("truncLen_f", tf), ("truncLen_r", tr)):
            if c in cols:
                f[cols[c]] = v
        if "note" in cols:
            f[cols["note"]] = ("strategyA_trunclen_ov%s" % r["overlap"]) if ok else "strategyA_keep_untruncated"
        changed.append((p, tf, tr, r["insert_len"], r["overlap"] if ok else "-"))
    out.append(f)
# append rows for plan projects that the table does not have yet (e.g. a project whose primers
# were only detected by the new screen and which therefore had no row before)
have = set(x[cols["project"]] for x in out[1:] if len(x) > cols["project"])
for p, r in planmap.items():
    if p in have:
        continue
    ok = r["note"].strip() == "ok"
    row = [""] * len(header)
    row[cols["project"]] = p
    if "primer_f" in cols: row[cols["primer_f"]] = (psrc.get(p, {}).get("primer_f") or "") if psrc else ""
    if "primer_r" in cols: row[cols["primer_r"]] = (psrc.get(p, {}).get("primer_r") or "") if psrc else ""
    if "truncLen_f" in cols: row[cols["truncLen_f"]] = r["tf"] if ok else "0"
    if "truncLen_r" in cols: row[cols["truncLen_r"]] = r["tr"] if ok else "0"
    if "note" in cols: row[cols["note"]] = ("strategyA_trunclen_ov%s" % r["overlap"]) if ok else "strategyA_keep_untruncated"
    out.append(row)
    changed.append((p, row[cols["truncLen_f"]], row[cols["truncLen_r"]], r["insert_len"], r["overlap"] if ok else "-"))

text = "\r\n".join(",".join(x) for x in out) + "\r\n"
bad = [ch for ch in text if ord(ch) > 127]
if bad:
    sys.exit("non-ASCII content found; refusing to write")
io.open(table, "w", encoding="utf-8-sig", newline="").write(text)
for c in changed:
    print("  %-14s truncLen=(%s,%s)  insert=%s  overlap=%s" % c)
print("wrote", table, os.path.getsize(table), "bytes")

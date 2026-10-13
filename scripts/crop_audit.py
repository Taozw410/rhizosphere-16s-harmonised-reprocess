# -*- coding: utf-8 -*-
"""Audit crop composition of processed / pending / candidate project pools."""
import csv, io, os, glob
from collections import Counter

R = r"D:\CodexWork\mic-re"


def load_crop_map():
    m = {}
    for f in sorted(glob.glob(os.path.join(R, "meta", "*.csv"))):
        try:
            rows = list(csv.DictReader(io.open(f, encoding="utf-8-sig")))
        except Exception:
            continue
        for r in rows:
            k = (r.get("study") or r.get("project") or "").strip()
            c = (r.get("crop") or "").strip().lower()
            if k and c in ("maize", "rice", "wheat"):
                m.setdefault(k, c)
    return m


crop = load_crop_map()

qc = list(csv.DictReader(io.open(os.path.join(R, "out50k", "marker_gene_qc.csv"), encoding="utf-8-sig")))
ok = [r["project"] for r in qc if r["verdict"] == "OK_16S"]
mixed = [r["project"] for r in qc if r["verdict"] == "CHECK_MIXED"]
excl = [r["project"] for r in qc if r["verdict"] not in ("OK_16S", "CHECK_MIXED")]

pend_short = ["PRJNA596848", "PRJNA996548", "PRJNA1476760", "PRJNA816337",
              "PRJNA665860", "PRJNA1282585", "PRJNA1334676", "PRJNA601112", "PRJNA753114"]
pend_long = ["PRJNA1157226", "PRJNA1225875", "PRJNA1261944", "PRJNA509515",
             "PRJNA520770", "PRJNA997935", "PRJNA607209"]

damaged = {"PRJNA430104"}

groups = {
    "usable_done": [p for p in ok if p not in damaged],
    "damaged_done": [p for p in ok if p in damaged],
    "mixed_done": mixed,
    "pending_short": pend_short,
    "pending_long": pend_long,
}

rows = []
for g, ps in groups.items():
    for p in ps:
        rows.append(dict(project=p, group=g, crop=crop.get(p, "?")))

with io.open(os.path.join(R, "out50k", "crop_audit.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["project", "group", "crop"])
    w.writeheader()
    w.writerows(rows)

print("%-14s %6s %6s %6s %6s" % ("group", "all", "maize", "rice", "wheat"))
for g, ps in groups.items():
    c = Counter(crop.get(p, "?") for p in ps)
    print("%-14s %6d %6d %6d %6d   unknown=%d" % (g, len(ps), c["maize"], c["rice"], c["wheat"], c["?"]))

pot = groups["usable_done"] + groups["pending_short"] + groups["pending_long"] + groups["mixed_done"]
c = Counter(crop.get(p, "?") for p in pot)
print()
print("potential total = %d  maize=%d rice=%d wheat=%d unknown=%d" % (len(pot), c["maize"], c["rice"], c["wheat"], c["?"]))
print("gap to 15/crop = maize %d rice %d wheat %d" % (15 - c["maize"], 15 - c["rice"], 15 - c["wheat"]))

used = set(ok) | set(mixed) | set(excl) | set(pend_short) | set(pend_long)
for pool in ("batch6_candidates.csv", "batch5_candidates.csv", "lightweight_60_final.csv",
             "eligible.csv", "candidates.csv"):
    path = os.path.join(R, "meta", pool)
    if not os.path.exists(path):
        continue
    try:
        rs = list(csv.DictReader(io.open(path, encoding="utf-8-sig")))
    except Exception:
        continue
    fresh = [r for r in rs if (r.get("study") or r.get("project") or "").strip() not in used]
    cc = Counter(((r.get("crop") or "?").strip().lower()) for r in fresh)
    print("%-28s rows=%-5d fresh=%-5d %s" % (pool, len(rs), len(fresh), dict(cc)))

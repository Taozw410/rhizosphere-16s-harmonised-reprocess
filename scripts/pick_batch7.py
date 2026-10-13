# -*- coding: utf-8 -*-
"""Pick the new crop-balanced batch from the eligible pool.

used = filesystem evidence (downloaded/processed dirs) + master_inventory + pending + known-bad.
Excludes fungal / long-read / metagenome titles and collapses same-trial series.
"""
import csv, io, os, glob, re
from collections import Counter

R = r"D:\CodexWork\mic-re"


def read_rows(path):
    try:
        return list(csv.DictReader(io.open(path, encoding="utf-8-sig")))
    except Exception:
        return []


used = set()
for r in read_rows(os.path.join(R, "out50k", "master_inventory.csv")):
    if r.get("project"):
        used.add(r["project"].strip())
for d in ["out50k", "out_stage2", "out"]:
    p = os.path.join(R, d)
    if os.path.isdir(p):
        for n in os.listdir(p):
            if re.match(r"^(PRJ|DRR|ERR|SRR)", n):
                used.add(n)
for d in ["raw", "raw50k", "raw50k_longread", "raw50k_w1done", "raw50k_excluded", "raw50k_excluded_its"]:
    p = os.path.join(R, "data", d)
    if os.path.isdir(p):
        for n in os.listdir(p):
            if os.path.isdir(os.path.join(p, n)):
                used.add(n)
rec = os.path.join(R, "_recycle")
if os.path.isdir(rec):
    for n in os.listdir(rec):
        for tok in re.findall(r"PRJ[A-Z]*\d+", n):
            used.add(tok)

known_bad = set()
for r in read_rows(os.path.join(R, "out50k", "rescreen_strict.csv")):
    if (r.get("verdict") or "").startswith("EXCLUDE"):
        known_bad.add(r["project"])
for r in read_rows(os.path.join(R, "out50k", "marker_gene_qc.csv")):
    if (r.get("verdict") or "").startswith("EXCLUDE"):
        known_bad.add(r["project"])
for r in read_rows(os.path.join(R, "meta", "batch3_candidates.csv")):
    pass

pending = ["PRJNA596848", "PRJNA996548", "PRJNA1476760", "PRJNA816337", "PRJNA665860",
           "PRJNA1282585", "PRJNA1334676", "PRJNA601112", "PRJNA753114", "PRJNA1157226",
           "PRJNA1225875", "PRJNA1261944", "PRJNA509515", "PRJNA520770", "PRJNA997935",
           "PRJNA607209"]
used |= set(pending)

b6 = set()
for r in read_rows(os.path.join(R, "meta", "batch6_candidates.csv")):
    b6.add(r["study"])
for r in read_rows(os.path.join(R, "meta", "batch5_candidates.csv")):
    b6.add(r["study"])


def score(r, is_b6):
    s = 0
    n = int(r.get("ill_amp") or 0)
    mr = int(r.get("med_reads") or 0)
    if r.get("funcmarker") == "0":
        s += 40
    if r.get("nonplant") == "0":
        s += 20
    if 20 <= n <= 200:
        s += 30
    elif 200 < n <= 400:
        s += 22
    elif n > 400:
        s += 10
    elif n >= 15:
        s += 12
    if mr >= 50000:
        s += 20
    elif mr >= 20000:
        s += 12
    t = (r.get("title") or "").lower()
    if "rhizosphere" in t or "rhizoplane" in t or "root" in t:
        s += 15
    if "16s" in t or "bacteri" in t or "prokaryot" in t or "microbi" in t:
        s += 15
    if ("its" in t or "fungal" in t or "fungi" in t) and "16s" not in t:
        s -= 60
    if "third-generation" in t or "nanopore" in t or "pacbio" in t or "long-read" in t:
        s -= 40
    if "metagenom" in t or "metatranscriptom" in t:
        s -= 60
    for w in ("pepper", "soybean", "soya", "tomato", "potato", "peanut", "cotton",
              "barley", "banana", "sugarcane", "citrus", "apple", "tea ", "cucumber", "tobacco"):
        if w in t:
            s -= 80
    if is_b6:
        s += 12
    return s


elig = read_rows(os.path.join(R, "meta", "eligible.csv"))
fresh = []
for r in elig:
    s = (r.get("study") or "").strip()
    if not s or s in used or s in known_bad:
        continue
    r["is_b6"] = s in b6
    r["score2"] = score(r, r["is_b6"])
    fresh.append(r)
fresh.sort(key=lambda r: (-r["score2"], r.get("study", "")))

best = {}
order = []
for r in fresh:
    sig = ((r.get("center") or "").strip().lower()[:40], (r.get("title") or "").strip().lower()[:45])
    if sig in best:
        continue
    best[sig] = r
    order.append(r)
fresh = order

want = {"wheat": 9, "maize": 6, "rice": 5}
reserve = {"wheat": 4, "maize": 4, "rice": 3}

sel = []
for crop in ("wheat", "maize", "rice"):
    pool = [r for r in fresh if r.get("crop") == crop]
    for r in pool[:want[crop]]:
        sel.append(("primary", crop, r["study"], r["score2"], r.get("ill_amp"), r.get("med_reads"), "Y" if r["is_b6"] else "n", r.get("title")))
    for r in pool[want[crop]:want[crop] + reserve[crop]]:
        sel.append(("reserve", crop, r["study"], r["score2"], r.get("ill_amp"), r.get("med_reads"), "Y" if r["is_b6"] else "n", r.get("title")))

out = os.path.join(R, "meta", "batch7_selection.csv")
with io.open(out, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["group", "crop", "study", "score", "ill_amp", "med_reads", "probed", "title"])
    for row in sel:
        w.writerow(list(row[:7]) + [(row[7] or "")[:120]])

print("used=%d known_bad=%d fresh=%d %s" % (len(used), len(known_bad), len(fresh), dict(Counter(r["crop"] for r in fresh))))
print()
for g in ("primary", "reserve"):
    sub = [x for x in sel if x[0] == g]
    print("== %s n=%d %s" % (g, len(sub), dict(Counter(x[1] for x in sub))))
    for row in sub:
        print("   %-7s %-14s score=%-4s n=%-5s reads=%-7s probed=%s %s" % (row[1], row[2], row[3], row[4], row[5], row[6], (row[7] or "")[:60]))
print()
print("wrote", out)

# -*- coding: utf-8 -*-
"""Choose the final new-batch list from ALL probes that passed the gates (not just the pre-screen shortlist)."""
import csv, io, os
from collections import Counter

R = r"D:\CodexWork\mic-re"


def rows(path, enc="utf-8-sig"):
    return list(csv.DictReader(io.open(path, encoding=enc)))


meta = {}
for f in ("batch6_candidates.csv", "batch7_selection.csv"):
    p = os.path.join(R, "meta", f)
    if not os.path.exists(p):
        continue
    for r in rows(p):
        k = r["study"]
        meta.setdefault(k, dict(crop=r.get("crop", "?"), samples=r.get("ill_amp", ""),
                                title=(r.get("title") or ""), center=(r.get("center") or ""),
                                score=int(r["score"]) if (r.get("score") or "").isdigit() else 0,
                                group=r.get("group", "")))

picks = rows(os.path.join(R, "meta", "batch7_selection.csv"))
prio = {}
for r in picks:
    prio[r["study"]] = 0 if r["group"] == "primary" else 1

known_bad = set()
for f in ("rescreen_strict.csv", "marker_gene_qc.csv"):
    p2 = os.path.join(R, "out50k", f)
    if not os.path.exists(p2):
        continue
    for r in rows(p2):
        if (r.get("verdict") or "").startswith("EXCLUDE"):
            known_bad.add(r["project"])

BAD_TITLE = ("fungi", "fungal", "its ", "buckwheat", "metagenom", "metatranscriptom",
             "third-generation", "nanopore", "pacbio", "pepper", "soybean", "tomato", "potato",
             "peanut", "cotton", "barley", "banana", "sugarcane", "citrus", "apple", "tea ")

scr = rows(os.path.join(R, "out50k", "batch7_screening.csv"))
scr = [r for r in scr if r["project"] not in known_bad]
ok = [r for r in scr if r["verdict"] == "OK"]
ma = {r["project"]: meta.get(r["project"], {}) for r in ok}
ok = [r for r in ok if not any(w in (ma[r["project"]].get("title", "") or "").lower() for w in BAD_TITLE)]
print("dropped by known-bad / off-target title; remaining OK:", len(ok))
print("OK probes:", len(ok), dict(Counter(r["crop"] for r in ok)))

cands = []
for r in ok:
    p = r["project"]
    m = meta.get(p, {})
    cands.append(dict(project=p, crop=r["crop"], prio=prio.get(p, 2),
                      score=m.get("score", 0), title=m.get("title", ""), center=m.get("center", ""),
                      samples=r.get("samples", ""), mode=r.get("mode", ""), rlen=r.get("rlen", ""),
                      qual=r.get("qual", ""), s16=r.get("S16", ""), ovl=r.get("ovl", "")))
cands.sort(key=lambda c: (c["prio"], -c["score"], c["project"]))

# collapse same-trial series (same centre + near-identical title)
seen = set(); uniq = []
for c in cands:
    sig = (c["center"].strip().lower()[:40], c["title"].strip().lower()[:45])
    if sig in seen:
        print("  twin skipped:", c["project"], "|", c["title"][:50])
        continue
    seen.add(sig); uniq.append(c)

target = {"wheat": 10, "rice": 8, "maize": 6}
total_cap = 20
chosen = []
for crop in ("wheat", "rice", "maize"):
    got = [c for c in uniq if c["crop"] == crop][:target[crop]]
    chosen += got
if len(chosen) < total_cap:
    for c in uniq:
        if len(chosen) >= total_cap:
            break
        if c not in chosen:
            chosen.append(c)
chosen = chosen[:total_cap]
chosen.sort(key=lambda c: (c["crop"], c["prio"], -c["score"]))

with io.open(os.path.join(R, "meta", "batch7_selected.txt"), "w", encoding="utf-8-sig", newline="") as f:
    for c in chosen:
        f.write(c["project"] + "\n")
with io.open(os.path.join(R, "meta", "batch7_final_selection.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["crop", "project", "priority", "score", "mode", "rlen", "qual", "S16", "ovl", "samples", "title"])
    for c in chosen:
        w.writerow([c["crop"], c["project"], c["prio"], c["score"], c["mode"], c["rlen"], c["qual"],
                    c["s16"], c["ovl"], c["samples"], c["title"][:90]])

print("\nselected =", len(chosen), dict(Counter(c["crop"] for c in chosen)))
for c in chosen:
    print("  %-6s %-14s prio=%d rlen=%-4s qual=%-3s S16=%-3s ovl=%-5s n=%-4s %s" %
          (c["crop"], c["project"], c["prio"], c["rlen"], c["qual"], c["s16"], c["ovl"], c["samples"], c["title"][:52]))
left = [c for c in uniq if c not in chosen]
print("\nstill available (backfill):", ", ".join(sorted(c["project"] for c in left)) or "(none)")

# -*- coding: utf-8 -*-
"""Strict probe screening (generalised).

usage: screen_batch7.py [pool1,pool2,...] [output_csv]
default: probe5,probe7 -> out50k/batch7_screening.csv

Gates (AGENTS.md 10.3 / 10.6):
  quality charset  < 6 distinct chars     -> EXCLUDE_QUALITY
  adapter signature > 20%                 -> EXCLUDE_ADAPTER
  ITS signature > 20% and > 16S signature -> EXCLUDE_ITS
  functional-gene signature > 15%         -> EXCLUDE_FUNC
  median read length < 100                -> EXCLUDE_SHORT
  single-end read length > 380            -> CHECK_LONGREAD
  paired overlap fraction < 15%           -> EXCLUDE_BROKEN_PAIR
"""
import csv, gzip, glob, os, io, sys
from collections import Counter

R = r"D:\CodexWork\mic-re"
pools = (sys.argv[1].split(",") if len(sys.argv) > 1 and sys.argv[1].strip() else ["probe5", "probe7"])
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(R, "out50k", "batch7_screening.csv")

SIG = {
    "ITS": ["CTTGGTCATTTAGAGGAAGTAA", "TCCGTAGGTGAACCTGCGG", "GAATCATCGAATCTTTGAACGC",
            "GCGTTCAAAGATTCGATGATTC", "GCTGCGTTCTTCATCGATGC", "CGATGAAGAACGCAGCG",
            "CCGCTTATTGATATGCTTAAG"],
    "16S": ["GTGCCAGCMGCCGCGGTAA", "GTGYCAGCMGCCGCGGTAA", "CCTACGGGNGGCWGCAG",
            "ACTCCTACGGGAGGCAGCAG", "ATTAGAWACCCBDGTAGTCC", "GACTACHVGGGTATCTAATCC",
            "GGACTACNVGGGTWTCTAAT", "TACGGGAGGCAGCAG"],
    "adapter": ["CTACACGACGCTCTTCCGATCT", "AATGATACGGCGACCACCGAGATCTACAC",
                "CAAGCAGAAGACGGCATACGAGAT"],
    "func": ["GGGGTTTCTACTGGTGGT", "GCIWTYTAYGGIAARGGIGG"],
}

_comp = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(_comp)[::-1]


def overlap_frac(f1, f2, n=150, min_ov=25, max_mm=3):
    a, b = [], []
    try:
        with gzip.open(f1, "rt", errors="replace") as fh:
            for i, l in enumerate(fh):
                if i % 4 == 1:
                    a.append(l.strip().upper())
                if len(a) >= n:
                    break
        with gzip.open(f2, "rt", errors="replace") as fh:
            for i, l in enumerate(fh):
                if i % 4 == 1:
                    b.append(l.strip().upper())
                if len(b) >= n:
                    break
    except Exception:
        return None
    hit = 0
    for x, y in zip(a, b):
        y = rc(y)
        for k in range(0, max(0, len(x) - min_ov) + 1):
            ov = min(len(x) - k, len(y))
            if ov < min_ov:
                break
            if sum(1 for p, q in zip(x[k:k + ov], y[:ov]) if p != q) <= max_mm:
                hit += 1
                break
    return 100.0 * hit / max(1, len(a))


def load_meta():
    m = {}
    for f in ("batch6_candidates.csv", "batch7_selection.csv", "batch8_candidates.csv"):
        p = os.path.join(R, "meta", f)
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(io.open(p, encoding="utf-8-sig")):
            k = r["study"]
            if k not in m:
                m[k] = dict(crop=r.get("crop", "?"), samples=r.get("ill_amp", ""),
                            title=(r.get("title") or "")[:70], group=r.get("group", ""))
    return m


meta = load_meta()
rows = []
for base_name in pools:
    base = os.path.join(R, "data", base_name)
    if not os.path.isdir(base):
        continue
    for proj in sorted(os.listdir(base)):
        d = os.path.join(base, proj)
        if not os.path.isdir(d):
            continue
        fs = sorted(glob.glob(os.path.join(d, "*.fastq.gz")))
        fs = [f for f in fs if os.path.getsize(f) > 1000]
        m = meta.get(proj, {})
        if not fs:
            rows.append(dict(project=proj, pool=base_name, crop=m.get("crop", "?"), mode="", qual=0,
                             rlen=0, ITS=0, S16=0, adpt=0, func=0, ovl="", samples=m.get("samples", ""),
                             verdict="NO_DATA"))
            continue
        paired = any(x.endswith("_1.fastq.gz") for x in fs)
        f1 = ([x for x in fs if x.endswith("_1.fastq.gz")] or fs)[0]
        q = set(); lens = {}; n = 0; hits = Counter()
        try:
            with gzip.open(f1, "rt", errors="replace") as fh:
                for i, line in enumerate(fh):
                    r = i % 4
                    if r == 1:
                        s = line.strip().upper(); n += 1
                        lens[len(s)] = lens.get(len(s), 0) + 1
                        for k, ps in SIG.items():
                            if any(x in s for x in ps):
                                hits[k] += 1
                    elif r == 3:
                        q.update(line.strip())
                    if n >= 400:
                        break
        except Exception:
            rows.append(dict(project=proj, pool=base_name, crop=m.get("crop", "?"), mode="?", qual=0,
                             rlen=0, ITS=0, S16=0, adpt=0, func=0, ovl="", samples=m.get("samples", ""),
                             verdict="READ_ERROR"))
            continue
        med = sorted(lens.items(), key=lambda x: -x[1])[0][0] if lens else 0
        its = 100.0 * hits["ITS"] / max(1, n); s16 = 100.0 * hits["16S"] / max(1, n)
        ad = 100.0 * hits["adapter"] / max(1, n); fn = 100.0 * hits["func"] / max(1, n)
        ovl = ""
        if paired:
            f2 = [x for x in fs if x.endswith("_2.fastq.gz")]
            if f2:
                o = overlap_frac(f1, f2[0])
                ovl = round(o, 1) if o is not None else None
        if len(q) < 6:
            v = "EXCLUDE_QUALITY"
        elif ad > 20:
            v = "EXCLUDE_ADAPTER"
        elif its > 20 and its > s16:
            v = "EXCLUDE_ITS"
        elif fn > 15:
            v = "EXCLUDE_FUNC"
        elif med < 100:
            v = "EXCLUDE_SHORT"
        elif (not paired) and med > 380:
            v = "CHECK_LONGREAD"
        elif paired and ovl not in ("", None) and ovl < 15:
            v = "EXCLUDE_BROKEN_PAIR"
        else:
            v = "OK"
        rows.append(dict(project=proj, pool=base_name, crop=m.get("crop", "?"),
                         mode=("paired" if paired else "single"), qual=len(q), rlen=med,
                         ITS=round(its), S16=round(s16), adpt=round(ad), func=round(fn), ovl=ovl,
                         samples=m.get("samples", ""), verdict=v))

rows.sort(key=lambda r: (r["verdict"], r["crop"], r["project"]))
with io.open(OUT, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["project", "pool", "crop", "mode", "qual", "rlen", "ITS", "S16",
                                      "adpt", "func", "ovl", "samples", "verdict"])
    w.writeheader(); w.writerows(rows)

print("pools:", ",".join(pools), " scanned=%d" % len(rows), "verdicts:", dict(Counter(r["verdict"] for r in rows)))
print()
print("%-14s %-7s %-7s %4s %5s %5s %5s %6s %-18s" % ("project", "crop", "mode", "qual", "rlen", "S16%", "ITS%", "ovl%", "verdict"))
for r in rows:
    print("%-14s %-7s %-7s %4s %5s %5s %5s %6s %-18s" % (r["project"], r["crop"], r["mode"], r["qual"],
                                                         r["rlen"], r["S16"], r["ITS"], r["ovl"], r["verdict"]))
ok = [r for r in rows if r["verdict"] == "OK"]
print("\nOK =", len(ok), dict(Counter(r["crop"] for r in ok)), "->", OUT)

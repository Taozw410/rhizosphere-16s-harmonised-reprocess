# -*- coding: utf-8 -*-
"""Gzip integrity check for downloaded FASTQ (finds truncated / CRC-broken files).

Usage: python verify_gzip_raw.py <raw_root> <p1,p2,...|ALL> [--workers 2] [--out report.csv]

Checks per file: gzip stream readable to EOF (CRC verified by Python's gzip at EOF),
size, zero-byte files, and unmatched _1/_2 pairs. Writes a CSV report and prints a summary.
"""
import os, sys, gzip, glob, io, csv, argparse
from concurrent.futures import ProcessPoolExecutor

def check(path):
    sz = os.path.getsize(path)
    if sz == 0:
        return (path, sz, "EMPTY", "")
    try:
        with gzip.open(path, "rb") as fh:
            n = 0
            while True:
                b = fh.read(1 << 22)
                if not b:
                    break
                n += len(b)
        return (path, sz, "OK", str(n))
    except Exception as e:
        return (path, sz, "BAD_GZIP", "%s: %s" % (type(e).__name__, e))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_root")
    ap.add_argument("projects")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    root = a.raw_root
    if a.projects.upper() == "ALL":
        projs = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))
    else:
        projs = [p.strip() for p in a.projects.split(",") if p.strip()]
    files = []
    for p in projs:
        files += sorted(glob.glob(os.path.join(root, p, "*.fastq.gz")))
    print("files to check:", len(files), " workers:", a.workers)
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for i, r in enumerate(ex.map(check, files, chunksize=4), 1):
            rel = os.path.relpath(r[0], root)
            rows.append((os.path.dirname(rel), os.path.basename(rel), r[1], r[2], r[3]))
            if r[2] not in ("OK",):
                print("  %-14s %-16s %s" % (r[2], os.path.basename(r[0]), r[3][:90]))
            if i % 200 == 0:
                print("   ...%d/%d" % (i, len(files)))
    ok = [r for r in rows if r[3] == "OK"]
    bad = [r for r in rows if r[3] != "OK"]
    print("OK:", len(ok), " BAD:", len(bad))
    for r in bad:
        print("   ", r)
    # unmatched pairs
    from collections import defaultdict
    seen = defaultdict(set)
    for proj, base, sz, st, info in rows:
        if base.endswith("_1.fastq.gz"):
            seen[(proj, base[:-len("_1.fastq.gz")])].add("1")
        elif base.endswith("_2.fastq.gz"):
            seen[(proj, base[:-len("_2.fastq.gz")])].add("2")
    un = [k for k, v in seen.items() if len(v) == 1]
    print("unmatched single-side runs:", len(un), un[:10])
    if a.out:
        with io.open(a.out, "w", encoding="utf-8-sig", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["project", "file", "bytes", "status", "info"])
            for r in rows:
                w.writerow(r)
        print("wrote", a.out)
    if bad:
        print("RESULT: BAD=%d  (repair before running downstream steps)" % len(bad))
        sys.exit(2)
    print("RESULT: OK  bad=0")

if __name__ == "__main__":
    main()

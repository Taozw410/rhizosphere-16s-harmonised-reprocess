# -*- coding: utf-8 -*-
"""After a subsample run, check that every sample of every project got its files in raw50k.

Usage: python check_subsample_counts.py <raw_root> <dst_root> <p1,p2,...> [--out report.csv]
Paired project: expects 2*dst files per sample that has both _1 and _2 in raw.
Single-end project: expects 1 dst file per raw file.
"""
import os, sys, glob, io, csv, argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_root"); ap.add_argument("dst_root"); ap.add_argument("projects")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    rows = []
    bad = 0
    for p in [x.strip() for x in a.projects.split(",") if x.strip()]:
        rfiles = sorted(os.path.basename(f) for f in glob.glob(os.path.join(a.raw_root, p, "*.fastq.gz")))
        dfiles = set(os.path.basename(f) for f in glob.glob(os.path.join(a.dst_root, p, "*.fastq.gz")))
        paired = any(f.endswith("_1.fastq.gz") for f in rfiles)
        if paired:
            samples = sorted(set(f[:-len("_1.fastq.gz")] for f in rfiles if f.endswith("_1.fastq.gz")))
            samples = [s for s in samples if (s + "_2.fastq.gz") in rfiles]
            need = [s + suf for s in samples for suf in ("_1.fastq.gz", "_2.fastq.gz")]
        else:
            need = rfiles
        miss = [f for f in need if f not in dfiles]
        extra = sorted(dfiles - set(need))
        st = "OK" if not miss else "MISSING %d" % len(miss)
        if miss:
            bad += 1
        rows.append((p, "paired" if paired else "single", len(need), len(dfiles), st,
                     "; ".join(miss[:5]) + (" ..." if len(miss) > 5 else ""), "; ".join(extra[:5])))
        print("%-14s %-7s need=%-4d have=%-4d %-12s %s %s" % (rows[-1][0], rows[-1][1], rows[-1][2],
              rows[-1][3], rows[-1][4], rows[-1][5], ("EXTRA:" + rows[-1][6]) if rows[-1][6] else ""))
    print("projects with missing files:", bad)
    if a.out:
        with io.open(a.out, "w", encoding="utf-8-sig", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["project", "mode", "need", "have", "status", "missing", "extra"])
            for r in rows:
                w.writerow(r)
        print("wrote", a.out)
    if bad:
        sys.exit(3)

if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Pick one representative sequence per section-17 feature (the most abundant ASV mapping to it)
and write it as FASTA, so a phylogeny can be built in the SAME feature space used everywhere else.

usage: python feature_reps.py
out:   out_analysis/feature_reps.fa  (+ feature_map already has feature -> level/reads)
"""
import csv, io, os, sys, collections
R = r"D:\CodexWork\mic-re"
sys.path.insert(0, os.path.join(R, "scripts"))
import tax_bins as TB
OUT = os.path.join(R, "out_analysis")


def main():
    feats_keep = set()
    for r in csv.DictReader(io.open(os.path.join(OUT, "feature_map.csv"), encoding="utf-8")):
        feats_keep.add(r["feature"])
    best = {}          # feature -> (reads, seq)
    sel = [r for r in csv.DictReader(io.open(os.path.join(R, "meta", "final45_selection.csv"), encoding="utf-8-sig"))
           if (r.get("rank") or "").strip()]
    for r in sel:
        pdir = os.path.join(R, r["base"], r["project"])
        reads = TB.load_reads(pdir)
        tax = TB.load_tax(pdir)
        for asv, n in reads.items():
            t = tax.get(asv)
            if t is None:
                continue
            name, level = TB.tax_label(t)
            if level == "organelle" or name not in feats_keep:
                continue
            if name not in best or n > best[name][0]:
                best[name] = (n, asv)
    with io.open(os.path.join(OUT, "feature_reps.fa"), "w", encoding="ascii", newline="\n") as fh:
        for f in sorted(feats_keep):
            if f in best:
                fh.write(">%s\n%s\n" % (f.replace(" ", "_"), best[f][1]))
    missing = [f for f in feats_keep if f not in best]
    print("features: %d ; representatives written: %d ; missing: %d" % (len(feats_keep), len(best), len(missing)))
    if missing:
        print("  missing:", ", ".join(missing[:10]))


if __name__ == "__main__":
    main()

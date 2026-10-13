# -*- coding: utf-8 -*-
"""Dominant phylum per section-17 feature (read-weighted across the 45 projects).

Needed for the life-history (copiotroph vs oligotroph) enrichment test: our feature space keeps only
genus/family labels, so the phylum has to be recovered from asv_taxonomy.csv.

usage: python feature_phylum.py  -> out_analysis/feature_phylum.csv
"""
import csv, io, os, sys, collections
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
sys.path.insert(0, os.path.join(R, "scripts"))
import tax_bins as TB

META = ["project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
        "genus_reads","family_reads","unassigned_reads","n_asv"]


def main():
    feats = [r["feature"] for r in csv.DictReader(io.open(os.path.join(OUT, "feature_map.csv"), encoding="utf-8"))]
    featset = set(feats)
    tot = collections.Counter()                      # (feature, phylum) -> reads
    sel = [r for r in csv.DictReader(io.open(os.path.join(R, "meta", "final45_selection.csv"), encoding="utf-8-sig"))
           if (r.get("rank") or "").strip()]
    for r in sel:
        pdir = os.path.join(R, r["base"], r["project"])
        reads = TB.load_reads(pdir); tax = TB.load_tax(pdir)
        for asv, n in reads.items():
            t = tax.get(asv)
            if t is None:
                continue
            name, level = TB.tax_label(t)
            if level == "organelle":
                continue
            f = name if level != "unassigned" else "unassigned"
            if f not in featset:
                continue
            ph = (t.get("Phylum") or "").strip()
            if not ph or ph.upper() == "NA":
                ph = "unassigned_phylum"
            tot[(f, ph)] += n
    best, ph_tot = {}, collections.Counter()
    for (f, ph), n in tot.items():
        ph_tot[ph] += n
        if f not in best or n > best[f][1]:
            best[f] = (ph, n)
    rows = []
    for f in feats:
        ph, n = best.get(f, ("unassigned_phylum", 0))
        rows.append(dict(feature=f, phylum=ph, reads_in_dominant=int(n), reads_total=int(sum(v for (ff, _), v in tot.items() if ff == f))))
    with io.open(os.path.join(OUT, "feature_phylum.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("features %d -> out_analysis/feature_phylum.csv" % len(rows))
    print("top phyla by reads:", ", ".join("%s=%d" % (p, n) for p, n in ph_tot.most_common(12)))


if __name__ == "__main__":
    main()

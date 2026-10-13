# -*- coding: utf-8 -*-
"""Build the classical table that the official FAPROTAX collapse_table.py expects:
first column = taxon name (feature label with the g__/f__ prefix stripped), then one column per sample.
Rows with no taxon identity ('unassigned') are dropped (FAPROTAX cannot map them).
usage: python faprotax_prep.py  -> out_analysis/faprotax_input.tsv
"""
import csv, io, os
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
META = ["project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
        "genus_reads","family_reads","unassigned_reads","n_asv"]


def main():
    rd = csv.reader(io.open(os.path.join(OUT, "feature_matrix.csv"), encoding="utf-8"))
    hdr = next(rd); feats = [c for c in hdr if c not in META]; gi = {c: i for i, c in enumerate(hdr)}
    rows = [r for r in rd if r and int(float(r[gi["bacteria_reads"]] or 0)) >= 1000]
    names = ["%s|%s" % (r[gi["project"]], r[gi["sample"]]) for r in rows]
    keep = [(f, f[3:]) for f in feats if f.startswith(("g__", "f__"))]
    with io.open(os.path.join(OUT, "faprotax_input.tsv"), "w", encoding="ascii", newline="\n") as fh:
        fh.write("taxon\t" + "\t".join(names) + "\n")
        for f, tax in keep:
            j = gi[f]
            fh.write(tax + "\t" + "\t".join(str(int(float(r[j] or 0))) for r in rows) + "\n")
    print("wrote out_analysis/faprotax_input.tsv : %d taxa x %d samples" % (len(keep), len(rows)))


if __name__ == "__main__":
    main()

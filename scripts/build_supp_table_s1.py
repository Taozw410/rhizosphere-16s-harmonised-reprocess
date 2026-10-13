"""Build supplementary Table S1 (the 45 selected projects and their processing parameters).

Source: meta/final45_selection.csv (selection record) + out_analysis/study_covariates.csv
(amplicon-length mode and unassigned fraction). Accessions are public identifiers and are
safe to publish; internal IDs are not: this file therefore contains no internal codes.

Run: D:\\bioenv\\envs\\riz16s\\python.exe scripts\\build_supp_table_s1.py
Out: supplementary/TableS1_projects.csv
"""

import csv
import io
import os
import sys

SEL = os.path.join("meta", "final45_selection.csv")
COV = os.path.join("out_analysis", "study_covariates.csv")
OUT = os.path.join("supplementary", "TableS1_projects.csv")

COLS = [
    ("project", "project"),
    ("crop", "crop"),
    ("centre", "institute"),
    ("n_samples", "n_samples"),
    ("mode", "sequencing_mode"),
    ("qual", "quality_encoding_levels"),
    ("asv", "n_asv"),
    ("fpass", "filter_pass_rate"),
    ("merge", "merge_rate"),
    ("keep", "retention_rate"),
    ("phylum_r", "phylum_assignment_rate"),
    ("family_r", "family_assignment_rate"),
    ("genus_r", "genus_assignment_rate"),
    ("rank", "selection_rank_within_crop"),
    ("marker", "marker_gene_label"),
    ("sens_flags", "sensitivity_flags"),
]


def read_csv(path):
    with io.open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    sel = [r for r in read_csv(SEL)
           if not (r.get("excl") or "").strip() and (r.get("rank") or "").strip().isdigit()
           and int(r["rank"]) <= 15]
    cov = {r["project"]: r for r in read_csv(COV)}
    # Platform clusters are published as stable but opaque ids (PC01, PC02, ...); the
    # underlying key is public metadata (institute + study title) and is kept out of the
    # supplementary file to stay readable.
    cluster_ids, seen = {}, []
    for r in sel:
        key = (r.get("cluster") or "").strip()
        if key and key not in seen:
            seen.append(key)
    for i, key in enumerate(seen, 1):
        cluster_ids[key] = "PC%02d" % i
    rows = []
    for r in sel:
        rec = {}
        for src, dst in COLS:
            rec[dst] = (r.get(src) or "").strip()
        c = cov.get(r["project"], {})
        mode_len = c.get("asv_len_mode") or ""
        rec["amplicon_length_mode_bp"] = mode_len
        if mode_len:
            bp = int(float(mode_len))
            rec["amplicon_region_class"] = "V4-like (<=320 bp)" if bp <= 320 else "V3-V4-like (>320 bp)"
        else:
            rec["amplicon_region_class"] = ""
        rec["unassigned_read_fraction_pct"] = c.get("unassigned_pct", "")
        rec["median_reads_per_sample"] = c.get("depth_median", "")
        rec["platform_cluster_id"] = cluster_ids.get((r.get("cluster") or "").strip(), "")
        rows.append(rec)
    order = {"maize": 0, "rice": 1, "wheat": 2}
    rows.sort(key=lambda r: (order.get(r["crop"], 9), int(r["selection_rank_within_crop"] or 99)))
    header = [dst for _, dst in COLS] + [
        "platform_cluster_id", "amplicon_length_mode_bp", "amplicon_region_class",
        "unassigned_read_fraction_pct", "median_reads_per_sample"]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with io.open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, lineterminator="\n")
        w.writeheader()
        for rec in rows:
            w.writerow(rec)
    print("wrote %s: %d projects (%s)" % (
        OUT, len(rows),
        ", ".join("%s %d" % (c, sum(1 for r in rows if r["crop"] == c)) for c in ("maize", "rice", "wheat"))))
    return 0


if __name__ == "__main__":
    sys.exit(main())

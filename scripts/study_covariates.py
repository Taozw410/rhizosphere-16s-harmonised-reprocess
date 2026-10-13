# -*- coding: utf-8 -*-
"""Build a study-level covariate table for the variance-decomposition analysis (one row per project).

Content: crop; technical covariates that plausibly differ between studies and could masquerade as a
crop effect - amplicon length (ASV length mode, a proxy for the sequenced region), per-sample
sequencing depth, number of ASVs, number of samples, unassigned (resolution-limited) read share,
and the mean within-project sample richness.

usage: python study_covariates.py  -> out_analysis/study_covariates.csv
"""
import csv, io, os, collections, statistics
import numpy as np
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
META = ["project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
        "genus_reads","family_reads","unassigned_reads","n_asv"]


def main():
    rd = csv.reader(io.open(os.path.join(OUT, "feature_matrix.csv"), encoding="utf-8"))
    hdr = next(rd); gi = {c: i for i, c in enumerate(hdr)}
    per = collections.defaultdict(lambda: dict(crop=None, n=0, bact=[], unas=[], nasv=[], gen=[], fam=[]))
    for r in rd:
        if not r:
            continue
        b = float(r[gi["bacteria_reads"]] or 0)
        if b < 1000:
            continue
        d = per[r[gi["project"]]]
        d["crop"] = r[gi["crop"]]; d["n"] += 1
        d["bact"].append(b)
        d["unas"].append(float(r[gi["unassigned_reads"]] or 0) / b * 100)
        d["nasv"].append(float(r[gi["n_asv"]] or 0))
        d["gen"].append(float(r[gi["genus_reads"]] or 0) / b * 100)
        d["fam"].append(float(r[gi["family_reads"]] or 0) / b * 100)
    rows = []
    for p, d in sorted(per.items()):
        # amplicon length proxy = ASV length mode in this project's table
        base = None
        for b in ("out_stage2", "out50k"):
            if os.path.exists(os.path.join(R, b, p, "asv_table.csv")):
                base = b; break
        mode = 0
        if base:
            cnt = collections.Counter()
            with io.open(os.path.join(R, base, p, "asv_table.csv"), encoding="utf-8-sig") as fh:
                rr = csv.reader(fh); next(rr, None)
                for row in rr:
                    if row:
                        cnt[len(row[0].strip().strip('"'))] += 1
            if cnt:
                mode = cnt.most_common(1)[0][0]
        rows.append(dict(project=p, crop=d["crop"], n_samples=d["n"],
                         depth_median=int(statistics.median(d["bact"])),
                         depth_mean=int(statistics.mean(d["bact"])),
                         asv_len_mode=mode, n_asv_median=int(statistics.median(d["nasv"])),
                         unassigned_pct=round(statistics.mean(d["unas"]), 2),
                         genus_pct=round(statistics.mean(d["gen"]), 2),
                         family_pct=round(statistics.mean(d["fam"]), 2)))
    with io.open(os.path.join(OUT, "study_covariates.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("projects %d -> out_analysis/study_covariates.csv" % len(rows))
    by = collections.defaultdict(list)
    for r in rows:
        by[r["crop"]].append(r)
    for c, rs in by.items():
        print("  %-6s n=%2d  amplicon_len_median=%.0f  depth_median=%.0f  unassigned=%.1f%%  genus=%.1f%%" % (
            c, len(rs), statistics.median(r["asv_len_mode"] for r in rs),
            statistics.median(r["depth_median"] for r in rs),
            statistics.mean(r["unassigned_pct"] for r in rs), statistics.mean(r["genus_pct"] for r in rs)))


if __name__ == "__main__":
    main()

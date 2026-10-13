"""Mechanical consistency check: every headline number in manuscript/*.md must exist in its designated source file.

Each check is (literal token, source file). Tokens are written exactly as they appear in the manuscript draft,
so a failure means either the draft or the source has drifted.

usage: python scripts/manuscript_number_check.py
"""
import io, os, sys

R = r"D:\CodexWork\mic-re"
OUT = os.path.join(R, "out_analysis")

CHECKS = [
    # data base
    ("1219", os.path.join(OUT, "feature_space_report.md")),
    ("2025", os.path.join(OUT, "feature_space_report.md")),
    ("64.0%", os.path.join(OUT, "feature_space_report.md")),
    ("21.4%", os.path.join(OUT, "feature_space_report.md")),
    ("26.53", os.path.join(OUT, "feature_space_report.md")),
    ("483", os.path.join(OUT, "genera_consistency_report.md")),
    ("48.2%", os.path.join(OUT, "genera_consistency_report.md")),
    # alpha
    ("7031", os.path.join(OUT, "alpha_diversity_report.md")),
    ("0.9993", os.path.join(OUT, "alpha_diversity_report.md")),
    ("167.0", os.path.join(OUT, "alpha_test_report.md")),
    ("5.52", os.path.join(OUT, "alpha_test_report.md")),
    ("13 / 45", os.path.join(OUT, "alpha_test_report.md")),
    # beta
    ("0.1329", os.path.join(OUT, "beta_all_report.md")),
    ("0.1371", os.path.join(OUT, "beta_all_report.md")),
    ("0.0936", os.path.join(OUT, "beta_all_report.md")),
    ("0.1087", os.path.join(OUT, "beta_all_report.md")),
    ("10.024", os.path.join(OUT, "dispersion_check.md")),
    ("0.0741", os.path.join(OUT, "ordination.md")),
    ("78.7", os.path.join(OUT, "beta_decomposition.md")),
    ("0.760", os.path.join(OUT, "beta_decomposition.md")),
    # variance decomposition
    ("0.6246", os.path.join(OUT, "variance_decomposition.md")),
    ("0.5898", os.path.join(OUT, "variance_decomposition.md")),
    ("0.157", os.path.join(OUT, "variance_decomposition.md")),
    ("0.752", os.path.join(OUT, "variance_decomposition.md")),
    # core / resolution / region
    ("universal core | 1", os.path.join(OUT, "core_edge_report.md")),
    ("22.7%", os.path.join(OUT, "core_resolution_matched.md")),
    ("PRJNA554261", os.path.join(OUT, "core_region_matched.md")),
    ("Haliangium", os.path.join(OUT, "core_edge_report.md")),
    # ML
    ("0.895", os.path.join(OUT, "ml_report2.md")),
    ("0.862", os.path.join(OUT, "ml_report2.md")),
    ("0.105", os.path.join(OUT, "ml_report2.md")),
    ("0.143", os.path.join(OUT, "ml_report2.md")),
    ("0.454", os.path.join(OUT, "ml_report2.md")),
    ("0.834", os.path.join(OUT, "ml_report2.md")),
    ("0.889", os.path.join(OUT, "ml_report2.md")),
    # differential abundance
    ("217", os.path.join(OUT, "da_report.md")),
    ("170", os.path.join(OUT, "da_report.md")),
    ("209", os.path.join(OUT, "da_report.md")),
    ("1501", os.path.join(OUT, "da_report.md")),
    ("59.0%", os.path.join(OUT, "sensitivity_report.md")),
    # functional
    ("0.2417", os.path.join(OUT, "functional_faprotax_table.csv")),
    ("16.330", os.path.join(OUT, "functional_faprotax_table.csv")),
    ("0.1640", os.path.join(OUT, "functional_faprotax_table.csv")),
    ("0.0966", os.path.join(OUT, "functional_faprotax_table.csv")),
    ("1.731", os.path.join(OUT, "functional_guilds.md")),
    ("0.191", os.path.join(OUT, "functional_guilds.md")),
    ("0.005", os.path.join(OUT, "water_regime_analysis.md")),
    # traits
    ("0.183", os.path.join(OUT, "trait_enrichment.md")),
    ("60.5", os.path.join(OUT, "trait_enrichment.md")),
    ("0.124", os.path.join(OUT, "trait_enrichment.md")),
    ("0.225", os.path.join(OUT, "trait_enrichment.md")),
    ("19.3", os.path.join(OUT, "trait_enrichment.md")),
    ("18.3", os.path.join(OUT, "trait_enrichment.md")),
    ("56.5", os.path.join(OUT, "trait_enrichment.md")),
    ("64.4", os.path.join(OUT, "trait_enrichment.md")),
    # selection / power
    ("1.0", os.path.join(OUT, "water_regime_analysis.md")),
    ("| 15 | 15 | 15 | 45 |", os.path.join(R, "out50k", "选集报告.md")),
    ("0.874", os.path.join(R, "out50k", "power_simulation_report_final45.md")),
    ("0.202", os.path.join(R, "out50k", "power_simulation_report_final45.md")),
    ("1218", os.path.join(R, "out50k", "power_simulation_report_final45.md")),
    # introduction / discussion, added 2026-10-13 (screening + amplicon-length span)
    ("22,317", os.path.join(R, "W1流程验证报告_20261006.md")),
    ("5,000–15,000", os.path.join(R, "W1流程验证报告_20261006.md")),
    ("88%–92%", os.path.join(R, "W1流程验证报告_20261006.md")),
    ("249", os.path.join(OUT, "study_covariates.csv")),
    ("482", os.path.join(OUT, "study_covariates.csv")),
]

# Computed checks: the number as written in the manuscript must equal the value
# computed from the source file (row counts and label tallies cannot be checked
# as literal tokens).
COUNT_CHECKS = [
    ("759", os.path.join(R, "meta", "ena_studies_maize_rice_wheat.csv"), "rows"),
    ("39,654", os.path.join(R, "meta", "ena_runs_maize_rice_wheat.csv"), "rows"),
    ("585", os.path.join(R, "meta", "study_frame.csv"), "rows"),
    ("121", os.path.join(R, "meta", "marker_class.csv"), "label:16S"),
    ("12", os.path.join(R, "meta", "marker_class.csv"), "label:ITS/18S"),
    ("11", os.path.join(R, "meta", "marker_class.csv"), "label:functional"),
    ("15", os.path.join(R, "meta", "marker_class.csv"), "label:mixed"),
    ("426", os.path.join(R, "meta", "marker_class.csv"), "label:unknown"),
]


def count_value(path, spec):
    import csv as _csv
    with io.open(path, encoding="utf-8-sig", errors="replace", newline="") as fh:
        rows = list(_csv.DictReader(fh))
    if spec == "rows":
        return len(rows)
    if spec.startswith("label:"):
        want = spec.split(":", 1)[1]
        return sum(1 for r in rows if (r.get("label") or "").strip() == want)
    raise ValueError("unknown count spec: %s" % spec)


def main():
    bad = 0
    for tok, src in CHECKS:
        if not os.path.exists(src):
            print("MISSING SOURCE :: %s" % src); bad += 1; continue
        txt = io.open(src, encoding="utf-8", errors="replace").read()
        ok = tok in txt
        if not ok:
            bad += 1
        print("%s  %-22s %s" % ("ok  " if ok else "FAIL", tok, os.path.basename(src)))
    for tok, src, spec in COUNT_CHECKS:
        if not os.path.exists(src):
            print("MISSING SOURCE :: %s" % src); bad += 1; continue
        want = tok.replace(",", "")
        got = str(count_value(src, spec))
        ok = want == got
        if not ok:
            bad += 1
        print("%s  %-22s %s (%s = %s)" % ("ok  " if ok else "FAIL", tok,
                                          os.path.basename(src), spec, got))
    print("\n%d checks, %d failures" % (len(CHECKS) + len(COUNT_CHECKS), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

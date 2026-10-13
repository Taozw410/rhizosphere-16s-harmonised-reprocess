"""Literature discovery driver for the manuscript Introduction / Discussion refs.

Reuses the nature-academic-search no-MCP fallback module (OpenAlex, stdlib only).
Writes one raw JSON per topic plus a compact markdown index for curation.

Run:  D:\\bioenv\\envs\\riz16s\\python.exe scripts\\refs_search.py
Note: needs outbound HTTPS (no proxy; the local 127.0.0.1:7888 proxy is down).
"""

import importlib.util
import json
import os
import sys

SKILL = r"C:\Users\hthome\.codex\skills\nature-academic-search\scripts\academic_search.py"
OUTDIR = os.path.join("manuscript", "refs")
RAWDIR = os.path.join(OUTDIR, "raw")

# (topic_id, query, per-topic result limit)
QUERIES = [
    ("t01_crop_genotype_effect", "rhizosphere microbiome crop genotype effect maize rice wheat", 8),
    ("t02_cross_site_reproducibility", "cross-site reproducibility soil microbial community biogeography", 8),
    ("t03_core_microbiome_scale", "core microbiome definition scale dependency rhizosphere soil", 8),
    ("t04_public_data_quality", "public sequencing data quality metadata reanalysis amplicon", 8),
    ("t05_amplicon_region_effect", "16S rRNA hypervariable region choice primer bias community comparison", 8),
    ("t06_dada2_parameters", "DADA2 amplicon sequence variant parameters truncation merging quality filtering", 8),
    ("t07_taxonomy_assigner_bias", "taxonomic assignment 16S SILVA naive Bayes IDTAXA accuracy comparison", 8),
    ("t08_electron_acceptor_paddy", "paddy soil electron acceptor sequence iron sulfate reduction methanogenesis", 8),
    ("t09_statistical_unit_multistudy", "random effect multi-study microbiome meta-analysis statistical unit", 8),
    ("t10_methanotroph_paddy", "methanotroph paddy soil rice rhizosphere methane oxidation", 8),
    ("t11_geobacter_iron_paddy", "Geobacter iron reduction paddy soil rice rhizosphere", 8),
    ("t12_flooding_aerobic_suppression", "flooding waterlogging soil microbial community aerobic heterotroph", 8),
    ("t13_paddy_nitrogen_cycling", "paddy soil nitrogen cycling nitrate respiration microbial community", 8),
    ("t14_functional_prediction_limits", "FAPROTAX 16S functional prediction limitations metagenome comparison", 8),
    ("t15_water_regime_rice", "alternate wetting and drying water management rice paddy soil microbiome", 8),
    ("t16_upland_rice", "upland rice aerobic soil rhizosphere microbial community", 8),
    ("t17_soil_biogeography", "soil microbial community distance decay biogeography continental scale", 8),
    ("t18_plant_recruitment", "plant genotype recruitment rhizosphere microbiome host selection", 8),
    ("t19_meta_analysis_rhizosphere", "meta-analysis rhizosphere microbiome crop type global patterns", 8),
    ("t20_16s_copy_number_bias", "16S rRNA gene copy number bias relative abundance estimation", 8),
    ("t21_rarefaction_depth", "sequencing depth rarefaction effect alpha diversity estimates", 8),
    ("t22_dispersion_confounding", "multivariate dispersion betadisper confounding community composition test", 8),
    ("t23_ml_crossstudy_classifier", "machine learning classification microbiome cross-study validation transferability", 8),
    ("t24_decontam_controls", "negative control decontamination low biomass amplicon sequencing", 8),
]

# Second battery: topics the first pass covered weakly (landmark methods papers,
# core-microbiome definition, study-level confounding, split-plot designs).
QUERIES_B2 = [
    ("u01_core_microbiome_definition", "defining quantifying core microbiome membership prevalence threshold", 8),
    ("u02_dada2_method", "DADA2 high resolution sample inference amplicon error model denoising", 8),
    ("u03_qiime2", "QIIME 2 reproducible interactive scalable extensible microbiome data science", 8),
    ("u04_decontam", "statistical identification contaminating sequences marker-gene sequencing decontam", 8),
    ("u05_faprotax", "FAPROTAX functional annotation prokaryotic taxa prediction cultivation-independent", 8),
    ("u06_beta_partitioning", "beta diversity partitioning turnover nestedness dissimilarity components", 8),
    ("u07_permanova_permdisp", "permutational multivariate analysis of variance dispersion homogeneity groups", 8),
    ("u08_compositional", "compositional data analysis microbiome differential abundance normalization", 8),
    ("u09_unifrac", "UniFrac phylogenetic method comparing microbial communities beta diversity", 8),
    ("u10_study_level_confounding", "study-level confounding batch effect reproducibility microbiome data integration", 8),
    ("u11_soil_continental_survey", "continental scale soil bacterial community survey land use cropland", 8),
    ("u12_rhizo_compartment_common_garden", "common garden experiment rhizosphere microbiome plant species soil origin", 8),
]

# Third battery: remaining gaps (core-microbiome definition, read merging,
# paddy sulfur/methylotrophy, method-driven between-study differences).
QUERIES_B3 = [
    ("v01_core_microbiome_challenges", "core microbiome definition challenges prospects membership", 8),
    ("v02_crop_species_field_survey", "rhizosphere microbiome differences between crop species large scale field survey", 8),
    ("v03_read_merging_overlap", "paired-end reads merging overlap amplicon length ASV recovery 16S", 8),
    ("v04_methylotrophy_paddy", "methylotrophic bacteria methanol utilization paddy soil rice", 8),
    ("v05_sulfate_reduction_paddy", "sulfate reducing bacteria paddy soil rice sulfur cycling", 8),
    ("v06_method_driven_differences", "methodological variation between studies drives microbial community differences", 8),
    ("v07_beneficial_rhizobacteria_yield", "rhizosphere microbiome field trial crop yield nitrogen fertilizer bacteria", 8),
    ("v08_environmental_controls_rhizosphere", "which factors shape rhizosphere microbiome assembly soil climate host", 8),
]

BATCHES = {1: QUERIES, 2: QUERIES_B2, 3: QUERIES_B3}


def load_module():
    spec = importlib.util.spec_from_file_location("academic_search", SKILL)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    batch = 1
    for i in (2, 3):
        if f"--batch{i}" in sys.argv:
            batch = i
    if "--batch" in sys.argv:
        batch = int(sys.argv[sys.argv.index("--batch") + 1])
    mod = load_module()
    queries = BATCHES[batch]
    rawdir = RAWDIR if batch == 1 else RAWDIR + str(batch)
    os.makedirs(rawdir, exist_ok=True)
    index = [f"# Raw OpenAlex hits, batch {batch} (for curation only)", ""]
    failures = []
    for topic, query, limit in queries:
        try:
            hits = mod.search(query=query, limit=limit)
        except Exception as exc:  # noqa: BLE001 - report and continue
            failures.append((topic, repr(exc)))
            print(f"FAIL {topic}: {exc!r}", file=sys.stderr)
            continue
        with open(os.path.join(rawdir, topic + ".json"), "w", encoding="utf-8") as fh:
            json.dump({"topic": topic, "query": query, "hits": hits}, fh, ensure_ascii=False, indent=1)
        index.append(f"## {topic}  |  `{query}`")
        index.append("")
        index.append("| # | title | year | journal | cites | doi |")
        index.append("| --- | --- | --- | --- | --- | --- |")
        for i, h in enumerate(hits, 1):
            title = (h.get("title") or "").replace("|", "/")
            journal = (h.get("journal") or "").replace("|", "/")
            index.append(f"| {i} | {title} | {h.get('year')} | {journal} | {h.get('cited_by_count')} | {h.get('doi')} |")
        index.append("")
        print(f"OK {topic}: {len(hits)} hits")
    idx_name = "_raw_index.md" if batch == 1 else f"_raw_index_batch{batch}.md"
    with open(os.path.join(OUTDIR, idx_name), "w", encoding="utf-8") as fh:
        fh.write("\n".join(index))
    if failures:
        print("FAILURES: " + "; ".join(f"{t}: {e}" for t, e in failures))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

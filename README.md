# Harmonised reprocessing of 45 public rhizosphere 16S projects

Code and derived data for the manuscript *Harmonised reprocessing of public 16S data reveals crop-level structuring
of root-associated bacterial communities along a paddy–anaerobe axis*.

## What is here

Uniform reprocessing of 45 publicly available rhizosphere 16S projects (15 each for maize, rice and wheat; 1,219
samples) into one shared 2,025-feature space, together with every statistical analysis behind the manuscript.

## What is not here

Raw sequence reads are **not** redistributed. They are public in SRA/ENA; the accessions and per-project processing
parameters for all 45 projects are in `supplementary/TableS1_projects.csv`.

## Layout

| Path | Contents |
| --- | --- |
| `scripts/` | Processing and analysis code (R and Python): pre-download probe screening, DADA2 runs, taxonomic gate, feature-space build, diversity, core–edge, classification, differential abundance, functional layer, variance partitioning, figures |
| `meta/` | Candidate frame, marker-gene assignments, per-project primer/truncation parameter tables, final selection record (plain-ASCII tables where they are read by R) |
| `out_analysis/` | Derived data and reports: feature space, distance matrices, per-module results tables and JSON summaries, figure source data |
| `out_analysis/figs/` | Manuscript figures (PNG) |
| `selection/` | Selection report, power simulation for the leave-one-study-out design, and the two-method taxonomy comparison table |
| `supplementary/` | Table S1: the 45 projects with accessions and processing parameters |
| `preregistration/` | Original (Chinese) analysis plan, frozen before data processing; kept unchanged as the pre-registration record |
| `MANIFEST.md` | Every staged file with its size and a SHA-256 prefix |

## Reproducing the analysis

1. Environment: R 4.5.3 (dada2 1.38.0, cutadapt 5.2, vegan 2.7.6, phyloseq, DECIPHER 3.6.0, ALDEx2, ANCOMBC,
zCompositions, iNEXT, ggplot2) and Python 3.12 (numpy, pandas, matplotlib).
2. Screening first: run the probe checks (`screen_quality`, `screen_marker`, `screen_probe`, `overlap_scan`) before
downloading anything. Marker-gene mismatches and read-length overruns only become visible later.
3. Processing: `run_dada2_parallel.R` per project, with truncation lengths taken from each project's measured
amplicon length (see the parameter tables in `meta/`), never from defaults.
4. Feature space and analyses: `build_feature_space.py`, then the module scripts; the manuscript's numbers are
reproduced by the tables in `out_analysis/`.

Templates and scripts are provided as they were used; some are batch drivers specific to this run and are kept for
provenance rather than as a supported pipeline.

## Licences

Code: MIT (see `LICENSE`). Derived data: CC BY 4.0 (see `LICENSE-DATA.md`).

## Citation

Please cite the manuscript and the archived release DOI (Zenodo) once assigned.

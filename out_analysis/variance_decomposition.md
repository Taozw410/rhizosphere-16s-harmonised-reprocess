# Variance decomposition: crop vs study vs study-level technical covariates

Sample level: 1192 samples from 45 projects. Study level: 45 project means.

## 1. Sample level (context: how much variation sits between studies)

| model | R2 | F | p |
| --- | --- | --- | --- |
| Bray-Curtis ~ project | 0.6246 | 43.38 | <0.001 |
| Bray-Curtis ~ crop | 0.0777 | 50.11 | <0.001 |
| Aitchison ~ project | 0.5898 | 37.48 | <0.001 |
| Aitchison ~ crop | 0.0696 | 44.47 | <0.001 |

## 2. Study level: crop alone vs technical covariates alone

| distance | model | R2 | F | p |
| --- | --- | --- | --- | --- |
| Bray-Curtis | crop | 0.1329 | 3.22 | <0.001 |
| Bray-Curtis | technical (all) | 0.3005 | 2.72 | <0.001 |
| Bray-Curtis | technical (no resolution) | 0.1411 | 1.64 | 0.002 |
| Aitchison | crop | 0.1371 | 3.34 | <0.001 |
| Aitchison | technical (all) | 0.2193 | 1.78 | <0.001 |
| Aitchison | technical (no resolution) | 0.1413 | 1.65 | <0.001 |

## 3. Variance partitioning (db-RDA): crop [a] vs technical [c], shared [b]

| distance | X2 = technical (all) | crop only [a] | shared [b] | technical only [c] | residual [d] |
| --- | --- | --- | --- | --- | --- |
| Bray-Curtis | yes | 0.058 | 0.157 | 0.033 | 0.752 |
| Aitchison | yes | 0.080 | 0.080 | 0.016 | 0.824 |
| Bray-Curtis | no resolution | 0.083 | 0.047 | 0.008 | 0.862 |
| Aitchison | no resolution | 0.098 | 0.057 | -0.002 | 0.847 |

> Read-out: [a] = crop variance NOT shared with the technical covariates;
> [b] = variance crop and technique explain jointly (inseparable in this design);
> [c] = technical-only.  If [b] dominates, the crop effect cannot be separated from method/resolution.
> The sample-level row `~ project` shows how much variation is between studies at all.

# Beta diversity over four distances (final 45 projects; pre-registration section 11)

Samples: 1192 (depth >= 1000 bacteria reads). Features: 2025. PERMANOVA = adonis2.
UniFrac computed in the section-17 feature space on an NJ tree of feature representatives
(midpoint-rooted in Python; 43 negative NJ branch lengths clamped to 0).

## 1. PERMANOVA

| distance | level | R2 | F | p |
| --- | --- | --- | --- | --- |
| Bray-Curtis | sample | 0.0777 | 50.11 | <0.001 |
| Bray-Curtis | project | 0.1329 | 3.22 | <0.001 |
| Aitchison | sample | 0.0696 | 44.47 | <0.001 |
| Aitchison | project | 0.1371 | 3.34 | <0.001 |
| unweighted UniFrac | sample | 0.0712 | 45.55 | <0.001 |
| unweighted UniFrac | project | 0.0936 | 2.17 | <0.001 |
| weighted UniFrac | sample | 0.0866 | 56.38 | <0.001 |
| weighted UniFrac | project | 0.1087 | 2.56 | <0.001 |

## 2. Dispersion check (betadisper, permutest 999)

| distance | F | p | mean distance to centroid (maize/rice/wheat) |
| --- | --- | --- | --- |
| Bray-Curtis | 10.024 | 0.001 | 0.5060 / 0.5070 / 0.5356 |
| Aitchison | 3.873 | 0.019 | 36.0579 / 34.5939 / 34.6509 |
| unweighted UniFrac | 69.382 | 0.001 | 0.3909 / 0.4802 / 0.4347 |
| weighted UniFrac | 9.771 | 0.001 | 0.2674 / 0.2856 / 0.2978 |

## 3. Leave-one-project-out (sample level, 99 perms)

| distance | R2 range | runs p<0.05 | runs |
| --- | --- | --- | --- |
| Bray-Curtis | 0.0748 - 0.0855 | 45 / 45 | 45 |
| Aitchison | 0.0665 - 0.0785 | 45 / 45 | 45 |
| unweighted UniFrac | 0.0688 - 0.0766 | 45 / 45 | 45 |
| weighted UniFrac | 0.0833 - 0.0940 | 45 / 45 | 45 |

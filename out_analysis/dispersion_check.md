# Dispersion check (betadisper) - section 11 requirement

Homogeneous dispersion is required for a PERMANOVA R2 to be interpreted as a location (composition) effect.
`permutest` tests the crop factor for differences in within-group dispersion.

| distance | F | p (permutest, 999) | mean distance to centroid (maize / rice / wheat) |
| --- | --- | --- | --- |
| Bray-Curtis | 10.024 | 0.001 | 0.5060 / 0.5070 / 0.5356 |
| Aitchison | 3.873 | 0.014 | 36.0579 / 34.5939 / 34.6509 |

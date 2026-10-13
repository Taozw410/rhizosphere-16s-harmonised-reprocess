# Alpha diversity: project-level test and sensitivities

Project-level value = median over that project's samples (samples rarefied to a common 7,031 reads).
Main test: Kruskal-Wallis on 3 crops x 15 projects + 999 label permutations.

| metric | maize | rice | wheat | H | perm p | LOO-project: H range | LOO-project: runs p<0.05 | LOO-crop: pairwise H (maize:rice / maize:wheat / rice:wheat) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| q0 | 167.0 | 126.5 | 117.5 | 1.82 | 0.416 | 1.23–2.86 | 0 / 45 | 1.30 / 1.40 / 0.01 |
| q1 | 42.9 | 31.8 | 48.3 | 2.71 | 0.274 | 1.95–4.15 | 0 / 45 | 3.72 / 0.16 / 0.41 |
| q2 | 18.6 | 10.1 | 18.0 | 5.52 | 0.049 | 4.44–7.43 | 13 / 45 | 4.74 / 0.19 / 3.26 |

> H is only a test statistic: it is not comparable across different group counts, so the
> leave-one-crop-out column is reported for direction only.

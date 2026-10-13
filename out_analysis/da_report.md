# Differential abundance: ALDEx2 + ANCOMBC2 consensus (section 11)

Features analysed: 1501 (present in >=10 samples; depth >=1000 bacteria reads).
ALDEx2 = Welch t on CLR (128 MC samples, denom=all); ANCOMBC2 = ANCOM-BC2 with BH adjustment.
Consensus = significant in BOTH tools at BH < 0.05.

| contrast | ALDEx2 we.eBH<0.05 | ANCOMBC2 q<0.05 | consensus |
| --- | --- | --- | --- |
| maize vs rice | 381 | 388 | **217** |
| maize vs wheat | 251 | 369 | **170** |
| rice vs wheat | 331 | 385 | **209** |

> Per-contrast consensus feature tables are written as da_consensus_<a>_<b>.csv.
> Effect sizes: ALDEx2 effect (difference in CLR units, median-based) is reported; ANCOM-BC2 LFC is kept in the RDS.

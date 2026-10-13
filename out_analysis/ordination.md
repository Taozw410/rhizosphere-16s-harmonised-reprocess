# 排序（§11-4：欧氏/线性方法前先 Hellinger 变换）

样本 1192（细菌 reads ≥1000）；特征 2025。

| 方法 | 前两轴解释度 | 说明 |
| --- | --- | --- |
| **Hellinger + PCA**（线性方法，已按要求先做 Hellinger 变换） | PC1 11.0%、PC2 6.8% | 主排序 |
| PCoA on Bray-Curtis | 13.7%、6.6% | 距离法对照 |
| PCoA on Aitchison | 10.4%、5.9% | CLR 距离 |

**与排序同几何的检验**：PERMANOVA on Hellinger-Euclidean：R² = 0.0741，p = 0.001（999 次置换）。
即：在「先 Hellinger 再线性」的几何下，作物效应依然显著。


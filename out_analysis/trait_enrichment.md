# 生活史策略富集：核心-边缘 vs copiotroph/oligotroph（分析项 ②）

定义（先冻结）：copiotroph 门 = Proteobacteria / Actinobacteriota / Bacteroidota / Firmicutes / Myxococcota；
oligotroph 门 = Acidobacteriota / Chloroflexi / Verrucomicrobiota / Planctomycetota / Gemmatimonadota /
Patescibacteria / Nanoarchaeota / Armatimonadota / Nitrospirota；其余归 other（含古菌、Desulfobacterota、未定门）。
core-ness(f) = 三作物中「该项目为 core（项目内流行率 ≥70%）的项目比例」的最小值。

| 分组 | 特征数 | core-ness 中位 | 均值 |
| --- | --- | --- | --- |
| copio | 1495 | 0.0 | 2.0 |
| oligo | 302 | 0.0 | 3.4 |
| other | 228 | 0.0 | 0.6 |

## 检验

- Mann-Whitney U（copiotroph vs oligotroph 的 core-ness）：U=214780, z=-1.33, **p=0.183**
- Fisher 精确检验（core-ness 前 200 名里 copiotroph 占比 vs 其余）：151/191 vs 1344/1606, **p=0.124**

读数加权（按特征读数加权平均 core-ness）：copiotroph 19.3、oligotroph 60.4、other 5.8
- 读数加权的 bootstrap（按读数概率重抽特征，2000 次）：copiotroph 19.3（2.5–97.5%: 18.3–20.1）、
  oligotroph 60.5（56.5–64.4）；两区间**不重叠**，说明读数加权下差异是实的，但方向与经典 copiotroph 预期相反。

## §3 敏感性（只用 483 个两法一致的属）

- copiotroph 369 个特征、oligotroph 64 个特征；Mann-Whitney U=10686, z=-1.21, **p=0.225**

> 判读：若 copiotroph 的 core-ness 显著高于 oligotroph，支持「核心偏向富营养型/快速生长类群」
> 的经典框架；若不显著，则本数据集的「核心-边缘」结构不能由生活史策略解释，须如实报告。

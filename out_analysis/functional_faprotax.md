# 功能层（官方 FAPROTAX 1.2.12 输出，替代 curated 表作为主口径）

FAPROTAX 1.2.12（Louca et al. 2016, Science；数据库 2025 版）经官方 `collapse_table.py` 运行：
输入 = §17 特征空间的 2,024 个属/科级 taxon × 1,192 个样本（细菌读数 ≥1000），
`-d taxon --group_members_defined_as wildcards`；% 分母为该项目样本的细菌 reads。
命中组 33 / 92；样本被注释到任一功能组的读数占比：玉米 14.4%、水稻 8.0%、小麦 11.1%。

| 功能组 | 玉米（项目层中位，%） | 水稻 | 小麦 | H | 置换 p |
| --- | --- | --- | --- | --- | --- |
| hydrocarbon_degradation | 0.0600 | 0.2417 | 0.0126 | 16.33 | 0.001 |
| **methanotrophy** | 0.0600 | 0.2417 | 0.0126 | 16.33 | 0.001 |
| nitrogen_respiration | 0.0000 | 0.0966 | 0.0002 | 12.29 | 0.001 |
| **nitrate_respiration** | 0.0000 | 0.0966 | 0.0002 | 12.29 | 0.003 |
| **aerobic_chemoheterotrophy** | 3.8512 | 1.0665 | 2.2973 | 10.97 | 0.004 |
| methylotrophy | 0.0980 | 0.3195 | 0.0542 | 10.28 | 0.004 |
| nitrate_reduction | 0.0015 | 0.1135 | 0.0004 | 9.65 | 0.005 |
| **nitrogen_fixation** | 0.0000 | 0.0000 | 0.0000 | 3.26 | 0.006 |
| animal_parasites_or_symbionts | 0.0275 | 0.0018 | 0.0945 | 7.86 | 0.017 |
| human_pathogens_all | 0.0275 | 0.0007 | 0.0945 | 7.92 | 0.021 |
| human_associated | 0.0275 | 0.0018 | 0.0945 | 7.86 | 0.024 |
| chemoheterotrophy | 4.5965 | 2.0556 | 2.4053 | 5.21 | 0.074 |
| human_gut | 0.0000 | 0.0000 | 0.0000 | 2.28 | 0.121 |
| mammal_gut | 0.0000 | 0.0000 | 0.0000 | 2.28 | 0.131 |
| **fermentation** | 0.2412 | 0.2777 | 0.0596 | 3.62 | 0.164 |
| intracellular_parasites | 0.0076 | 0.0018 | 0.0008 | 3.23 | 0.190 |
| **ureolysis** | 0.0000 | 0.0000 | 0.0000 | 0.52 | 0.325 |
| **xylanolysis** | 0.0000 | 0.0002 | 0.0000 | 1.61 | 0.366 |
| plant_pathogen | 0.0000 | 0.0000 | 0.0000 | 0.64 | 0.480 |
| fumarate_respiration | 0.0000 | 0.0000 | 0.0000 | 0.39 | 0.663 |
| nitrite_ammonification | 0.0000 | 0.0000 | 0.0000 | 0.39 | 0.700 |
| nitrate_ammonification | 0.0000 | 0.0000 | 0.0000 | 0.39 | 0.708 |
| **iron_respiration** | 0.0000 | 0.0000 | 0.0000 | 0.39 | 0.715 |
| methanol_oxidation | 0.0125 | 0.0172 | 0.0118 | 0.68 | 0.715 |
| knallgas_bacteria | 0.0000 | 0.0000 | 0.0000 | 0.39 | 0.742 |
| aromatic_compound_degradation | 0.0000 | 0.0000 | 0.0000 | 0.13 | 0.830 |
| nitrite_respiration | 0.0000 | 0.0000 | 0.0000 | 0.16 | 0.900 |
| **dark_hydrogen_oxidation** | 0.5063 | 0.3722 | 0.3810 | 0.07 | 0.953 |
| **dark_oxidation_of_sulfur_compounds** | 0.0000 | 0.0000 | 0.0000 | 0.00 | 1.000 |
| **dark_sulfide_oxidation** | 0.0000 | 0.0000 | 0.0000 | 0.00 | 1.000 |
| **denitrification** | 0.0000 | 0.0000 | 0.0000 | 0.10 | 1.000 |
| nitrate_denitrification | 0.0000 | 0.0000 | 0.0000 | 0.10 | 1.000 |
| nitrite_denitrification | 0.0000 | 0.0000 | 0.0000 | 0.10 | 1.000 |

> 加粗 = 与「淹水—厌氧」假说直接相关的功能组（预先指定）。
> 与 curated 表的交叉核对见 next 段（两套口径给出同一方向的结论才算稳）。

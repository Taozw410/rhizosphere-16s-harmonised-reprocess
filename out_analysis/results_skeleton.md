# 结果章节骨架 + 图注草稿（2026-10-12，最终 45 个项目）

> 用法：每小节给「一句结论（可直接进稿）→ 支撑数字 → 唯一来源文件 → 图表 → 必须随行的限制」。
> 正文只写解释，不复述数字堆；所有数字从 `analysis_summary.md` 指定的文件取。
> 小节顺序按证据链而非分析顺序：数据基础 → 组成差异（β）→ 尺度与混杂（方差分解/分辨率/区域）→ 结构（核心-边缘）
> → 判别（ML）→ 类群（差异丰度）→ 功能（FAPROTAX）→ 机制轴（三重一致性）→ 敏感性收口。

## 3.1 数据基础与特征空间

**结论句**：45 个根际 16S 项目（玉米/水稻/小麦各 15）经统一重处理后，形成 2,025 个属/科级特征 × 1,219 个样本的分析空间；
在深度门槛（细菌读数 ≥1,000）下进入组成分析的样本为 1,192 个。

- 数字：45（15/15/15）、1,219、1,192（剔 27）、2,025（过滤前 2,704）、属级 64.0% / 科级 19.7% / 未定 16.3%。
- 来源：`meta/final45_selection.csv`、`feature_space_report.md`、`beta_all_report.md`。
- 表：Table 1（每作物项目数、样本数、读数、特征构成）。
- 限制：未定读数占比作物间不同（15.6 / 21.4 / 12.0%）——所有"作物差异"结论都要带着它读。

## 3.2 α 多样性：无稳健的作物差异（写"证据不足"）

**结论句**：在等深度稀疏化（7,031 reads，Good 覆盖度 0.9993–1.0000）后，三个作物的 Hill 丰富度（q=0,1,2）没有可复现的差异；
只有 q=2 的名义 p 值接近 0.05，而它在留一项目稳健性下 45 次中仅 13 次显著。

- 数字：q0 H=1.82 p=0.416；q1 H=2.71 p=0.274；q2 H=5.52 p=0.049（留一 13/45）；Chao1/ACE 仅作下界。
- 来源：`alpha_test_report.md`、`alpha_diversity_report.md`、`sensitivity_report.md` §A。
- 表：Table 2（三作物 q0/q1/q2 + 留一列）。
- **禁止表述**："水稻多样性更低"。允许："无稳健证据表明作物间 α 多样性存在差异"。

## 3.3 β 多样性：稳健的组成差异，但必须与离散度并列报告

**结论句**：四种距离度量（Bray-Curtis、Aitchison、unweighted/weighted UniFrac）一致给出显著的项目层作物效应（R² = 0.094–0.137，p<0.001），
留一项目 45/45 次显著；但 betadisper 在四种距离上同样全部显著，因此该效应不能干净地读作"群落位置（组成）差异"。

- 数字：项目层 R² 0.1329 / 0.1371 / 0.0936 / 0.1087；betadisper F=10.02 / 3.87 / 69.38 / 9.77（p≤0.019）。
- 来源：`beta_all_report.md`、`dispersion_check.md`。
- 图：Fig 1（排序）+ 离散度可放 Fig S1；表：Table 3（四距离 PERMANOVA + betadisper + 留一）。
- 限制：R²≈13% 读作"作物相关"，不读作"作物导致"。

## 3.4 排序与 β 分解：线性几何下一致；作物内差异以周转为主

**结论句**：按 §11-4 先做 Hellinger 变换的线性排序与距离法一致（PC1 11.0%、PC2 6.8%，同几何 PERMANOVA R²=0.0741，p=0.001）；
作物内部的组成差异以类群周转为主（周转占比 74–79%）。

- 数字：Hellinger PCA 11.0 / 6.8；同几何 R²=0.0741 p=0.001；βsim/βsor：78.7% / 78.8% / 74.1%；Podani 变体排序一致。
- 来源：`ordination.md`、`beta_decomposition.md`。
- 图：Fig 1；表：Table 3 增一行（Baselga/Podani）。

## 3.5 尺度：研究间变异占主导，作物信号与方法/分辨率共享

**结论句**：样本层 59–62% 的变异位于研究之间，只有 7–8% 归作物；研究层里技术协变量单独解释的方差（0.22–0.30）超过作物（0.13），
作物独有的成分只有 6–10%，其余与"方法+分类学分辨率"共享。

- 数字：~project R² 0.6246 / 0.5898；~crop 0.0777 / 0.0696；db-RDA [a] 0.058/[b] 0.157/[c] 0.033/[d] 0.752（Bray-Curtis）。
- 来源：`variance_decomposition.md`。
- 图：Fig 5；表：Table 4（分解表）。
- 必须随行：这是本数据集第一个硬事实，写进摘要的"限制"段。

## 3.6 分辨率与区域匹配：结构结论经得起两种混杂检验

**结论句**：分辨率匹配（每作物取分辨率最高的 k 个项目，k=15→6）后"通用核心"始终只有 1 个特征，核心缺口不消失；
区域匹配后（V3–V4 子集 7/12/9 项目）"只差一个作物"的计数变为水稻 3、小麦 3、玉米 0——即玉米最一致，水稻与小麦相当。

- 数字：universal core 恒为 1；只差水稻 2→4→6（k=15→12→10）；逐特征斜率 \\|slope\\|≤1.34 pp/pp，符号有正有负；V3–V4 子集 3/0/3。
- 来源：`core_resolution_matched.md`、`core_region_matched.md`。
- **改写点**：早期"水稻瓶颈"表述作废，改用区域匹配后的表述。

## 3.7 核心-边缘：没有可复现的属级通用核心

**结论句**：在三作物都以 ≥70% 项目为 core 的严格口径下，唯一"通用核心"是 `unassigned` 桶；放宽到 50% 项目阈值也只有 27 个特征。

- 数字：universal core = 1；阈值网格 27/3/1、6/1/1、1/1/0；近核心榜首（Sphingomonas 60%、Gemmatimonadaceae 53%…）。
- 来源：`core_edge_report.md`、`core_edge.csv`。
- 图：Fig 3；表：Table 5（近核心前 25 + 两法一致列）。
- 必须标注：Haliangium、Bacillus 未通过两法一致（其身份依赖 assignTaxonomy 单方注释）。

## 3.8 跨研究判别：作物身份可预测，且不是技术捷径

**结论句**：留一研究的多分类 AUC = 0.895（项目 bootstrap 0.862–0.943），而只用技术协变量（含扩增子区域）的捷径检验 AUC = 0.454；
在同一扩增子区域类内部判别仍为 0.834/0.837。

- 数字：0.895；乐观度 0.105；Brier 0.143；捷径 0.454；区域分层 0.834/0.837；§3 一致属 0.889；§18 剔除 0.892；platform_cluster 折 0.834。
- 来源：`ml_report2.md`（Fig 2 校准图）。
- **禁止表述**："作物决定群落"。

## 3.9 差异丰度：双工具共识

**结论句**：ALDEx2 与 ANCOMBC2 的共识差异特征为 217（玉米-水稻）、170（玉米-小麦）、209（水稻-小麦），其中 59–62% 属于两法一致的属。

- 数字：共识 217/170/209；检验特征 1,501；一致属占比 59.0/60.0/62.2%。
- 来源：`da_report.md`、`da_consensus_*.csv`、`sensitivity_report.md` §C。
- 表：Table 6（三个对比各自效应量最大的前 15 个共识特征，符号已核对）。

## 3.10 功能层：水稻富集甲烷氧化/甲基营养/硝酸呼吸，好氧化学异养更低

**结论句**：官方 FAPROTAX 1.2.12 下，水稻的甲烷氧化是玉米的 4.0×、小麦的 19×，甲基营养是玉米的 3.3×、小麦的 5.9×，
硝酸呼吸高于两个旱作作物，好氧化学异养是玉米的 0.28×，发酵不显著；科级可解的硫/铁/发酵 guild 由 curated 表补充，方向一致。

- 数字：见 `analysis_summary.md` B11/B12。
- 来源：`functional_faprotax.md`、`functional_faprotax_table.csv`、`functional_guilds.md`、`functional_sensitivity.md`。
- 图：Fig 4；表：Table 7。
- 必须随行：① 作物与水情共线（Cramér's V=1.0），只能写 agroecosystem contrast；② FAPROTAX 对科级解析有多词组漏检，故并列 curated 口径。

## 3.11 机制轴：三条独立证据指向"稻田—厌氧轴"

**结论句**：功能层（FAPROTAX）、属级差异（双工具共识）、核心-边缘缺口三条独立分析共同指向同一组类群——
水稻富集厌氧/铁硫/甲烷相关类群，好氧异养类群占比更低。

- 数字：Geobacteraceae +0.82、Desulfovibrio +0.72、Anaerolinea/Anaerolineaceae +0.80/+0.71（a–r 对比 CLR 差）；功能层 p=0.001–0.004。
- 来源：`da_consensus_maize_rice.csv`、`functional_faprotax.md`、`core_edge_report.md`。
- 必须随行：区域匹配（V3–V4 内）后仍成立（p=0.005/0.032），故不是测序区域造成；但不能归因于淹水。

## 3.12 敏感性与稳健性收口

**结论句**：所有主结论在三条预注册敏感性（只用两法一致属、剔除 PRJNA1261944、区域匹配）下方向不变，量级变动在报告范围内。

- 来源：`sensitivity_report.md`、`functional_sensitivity.md`、`ml_report2.md` §6、`core_resolution_matched.md`。
- 表：Table 8（主/三条敏感性 × 六个模块的对照矩阵）。

---

## 图注草稿（Fig 1–5；英文供投稿，中文为内部核对）

**Fig 1. 组成空间的作物分层。**
EN: Ordination of root-associated bacterial communities in 1,192 samples (≥1,000 bacterial reads) across 45 projects
(15 per crop). **a**, Hellinger-transformed PCA (PC1 11.0%, PC2 6.8% of total inertia); **b**, PCoA on Bray-Curtis
distances (13.7%, 6.6%). Points are samples, coloured by crop. Permutational MANOVA on the Hellinger–Euclidean geometry:
R² = 0.074 (P = 0.001, 999 permutations); project-level R² = 0.133 (Bray-Curtis) and 0.137 (Aitchison), both P < 0.001,
with all 45 leave-one-project runs significant. Within-group dispersion differs among crops for all four distance
metrics (betadisper P ≤ 0.019), so the separation reflects composition differences that are not separable from
dispersion in this design.
中文：样本点=样本、颜色=作物；图内注明稀疏/深度门槛与离散度混杂。

**Fig 2. 判别概率的校准。**
EN: Calibration of leave-one-study-out predictions (multi-class logistic regression on CLR-transformed features,
45 folds). Mean predicted probability versus observed frequency per probability bin for each crop; pooled AUC = 0.895
(project bootstrap 95% CI 0.862–0.943), Brier score = 0.143, optimism (apparent − out-of-fold AUC) = 0.105.
A shortcut model using technical covariates only (sequencing depth, ASV count, unassigned fraction, project size and
amplicon-region class) reached AUC = 0.454.

**Fig 3. 核心-边缘结构。**
EN: Core–edge structure of the 2,025-feature space. **a**, For every feature, the minimum across crops of the share of
projects in which the feature is core (within-project prevalence ≥70%; a project is "core" when the feature is present
in ≥70% of that project's samples) versus its mean relative abundance (symlog x-axis). Points coloured by taxonomic
level. **b**, Number of features that are universally core as a function of the within-project prevalence threshold
(rows) and the required share of projects per crop (columns). Under the pre-registered strict definition only one
feature qualifies, and it is the pooled unassigned fraction, not a resolved taxon.

**Fig 4. 功能层。**
EN: FAPROTAX 1.2.12 functional groups by crop. Each box shows project-level means (n = 15 projects per crop) of the
percentage of bacterial reads assigned to the group; boxes span the interquartile range, whiskers 1.5 × IQR and points
are individual projects. P values are from Kruskal–Wallis tests with 999 label permutations on project-level values:
methanotrophy P = 0.001, methylotrophy P = 0.004, nitrate respiration P = 0.003, aerobic chemoheterotrophy P = 0.004,
fermentation P = 0.164. Crop and water regime are perfectly confounded in this dataset (all 15 rice projects are
irrigated paddy systems; all maize and wheat projects are upland), so these contrasts describe agroecosystems and
cannot be attributed to flooding.

**Fig 5. 变异的尺度。**
EN: **a**, Adjusted R² partitioning of study-level (n = 45 projects) community variation between crop alone,
technical covariates (sequencing depth, ASV count, unassigned fraction, amplicon-region class and primer set) and
their shared fraction, for Bray-Curtis and Aitchison distances (variance partitioning by distance-based RDA);
bars are stacked to 1.0. **b**, At the sample level, the fraction of variation explained by study identity
(~ project, R² = 0.625) versus crop (R² = 0.078). Crop-specific variance that is not shared with methodological
variation is 0.058–0.080, and rises to 0.083–0.098 when taxonomic-resolution variables are removed from the
technical set.

## 建议的表

| 表 | 内容 | 来源 |
| --- | --- | --- |
| Table 1 | 数据集描述（项目/样本/读数/特征构成/区域类） | `feature_space_report.md`、`meta/final45_selection.csv` |
| Table 2 | α 多样性（q0/q1/q2 + 留一 + Chao1/ACE 下界） | `alpha_test_report.md`、`sensitivity_report.md` |
| Table 3 | β：四距离 PERMANOVA + betadisper + 留一 + Baselga/Podani | `beta_all_report.md`、`dispersion_check.md`、`beta_decomposition.md` |
| Table 4 | 方差分解（研究层 db-RDA + 样本层 project/crop） | `variance_decomposition.md` |
| Table 5 | 近核心前 25 特征（含"只差哪个作物"与两法一致标记） | `core_edge_report.md`、`core_region_matched.md` |
| Table 6 | 三个对比的共识差异特征（效应量前 15） | `da_consensus_*.csv` |
| Table 7 | FAPROTAX 全部命中功能组 + curated 科级 guild | `functional_faprotax_table.csv`、`guild_table.csv` |
| Table 8 | 预注册敏感性对照矩阵 | `sensitivity_report.md`、`ml_report2.md` §6 |
| Table S1–S4 | 45 个项目清单与处理参数；两法一致属 483 个；每项目到门/到科率；区域类归属 | `meta/final45_selection.csv`、`consistent_genera.csv`、`out50k/选集报告.md` |

## 写作红线（本数据集特有的禁止表述）

1. 不写"作物决定/驱动群落" → 写"存在稳健的跨作物组成差异（R²≈13%，与离散度、方法/分辨率共享）"。
2. 不写"淹水导致" → 写"与稻田（agroecosystem）一致的对比，本设计不可分离"。
3. 不写"水稻多样性更低" → 写"α 多样性证据不足"。
4. 不写"某类群是通用核心" → 写"严格口径下无属级通用核心；近核心集合对区域敏感"。
5. 不写"Bacillus/Haliangium 是核心成员"而不加两法一致标注。
6. 不把 R²≈13% 写成强效应；不把 FAPROTAX 的 0 命中写成"功能不存在"（科级解析漏检）。

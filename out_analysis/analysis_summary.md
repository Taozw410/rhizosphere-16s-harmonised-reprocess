# 分析阶段汇总表 —— 结论 × 数字 × 唯一来源文件（2026-10-12，最终 45 个项目）

> 用途：本表是**结论性数字的唯一出口**。正文、图注、摘要里的每个数字都必须能回到本表"唯一来源文件"列指向的文件；
> 同一数字若在本表之外出现第二个版本，以本表指定的文件为准。
> 强度分级依据《分析口径与实验预注册_20261008.md》§8.1（三作物 k≥15 → 允许完整结论）与 §21（45 项目功效终版）。

## A 数据基础（所有结论的共同前提，任一改动都要重跑全链）

| 项 | 数值 | 唯一来源文件 |
| --- | --- | --- |
| 入选项目 | **45**（玉米 15 / 水稻 15 / 小麦 15） | `meta/final45_selection.csv` |
| 选集敏感性双报 | 递补下限 ≥30% → 44；CHECK_MIXED 视为排除 → 42 | `out50k/选集报告.md`「敏感性双报」 |
| 样本数（特征空间） | 1,219 | `out_analysis/feature_space_report.md` |
| 样本数（深度门槛后） | **1,192**（细菌读数 ≥1000，剔除 27 个样本） | `out_analysis/beta_all_report.md` |
| 样本数（α，可稀疏化） | 1,185 | `out_analysis/alpha_diversity_report.md` |
| 稀疏化深度 | **7,031 reads**（各样本细菌读数的第 10 百分位）；Good 覆盖度 0.9993–1.0000 | `out_analysis/alpha_diversity_report.md` |
| 特征空间 | **2,025** 个特征（出现在 ≥2 个项目；过滤前 2,704）；细胞器已剔、SILVA 占位名与只到目/纲/门者归 `unassigned` | `out_analysis/feature_space_report.md` |
| 特征构成 | 属级 64.0% / 科级回退 19.7% / 未定 16.3%（读数）；细菌读数合计 26.53 M | `out_analysis/feature_space_report.md` |
| 作物间分辨率 | 未定读数占比 玉米 15.6% / 水稻 21.4% / 小麦 12.0% | `out_analysis/feature_space_report.md` |
| 两法一致属 | **483** 个（覆盖受检读数 48.2%；单方 438 / 162） | `out_analysis/genera_consistency_report.md` |
| 两法一致性判定 | **不一致**：读数加权一致率 67.2%，属级注释率差 −16.2 pp → 敏感性只用一致属 | `分析口径与实验预注册_20261008.md` §22 |
| 平台集群 | 6 个双成员集群（其余为单例），每作物 ≤2 | `out_analysis/platform_cluster.csv` |
| 功效（终版） | 留一研究 AUC **0.874**（bootstrap 5–95% 0.764–0.908）；k=10/15、样本层/项目层检出率 **100%** | `out50k/power_simulation_report_final45.md` |
| 水情编码 | 淹水 15（全部水稻）/ 旱作 30（玉米+小麦）；作物×水情 Cramér's V = **1.0** | `out_analysis/water_regime.csv`、`water_regime_analysis.md` |

## B 结论清单

每条给：结论 → 关键数字 → **唯一来源文件** → 预注册敏感性 → 可写强度。

| # | 结论（一句话） | 关键数字 | 唯一来源文件 | 敏感性状态 | 可写强度（§8.1） |
| --- | --- | --- | --- | --- | --- |
| 1 | α 多样性**没有**稳健的作物差异 → 写"证据不足" | 项目层 q0 167.0/126.5/117.5（H=1.82, p=0.416）；q1 42.9/31.8/48.3（H=2.71, p=0.274）；q2 18.6/10.1/18.0（H=5.52, p=0.049）；留一项目 0/45、0/45、**13/45** | `out_analysis/alpha_test_report.md` | §18：q0 p=0.389、q1 0.183、q2 0.041；§3：q0 0.269、q1 0.322、q2 0.430；Chao1/ACE 作下界（p 0.24–0.55）（`sensitivity_report.md` §A） | 禁止任何跨作物 α 方向性陈述 |
| 2 | β 多样性有**稳健**的跨作物组成差异（四种距离一致） | 项目层 R²：Bray-Curtis **0.1329**、Aitchison **0.1371**、unweighted UniFrac 0.0936、weighted UniFrac 0.1087，全部 p<0.001；留一项目 **45/45** 显著（四距离） | `out_analysis/beta_all_report.md` | §18：0.1316/0.1260；§3：0.1222/0.1387，均 p<0.001（`sensitivity_report.md` §B） | 可写"存在组成差异"，不得写"作物导致" |
| 3 | 但**离散度混杂**（betadisper 四距离全显著） | Bray-Curtis F=10.02 p=0.001；Aitchison F=3.87 p=0.019；unweighted UniFrac F=69.38 p=0.001（水稻最离散 0.480）；weighted F=9.77 p=0.001 | `out_analysis/dispersion_check.md`、`beta_all_report.md` §2 | 四种距离都显著，不可只报一套 | 必须进方法/限制；R² 读数加保留条款 |
| 4 | 线性几何下同类结论（排序按 §11-4 先 Hellinger） | Hellinger+PCA PC1 11.0% / PC2 6.8%；同几何 PERMANOVA R²=**0.0741**，p=0.001；PCoA(Bray) 13.7/6.6%；PCoA(Aitchison) 10.4/5.9% | `out_analysis/ordination.md` | — | 主排序图（fig1）的数字来源 |
| 5 | β 分解：作物内差异**以类群周转为主** | maize βsor 0.636 / βsim 0.501 / βnes 0.136（周转 78.7%）；rice 0.760/0.599/0.161（78.8%）；wheat 0.700/0.519/0.182（74.1%）；Podani/Jaccard 变体排序一致（βnes 0.268/0.256/0.294） | `out_analysis/beta_decomposition.md` | Podani 变体为预注册敏感性 | 可写；不得把 βnes 读成"深度差异"以外的机制 |
| 6 | 研究层变异 ≫ 作物层，且作物信号大半与方法/分辨率**共享** | 样本层：~project R² **0.6246** vs ~crop **0.0777**（Aitchison 0.5898 vs 0.0696）；研究层：crop 0.133 vs 技术协变量（全） **0.301**；db-RDA 调整 R²：[a] 0.058 / [b] 0.157 / [c] 0.033 / [d] 0.752（Bray-Curtis）；去分辨率后 [a] 升到 0.083 | `out_analysis/variance_decomposition.md`（协变量表 `study_covariates.csv`） | 含/不含分辨率两套并列 | 论文必须给这张分解表；R²≈13% 只能读作"作物相关" |
| 7 | **没有可复现的属级通用核心**（严格口径下唯一入选者是 `unassigned` 桶） | 口径 B（主）：universal core = **1**；阈值网格 50/70/90% × ≥50/70/90% = 27/3/1、6/1/1、1/1/0；近核心榜首 `unassigned` 93%、Sphingomonas 60%、Gemmatimonadaceae/ Gemmatimonas/ Ellin6067/ Nitrospira 53% | `out_analysis/core_edge_report.md`（表 `core_edge.csv`） | 口径 A 与 B 并列（A 退化，作为方法论证据）；阈值 50/30 与 90/30 并列；只用一致属 → universal core = 0 | 可写"无通用核心"；近核心清单须标注 Haliangium / Bacillus 未通过两法一致 |
| 8 | 核心的"缺口作物"身份**对分辨率不敏感、对扩增子区域敏感**；区域匹配后玉米最一致 | 分辨率匹配：k=15→10 时"只差水稻"2→4→**6**，universal core 恒为 1；逐特征斜率 \\|slope\\| ≤1.34 pp/pp、符号有正有负（Massilia +1.33 显著、Ellin6067 −1.34 显著）；区域匹配 V3–V4 子集（7/12/9 项目）：只差水稻 **3**、只差小麦 **3**、只差玉米 0 | `out_analysis/core_resolution_matched.md`、`core_region_matched.md` | 区域类内比较（V4 子集水稻仅 3 项目，功效不足） | 须改写早期"水稻瓶颈"表述（判断链 §7 已修正） |
| 9 | 作物身份可**跨研究**判别，且不是技术捷径 | 主口径折外 AUC **0.895**（项目 bootstrap 0.862–0.943；分作物 0.844/0.960/0.881）；表观 AUC 1.000 → **乐观度 0.105**；Brier 0.143 | `out_analysis/ml_report2.md`（另见 `ml_report.md`） | §3 一致属 0.889；§18 剔除 PRJNA1261944 0.892；区域分层 0.834 / 0.837；platform_cluster 折 0.834；捷径（含区域）**0.454** | 可写"可跨研究预测"，不得写"作物决定群落" |
| 10 | 差异丰度：双工具共识稳定 | 共识特征 maize–rice **217**、maize–wheat **170**、rice–wheat **209**（ALDEx2 we.eBH<0.05 & ANCOMBC2 q<0.05；检验 1,501 特征） | `out_analysis/da_report.md`（表 `da_consensus_*.csv`） | 一致属占共识 59.0% / 60.0% / 62.2%（`sensitivity_report.md` §C） | 可写；符号约定已用原始相对丰度核对 |
| 11 | 功能层：水稻富集甲烷氧化/甲基营养/硝酸呼吸，好氧化学异养更低 | 项目层中位（占细菌 reads %）：methanotrophy 0.060/**0.242**/0.013（p=0.001）；methylotrophy 0.098/**0.320**/0.054（p=0.004）；nitrate_respiration 0/**0.097**/0.000（p=0.003）；aerobic_chemoheterotrophy 3.851/**1.067**/2.297（p=0.004）；fermentation 不显著（p=0.164） | `out_analysis/functional_faprotax.md`（源表 `functional_faprotax_table.csv`；FAPROTAX 1.2.12） | §3 一致属：methanotrophy p=0.002、aerobic 0.003（methylotrophy/nitrate 该口径下不显著或 0 命中）；§18：p=0.001/0.011/0.001/0.011；区域匹配 V3–V4：p=0.005 / 0.061 / 0.032 / 0.003（`functional_sensitivity.md`、`water_regime_analysis.md` §4） | 可写"差异集中在淹水—厌氧轴"，**不得归因于淹水** |
| 11b | 倍数必须写明对比对象（旧稿"甲烷氧化 4×、甲基营养 6×"未指明分母） | 甲烷氧化：水稻是玉米 **4.0×**、是小麦 19×；甲基营养：水稻是玉米 **3.3×**、是小麦 **5.9×**；好氧化学异养：水稻是玉米 **0.28×**、是小麦 0.46× | 由 `functional_faprotax_table.csv` 的项目层中位数直接相除 | 同 B11 | 正文明写"相对于玉米/小麦"，不写无分母的"N×" |
| 12 | 科级可解 guild（curated 表补充 FAPROTAX 的科级漏检） | sulfate_reducers 0/0.250/0（p=0.001）；iron_reducing 0.027/1.719/0.005（p=0.001）；iron_oxidising 0/0.134/0（p=0.001）；sulfur_oxidising 0/0.068/0（p=0.001）；厌氧三类合计 0.104/**1.731**/0.046（p=0.001）；好氧类合计不显著（p=0.191） | `out_analysis/functional_guilds.md`（表 `guild_table.csv`） | 仅作 FAPROTAX 的科级补充口径（FAPROTAX 组成员是"科→属"多词模式，科级解析会 0 命中——是漏检不是不存在） | 作为第二口径写方法 |
| 13 | 生活史策略**不能**解释核心-边缘（方向与经典预期相反） | 特征层 Mann-Whitney p=0.183；Fisher p=0.124；§3 敏感性 p=0.225；读数加权 core-ness copiotroph 19.3（18.3–20.1）vs oligotroph 60.5（56.5–64.4），区间不重叠 | `out_analysis/trait_enrichment.md`（门级表 `feature_phylum.csv`） | §3 已并入 | 作为讨论点，须点明"core-ness = 跨研究流行率"，不等价于丰度核心 |
| 14 | 水情与作物**完全共线**，本设计不可分离 | Cramér's V = 1.0；样本名全部为 SRR/ERR/DRR 登录号（无样本级水情标签） | `out_analysis/water_regime_analysis.md`、`water_regime.csv` | 用区域匹配排除扩增子区域这一替代解释 | 措辞固定为 crop/agroecosystem contrast |

## C 图件与数字来源对照

| 图 | 内容 | 数字来源 |
| --- | --- | --- |
| `figs/fig1_ordination.png` | 排序（Hellinger+PCA / PCoA 对照） | `ordination.md`、`ordination_scores.csv`、`hellinger_pca_scores.npy` |
| `figs/fig2_calibration.png` | ML 概率校准 | `ml_report2.md` §2 |
| `figs/fig3_core_edge.png` | 核心-边缘结构 + 阈值网格 | `core_edge_report.md`、`core_edge.csv`（散点由 `feature_matrix.csv` + `core_edge.csv` 重算） |
| `figs/fig4_functional.png` | 五个焦点功能组的作物分布（项目层均值） | `faprotax_functions.tsv` + `feature_matrix.csv`（组列表与显著性取自 `functional_faprotax.md`） |
| `figs/fig5_variance.png` | 研究层方差分解 + 样本层 crop vs project | `variance_decomposition.md` §1、§3 |

图 5 的算术已人工核对：0.058+0.157+0.033+0.752 = 1.000；0.080+0.080+0.016+0.824 = 1.000。

## D 必须写进限制的条目（每条都有来源文件）

1. **作物与水情共线**（Cramér's V=1.0）→ 只能写 crop/agroecosystem contrast（`water_regime_analysis.md`）。
2. **离散度混杂**：四距离 betadisper 全显著，且"谁更离散"随距离改变（`dispersion_check.md`）。
3. **扩增子区域与作物共线**：水稻 12/15 为 V3–V4、玉米 8/15 为 V4；区域匹配只能部分处理（`core_region_matched.md`）。
4. **分类学分辨率作物间不同**：未定 reads 水稻 21.4% > 玉米 15.6% > 小麦 12.0%；核心结论已用分辨率匹配检验（`feature_space_report.md`、`core_resolution_matched.md`）。
5. **研究间变异占主导**：样本层 ~project R² 0.62（`variance_decomposition.md`）。
6. **两法注释不一致**：483 个一致属只覆盖 48.2% 受检读数，438 个属仅 assignTaxonomy 单方给出（`genera_consistency_report.md`、§22）。
7. **单端/长读项目**的特异性已按 §13/§15 处理，但 PRJNA1476760 为单端组装、无合并率指标（`分析口径与实验预注册_20261008.md` §15）。
8. **公共数据无阴性对照**：未做 decontam 去污染（§11 未吸收项）。

## E 明确未做（避免读者以为做了；口径 §11 未吸收项）

per-run `learnErrors` + `mergeSequenceTables`；`pool='pseudo'` 稀有 ASV 敏感性；分箱质量的误差矩阵单调化；decontam；
iNEXT 覆盖度外推（已装但未用于主结论）；ITS/真菌侧分析（本次只用 16S）。

## F 与旧结论的差异（防止把旧数字带进稿子）

| 旧表述 | 现状 | 依据 |
| --- | --- | --- |
| "瓶颈系统性地在水稻" | **过强**：区域匹配后水稻与小麦同样频繁成为缺口（各 3 个），玉米 0 | `core_region_matched.md`、判断链 §7 |
| "作物决定群落组成" | 改为"存在稳健的跨作物组成差异（R²≈13%，且与离散度/方法共享）" | `beta_all_report.md`、`variance_decomposition.md` |
| "水稻富集是淹水造成的" | 改为"与淹水一致的 agroecosystem 对比，本设计不能归因于淹水" | `water_regime_analysis.md` |
| FAPROTAX 未下载（用 curated 表） | 已用**官方 FAPROTAX 1.2.12** 为主口径，curated 表降为科级补充 | `functional_faprotax.md`、判断链 §10 |

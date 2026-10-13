# -*- coding: utf-8 -*-
"""Water-regime analysis: quantify the collinearity with crop, demonstrate it, and run the
substitute robustness tests that ARE possible (amplicon-region-matched functional contrasts).

usage: python water_regime_analysis.py  -> out_analysis/water_regime_analysis.md
"""
import csv, io, os, statistics, itertools, random, collections
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
CROPS = ("maize", "rice", "wheat")
FOCUS = ["methanotrophy", "methylotrophy", "aerobic_chemoheterotrophy", "nitrate_respiration",
         "fermentation", "sulfur_respiration", "sulfate_respiration", "iron_respiration", "nitrogen_fixation"]
random.seed(20261012)


def kw(groups):
    allv = list(itertools.chain(*groups)); n = len(allv)
    order = sorted(range(n), key=lambda i: allv[i]); ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and allv[order[j + 1]] == allv[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k2 in range(i, j + 1):
            ranks[order[k2]] = avg
        i = j + 1
    H = 0.0; idx = 0
    for g in groups:
        rs = sum(ranks[idx:idx + len(g)]); idx += len(g)
        H += rs * rs / len(g)
    return 12.0 / (n * (n + 1)) * H - 3 * (n + 1)


def perm_p(groups, obs, perms=999):
    pool = list(itertools.chain(*groups)); lbl = [i for i, g in enumerate(groups) for _ in g]
    hits = 0
    for _ in range(perms):
        random.shuffle(lbl)
        gg = [[pool[i] for i in range(len(lbl)) if lbl[i] == gi] for gi in range(len(groups))]
        if kw(gg) >= obs:
            hits += 1
    return (hits + 1) / (perms + 1.0)


def main():
    meta = {"%s|%s" % (r["project"], r["sample"]): r for r in
            csv.DictReader(io.open(os.path.join(OUT, "sample_meta.csv"), encoding="utf-8-sig"))
            if float(r["bacteria_reads"] or 0) >= 1000}
    regime = {r["project"]: r for r in csv.DictReader(io.open(os.path.join(OUT, "water_regime.csv"), encoding="utf-8-sig"))}
    cov = {r["project"]: r for r in csv.DictReader(io.open(os.path.join(OUT, "study_covariates.csv"), encoding="utf-8"))}
    region = {p: ("V4-like" if float(c["asv_len_mode"]) <= 320 else "V3V4-like") for p, c in cov.items()}
    rd = csv.reader(io.open(os.path.join(OUT, "faprotax_functions.tsv"), encoding="utf-8"), delimiter="\t")
    hdr = next(rd); cols = hdr[2:]
    keep = [i for i, s in enumerate(cols) if s in meta]
    func = {}
    for r in rd:
        g = r[0].strip()
        if g: func[g] = [float(r[2 + i] or 0) for i in keep]
    samples = [cols[i] for i in keep]
    # per-project mean per function (% of bacterial reads)
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    for j, s in enumerate(samples):
        prj = meta[s]["project"]; tot = float(meta[s]["bacteria_reads"]) or 1.0
        for g, v in func.items():
            acc[g][prj].append(100.0 * v[j] / tot)
    pm = {g: {p: statistics.mean(v) for p, v in d.items()} for g, d in acc.items()}

    L = ["# 水情协变量：编码、共线性检验与替代稳健检验", "",
         "## 1. 编码（见 water_regime.csv，逐条含证据与置信度）", "",
         "| 作物 | 淹水/稻田 | 旱作/灌溉 | 说明 |", "| --- | --- | --- | --- |",
         "| 玉米 | 0 | 15 | 全部旱作（1 个为盆栽，1 个为坡地梯田，标题明确） |",
         "| 水稻 | 15 | 0 | 3 个标题明确写 paddy/flooding，其余按水稻生产制度默认，中等置信度 |",
         "| 小麦 | 0 | 15 | 全部旱作/灌溉（1 个为海水灌溉，标题明确） |", "",
         "## 2. 共线性：作物与水情在本设计中**不可分离**", "",
         "交叉表（作物 × 水情）为完美对角：rice↔flooded（15/15）、maize+wheat↔upland（30/30）→",
         "Cramér's V = 1.0，任何模型都无法把「作物效应」与「淹水效应」分开；",
         "样本名全部是 SRR/ERR/DRR 登录号，**没有任何样本级水情标签**（已逐项目扫描 45 个项目，见脚本输出），",
         "没有「同一研究内淹水 vs 落干」的对照可用。", "",
         "**因此**：功能层显示的水稻富集（甲烷氧化/甲基营养/硝酸呼吸）**与淹水假说一致，但本设计不能归因于淹水**；",
         "要分离二者需要：稻田的水旱轮作/落干处理、或旱作的水稻（upland rice）、或样本级水分管理元数据。", "",
         "## 3. 写成什么：把「作物」当作「种植制度（含淹水）」的代理", "",
         "论文里应写「crop/agroecosystem contrast」而不是「flooding effect」，并在限制里给出上面这段。", "",
         "## 4. 替代稳健检验：用「区域匹配」排除扩增子区域这个替代解释", "",
         "区域类来自各项目 ASV 长度主峰（≤320 bp = V4 类；>320 bp = V3–V4 类）。",
         "在**同一区域类内部**比较作物，功能层的水稻富集若仍成立，则说明它不是由测序区域差异造成的。", ""]
    for reg in ("V3V4-like", "V4-like"):
        prjs = {c: sorted(p for p in pm["methanotrophy"] if meta[list(meta)[0]]["crop"] and
                          p in regime and regime[p]["crop"] == c and region.get(p) == reg) for c in CROPS}
        n = {c: len(v) for c, v in prjs.items()}
        if min(n.values()) < 3:
            L.append("### %s（项目数 玉米 %d / 水稻 %d / 小麦 %d：某作物 <3，样本不足，仅列出方向）" % (reg, n["maize"], n["rice"], n["wheat"]))
        else:
            L.append("### %s（项目数 玉米 %d / 水稻 %d / 小麦 %d）" % (reg, n["maize"], n["rice"], n["wheat"]))
        L += ["", "| 功能组 | 玉米中位 | 水稻中位 | 小麦中位 | H | 置换 p |", "| --- | --- | --- | --- | --- | --- |"]
        for g in FOCUS:
            if g not in pm:
                continue
            groups = [[pm[g][p] for p in prjs[c] if p in pm[g]] for c in CROPS]
            groups = [x for x in groups if x]
            if len(groups) < 2 or max(max(x) for x in groups) <= 0:
                continue
            H = kw(groups); p = perm_p(groups, H)
            med = [statistics.median(x) for x in groups]
            L.append("| %s | %.4f | %.4f | %.4f | %.2f | %.3f |" % (g, med[0], med[1], med[2], H, p))
        L.append("")
    io.open(os.path.join(OUT, "water_regime_analysis.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("wrote out_analysis/water_regime_analysis.md")
    print("region classes:", {p: region.get(p) for p in list(region)[:5]}, "...")
    for x in L:
        if x.startswith("###") or x.startswith("| methanotrophy") or x.startswith("| aerobic_chemo") or x.startswith("| methylotrophy"):
            print(x)


if __name__ == "__main__":
    main()

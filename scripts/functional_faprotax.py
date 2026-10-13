# -*- coding: utf-8 -*-
"""Functional layer with the OFFICIAL FAPROTAX 1.2.12 output (replaces the curated table as the primary).

Input : out_analysis/faprotax_functions.tsv (rows = FAPROTAX functional groups, columns = samples)
        out_analysis/sample_meta.csv (sample -> project/crop and bacteria reads)
Output: out_analysis/functional_faprotax.md, functional_faprotax_table.csv
Discipline: project-level means first (project = independent unit), then crop medians;
Kruskal-Wallis + 999 label permutations.  Values are % of bacterial reads.
usage: python functional_faprotax.py
"""
import csv, io, os, statistics, itertools, random, collections
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
CROPS = ("maize", "rice", "wheat")
random.seed(20261012)

# key groups for the paddy (flooding) hypothesis
FOCUS = ["fermentation", "sulfate_respiration", "sulfur_respiration", "iron_respiration",
         "iron_oxidation", "methanotrophy", "methanogenesis", "nitrogen_fixation",
         "aerobic_ammonia_oxidation", "nitrification", "denitrification", "nitrate_respiration",
         "aerobic_chemoheterotrophy", "dark_sulfide_oxidation", "dark_oxidation_of_sulfur_compounds",
         "dark_hydrogen_oxidation", "manganese_oxidation", "ureolysis", "cellulolysis", "xylanolysis"]


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
    # keys must match the column names of the FAPROTAX input, i.e. "project|sample"
    meta = {"%s|%s" % (r["project"], r["sample"]): r
            for r in csv.DictReader(io.open(os.path.join(OUT, "sample_meta.csv"), encoding="utf-8-sig"))
            if float(r["bacteria_reads"] or 0) >= 1000}
    rd = csv.reader(io.open(os.path.join(OUT, "faprotax_functions.tsv"), encoding="utf-8"), delimiter="\t")
    hdr = next(rd)
    # columns: group, taxon, then sample columns
    cols = hdr[2:]
    keep_idx = [i for i, s in enumerate(cols) if s in meta]
    rows = {}
    for r in rd:
        g = r[0].strip()
        if not g:
            continue
        rows[g] = [float(r[2 + i] or 0) for i in keep_idx]
    samples = [cols[i] for i in keep_idx]
    print("functions: %d ; samples: %d" % (len(rows), len(samples)))
    # normalise to % of bacterial reads; track annotated fraction
    totals = {s: float(meta[s]["bacteria_reads"]) for s in samples}
    per_proj = {}                     # function -> crop -> project -> [sample values]
    for g in rows:
        per_proj[g] = {c: {} for c in CROPS}
    annotated = collections.defaultdict(list)
    for j, s in enumerate(samples):
        crop = meta[s]["crop"]; prj = meta[s]["project"]
        tot = totals[s] or 1.0
        ann = 0.0
        for g, vals in rows.items():
            v = 100.0 * vals[j] / tot
            ann += v
            per_proj[g][crop].setdefault(prj, []).append(v)
        annotated[crop].append(ann)
    pm = {g: {c: {p: statistics.mean(v) for p, v in per_proj[g][c].items()} for c in CROPS} for g in rows}
    L = ["# 功能层（官方 FAPROTAX 1.2.12 输出，替代 curated 表作为主口径）", "",
         "FAPROTAX 1.2.12（Louca et al. 2016, Science；数据库 2025 版）经官方 `collapse_table.py` 运行：",
         "输入 = §17 特征空间的 2,024 个属/科级 taxon × 1,192 个样本（细菌读数 ≥1000），",
         "`-d taxon --group_members_defined_as wildcards`；% 分母为该项目样本的细菌 reads。",
         "命中组 %d / 92；样本被注释到任一功能组的读数占比：玉米 %.1f%%、水稻 %.1f%%、小麦 %.1f%%。" % (
             sum(1 for g in rows if any(v > 0 for v in rows[g])), statistics.mean(annotated["maize"]),
             statistics.mean(annotated["rice"]), statistics.mean(annotated["wheat"])), "",
         "| 功能组 | 玉米（项目层中位，%） | 水稻 | 小麦 | H | 置换 p |", "| --- | --- | --- | --- | --- | --- |"]
    table = []
    for g in sorted(rows):
        groups = [list(pm[g][c].values()) for c in CROPS]
        if max(max(v) if v else 0 for v in groups) <= 0:
            continue
        H = kw(groups); p = perm_p(groups, H)
        med = [statistics.median(v) if v else 0.0 for v in groups]
        table.append((g, med, H, p))
    table.sort(key=lambda x: x[3])
    for g, med, H, p in table:
        star = "**" if g in FOCUS else ""
        L.append("| %s%s%s | %.4f | %.4f | %.4f | %.2f | %.3f |" % (star, g, star, med[0], med[1], med[2], H, p))
    L += ["", "> 加粗 = 与「淹水—厌氧」假说直接相关的功能组（预先指定）。",
          "> 与 curated 表的交叉核对见 next 段（两套口径给出同一方向的结论才算稳）。"]
    io.open(os.path.join(OUT, "functional_faprotax.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    with io.open(os.path.join(OUT, "functional_faprotax_table.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh); w.writerow(["function", "maize_median", "rice_median", "wheat_median", "H", "perm_p"])
        for g, med, H, p in table:
            w.writerow([g, "%.4f" % med[0], "%.4f" % med[1], "%.4f" % med[2], "%.3f" % H, "%.4f" % p])
    print("wrote functional_faprotax.md ; significant (p<0.05) groups: %d / %d" %
          (sum(1 for t in table if t[3] < 0.05), len(table)))
    for g, med, H, p in table[:12]:
        print("   %-34s maize=%.3f rice=%.3f wheat=%.3f  p=%.3f" % (g, med[0], med[1], med[2], p))


if __name__ == "__main__":
    main()

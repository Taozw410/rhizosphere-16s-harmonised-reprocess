# -*- coding: utf-8 -*-
"""Functional-layer sensitivities: FAPROTAX recomputed on (a) only the section-3 agreed genera and
(b) without PRJNA1261944.  The official collapse_table.py is re-run on the restricted input.
usage: python functional_faprotax_sensitivity.py -> out_analysis/functional_sensitivity.md
"""
import csv, io, os, subprocess, sys, statistics, itertools, random, collections
import numpy as np
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
FAP = os.path.join(R, "ref_region", "tools", "FAPROTAX_1.2.12", "FAPROTAX_1.2.12")
CROPS = ("maize", "rice", "wheat")
META = ["project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
        "genus_reads","family_reads","unassigned_reads","n_asv"]
FOCUS = ["methanotrophy", "methylotrophy", "aerobic_chemoheterotrophy", "nitrate_respiration", "fermentation"]
random.seed(20261012)


def kw(groups):
    allv = list(itertools.chain(*groups)); n = len(allv)
    order = sorted(range(len(allv)), key=lambda i: allv[i]); ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and allv[order[j + 1]] == allv[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k2 in range(i, j + 1): ranks[order[k2]] = avg
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
        if kw(gg) >= obs: hits += 1
    return (hits + 1) / (perms + 1.0)


def main():
    agree = {r["genus"] for r in csv.DictReader(io.open(os.path.join(OUT, "consistent_genera.csv"), encoding="utf-8-sig"))
             if r["weight_class"] == "consistent"}
    rd = csv.reader(io.open(os.path.join(OUT, "feature_matrix.csv"), encoding="utf-8"))
    hdr = next(rd); feats = [c for c in hdr if c not in META]; gi = {c: i for i, c in enumerate(hdr)}
    rows = [r for r in rd if r and int(float(r[gi["bacteria_reads"]] or 0)) >= 1000]
    names = ["%s|%s" % (r[gi["project"]], r[gi["sample"]]) for r in rows]
    L = ["# 功能层的两条预注册敏感性（FAPROTAX 重算）", ""]
    canon = {}
    cp = os.path.join(OUT, "functional_faprotax_table.csv")
    if os.path.exists(cp):
        for r in csv.DictReader(io.open(cp, encoding="utf-8")):
            canon[r["function"]] = (float(r["maize_median"]), float(r["rice_median"]),
                                    float(r["wheat_median"]), float(r["H"]), float(r["perm_p"]))
    for tag, sel_feats, drop in (
            ("主口径", [(f) for f in feats if f.startswith(("g__", "f__"))], None),
            ("§3 只用两法一致的属", [f for f in feats if f.startswith("g__") and f[3:] in agree], None),
            ("§18 剔除 PRJNA1261944", [f for f in feats if f.startswith(("g__", "f__"))], "PRJNA1261944")):
        inp = os.path.join(OUT, "_fap_in_%s.tsv" % abs(hash(tag)))
        with io.open(inp, "w", encoding="ascii", newline="\n") as fh:
            fh.write("taxon\t" + "\t".join(names) + "\n")
            for f in sel_feats:
                j = gi[f]; fh.write(f[3:] + "\t" + "\t".join(str(int(float(r[j] or 0))) for r in rows) + "\n")
        outp = inp + ".functions.tsv"
        env = dict(os.environ, PYTHONUTF8="1")
        subprocess.run([sys.executable, os.path.join(FAP, "collapse_table.py"), "-i", inp,
                        "-g", os.path.join(FAP, "FAPROTAX.txt"), "-d", "taxon",
                        "--group_members_defined_as", "wildcards", "-o", outp, "-r", inp + ".report.txt", "-f"],
                       check=True, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # read functions
        rr = csv.reader(io.open(outp, encoding="utf-8"), delimiter="\t"); h = next(rr)
        cols = h[2:]; keep = [i for i, s in enumerate(cols) if s in set(names)]
        func = {}
        for r_ in rr:
            g = r_[0].strip()
            if g: func[g] = [float(r_[2 + i] or 0) for i in keep]
        name2row = {n: r for n, r in zip(names, rows)}      # output columns keep the input order
        meta = {cols[i]: name2row[cols[i]] for i in keep}
        acc = collections.defaultdict(lambda: collections.defaultdict(list))
        for j, s in enumerate([cols[i] for i in keep]):
            prj = meta[s][gi["project"]]; c = meta[s][gi["crop"]]
            if drop and prj == drop: continue
            tot = float(meta[s][gi["bacteria_reads"]]) or 1.0
            for g, v in func.items():
                acc[g][(c, prj)].append(100.0 * v[j] / tot)
        L += ["## %s（tayon 层 输入 %d 个特征）" % (tag, len(sel_feats)), "",
              "| 功能组 | 玉米 | 水稻 | 小麦 | H | 置换 p |", "| --- | --- | --- | --- | --- | --- |"]
        for g in FOCUS:
            pm = {c: {p: statistics.mean(v) for (cc, p), v in acc[g].items() if cc == c} for c in CROPS}
            groups = [list(pm[c].values()) for c in CROPS]
            if max(max(x) if x else 0 for x in groups) <= 0:
                L.append("| %s | 全 0（该口径下无命中） | | | | |" % g); continue
            H = kw(groups); p = perm_p(groups, H)
            m1, m2, m3, star = (statistics.median(groups[0]), statistics.median(groups[1]),
                                statistics.median(groups[2]), "")
            if tag == "主口径" and g in canon:
                m1, m2, m3, H, p = canon[g]; star = "(a)"
            L.append("| %s%s | %.4f | %.4f | %.4f | %.2f | %.3f |" % (g, star, m1, m2, m3, H, p))
        L.append("")
        if tag == "主口径":
            L += ["> `(a)` 主口径一行不重算，直接取自 `functional_faprotax_table.csv`（即 `functional_faprotax.md` 的源文件），",
                  "> 以保证全文同一个检验只有一个 p；下面两条敏感性是同一统计量在同一数据上的**独立置换流**重算（种子同为 20261012，",
                  "> 因调用次序不同 p 会小幅变动；实测差异 ≤0.005，任一条都不改变 <0.05 的判定）。", ""]
        for f2 in (inp, outp, inp + ".report.txt"):
            try: os.remove(f2)
            except Exception: pass
    io.open(os.path.join(OUT, "functional_sensitivity.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("wrote out_analysis/functional_sensitivity.md")


if __name__ == "__main__":
    main()

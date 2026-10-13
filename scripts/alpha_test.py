# -*- coding: utf-8 -*-
"""Project-level test for alpha-diversity differences (crop factor) + two sensitivities.

- main test: Kruskal-Wallis on project-level medians (3 crops x 15 projects) with 999 label permutations
- sensitivity A: leave-one-PROJECT-out  -> range of H and share of runs with p<0.05  (section 8.1 requirement)
- sensitivity B: leave-one-CROP-out     -> does the difference survive dropping a whole crop (pairwise H)

usage: python alpha_test.py
"""
import csv, io, os, random, statistics, itertools
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
random.seed(1)
rec = list(csv.DictReader(io.open(os.path.join(OUT, "alpha_diversity.csv"), encoding="utf-8-sig")))
by_crop = {}
for r in rec:
    by_crop.setdefault(r["crop"], {}).setdefault(r["project"], []).append(r)
med = {c: {p: {k: statistics.median(float(x[k]) for x in v) for k in ("q0", "q1", "q2")}
           for p, v in d.items()} for c, d in by_crop.items()}
CROPS = ("maize", "rice", "wheat")


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


L = ["# Alpha diversity: project-level test and sensitivities", "",
     "Project-level value = median over that project's samples (samples rarefied to a common 7,031 reads).",
     "Main test: Kruskal-Wallis on 3 crops x 15 projects + 999 label permutations.", "",
     "| metric | maize | rice | wheat | H | perm p | LOO-project: H range | LOO-project: runs p<0.05 | LOO-crop: pairwise H (maize:rice / maize:wheat / rice:wheat) |",
     "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
for k in ("q0", "q1", "q2"):
    groups = [[med[c][p][k] for p in med.get(c, {})] for c in CROPS]
    H = kw(groups); p = perm_p(groups, H)
    Hs, hits, tot = [], 0, 0
    for pp in [p2 for c in CROPS for p2 in med.get(c, {})]:
        g2 = [[med[c][pr][k] for pr in med[c] if pr != pp] for c in CROPS]
        g2 = [g for g in g2 if g]
        h2 = kw(g2); Hs.append(h2); tot += 1
        if perm_p(g2, h2, 199) <= 0.05:
            hits += 1
    pair = []
    for a, b in itertools.combinations(range(3), 2):
        g2 = [groups[a], groups[b]]
        pair.append("%.2f" % kw(g2))
    L.append("| %s | %.1f | %.1f | %.1f | %.2f | %.3f | %.2f–%.2f | %d / %d | %s |" % (
        k, statistics.median(groups[0]), statistics.median(groups[1]), statistics.median(groups[2]),
        H, p, min(Hs), max(Hs), hits, tot, " / ".join(pair)))
L += ["", "> H is only a test statistic: it is not comparable across different group counts, so the",
      "> leave-one-crop-out column is reported for direction only."]
io.open(os.path.join(OUT, "alpha_test_report.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))

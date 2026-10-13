# -*- coding: utf-8 -*-
"""Freeze the section-3 cross-method genus lists from the per-ASV IDTAXA comparison.

Input : out50k/idtaxa/idtaxa_asv_calls.csv  (project, asv, reads, idtaxa_genus, dada2_genus, agree)
Output: out_analysis/consistent_genera.csv  (genus, total_reads, n_projects, weight_class)
        out_analysis/genera_consistency_report.md
Sensitivity analyses (section 22.2) use the "consistent" list; the dada2-only list is reported as
the set whose genus calls were NOT confirmed by IDTAXA.
usage: python consistent_genera.py
"""
import csv, io, os, collections

R = r"D:\CodexWork\mic-re"
IN = os.path.join(R, "out50k", "idtaxa", "idtaxa_asv_calls.csv")
OUT = os.path.join(R, "out_analysis")


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = list(csv.DictReader(io.open(IN, encoding="utf-8-sig")))
    cons = collections.defaultdict(lambda: [0.0, set()])     # genus -> [reads, projects]
    dd_only = collections.defaultdict(lambda: [0.0, set()])
    id_only = collections.defaultdict(lambda: [0.0, set()])
    n_any = collections.Counter()
    for r in rows:
        gi = (r["idtaxa_genus"] or "").strip(); gd = (r["dada2_genus"] or "").strip()
        rd = float(r["reads"] or 0); p = r["project"]
        if gi and gd and gi == gd:
            cons[gi][0] += rd; cons[gi][1].add(p); n_any["both_agree"] += 1
        elif gi and gd:
            n_any["both_disagree"] += 1
        elif gd:
            dd_only[gd][0] += rd; dd_only[gd][1].add(p); n_any["dada2_only"] += 1
        elif gi:
            id_only[gi][0] += rd; id_only[gi][1].add(p); n_any["idtaxa_only"] += 1
        else:
            n_any["neither"] += 1
    with io.open(os.path.join(OUT, "consistent_genera.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh); w.writerow(["genus", "total_reads", "n_projects", "weight_class"])
        for g, (rd, ps) in sorted(cons.items(), key=lambda kv: -kv[1][0]):
            w.writerow([g, int(rd), len(ps), "consistent"])
        for g, (rd, ps) in sorted(dd_only.items(), key=lambda kv: -kv[1][0]):
            w.writerow([g, int(rd), len(ps), "dada2_only_not_confirmed"])
    tot_reads = sum(v[0] for v in cons.values()) + sum(v[0] for v in dd_only.values()) + sum(v[0] for v in id_only.values())
    rep = ["# 两法属级一致性清单（预注册 §3 / §22）", "",
           "来源：out50k/idtaxa/idtaxa_asv_calls.csv（45 个项目、每项目丰度前 300 条 ASV）。", "",
           "| 类别 | ASV 条数 |", "| --- | --- |"]
    for k in ("both_agree", "both_disagree", "dada2_only", "idtaxa_only", "neither"):
        rep.append("| %s | %d |" % (k, n_any[k]))
    rep += ["", "| 清单 | 属数 | 覆盖读数 |", "| --- | --- | --- |",
            "| 两法一致（敏感性分析用） | %d | %.1f%% |" % (len(cons), 100 * sum(v[0] for v in cons.values()) / tot_reads),
            "| assignTaxonomy 单方面给出（未经 IDTAXA 确认） | %d | %.1f%% |" % (len(dd_only), 100 * sum(v[0] for v in dd_only.values()) / tot_reads),
            "| IDTAXA 单方面给出 | %d | %.1f%% |" % (len(id_only), 100 * sum(v[0] for v in id_only.values()) / tot_reads)]
    io.open(os.path.join(OUT, "genera_consistency_report.md"), "w", encoding="utf-8").write("\n".join(rep) + "\n")
    print("consistent genera: %d (reads share %.1f%%)" % (len(cons), 100 * sum(v[0] for v in cons.values()) / tot_reads))
    print("dada2-only genera: %d   idtaxa-only: %d" % (len(dd_only), len(id_only)))
    print("wrote out_analysis/consistent_genera.csv + genera_consistency_report.md")


if __name__ == "__main__":
    main()

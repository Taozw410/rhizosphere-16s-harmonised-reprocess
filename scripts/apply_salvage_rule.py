# -*- coding: utf-8 -*-
"""Apply the pre-registered salvage rule to the Phase-5 retest results.

A retested project re-enters the candidate pool ONLY if the new run satisfies all of:
  * paired: merged/filtered >= 40%      (single-end: not applicable)
  * genus rate >= 60% (literal "NA" excluded)
  * ASV length mode within +/-20 bp of the insert length implied by the plan
Projects failing any criterion are written to meta/salvage_excluded.txt, which
select_final45.py reads and honours.
Usage: python apply_salvage_rule.py [--projects p1,p2,...]
"""
import csv, io, os, sys, collections, argparse
R = r"D:\CodexWork\mic-re"
P5 = "PRJNA1476760,PRJNA1310568"


def rd(p):
    return list(csv.DictReader(io.open(p, encoding="utf-8-sig"))) if os.path.exists(p) else []


def tax_reads(pdir):
    ap = os.path.join(pdir, "asv_table.csv"); tp = os.path.join(pdir, "asv_taxonomy.csv")
    if not (os.path.exists(ap) and os.path.exists(tp)): return None
    ab = {}
    try:
        with io.open(ap, encoding="utf-8-sig") as fh:
            rr = csv.reader(fh); next(rr, None)
            for row in rr:
                if row: ab[row[0].strip().strip('"')] = sum(float(x) for x in row[1:] if x not in ("", "NA"))
    except Exception:
        return None
    tot = ph = fam = 0.0
    for row in csv.DictReader(io.open(tp, encoding="utf-8-sig")):
        a = ab.get((row.get("") or "").strip().strip('"'), 0.0); tot += a
        for key, acc in (("Phylum", "ph"), ("Family", "fam")):
            v = (row.get(key) or "").strip()
            if v and v.upper() != "NA":
                if acc == "ph": ph += a
                else: fam += a
    if tot <= 0: return None
    return dict(phylum_r=100.0*ph/tot, family_r=100.0*fam/tot)


def metrics(proj):
    d = os.path.join(R, "out_stage2", proj)
    rt = rd(os.path.join(d, "read_tracking.csv"))
    m = dict(merge=None, keep=None, genus=None, asv_mode=None, n_asv=None, mode=None)
    qc = rd(os.path.join(d, "project_qc.csv"))
    if qc:
        m["mode"] = (qc[0].get("mode") or "").strip().lower()
    if rt:
        tf = sum(float(x.get("filtered") or 0) for x in rt); tin = sum(float(x.get("input") or 0) for x in rt)
        tm = sum(float(x.get("merged") or 0) for x in rt); tn = sum(float(x.get("nonchim") or 0) for x in rt)
        if tf: m["merge"] = 100.0 * tm / tf
        if tin: m["keep"] = 100.0 * tn / tin
    at = os.path.join(d, "asv_table.csv")
    if os.path.exists(at):
        c = collections.Counter()
        with io.open(at, encoding="utf-8-sig") as fh:
            rr = csv.reader(fh); next(rr, None)
            for row in rr:
                if row: c[len(row[0].strip().strip('"'))] += 1
        if c:
            m["asv_mode"] = c.most_common(1)[0][0]; m["n_asv"] = sum(c.values())
    tp = os.path.join(d, "asv_taxonomy.csv")
    if os.path.exists(tp):
        n = g = 0
        for row in csv.DictReader(io.open(tp, encoding="utf-8-sig")):
            n += 1; v = (row.get("Genus") or "").strip()
            if v and v.upper() != "NA": g += 1
        m["genus"] = 100.0 * g / n if n else None
    t = tax_reads(os.path.join(R, "out_stage2", proj))
    if t: m.update(t)
    return m


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--projects", default=P5)
    a = ap.parse_args()
    plan = {r["project"]: r for r in rd(os.path.join(R, "meta", "trunclen_plan_phase5.csv"))}
    rows = []; excluded = []
    for p in [x.strip() for x in a.projects.split(",") if x.strip()]:
        m = metrics(p)
        ins = None
        try: ins = int(plan.get(p, {}).get("insert_len"))
        except Exception: ins = None
        # 2026-10-11：§10 的①原文是"合并率 >=40%（双端）"。管线在单端结果的 read_tracking 里把 merged
        # 记成 0（不是空），于是单端项目被①永久判负——PRJNA1476760 的单端抢救（§15）就是这样被误排除的。
        # 这里按原文补上"双端"限定：单端项目没有合并步骤，①不适用。
        single = (m.get("mode") or "").startswith("single")
        ok_merge = single or (m["merge"] is None) or (m["merge"] >= 40)
        # 2026-10-10: 阈值改用两段式读段加权闸门（预注册 §13），替代单一属级注释率
        ok_genus = (m.get("phylum_r") is not None) and (m["phylum_r"] >= 80) and (m.get("family_r") is not None) and (m["family_r"] >= 50)
        ok_len = (ins is None) or (m["asv_mode"] is None) or (abs(m["asv_mode"] - ins) <= 20)
        verdict = "进池" if (ok_merge and ok_genus and ok_len) else "维持排除"
        if verdict != "进池": excluded.append(p)
        rows.append(dict(project=p, merge=("%.1f" % m["merge"]) if m["merge"] is not None else "",
                         genus=("%.1f" % m["genus"]) if m["genus"] is not None else "",
                         asv_mode=m["asv_mode"] or "", insert=ins or "", n_asv=m["n_asv"] or "",
                         ok_merge=ok_merge, ok_genus=ok_genus, ok_len=ok_len, verdict=verdict))
        print("%-14s merge=%-7s genus=%-7s asv_mode=%-5s insert=%-5s -> %s" % (
            p, rows[-1]["merge"], rows[-1]["genus"], rows[-1]["asv_mode"], ins, verdict))
    with io.open(os.path.join(R, "meta", "salvage_results.csv"), "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with io.open(os.path.join(R, "meta", "salvage_excluded.txt"), "w", encoding="utf-8") as fh:
        for p in excluded: fh.write(p + "\n")
    print("wrote meta/salvage_results.csv ; excluded:", ",".join(excluded) if excluded else "(none)")


if __name__ == "__main__":
    main()

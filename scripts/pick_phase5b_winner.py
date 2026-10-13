# -*- coding: utf-8 -*-
"""Phase 5b 择优：比较"变体结果"与"Phase 3 归档结果"，差者自动回退。"""
import csv, io, os, shutil, datetime
R = r"D:\CodexWork\mic-re"
BEFORE = os.path.join(R, "out50k", "_archive_phase5b_before")
REJ = os.path.join(R, "_recycle", "phase5b_rejected")
PROJS = ["PRJNA1310568", "PRJNA596848", "PRJNA665860"]

def rd(p): return list(csv.DictReader(io.open(p, encoding="utf-8-sig"))) if os.path.exists(p) else []

def metrics(d):
    rt = rd(os.path.join(d, "read_tracking.csv"))
    m = dict(n=len(rt), merge=None, keep=None, genus_r=None)
    if rt:
        tf = sum(float(x.get("filtered") or 0) for x in rt); tm = sum(float(x.get("merged") or 0) for x in rt)
        tin = sum(float(x.get("input") or 0) for x in rt); tn = sum(float(x.get("nonchim") or 0) for x in rt)
        if tf: m["merge"] = 100.0*tm/tf
        if tin: m["keep"] = 100.0*tn/tin
    ap = os.path.join(d, "asv_table.csv"); tp = os.path.join(d, "asv_taxonomy.csv")
    if os.path.exists(ap) and os.path.exists(tp):
        ab = {}
        with io.open(ap, encoding="utf-8-sig") as fh:
            rr = csv.reader(fh); next(rr, None)
            for row in rr:
                if row: ab[row[0].strip().strip('"')] = sum(float(x) for x in row[1:] if x not in ("", "NA"))
        tot = g = 0.0
        for row in csv.DictReader(io.open(tp, encoding="utf-8-sig")):
            a = ab.get((row.get("") or "").strip().strip('"'), 0.0); tot += a
            v = (row.get("Genus") or "").strip()
            if v and v.upper() != "NA": g += a
        if tot: m["genus_r"] = 100.0*g/tot
    return m

rows = []
for p in PROJS:
    b = metrics(os.path.join(BEFORE, p)); a = metrics(os.path.join(R, "out_stage2", p))
    ok = (a["n"] >= b["n"]) and ((a["keep"] or 0) >= (b["keep"] or 0) - 1.0) and (a["merge"] is None or a["merge"] >= 40)
    bdir = os.path.join(BEFORE, p)
    if not ok and not os.path.isdir(bdir):
        # 2026-10-10: 参照版本已被上一次运行搬走时，绝不能再把变体搬走——否则项目会从 out_stage2 消失，
        # 而 shutil.move 会抛 FileNotFoundError 中断本脚本。此时保留变体并打印，交人工判断。
        print("!! %s: no reference version in %s -- keep variant in place, needs manual review" % (p, BEFORE))
        verdict = "保留变体（无参照版本，未自动回退）"
    else:
        verdict = "保留变体" if ok else "回退到 Phase 3 版本"
        if not ok:
            os.makedirs(REJ, exist_ok=True)
            dst = os.path.join(REJ, "%s_%s" % (p, datetime.datetime.now().strftime("%H%M%S")))
            shutil.move(os.path.join(R, "out_stage2", p), dst)
            shutil.move(bdir, os.path.join(R, "out_stage2", p))
    rows.append(dict(project=p, before_samples=b["n"], after_samples=a["n"],
                     before_keep=("%.1f" % b["keep"]) if b["keep"] is not None else "",
                     after_keep=("%.1f" % a["keep"]) if a["keep"] is not None else "",
                     before_merge=("%.1f" % b["merge"]) if b["merge"] is not None else "",
                     after_merge=("%.1f" % a["merge"]) if a["merge"] is not None else "",
                     after_genus_r=("%.1f" % a["genus_r"]) if a["genus_r"] is not None else "", verdict=verdict))
    print("%-14s 样本 %d->%d  保留率 %s->%s  合并率 %s->%s  -> %s" % (
        p, b["n"], a["n"], rows[-1]["before_keep"], rows[-1]["after_keep"],
        rows[-1]["before_merge"], rows[-1]["after_merge"], verdict))
with io.open(os.path.join(R, "out50k", "phase5b_decision.csv"), "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("wrote out50k/phase5b_decision.csv")

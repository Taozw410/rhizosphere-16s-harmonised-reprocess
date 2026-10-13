# -*- coding: utf-8 -*-
"""Quarantine projects that the gates have settled as unusable.

Nothing is deleted: every path is MOVED into _recycle/unusable_<date>/<project>/<slot>/.
Selection rule (settled evidence only):
  * gate tier 3 (classify_w1)                                  -> unusable
  * marker_gene_qc verdict EXCLUDE_MARKER                      -> unusable
  * finished run (asv_taxonomy + read_tracking present) with genus rate < 20% (literal "NA" excluded)
  * explicit: placeholder-quality projects recorded in the gate output
Dry run by default; pass --apply to move.

Usage: python quarantine_unusable.py [--apply] [--date 20261009]
"""
import csv, io, os, re, shutil, argparse, datetime

R = r"D:\CodexWork\mic-re"
RE_PRJ = re.compile(r"^PRJ[A-Z]*\d+$")


def rd(p):
    if not os.path.exists(p):
        return []
    try:
        return list(csv.DictReader(io.open(p, encoding="utf-8-sig")))
    except Exception:
        return []


def genus_rate(pdir):
    tax = os.path.join(pdir, "asv_taxonomy.csv")
    if not os.path.exists(tax) or not os.path.exists(os.path.join(pdir, "read_tracking.csv")):
        return None
    n = g = 0
    for row in csv.DictReader(io.open(tax, encoding="utf-8-sig")):
        n += 1
        v = (row.get("Genus") or "").strip()
        if v and v.upper() != "NA":
            g += 1
    return (100.0 * g / n) if n else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--date", default=datetime.datetime.now().strftime("%Y%m%d"))
    a = ap.parse_args()

    reasons = {}

    def mark(p, why):
        if p and RE_PRJ.match(p):
            reasons.setdefault(p, [])
            if why not in reasons[p]:
                reasons[p].append(why)

    # 1) gate tier 3
    for f in ("out_stage2", "out50k"):
        for r in rd(os.path.join(R, f, "W1_usability_tiers.csv")):
            if (r.get("tier") or "").strip() == "3":
                mark(r.get("project"), "闸门 tier3：%s" % (r.get("reason") or "")[:60])
    # 2) marker gene EXCLUDE_MARKER
    for r in rd(os.path.join(R, "out50k", "marker_gene_qc.csv")):
        if (r.get("verdict") or "").strip() == "EXCLUDE_MARKER":
            mark(r.get("project"), "标记基因闸门：非 16S（%s）" % (r.get("its") and ("ITS=%s" % r.get("its")) or r.get("other") or ""))
    # 3) finished runs whose genus rate is below the 20% marker-gate threshold
    for base in ("out_stage2", "out50k"):
        d0 = os.path.join(R, base)
        if not os.path.isdir(d0):
            continue
        for p in sorted(os.listdir(d0)):
            if not RE_PRJ.match(p):
                continue
            g = genus_rate(os.path.join(d0, p))
            if g is not None and g < 20:
                mark(p, "属级注释率 %.1f%%（<20%%，标记基因不符）" % g)

    # 4) probe-screening verdicts, but ONLY for projects that never got a processed verdict.
    # A processed verdict (marker gate / tier) always supersedes a probe-stage verdict: e.g.
    # PRJNA1479906 carries EXCLUDE_ITS in rescreen_strict.csv yet is OK_16S (genus 75.4%) after
    # processing, because the probe had read an ITS file of a mixed project.
    import glob as _glob
    processed = set()
    for r in rd(os.path.join(R, "out50k", "master_inventory.csv")):
        if (r.get("verdict") or "").strip() in ("OK_16S", "CHECK_MIXED"):
            processed.add(r.get("project"))
    for base in ("out_stage2", "out50k"):
        d0 = os.path.join(R, base)
        if os.path.isdir(d0):
            for p in os.listdir(d0):
                if os.path.exists(os.path.join(d0, p, "project_qc.csv")):
                    processed.add(p)
    for f in sorted(_glob.glob(os.path.join(R, "out50k", "*screening*.csv"))) + \
             [os.path.join(R, "out50k", "rescreen_strict.csv")]:
        for r in rd(f):
            v = (r.get("verdict") or "").strip()
            p = r.get("project")
            if v.startswith("EXCLUDE") and p not in processed:
                mark(p, "筛查判定 %s（%s）" % (v, os.path.basename(f)))

    slots = [("raw", os.path.join(R, "data", "raw")),
             ("raw50k", os.path.join(R, "data", "raw50k")),
             ("raw50k_longread", os.path.join(R, "data", "raw50k_longread")),
             ("raw50k_w1done", os.path.join(R, "data", "raw50k_w1done")),
             ("raw50k_excluded", os.path.join(R, "data", "raw50k_excluded")),
             ("raw50k_excluded_its", os.path.join(R, "data", "raw50k_excluded_its")),
             ("out_stage2", os.path.join(R, "out_stage2")),
             ("out50k", os.path.join(R, "out50k"))]

    dest_root = os.path.join(R, "_recycle", "unusable_%s" % a.date)
    if not os.path.abspath(dest_root).startswith(os.path.abspath(os.path.join(R, "_recycle"))):
        raise SystemExit("refuse: destination outside _recycle")

    rows = []
    for p in sorted(reasons):
        moved = []
        size = 0
        for slot, base in slots:
            src = os.path.join(base, p)
            if os.path.isdir(src):
                sz = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(src) for f in fs)
                size += sz
                moved.append((slot, src, sz))
        rows.append((p, "; ".join(reasons[p]), size, moved))

    print("%-14s %-9s %s" % ("project", "size", "reason"))
    for p, why, size, moved in rows:
        print("%-14s %8.1f MB  %s" % (p, size / 1e6, why))
    print("total projects: %d   total size: %.2f GB" % (len(rows), sum(r[2] for r in rows) / 1e9))

    manifest = os.path.join(R, "out50k", "unusable_quarantine_%s.csv" % a.date)
    with io.open(manifest, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["project", "reason", "size_bytes", "moved_paths", "new_path"])
        for p, why, size, moved in rows:
            w.writerow([p, why, size, "; ".join(m[1] for m in moved),
                        os.path.join(dest_root, p)])
    print("manifest:", manifest)

    if not a.apply:
        print("DRY RUN (pass --apply to move). Nothing was changed.")
        return

    log = io.open(os.path.join(R, "out50k", "unusable_quarantine_%s.log" % a.date), "w", encoding="utf-8")
    log.write("quarantine run %s\n" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for p, why, size, moved in rows:
        for slot, src, sz in moved:
            dst = os.path.join(dest_root, p, slot)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.move(src, dst)
            log.write("MOVED %-14s %-18s %-9.1f MB -> %s\n" % (p, slot, sz / 1e6, dst))
    log.close()
    print("moved %d projects into %s" % (len(rows), dest_root))


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Before/after comparison for a truncLen rerun (A/B), written as Markdown.

Usage: python trunclen_ab_report.py <before_root> <after_root> <p1,p2,...> [--out file.md]
Reads, per project and per side: asv_table.csv (rows = ASVs), read_tracking.csv, asv_taxonomy.csv.
Reports: n_asv, ASV length median/mode, total reads, filtered/input, merged/filtered, nonchim/input,
and the genus assignment rate with the literal "NA" excluded.
"""
import sys, os, csv, io, collections, statistics

def genus_rate(path):
    if not os.path.exists(path):
        return None
    n = g = 0
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            n += 1
            v = (row.get("Genus") or "").strip()
            if v and v.upper() != "NA":
                g += 1
    return (100.0 * g / n) if n else None

def describe(d):
    out = dict(n_asv=None, med=None, mode=None, total=None, samples=None,
               inp=None, filt=None, merged=None, nonch=None, pass_pct=None,
               merge_pct=None, retain_pct=None, genus=None)
    t = os.path.join(d, "asv_table.csv")
    if os.path.exists(t):
        lens = []
        tot = 0
        with open(t, encoding="utf-8-sig", newline="") as fh:
            r = csv.reader(fh)
            next(r, None)
            for row in r:
                if not row:
                    continue
                lens.append(len(row[0].strip().strip('"')))
                for x in row[1:]:
                    if x not in ("", "NA"):
                        try:
                            tot += int(float(x))
                        except ValueError:
                            pass
        if lens:
            out["n_asv"] = len(lens)
            out["med"] = statistics.median(lens)
            out["mode"] = collections.Counter(lens).most_common(1)[0][0]
            out["total"] = tot
    rt = os.path.join(d, "read_tracking.csv")
    if os.path.exists(rt):
        with open(rt, encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.DictReader(fh))
        out["samples"] = len(rows)
        out["inp"] = sum(int(float(x["input"])) for x in rows if x.get("input"))
        out["filt"] = sum(int(float(x["filtered"])) for x in rows if x.get("filtered"))
        out["merged"] = sum(int(float(x.get("merged") or x.get("denoised") or 0)) for x in rows)
        out["nonch"] = sum(int(float(x.get("nonchim") or 0)) for x in rows)
        if out["inp"]:
            out["pass_pct"] = 100.0 * out["filt"] / out["inp"]
        if out["filt"]:
            out["merge_pct"] = 100.0 * out["merged"] / out["filt"]
        if out["inp"]:
            out["retain_pct"] = 100.0 * out["nonch"] / out["inp"]
    out["genus"] = genus_rate(os.path.join(d, "asv_taxonomy.csv"))
    return out

def fmt(v, nd=1):
    if v is None:
        return "-"
    if isinstance(v, float):
        return ("%." + str(nd) + "f") % v
    return str(v)

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    outfile = ""
    if "--out" in sys.argv:
        outfile = sys.argv[sys.argv.index("--out") + 1]
    if len(args) < 3:
        sys.exit(__doc__)
    before_root, after_root, projs = args[0], args[1], [p.strip() for p in args[2].split(",") if p.strip()]
    buf = io.StringIO()
    buf.write("| project | side | n_asv | ASV_len_median | ASV_len_mode | total_reads | samples | filtered/input | merged/filtered | nonchim/input | genus_rate |\n")
    buf.write("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |\n")
    for p in projs:
        for side, root in (("before", before_root), ("after", after_root)):
            d = os.path.join(os.path.expanduser(root), p)
            s = describe(d)
            buf.write("| %s | %s | %s | %s | %s | %s | %s | %s%% | %s%% | %s%% | %s |\n" % (
                p, side, fmt(s["n_asv"]), fmt(s["med"]), fmt(s["mode"]), fmt(s["total"]), fmt(s["samples"]),
                fmt(s["pass_pct"]), fmt(s["merge_pct"]), fmt(s["retain_pct"]),
                ("-" if s["genus"] is None else "%.1f%%" % s["genus"])))
    text = buf.getvalue()
    sys.stdout.write(text)
    if outfile:
        with io.open(outfile, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        print("wrote " + outfile)

if __name__ == "__main__":
    main()

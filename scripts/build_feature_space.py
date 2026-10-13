# -*- coding: utf-8 -*-
"""Build the section-17 feature space for the final 45 projects (genus-primary + family fallback).

Rules (pre-registration section 17.1, implemented once in scripts/tax_bins.py and reused here):
  Genus valid (non-empty, not literal NA, not a SILVA placeholder)  -> feature "g__<genus>"
  Family valid                                                      -> feature "f__<family>"
  otherwise                                                         -> feature "unassigned" (single bucket)
  host organelles (mitochondria/chloroplast/plastid/host genus)      -> REMOVED from the totals

Outputs (out_analysis/):
  feature_matrix.csv  rows = samples ("project|sample"), cols = features, values = read counts
  sample_meta.csv     per-sample metadata incl. bacteria/organelle/unassigned reads
  feature_map.csv     feature -> level, total reads, number of projects present
  feature_space_report.md
usage: python build_feature_space.py
"""
import csv, io, os, sys, collections

R = r"D:\CodexWork\mic-re"
sys.path.insert(0, os.path.join(R, "scripts"))
import tax_bins as TB          # single definition of the section-17 label rule

OUT = os.path.join(R, "out_analysis")


def read_table(pdir):
    """-> (samples, {asv: [counts...]}) from asv_table.csv"""
    ap = os.path.join(pdir, "asv_table.csv")
    with io.open(ap, encoding="utf-8-sig") as fh:
        rr = csv.reader(fh)
        hdr = next(rr)
        samples = [x.strip().strip('"') for x in hdr[1:]]
        per = {}
        for row in rr:
            if not row:
                continue
            asv = row[0].strip().strip('"')
            per[asv] = [float(x) if x not in ("", "NA") else 0.0 for x in row[1:len(samples) + 1]]
    return samples, per


def main():
    os.makedirs(OUT, exist_ok=True)
    sel = [r for r in csv.DictReader(io.open(os.path.join(R, "meta", "final45_selection.csv"), encoding="utf-8-sig"))
           if (r.get("rank") or "").strip()]
    print("selected projects: %d" % len(sel))

    rows = []                    # sample rows: dict(meta + {feature: count})
    feat_tot = collections.Counter()
    feat_proj = collections.defaultdict(set)
    crop_fallback = collections.defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])   # bacteria, genus, family, unassigned
    n_flag = []

    for r in sel:
        proj, base, crop = r["project"], r["base"], r["crop"]
        pdir = os.path.join(R, base, proj)
        if not os.path.exists(os.path.join(pdir, "asv_table.csv")):
            n_flag.append("%s missing asv_table" % proj); continue
        samples, per = read_table(pdir)
        tax = TB.load_tax(pdir)
        # label every ASV once
        lab = {}
        for asv in per:
            t = tax.get(asv)
            if t is None:
                lab[asv] = ("unassigned", "unassigned")
            else:
                name, level = TB.tax_label(t)
                lab[asv] = (name if level != "organelle" else "ORGANELLE", level)
        # aggregate per sample
        acc = [collections.Counter() for _ in samples]
        org = [0.0] * len(samples)
        for asv, counts in per.items():
            name, level = lab[asv]
            for i, c in enumerate(counts):
                if not c:
                    continue
                if level == "organelle":
                    org[i] += c
                else:
                    acc[i][name] += c
                    feat_tot[name] += c
                    feat_proj[name].add(proj)
        for i, s in enumerate(samples):
            d = acc[i]
            bact = sum(d.values())
            g = sum(v for k, v in d.items() if k.startswith("g__"))
            f = sum(v for k, v in d.items() if k.startswith("f__"))
            u = d.get("unassigned", 0.0)
            a = crop_fallback[crop]
            a[0] += bact; a[1] += g; a[2] += f; a[3] += u
            row = dict(project=proj, base=base, crop=crop, sample=s,
                       all_reads=bact + org[i], bacteria_reads=bact, organelle_reads=org[i],
                       genus_reads=g, family_reads=f, unassigned_reads=u,
                       n_asv=int(sum(1 for asv in per if per[asv][i] > 0)))
            row.update({k: v for k, v in d.items() if v > 0})
            rows.append(row)

    # keep features present in >=2 projects (singletons cannot inform any cross-study structure)
    keep = {f for f, s in feat_proj.items() if len(s) >= 2}
    feats = sorted(keep, key=lambda f: -feat_tot[f])
    print("samples: %d   features: %d (kept; %d total before >=2-project filter)" % (len(rows), len(feats), len(feat_tot)))

    cols = ["project", "base", "crop", "sample", "all_reads", "bacteria_reads", "organelle_reads",
            "genus_reads", "family_reads", "unassigned_reads", "n_asv"] + feats
    with io.open(os.path.join(OUT, "feature_matrix.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in rows:
            # missing features must be written as an explicit 0 (an empty cell would be read as NA in R)
            w.writerow([(int(r[c]) if isinstance(r.get(c), float)
                         else (r[c] if r.get(c) is not None else 0)) for c in cols])

    with io.open(os.path.join(OUT, "sample_meta.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["project", "base", "crop", "sample", "all_reads", "bacteria_reads",
                                           "organelle_reads", "genus_reads", "family_reads", "unassigned_reads", "n_asv"])
        w.writeheader()
        for r in rows:
            w.writerow({k: (int(r[k]) if k.endswith("reads") or k == "n_asv" else r[k]) for k in w.fieldnames})

    with io.open(os.path.join(OUT, "feature_map.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["feature", "level", "total_reads", "n_projects"])
        for f in feats:
            w.writerow([f, "genus" if f.startswith("g__") else ("unassigned" if f == "unassigned" else "family"),
                        int(feat_tot[f]), len(feat_proj[f])])

    rep = ["# 分析用特征空间（预注册 §17，最终 45 个项目）", "",
           "样本 %d；保留特征 %d（出现在 ≥2 个项目里；过滤前共 %d）。" % (len(rows), len(feats), len(feat_tot)),
           "细胞器 reads 已剔除；SILVA 占位名（uncultured/metagenome…）与只到目/纲/门的 reads 归入 `unassigned` 单桶。", "",
           "| 作物 | 细菌读数(M) | 属级 | 科级回退 | 未定 |", "| --- | --- | --- | --- | --- |"]
    for crop in ("maize", "rice", "wheat"):
        a = crop_fallback[crop]
        if a[0] <= 0:
            continue
        rep.append("| %s | %.2f | %.1f%% | %.1f%% | %.1f%% |" % (crop, a[0] / 1e6, 100 * a[1] / a[0], 100 * a[2] / a[0], 100 * a[3] / a[0]))
    tot = [sum(crop_fallback[c][i] for c in crop_fallback) for i in range(4)]
    rep.append("| 全部 | %.2f | %.1f%% | %.1f%% | %.1f%% |" % (tot[0] / 1e6, 100 * tot[1] / tot[0], 100 * tot[2] / tot[0], 100 * tot[3] / tot[0]))
    rep += ["", "## 特征构成", "",
            "| 层级 | 特征数 | 读数(M) |", "| --- | --- | --- |"]
    for lvl in ("genus", "family", "unassigned"):
        fs = [f for f in feats if (f.startswith("g__") if lvl == "genus" else (f == "unassigned" if lvl == "unassigned" else f.startswith("f__")))]
        rep.append("| %s | %d | %.2f |" % (lvl, len(fs), sum(feat_tot[f] for f in fs) / 1e6))
    if n_flag:
        rep += ["", "## 异常", ""] + ["- " + x for x in n_flag]
    io.open(os.path.join(OUT, "feature_space_report.md"), "w", encoding="utf-8").write("\n".join(rep) + "\n")
    print("wrote out_analysis/{feature_matrix.csv,sample_meta.csv,feature_map.csv,feature_space_report.md}")
    for line in rep[3:14]:
        print(line)


if __name__ == "__main__":
    main()

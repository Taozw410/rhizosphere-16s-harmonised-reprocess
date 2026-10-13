# -*- coding: utf-8 -*-
"""预注册 §17 的执行件：属级为主 + 科级回退的特征标签，以及回退规模报告。

标签规则（先写死，见 分析口径与实验预注册_20261008.md §17.1）：
  宿主机体细胞器（任一层级含 mitochond/chloroplast/plastid、或 Family/Genus 以宿主植物属名开头）
        -> organelle   （从分析总体中剔除，单独报告占比）
  Genus 非空、非字面 "NA"、非占位名   -> g__<属>     （属级）
  Family 非空、非字面 "NA"、非占位名  -> f__<科>     （科级回退 bin）
  其余（只到目/纲/门，或名字是占位名） -> unassigned

占位名（不建 bin）：含 uncultured / unidentified / metagenome / unclassified /
environmental sample / ambiguous taxa。

用法：
  python tax_bins.py [--selection meta/final45_selection.csv] [--out out50k/tax_bin_report.md]

内置一致性自检：按 ASV 聚合出的到门/到科/到属率（以含细胞器的全部 reads 为分母）
必须与 select_final45.tax_reads() 一致（差 <=0.1 个百分点），否则报错退出。
"""
import argparse
import csv
import io
import os
import re
import sys

R = r"D:\CodexWork\mic-re"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

LEVELS = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
ORGANELLE_TOKENS = ("mitochond", "chloroplast", "plastid")
HOST_GENERA = ("Oryza", "Zea", "Triticum", "Arabidopsis", "Nicotiana", "Glycine", "Medicago")
PLACEHOLDER_TOKENS = ("uncultured", "unidentified", "metagenome", "unclassified",
                      "environmental sample", "ambiguous taxa")


def clean(v):
    v = (v or "").strip()
    return "" if (not v or v.upper() == "NA") else v


def is_organelle(tax):
    vals = [clean(tax.get(k)) for k in LEVELS]
    blob = " | ".join(v.lower() for v in vals if v)
    if any(tok in blob for tok in ORGANELLE_TOKENS):
        return True
    for v in vals:
        if v and v.split()[0] in HOST_GENERA:
            return True
    return False


def is_placeholder(name):
    low = name.lower()
    return any(tok in low for tok in PLACEHOLDER_TOKENS)


def tax_label(tax):
    """返回 (特征名, 层级)；层级为 organelle / genus / family / unassigned。"""
    if is_organelle(tax):
        return "", "organelle"
    g = clean(tax.get("Genus"))
    if g and not is_placeholder(g):
        return "g__" + g, "genus"
    f = clean(tax.get("Family"))
    if f and not is_placeholder(f):
        return "f__" + f, "family"
    return "unassigned", "unassigned"


def load_reads(pdir):
    ap = os.path.join(pdir, "asv_table.csv")
    out = {}
    with io.open(ap, encoding="utf-8-sig") as fh:
        rr = csv.reader(fh)
        next(rr, None)
        for row in rr:
            if not row:
                continue
            out[row[0].strip().strip('"')] = sum(float(x) for x in row[1:] if x not in ("", "NA"))
    return out


def load_tax(pdir):
    out = {}
    for row in csv.DictReader(io.open(os.path.join(pdir, "asv_taxonomy.csv"), encoding="utf-8-sig")):
        out[(row.get("") or "").strip().strip('"')] = row
    return out


def project_bins(pdir):
    reads = load_reads(pdir)
    tax = load_tax(pdir)
    allr = org = gen = fam = una = unap = 0.0
    dep = {}
    bins = {}
    for k, n in reads.items():
        t = tax.get(k)
        if t is None:
            continue
        allr += n
        for lvl in ("Phylum", "Family", "Genus"):
            if clean(t.get(lvl)):
                dep[lvl] = dep.get(lvl, 0.0) + n
        name, level = tax_label(t)
        if level == "organelle":
            org += n
        elif level == "genus":
            gen += n
        elif level == "family":
            fam += n
            bins[name] = bins.get(name, 0.0) + n
        else:
            una += n
            f = clean(t.get("Family"))
            if f and is_placeholder(f):
                unap += n
    if not allr:
        return None
    return dict(all_reads=allr, bacteria=allr - org, organelle=org, genus=gen, family=fam,
                unassigned=una, unassigned_placeholder=unap, depth=dep, bins=bins)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", default=os.path.join(R, "meta", "final45_selection.csv"))
    ap.add_argument("--out", default=os.path.join(R, "out50k", "tax_bin_report.md"))
    ap.add_argument("--skip-check", action="store_true")
    a = ap.parse_args()

    rows = [r for r in csv.DictReader(io.open(a.selection, encoding="utf-8-sig"))
            if (r.get("rank") or "").strip()]
    if not rows:
        print("no selected projects in", a.selection)
        return 1

    import select_final45 as S
    keys = ("all_reads", "bacteria", "organelle", "genus", "family", "unassigned", "unassigned_placeholder")
    per_crop = {}
    bins_all = {}
    bad = []
    n_ok = 0
    for r in rows:
        proj = r["project"]
        pdir = os.path.join(R, r["base"], proj)
        if not (os.path.exists(os.path.join(pdir, "asv_table.csv")) and
                os.path.exists(os.path.join(pdir, "asv_taxonomy.csv"))):
            print("skip (incomplete):", proj)
            continue
        b = project_bins(pdir)
        if not b:
            continue
        n_ok += 1
        if not a.skip_check:
            t = S.tax_reads(pdir) or {}
            for lvl, key in (("Phylum", "phylum_r"), ("Family", "family_r"), ("Genus", "genus_r")):
                mine = 100.0 * b["depth"].get(lvl, 0.0) / b["all_reads"]
                theirs = t.get(key)
                if theirs is not None and abs(mine - theirs) > 0.1:
                    bad.append("%s %s: mine=%.2f theirs=%.2f" % (proj, key, mine, theirs))
        c = per_crop.setdefault(r["crop"], dict((k, 0.0) for k in keys))
        for k in keys:
            c[k] += b[k]
        for name, n in b["bins"].items():
            bins_all[name] = bins_all.get(name, 0.0) + n

    if bad:
        print("CONSISTENCY CHECK FAILED (label rule disagrees with select_final45.tax_reads):")
        for x in bad[:20]:
            print("  ", x)
        return 2

    order = [c for c in ("maize", "rice", "wheat") if c in per_crop]
    tot = dict((k, sum(per_crop[c][k] for c in order)) for k in keys)
    rep = ["# 属级为主 + 科级回退：读数归属与回退规模（预注册 §17）", "",
           "入选项目 %d 个（实际读到结果 %d 个）；标签规则见脚本 docstring 与口径文件 §17.1（先写死，不得事后调整）。" % (len(rows), n_ok), "",
           "细胞器 reads 已从分析总体中剔除；下表百分比的分母是**细菌 reads**（全部 reads 减细胞器）。", "",
           "| 作物 | 细菌读数(M) | 属级 | 科级回退 | 未定 | 细胞器占全部 reads |",
           "| --- | --- | --- | --- | --- | --- |"]
    for c in order + (["全部"] if len(order) > 1 else []):
        d = dict((k, tot[k]) for k in keys) if c == "全部" else per_crop[c]
        rep.append("| %s | %.2f | %.1f%% | %.1f%% | %.1f%% | %.1f%% |" % (
            c, d["bacteria"] / 1e6, 100 * d["genus"] / d["bacteria"],
            100 * d["family"] / d["bacteria"], 100 * d["unassigned"] / d["bacteria"],
            100 * d["organelle"] / d["all_reads"]))
    rep += ["", "## 科级回退规模（按作物；这是结论强度的边界，必须进方法）", "",
            "| 作物 | 占该作物细菌读数 | 占全部细菌读数 |",
            "| --- | --- | --- |"]
    for c in order:
        d = per_crop[c]
        rep.append("| %s | %.1f%% | %.1f%% |" % (c, 100 * d["family"] / d["bacteria"],
                                                 100 * d["family"] / tot["bacteria"]))
    rep += ["", "## 未定类别的构成（透明化）", "",
            "| 作物 | 未定合计 | 其中占位名校名 | 其中只到目/纲/门 |",
            "| --- | --- | --- | --- |"]
    for c in order:
        d = per_crop[c]
        rep.append("| %s | %.1f%% | %.1f%% | %.1f%% |" % (
            c, 100 * d["unassigned"] / d["bacteria"],
            100 * d["unassigned_placeholder"] / d["bacteria"],
            100 * (d["unassigned"] - d["unassigned_placeholder"]) / d["bacteria"]))
    rep += ["", "## 回退读数最多的科（前 15，供正文/补充材料说明）", "",
            "| 科级 bin | 读数 | 占全部细菌读数 |", "| --- | --- | --- |"]
    for name, n in sorted(bins_all.items(), key=lambda kv: -kv[1])[:15]:
        rep.append("| %s | %.2f M | %.2f%% |" % (name, n / 1e6, 100 * n / tot["bacteria"]))
    io.open(a.out, "w", encoding="utf-8").write("\n".join(rep) + "\n")
    print("wrote", a.out, " (projects:", n_ok, ")")
    return 0


if __name__ == "__main__":
    sys.exit(main())

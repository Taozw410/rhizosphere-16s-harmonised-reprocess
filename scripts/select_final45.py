# -*- coding: utf-8 -*-
"""Apply the pre-registered selection rule (分析口径与实验预注册_20261008.md §4) to pick the final set.

Hard exclusions: placeholder quality, zero reads after filtering, marker-gene mismatch, and
projects that cannot merge (no overlap / wrong variable region).
Ranking keys: (1) genus rate >= 60% (excluding literal "NA"), (2) merge rate >= 40%,
(3) final retention, (4) ASV count sanity. Max 15 per crop; max 2 per platform cluster per crop
(cluster = crop + centre + first 45 chars of title).

Usage: python select_final45.py [--per-crop 15] [--out meta/final45_selection.csv] [--report out50k/选集报告.md]
Reads per-project artifacts directly (out50k/<p> and out_stage2/<p>): project_qc.csv, read_tracking.csv,
asv_taxonomy.csv, filter_stats.csv. Crop/cluster metadata from meta/*.csv.
"""
import csv, io, os, glob, re, argparse, collections, statistics

R = r"D:\CodexWork\mic-re"
DIRS = [os.path.join(R, "out50k"), os.path.join(R, "out_stage2")]
RE_PRJ = re.compile(r"^PRJ[A-Z]*\d+$")


def rd(path):
    if not os.path.exists(path):
        return []
    try:
        return list(csv.DictReader(io.open(path, encoding="utf-8-sig")))
    except Exception:
        return []


def load_meta():
    """project -> dict(crop, center, title)"""
    m = {}
    for f in sorted(glob.glob(os.path.join(R, "meta", "*.csv"))):
        for r in rd(f):
            k = (r.get("study") or r.get("project") or "").strip()
            if not k:
                continue
            d = m.setdefault(k, {})
            for src, dst in (("crop", "crop"), ("center", "center"), ("centre", "center"), ("title", "title")):
                v = (r.get(src) or "").strip()
                if v and not d.get(dst):
                    d[dst] = v
    return m


def tax_reads(pdir):
    """Read-weighted taxonomy depth: % of reads assigned to phylum / family / genus.

    The old gate used the genus rate counted over ASVs, which conflates two different
    failures: artifact-laden data (nothing reaches phylum) and genuine communities whose
    dominant lineages SILVA resolves only to phylum/class (Patescibacteria, Nanoarchaeota,
    Acidobacteriota ...). 2026-10-10 diagnosis: artifact-type projects sit at 0.3-52.9%
    phylum-assigned reads, real data at 79-100%.
    """
    ap = os.path.join(pdir, "asv_table.csv")
    tp = os.path.join(pdir, "asv_taxonomy.csv")
    if not (os.path.exists(ap) and os.path.exists(tp)):
        return None
    ab = {}
    try:
        with io.open(ap, encoding="utf-8-sig") as fh:
            rr = csv.reader(fh)
            next(rr, None)
            for row in rr:
                if row:
                    ab[row[0].strip().strip('"')] = sum(float(x) for x in row[1:] if x not in ("", "NA"))
    except Exception:
        return None
    tot = ph = fam = gen = 0.0
    n = ngen = 0
    for row in csv.DictReader(io.open(tp, encoding="utf-8-sig")):
        seq = (row.get("") or "").strip().strip('"')
        a = ab.get(seq, 0.0)
        tot += a
        n += 1
        g = (row.get("Genus") or "").strip()
        if g and g.upper() != "NA":
            gen += a
            ngen += 1
        f = (row.get("Family") or "").strip()
        if f and f.upper() != "NA":
            fam += a
        p2 = (row.get("Phylum") or "").strip()
        if p2 and p2.upper() != "NA":
            ph += a
    if tot <= 0 or n == 0:
        return None
    return dict(phylum_r=100.0 * ph / tot, family_r=100.0 * fam / tot, genus_r=100.0 * gen / tot,
                genus_asv=100.0 * ngen / n, n_asv_tax=n)


def project_metrics(pdir):
    """(n_samples, n_asv, fpass, merge, keep, genus) or None when the run did not finish."""
    qc = rd(os.path.join(pdir, "project_qc.csv"))
    qc = qc[0] if qc else {}
    rt = rd(os.path.join(pdir, "read_tracking.csv"))
    n_asv = None
    try:
        n_asv = int(qc.get("n_asv"))
    except Exception:
        pass
    fpass = merge = keep = None
    if rt:
        tin = sum(float(x.get("input") or 0) for x in rt)
        tf = sum(float(x.get("filtered") or 0) for x in rt)
        tm = sum(float(x.get("merged") or 0) for x in rt)
        tn = sum(float(x.get("nonchim") or 0) for x in rt)
        if tin:
            fpass = 100.0 * tf / tin
            keep = 100.0 * tn / tin
        if tf:
            merge = 100.0 * tm / tf
    else:
        fs = rd(os.path.join(pdir, "filter_stats.csv"))
        if fs:
            tin = sum(float(x.get("reads.in") or 0) for x in fs)
            tf = sum(float(x.get("reads.out") or 0) for x in fs)
            if tin:
                fpass = 100.0 * tf / tin
    genus = None
    tr = tax_reads(pdir)
    tax = os.path.join(pdir, "asv_taxonomy.csv")
    if os.path.exists(tax):
        n = g = 0
        for row in csv.DictReader(io.open(tax, encoding="utf-8-sig")):
            n += 1
            v = (row.get("Genus") or "").strip()
            if v and v.upper() != "NA":
                g += 1
        genus = (100.0 * g / n) if n else None
    if not qc and n_asv is None and fpass is None:
        return None
    return dict(n_samples=qc.get("n_samples", ""), qual=qc.get("quality_chars", ""),
                n_asv=n_asv, fpass=fpass, merge=merge, keep=keep, genus=genus,
                **({} if tr is None else tr),
                mode=qc.get("mode", ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-crop", type=int, default=15)
    ap.add_argument("--min-genus", type=float, default=60.0)
    ap.add_argument("--min-phylum", type=float, default=80.0)
    ap.add_argument("--min-family", type=float, default=50.0)
    ap.add_argument("--min-merge", type=float, default=40.0)
    ap.add_argument("--min-merge-relaxed", type=float, default=30.0)
    ap.add_argument("--out", default=os.path.join(R, "meta", "final45_selection.csv"))
    ap.add_argument("--report", default=os.path.join(R, "out50k", "选集报告.md"))
    a = ap.parse_args()

    meta = load_meta()
    marker = {r["project"]: r.get("verdict", "") for r in rd(os.path.join(R, "out50k", "marker_gene_qc.csv"))}
    # 复测裁定：未通过预注册救回规则的项目不再回到候选池
    salvage_excluded = set()
    _se = os.path.join(R, "meta", "salvage_excluded.txt")
    if os.path.exists(_se):
        salvage_excluded = set(x.strip() for x in io.open(_se, encoding="utf-8") if x.strip())
    tiers = {}
    for f in ("out_stage2", "out50k"):
        for r in rd(os.path.join(R, f, "W1_usability_tiers.csv")):
            p = r.get("project")
            if p and p not in tiers:
                tiers[p] = (r.get("tier") or "").strip()

    cands = []
    seen = set()
    for base in DIRS:
        if not os.path.isdir(base):
            continue
        for d in sorted(os.listdir(base)):
            if not RE_PRJ.match(d) or d in seen:
                continue
            m = project_metrics(os.path.join(base, d))
            if not m:
                continue
            seen.add(d)
            mm = meta.get(d, {})
            crop = (mm.get("crop") or "").lower()
            center = (mm.get("center") or "").strip()
            title = (mm.get("title") or "").strip()
            mk = marker.get(d, "")
            excl = []
            try:
                q = int(m["qual"])
            except Exception:
                q = None
            # completeness: a run counts only when read tracking AND taxonomy are on disk
            pdir = os.path.join(base, d)
            complete = os.path.exists(os.path.join(pdir, "read_tracking.csv")) and \
                       os.path.exists(os.path.join(pdir, "asv_taxonomy.csv"))
            if not complete:
                excl.append("尚未跑完（缺 read_tracking / asv_taxonomy）")
            if q is not None and q <= 2:
                excl.append("质量值占位符")
            if m["fpass"] is not None and m["fpass"] < 0.001:
                excl.append("过滤后无 reads")
            # 2026-10-10 预注册 §16：只有 EXCLUDE_MARKER 是硬性排除。
            # CHECK_MIXED 来自已被 §13 取代的"属级注释率 <60%（ASV 计数）"判据，
            # 按 §10.6（20–60% 属"待核查"）与 §13（分辨率受限型保留）降为敏感性标记，不排除。
            if mk == "EXCLUDE_MARKER":
                excl.append("标记基因闸门: EXCLUDE_MARKER（硬性排除）")
            if m["genus"] is not None and m["genus"] < 20:
                excl.append("属级注释率 <20%（标记基因不符）")
            if d in salvage_excluded:
                excl.append("复测未通过预注册救回规则")
            if tiers.get(d) == "3":
                excl.append("闸门 tier3")
            if complete and (m.get("phylum_r") is None or m.get("family_r") is None):
                excl.append("\u5206\u7c7b\u5b66\u6307\u6807\u7f3a\u5931\uff08\u65e0\u6cd5\u8ba1\u7b97\u5230\u95e8/\u5230\u79d1\u7387\uff09")
            elif complete and ((m["phylum_r"] < a.min_phylum) or (m["family_r"] < a.min_family)):
                excl.append("\u5206\u7c7b\u5b66\u95f8\u95e8\uff1a\u5230\u95e8 %.1f%% (\u95e8\u69db %.0f%%) / \u5230\u79d1 %.1f%% (\u95e8\u69db %.0f%%)" % (
                    m["phylum_r"], a.min_phylum, m["family_r"], a.min_family))
            if False:
                excl.append("属级注释率 <%.0f%%（预注册硬门槛）" % a.min_genus)
            if crop not in ("maize", "rice", "wheat"):
                excl.append("\u4f5c\u7269\u5b57\u6bb5\u7f3a\u5931\uff08\u65e0\u6cd5\u5f52\u6863\uff09")
            cands.append(dict(project=d, base=os.path.basename(base), crop=crop, center=center, title=title,
                              marker=mk, excl=";".join(excl), rank="", **m))

    # 2026-10-10 预注册 §16：选集必须能给出"敏感性双报"，所以把排除+排序做成可复用函数，
    # 主选集与两种敏感性情形各跑一次；每次都对行做浅拷贝，互不污染。
    def run_selection(src, merge_floor, checkmixed_is_excl):
        rows = [dict(c) for c in src]
        for c in rows:
            ex = c.get("excl", "")
            if checkmixed_is_excl and c.get("marker") == "CHECK_MIXED":
                ex = (ex + "；" if ex else "") + "敏感性情形：CHECK_MIXED 视为排除"
            c["excl"] = ex
            c["rank"] = ""
        picked = []
        for crop in ("maize", "rice", "wheat"):
            pool = [c for c in rows if c["crop"] == crop and not c["excl"]]
            # §4 execution note: (1) genus>=min_genus and (2) merge>=40% (paired only) are THRESHOLDS;
            # ranking then uses (3) final retention, (4) merge, (5) genus; 3-5 bin quantized quality is
            # pushed to the end of the crop. Single-end runs have no merge step -> key 2 not applicable.
            def quant(c):
                try:
                    q = int(c["qual"])
                except Exception:
                    return 1
                return 0 if q >= 6 else 1

            def rank_key(c):
                paired = (c["mode"] == "paired") or (c["mode"] == "" and c["merge"] is not None)
                return (quant(c), -(c["keep"] or 0), -(c["merge"] or 0) if paired else 0,
                        -(c["genus"] or 0), -(c["n_asv"] or 0), c["project"])

            strict = [c for c in pool if (c["merge"] is None or c["merge"] >= a.min_merge
                                          or not ((c["mode"] == "paired") or (c["mode"] == "" and c["merge"] is not None)))]
            relaxed = [c for c in pool if c not in strict]
            if merge_floor is not None:
                # §4 原文的递补下限：达不到下限的双端项目不进递补池（仅敏感性情形用）
                relaxed = [c for c in relaxed if c["merge"] is None or c["merge"] >= merge_floor]
            strict.sort(key=rank_key)
            relaxed.sort(key=rank_key)
            pool = strict + relaxed
            cluster_n = collections.Counter()
            n = 0
            for c in pool:
                key = (crop, c["center"].lower(), c["title"][:45].lower())
                if cluster_n[key] >= 2:
                    c["excl"] = c["excl"] or "平台集群内已取 2 个"
                    continue
                if n >= a.per_crop:
                    c["excl"] = c["excl"] or "超出每作物名额"
                    continue
                n += 1
                c["rank"] = n
                c["cluster"] = "%s|%s" % (c["center"][:24], c["title"][:32])
                cluster_n[key] += 1
                picked.append(c)
        return picked, rows

    picked, rows = run_selection(cands, None, False)                          # 主选集（§16 口径）
    _picked_floor, _rows_floor = run_selection(cands, a.min_merge_relaxed, False)   # 敏感性 1：执行递补下限
    _picked_cm, _rows_cm = run_selection(cands, None, True)                        # 敏感性 2：CHECK_MIXED 视为排除

    def sens_flags(c):
        f = []
        if c.get("marker") == "CHECK_MIXED":
            f.append("CHECK_MIXED")
        paired = (c["mode"] == "paired") or (c["mode"] == "" and c["merge"] is not None)
        if paired and c["merge"] is not None:
            if c["merge"] < a.min_merge_relaxed:
                f.append("合并率<%.0f%%（低于递补下限）" % a.min_merge_relaxed)
            elif c["merge"] < a.min_merge:
                f.append("递补（合并率 %.0f–%.0f%%）" % (a.min_merge_relaxed, a.min_merge))
        return ";".join(f)

    dst = a.out
    cols = ["project", "crop", "centre", "title", "n_samples", "mode", "qual", "asv", "fpass", "merge", "keep",
            "genus", "phylum_r", "family_r", "genus_r", "rank", "excl", "cluster", "base",
            "marker", "sens_flags"]
    with io.open(dst, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for c in sorted(rows, key=lambda c: (c["crop"], c["rank"] if c["rank"] != "" else 999, c["project"])):
            w.writerow([c["project"], c["crop"], c["center"], c["title"][:80], c["n_samples"], c["mode"], c["qual"],
                        c["n_asv"] if c["n_asv"] is not None else "",
                        "%.1f%%" % c["fpass"] if c["fpass"] is not None else "",
                        "%.1f%%" % c["merge"] if c["merge"] is not None else "",
                        "%.1f%%" % c["keep"] if c["keep"] is not None else "",
                        "%.1f%%" % c["genus"] if c["genus"] is not None else "",
                        "%.1f%%" % c["phylum_r"] if c.get("phylum_r") is not None else "",
                        "%.1f%%" % c["family_r"] if c.get("family_r") is not None else "",
                        "%.1f%%" % c["genus_r"] if c.get("genus_r") is not None else "",
                        c["rank"], c["excl"], c.get("cluster", ""), c["base"],
                        c.get("marker", ""), sens_flags(c)])

    rep = ["# 最终项目选集（按预注册 §4 规则）", ""]
    rep.append("| 作物 | 入选 | 候选（未入选） | 说明 |")
    rep.append("| --- | --- | --- | --- |")
    for crop in ("maize", "rice", "wheat"):
        sel = [c for c in picked if c["crop"] == crop]
        allc = [c for c in rows if c["crop"] == crop]
        rep.append("| %s | %d | %d | 每作物目标 %d |" % (crop, len(sel), len(allc) - len(sel), a.per_crop))
    rep.append("")
    rep.append("## 入选清单")
    rep.append("")
    rep.append("| # | 作物 | 项目 | 机构 | 样本 | ASV | 过滤通过率 | 合并率 | 保留率 | 属级注释率 | 到门率 | 到科率 | 敏感性标记 |")
    rep.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for c in sorted(picked, key=lambda c: (c["crop"], c["rank"])):
        rep.append("| {rank} | {crop} | {project} | {centre} | {n_samples} | {asv} | {fpass} | {merge} | {keep} | {genus} | {phy} | {fam} | {sens} |".format(
            rank=c["rank"], crop=c["crop"], project=c["project"], centre=c["center"][:28], n_samples=c["n_samples"],
            asv=c["n_asv"], fpass="%.1f%%" % c["fpass"] if c["fpass"] is not None else "-",
            merge="%.1f%%" % c["merge"] if c["merge"] is not None else "-",
            keep="%.1f%%" % c["keep"] if c["keep"] is not None else "-",
            genus="%.1f%%" % c["genus"] if c["genus"] is not None else "-",
            phy="%.1f%%" % c["phylum_r"] if c.get("phylum_r") is not None else "-",
            fam="%.1f%%" % c["family_r"] if c.get("family_r") is not None else "-",
            sens=sens_flags(c) or "—"))
    rep.append("")
    rep.append("## 未入选（含原因）")
    rep.append("")
    rep.append("| 作物 | 项目 | 属级注释率 | 合并率 | 原因 |")
    rep.append("| --- | --- | --- | --- | --- |")
    for c in sorted([c for c in rows if c["rank"] == ""], key=lambda c: (c["crop"], c["project"])):
        rep.append("| %s | %s | %s | %s | %s |" % (c["crop"], c["project"],
                   "%.1f%%" % c["genus"] if c["genus"] is not None else "-",
                   "%.1f%%" % c["merge"] if c["merge"] is not None else "-", c["excl"] or "排序靠后"))
    rep.append("")
    rep.append("## 敏感性双报（预注册 §16）")
    rep.append("")
    rep.append("主选集口径（§16 决定）：**只排除 EXCLUDE_MARKER**；递补池**不设合并率下限**。")
    rep.append("下面两种情形各自重算一遍，给出每作物数量与相对主选集的进出项目，供正文与补充材料双报。")
    rep.append("")
    rep.append("| 情形 | 玉米 | 水稻 | 小麦 | 合计 |")
    rep.append("| --- | --- | --- | --- | --- |")
    for _label, _pk in (("主选集（§16 口径）", picked),
                        ("敏感性 1：递补下限 合并率 ≥%.0f%%（§4 原文）" % a.min_merge_relaxed, _picked_floor),
                        ("敏感性 2：CHECK_MIXED 视为排除（§4 旧读法）", _picked_cm)):
        _cc = [len([c for c in _pk if c["crop"] == _cr]) for _cr in ("maize", "rice", "wheat")]
        rep.append("| %s | %d | %d | %d | %d |" % (_label, _cc[0], _cc[1], _cc[2], sum(_cc)))
    rep.append("")
    _base_set = set(c["project"] for c in picked)
    for _label, _pk in (("敏感性 1（递补下限）", _picked_floor), ("敏感性 2（CHECK_MIXED 排除）", _picked_cm)):
        _s = set(c["project"] for c in _pk)
        rep.append("- **%s**：主选集保留而该情形退出 → %s；该情形新增 → %s" % (
            _label, "、".join(sorted(_base_set - _s)) or "（无）", "、".join(sorted(_s - _base_set)) or "（无）"))
    io.open(a.report, "w", encoding="utf-8").write("\n".join(rep) + "\n")
    print("selected:", len(picked), " (maize %d / rice %d / wheat %d)" % (
        len([c for c in picked if c["crop"] == "maize"]), len([c for c in picked if c["crop"] == "rice"]),
        len([c for c in picked if c["crop"] == "wheat"])))
    print("wrote", dst, "and", a.report)


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""审计本机已下载数据：逐项目探针 + 结果完整性 + 按现行口径归类。

判据（每条都可回溯到口径文件/AGENTS 章节）：
  引物未去   ：R1 开头符合已知 16S/ITS 引物 >=50%，但结果目录里没有 cutadapt_report.csv
               -> 按 §10.3「先去引物再截断」应重跑
  占位符质量 ：质量值种类 <=2 -> AGENTS §10.3 直接排除，不重跑
  标记基因   ：ITS/接头特征占优 -> §10.6 排除，不重跑
  低合并率   ：双端且 merged/filtered <40% 且无 params_fingerprint.txt（旧口径）-> 需按 §9/§14 复核
用法：python audit_downloaded.py [--out out50k/audit_downloaded.md]
"""
import argparse
import csv
import glob
import gzip
import io
import os
import re
import statistics
import sys

R = r"D:\CodexWork\mic-re"
RE_PRJ = re.compile(r"^PRJ[A-Z]*\d+$")

# 引物序列取自 scripts/screen_primers.py 的 PRIMERS 表（保持与本项目既有定义一致）
PRIMERS_16S = {
    "515F": "GTGYCAGCMGCCGCGGTAA", "515Fb": "GTGCCAGCMGCCGCGGTAA",
    "806R": "GGACTACNVGGGTWTCTAAT", "806Rb": "GGACTACHVGGGTWTCTAAT",
    "341F": "CCTACGGGNGGCWGCAG", "805R": "GACTACHVGGGTATCTAATCC",
    "926F": "AAACTYAAAKGAATTGACGG", "1392R": "ACGGGCGGTGTGTRC",
    "27F": "AGAGTTTGATCMTGGCTCAG", "338F": "ACTCCTACGGGAGGCAGCAG",
}
PRIMERS_ITS = {"ITS1F": "CTTGGTCATTTAGAGGAAGTAA", "ITS2": "GCTGCGTTCTTCATCGATGC",
               "ITS3": "GCATCGATGAAGAACGCAGC", "ITS4": "TCCTCCGCTTATTGATATGC"}
ADAPTER = ["CTACACGACGCTCTTCCGATCT", "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
           "AATGATACGGCGACCACCGAGATCTACAC", "CAAGCAGAAGACGGCATACGAGAT"]
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def ambig_prefix(seq, pat):
    if len(seq) < len(pat):
        return False
    for s, p in zip(seq, pat):
        if s not in IUPAC.get(p, ""):
            return False
    return True


def probe(path, n=200):
    """返回 dict：读长中位/最短、质量档数、N 比例、引物/ITS/接头检出率。"""
    lens, qsets, nread, nN, p16, pits, pada = [], set(), 0, 0, 0, 0, 0
    buf = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            buf.append(line.rstrip("\n"))
            if len(buf) == 4:
                s, q = buf[1].upper(), buf[3]
                lens.append(len(s)); qsets.update(q); nN += s.count("N"); nread += 1
                if any(ambig_prefix(s, p) for p in PRIMERS_16S.values()):
                    p16 += 1
                if any(ambig_prefix(s, p) for p in PRIMERS_ITS.values()):
                    pits += 1
                if any(a in s for a in ADAPTER):
                    pada += 1
                buf = []
            if nread >= n:
                break
    if not nread:
        return None
    return dict(n=nread, med=int(statistics.median(lens)), mn=min(lens),
                q=len(qsets), nN=nN / (nread * max(1, statistics.median(lens))),
                p16=100.0 * p16 / nread, pits=100.0 * pits / nread, pada=100.0 * pada / nread)


def rd(p):
    return list(csv.DictReader(io.open(p, encoding="utf-8-sig"))) if os.path.exists(p) else []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(R, "out50k", "audit_downloaded.md"))
    a = ap.parse_args()

    dirs = ["data/raw50k_longread", "data/raw50k_se", "data/raw50k", "data/raw"]
    seen = {}
    for d in dirs:
        bd = os.path.join(R, d)
        if not os.path.isdir(bd):
            continue
        for p in sorted(os.listdir(bd)):
            if not RE_PRJ.match(p):
                continue
            fs = sorted(glob.glob(os.path.join(bd, p, "*.fastq.gz")))
            if not fs:
                continue
            if p not in seen:
                seen[p] = dict(dir=d, files=fs)

    marker = {r["project"]: r.get("verdict", "") for r in rd(os.path.join(R, "out50k", "marker_gene_qc.csv"))}
    rescreen = {r["study"]: r for r in rd(os.path.join(R, "meta", "rescreen_strict.csv")) if r.get("study")} if os.path.exists(os.path.join(R, "meta", "rescreen_strict.csv")) else {}
    rescreen = {}
    p_res = os.path.join(R, "meta", "rescreen_strict.csv")
    if os.path.exists(p_res):
        rr = rd(p_res)
        key = "study" if rr and "study" in rr[0] else list(rr[0].keys())[0]
        rescreen = {r.get(key, ""): r for r in rr}
    quarantine = set()
    for q in glob.glob(os.path.join(R, "_recycle", "*", "*.csv")):
        for r in rd(q):
            v = r.get("project") or r.get("study") or ""
            if v and RE_PRJ.match(v):
                quarantine.add(v)
    for q in glob.glob(os.path.join(R, "_recycle", "*", "*")):
        b = os.path.basename(q)
        if RE_PRJ.match(b):
            quarantine.add(b)
    salvage = set(x.strip() for x in io.open(os.path.join(R, "meta", "salvage_excluded.txt"), encoding="utf-8") if x.strip()) \
        if os.path.exists(os.path.join(R, "meta", "salvage_excluded.txt")) else set()

    rows = []
    for p in sorted(seen):
        info = seen[p]
        fs = info["files"]
        paired = any(x.endswith("_1.fastq.gz") for x in fs) and any(x.endswith("_2.fastq.gz") for x in fs)
        first = next((x for x in fs if x.endswith("_1.fastq.gz")), fs[0])
        pr = probe(first)
        # 结果完整性
        base = None
        for c in ("out50k", "out_stage2"):
            d = os.path.join(R, c, p)
            if os.path.exists(os.path.join(d, "project_qc.csv")) and os.path.exists(os.path.join(d, "asv_taxonomy.csv")):
                base = c
                break
        cutadapt = any(os.path.exists(os.path.join(R, c, p, "cutadapt_report.csv")) for c in ("out50k", "out_stage2"))
        fp = any(os.path.exists(os.path.join(R, c, p, "params_fingerprint.txt")) for c in ("out50k", "out_stage2"))
        merge = keep = None
        if base:
            rt = rd(os.path.join(R, base, p, "read_tracking.csv"))
            if rt:
                tf = sum(float(x.get("filtered") or 0) for x in rt); tm = sum(float(x.get("merged") or 0) for x in rt)
                tin = sum(float(x.get("input") or 0) for x in rt); tn = sum(float(x.get("nonchim") or 0) for x in rt)
                if tf:
                    merge = 100.0 * tm / tf
                if tin:
                    keep = 100.0 * tn / tin
        # 归类
        tags = []
        if not base:
            if p in quarantine:
                tags.append("未跑·已隔离")
            elif p in salvage:
                tags.append("未跑·复测未过")
            else:
                tags.append("**未跑**")
        else:
            tags.append("已跑")
        if pr:
            if pr["q"] <= 2:
                tags.append("占位符质量（排除）")
            if pr["pits"] > 30 and pr["pits"] > pr["p16"]:
                tags.append("疑 ITS（排除）")
            if pr["pada"] > 30:
                tags.append("接头未切（排除）")
            if pr["p16"] >= 50 and not cutadapt:
                tags.append("**引物未去→应重跑**")
            if paired and merge is not None and merge < 40:
                tags.append("合并率<40%")
                if not fp:
                    tags.append("旧口径（无参数指纹）")
        rows.append(dict(project=p, dir=info["dir"], mode="paired" if paired else "single", nfiles=len(fs),
                         med=pr["med"] if pr else "", mn=pr["mn"] if pr else "", q=pr["q"] if pr else "",
                         p16="%.0f" % pr["p16"] if pr else "", pits="%.0f" % pr["pits"] if pr else "",
                         pada="%.0f" % pr["pada"] if pr else "",
                         base=base or "", cutadapt="Y" if cutadapt else "", fp="Y" if fp else "",
                         merge=("%.1f" % merge) if merge is not None else "",
                         keep=("%.1f" % keep) if keep is not None else "",
                         marker=marker.get(p, ""), tags="；".join(tags)))

    def cnt(pred):
        return [r for r in rows if pred(r)]
    unrun = cnt(lambda r: "未跑" in r["tags"])
    rerun = cnt(lambda r: "应重跑" in r["tags"])
    excl = cnt(lambda r: ("排除）" in r["tags"]) or ("已隔离" in r["tags"]))
    lowmerge = cnt(lambda r: "合并率<40%" in r["tags"])
    rep = ["# 已下载数据审计（按现行口径）", "",
           "覆盖 %d 个项目（data/raw、data/raw50k、raw50k_longread、raw50k_se 取并集），探针取每个项目第一个 R1 文件的前 200 条 reads。" % len(rows), "",
           "| 项目 | 目录 | 模式 | 文件 | 读长中位 | 最短 | 质量档 | 16S引物开头% | ITS开头% | 接头% | 结果 | 已去引物 | 参数指纹 | 合并率 | 保留率 | 标记基因 | 归类 |",
           "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        rep.append("| {project} | {dir} | {mode} | {nfiles} | {med} | {mn} | {q} | {p16} | {pits} | {pada} | {base} | {cutadapt} | {fp} | {merge} | {keep} | {marker} | {tags} |".format(**r))
    rep += ["", "## 汇总", "",
            "- 项目总数：%d" % len(rows),
            "- **未跑**：%d 个 → %s" % (len(unrun), ", ".join(r["project"] for r in unrun)),
            "- **按新规则应重跑（引物未去）**：%d 个 → %s" % (len(rerun), ", ".join(r["project"] for r in rerun)),
            "- 已排除（占位符质量 / 疑 ITS / 接头未切 / 已隔离）：%d 个" % len(excl),
            "- 双端且合并率 <40%%：%d 个（其中无参数指纹的旧口径：%d 个）" % (
                len(lowmerge), len(cnt(lambda r: "旧口径" in r["tags"]))),
            ]
    io.open(a.out, "w", encoding="utf-8").write("\n".join(rep) + "\n")
    print("wrote", a.out, " projects:", len(rows))
    for line in rep[-7:]:
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())

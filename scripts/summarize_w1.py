import re
# 汇总：读取所有产出，生成统一 QC 报告
# 用法: summarize_w1.py [输出目录]   默认 out50k
import csv, os, gzip, glob, io, sys

OUT = sys.argv[1] if len(sys.argv) > 1 else r"D:\CodexWork\mic-re\out50k"
META = r"D:\CodexWork\mic-re\meta"
EXC = r"D:\CodexWork\mic-re\data\raw50k_excluded"

def rd(p):
    if not os.path.exists(p): return None
    with io.open(p, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

proj_meta = {}
for r in rd(os.path.join(META, "lightweight_60_final.csv")) or []:
    proj_meta[r["study"]] = r
for r in rd(os.path.join(META, "w1_replacement_candidates.csv")) or []:
    proj_meta.setdefault(r["study"], r)

primers = {r["project"]: r for r in (rd(os.path.join(META, "primer_table.csv")) or [])}

rows = []
for d in sorted(os.listdir(OUT)):
    full = os.path.join(OUT, d)
    if not os.path.isdir(full) or not re.match(r"^PRJ[A-Z]*\d+$", d): continue
    qc = rd(os.path.join(full, "project_qc.csv"))
    if not qc: continue
    qc = qc[0]
    tax = rd(os.path.join(full, "asv_taxonomy.csv"))
    genus_pct = ""
    if tax:
        n = len(tax)
        g = sum(1 for t in tax if (t.get("Genus") or "").strip() and (t.get("Genus") or "").strip().upper() != "NA")
        genus_pct = f"{100*g/n:.1f}%" if n else ""
    trim_pct = ""
    cr = rd(os.path.join(full, "cutadapt_report.csv"))
    if cr:
        vals = [float(x["pct"]) for x in cr if x.get("pct") not in (None, "", "NA")]
        if vals: trim_pct = f"{sum(vals)/len(vals):.1f}%"
    rt = rd(os.path.join(full, "read_tracking.csv")) or []
    keep = ""
    merge = ""
    if rt:
        try:
            tin = sum(int(float(x["input"])) for x in rt)
            tno = sum(int(float(x["nonchim"])) for x in rt)
            keep = f"{100*tno/tin:.1f}%" if tin else ""
            if "merged" in rt[0]:
                tm = sum(int(float(x["merged"])) for x in rt)
                merge = f"{100*tm/tin:.1f}%" if tin else ""
        except Exception: pass
    # 读长（取第一个 raw50k 文件）
    rlen = ""
    rd_dir = os.path.join(r"D:\CodexWork\mic-re\data\raw50k", d)
    try:
        fs = sorted(glob.glob(os.path.join(rd_dir, "*.fastq.gz")))
        if fs:
            with gzip.open(fs[0], "rt", errors="replace") as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 1:
                        rlen = str(len(line.strip())); break
    except Exception: pass
    m = proj_meta.get(d, {})
    rows.append(dict(
        project=d,
        centre=(m.get("center", "") or "")[:32],
        year=(m.get("date", "") or "")[:4],
        samples=qc.get("n_samples", ""),
        mode=qc.get("mode", ""),
        qual=qc.get("quality_chars", ""),
        placeholder=qc.get("placeholder_quality", ""),
        primer=(primers.get(d, {}).get("primer_f", "") or "")[:12],
        cutadapt=trim_pct,
        asv=qc.get("n_asv", ""),
        genus=genus_pct,
        rlen=rlen,
        merge=merge,
        reads_kept=keep,
    ))

excluded = []
if os.path.isdir(EXC):
    for d in sorted(os.listdir(EXC)):
        if os.path.isdir(os.path.join(EXC, d)):
            excluded.append(d)

lines = []
lines.append("# W1 汇总报告（自动生成）\n")
lines.append(f"生成时间：{__import__('datetime').datetime.now():%Y-%m-%d %H:%M:%S}\n")
lines.append(f"完成项目数：**{len(rows)}**\n")
if excluded:
    lines.append(f"已排除（占位符质量）：{', '.join(excluded)}\n")
lines.append("\n## 逐项目 QC\n")
hdr = ["项目", "机构", "年", "样本", "类型", "读长", "质量值种类", "占位符", "去引物%", "ASV数", "属级注释率", "合并率", "最终保留率"]
lines.append("| " + " | ".join(hdr) + " |")
lines.append("|" + "---|" * len(hdr))
for r in rows:
    lines.append("| {project} | {centre} | {year} | {samples} | {mode} | {rlen} | {qual} | {placeholder} | {cutadapt} | {asv} | {genus} | {merge} | {reads_kept} |".format(**r))

tot_s = sum(int(float(r["samples"])) for r in rows if str(r["samples"]).strip().isdigit())
tot_a = sum(int(float(r["asv"])) for r in rows if str(r["asv"]).strip().isdigit())
lines.append(f"\n**合计样本数：{tot_s}　合计 ASV：{tot_a}**\n")

with io.open(os.path.join(OUT, "W1_summary.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print("projects summarized:", len(rows))
print("total samples:", tot_s, " total ASVs:", tot_a)
print("written:", os.path.join(OUT, "W1_summary.md"))

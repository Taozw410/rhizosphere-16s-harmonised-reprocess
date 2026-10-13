import gzip, glob, os, io, csv
IUPAC = {"A":"A","C":"C","G":"G","T":"T","U":"T","R":"AG","Y":"CT","S":"GC","W":"AT","K":"GT","M":"AC","B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT"}
def match(seq, pat):
    if len(seq) < len(pat): return False
    return all(s in IUPAC.get(p,"") for s,p in zip(seq,pat))
PRIMERS = {"515F":"GTGYCAGCMGCCGCGGTAA","806R":"GGACTACNVGGGTWTCTAAT","341F":"CCTACGGGNGGCWGCAG",
           "805R":"GACTACHVGGGTATCTAATCC","338F":"ACTCCTACGGGAGGCAGCAG","1392R":"ACGGGCGGTGTGTRC",
           "27F":"AGAGTTTGATCMTGGCTCAG","ITS1F":"CTTGGTCATTTAGAGGAAGTAA"}
meta = {r["study"]: r for r in csv.DictReader(io.open(r"D:\CodexWork\mic-re\meta\stage2_projects.csv", encoding="utf-8-sig"))}
base = r"D:\CodexWork\mic-re\data\probe2"
rows=[]
for proj in sorted(os.listdir(base)):
    d = os.path.join(base, proj)
    if not os.path.isdir(d): continue
    fs = sorted(glob.glob(os.path.join(d,"*.fastq.gz")))
    if not fs: continue
    f1 = [f for f in fs if f.endswith("_1.fastq.gz")] or [fs[0]]
    f = f1[0]
    q=set(); lens={}; seqs=[]
    try:
        with gzip.open(f,"rt",errors="replace") as fh:
            for i,line in enumerate(fh):
                r=i%4
                if r==1:
                    s=line.strip().upper(); lens[len(s)]=lens.get(len(s),0)+1
                    if len(seqs)<300: seqs.append(s)
                elif r==3: q.update(line.strip())
                if i>=4000: break
    except Exception:
        q=set(); lens={1:1}; seqs=[]
    lo,hi=(min(q),max(q)) if q else (0,0)
    hits=[n for n,p in PRIMERS.items() if sum(1 for s in seqs if match(s,p))>=150]
    med=sorted(lens.items(),key=lambda x:-x[1])[0][0] if lens else 0
    marker="16S"
    for s in seqs[:50]:
        if "GAATCATCGAATCTTTGAACGC" in s[:60] or "TCCGTAGGTGAACCTGCGG" in s[:60]: marker="ITS?"; break
    tier = "EXCLUDE" if len(q)<=2 else ("CHECK" if len(q)<=5 else "OK")
    m = meta.get(proj, {})
    rows.append(dict(project=proj, crop=m.get("crop",""), tier_orig=m.get("tier",""),
                     centre=(m.get("center","") or "")[:28], qual=len(q), qmin=lo, qmax=hi,
                     rlen=med, marker=marker, primers=",".join(hits) if hits else "-",
                     verdict=tier, samples=m.get("ill_amp","")))
rows.sort(key=lambda r: ({"OK":0,"CHECK":1,"EXCLUDE":2}[r["verdict"]], r["crop"], r["project"]))
out = [f"# 第二阶段（46 个项目）探针筛查结果（2026-10-07）\n",
       "## 判定规则\n",
       "- **OK**：质量值 ≥6 档，可直接进入下载与处理",
       "- **CHECK**：质量值 3–5 档（被量化），可处理但必须做敏感性分析",
       "- **EXCLUDE**：质量值 ≤2 档（占位符），无法估计错误率，直接排除\n",
       "## 汇总\n"]
for v,lab in [("OK","可直接使用"),("CHECK","需加敏感性分析"),("EXCLUDE","排除")]:
    n=len([r for r in rows if r["verdict"]==v])
    out.append(f"- **{v}（{lab}）**：{n} 个项目")
out.append("")
for v,lab in [("OK","✅ 可直接使用"),("CHECK","⚠️ 需加敏感性分析"),("EXCLUDE","❌ 排除")]:
    sub=[r for r in rows if r["verdict"]==v]
    if not sub: continue
    out.append(f"\n## {lab}（{len(sub)} 个）\n")
    out.append("| 项目 | 作物 | 机构 | 质量档 | 读长 | 标记 | 引物 | 样本 |")
    out.append("|---|---|---|---|---|---|---|---|")
    for r in sub:
        out.append(f"| {r['project']} | {r['crop']} | {r['centre']} | {r['qual']} | {r['rlen']} | {r['marker']} | {r['primers']} | {r['samples']} |")
with io.open(r"D:\CodexWork\mic-re\out50k\stage2_screening.md","w",encoding="utf-8") as f:
    f.write("\n".join(out))
with io.open(r"D:\CodexWork\mic-re\out50k\stage2_screening.csv","w",encoding="utf-8-sig",newline="") as f:
    w=csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
c={v:len([r for r in rows if r["verdict"]==v]) for v in ("OK","CHECK","EXCLUDE")}
print("OK =",c["OK"]," CHECK =",c["CHECK"]," EXCLUDE =",c["EXCLUDE"])

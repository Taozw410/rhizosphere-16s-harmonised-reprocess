import gzip, glob, os
IUPAC = {"A":"A","C":"C","G":"G","T":"T","U":"T","R":"AG","Y":"CT","S":"GC","W":"AT","K":"GT","M":"AC","B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT"}
def match(seq, pat):
    if len(seq) < len(pat): return False
    return all(s in IUPAC.get(p,"") for s,p in zip(seq,pat))
PRIMERS = {"515F":"GTGYCAGCMGCCGCGGTAA","806R":"GGACTACNVGGGTWTCTAAT","341F":"CCTACGGGNGGCWGCAG",
           "805R":"GACTACHVGGGTATCTAATCC","338F":"ACTCCTACGGGAGGCAGCAG","1392R":"ACGGGCGGTGTGTRC",
           "27F":"AGAGTTTGATCMTGGCTCAG","ITS1F":"CTTGGTCATTTAGAGGAAGTAA"}
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
        rows.append((proj,"READ_ERR",0,0,"","")); continue
    lo,hi=(min(q),max(q)) if q else (0,0)
    placeholder = len(q)<=2
    hits=[n for n,p in PRIMERS.items() if sum(1 for s in seqs if match(s,p))>=150]
    med=sorted(lens.items(),key=lambda x:-x[1])[0][0] if lens else 0
    # 标记基因初判
    marker="16S"
    for s in seqs[:50]:
        if "GAATCATCGAATCTTTGAACGC" in s[:60] or "TCCGTAGGTGAACCTGCGG" in s[:60]: marker="ITS?"; break
    rows.append((proj, len(q), med, med, marker if marker!="16S" else "16S", ",".join(hits) if hits else "-"))
print(f"{'project':<15}{'qual':>5}{'len':>6}  {'marker':<8} primers")
for r in rows:
    print(f"{r[0]:<15}{r[1]:>5}{r[2]:>6}  {r[4]:<8} {r[5]}")
print()
print("占位符(需排除):", [r[0] for r in rows if isinstance(r[1],int) and r[1]<=2])
print("疑似 ITS 混杂:", [r[0] for r in rows if r[4]=="ITS?"])

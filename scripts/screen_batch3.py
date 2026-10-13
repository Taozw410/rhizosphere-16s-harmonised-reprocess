# 批次3 综合筛查：标记基因（全域双向）+ 质量编码 + 引物 + 读长
import csv, gzip, glob, os, io
R = r"D:\CodexWork\mic-re"
SIG = {
 "ITS": ["CTTGGTCATTTAGAGGAAGTAA","TCCGTAGGTGAACCTGCGG","GAATCATCGAATCTTTGAACGC",
         "GCGTTCAAAGATTCGATGATTC","GCTGCGTTCTTCATCGATGC","CGATGAAGAACGCAGCG","CCGCTTATTGATATGCTTAAG"],
 "16S": ["GTGCCAGCMGCCGCGGTAA","GTGYCAGCMGCCGCGGTAA","CCTACGGGNGGCWGCAG","ACTCCTACGGGAGGCAGCAG",
         "ATTAGAWACCCBDGTAGTCC","GACTACHVGGGTATCTAATCC","GGACTACNVGGGTWTCTAAT","TACGGGAGGCAGCAG"],
 "adapter": ["CTACACGACGCTCTTCCGATCT","AATGATACGGCGACCACCGAGATCTACAC","CAAGCAGAAGACGGCATACGAGAT"],
 "amoA": ["GGGGTTTCTACTGGTGGT"], "nifH": ["GCIWTYTAYGGIAARGGIGG"],
}
PRIMERS = {"515F":"GTGYCAGCMGCCGCGGTAA","806R":"GGACTACNVGGGTWTCTAAT","341F":"CCTACGGGNGGCWGCAG",
           "805R":"GACTACHVGGGTATCTAATCC","338F":"ACTCCTACGGGAGGCAGCAG","ITS1F":"CTTGGTCATTTAGAGGAAGTAA"}
meta = {r["study"]: r for r in csv.DictReader(io.open(os.path.join(R,"meta","batch3_candidates.csv"), encoding="utf-8-sig"))}
base = os.path.join(R,"data","probe4")
rows=[]
for proj in sorted(os.listdir(base)):
    d=os.path.join(base,proj)
    if not os.path.isdir(d): continue
    fs=sorted(glob.glob(os.path.join(d,"*.fastq.gz")))
    if not fs: continue
    paired = any(x.endswith("_1.fastq.gz") for x in fs)
    f=([x for x in fs if x.endswith("_1.fastq.gz")] or fs)[0]
    q=set(); lens={}; n=0; hits={k:0 for k in SIG}
    try:
        with gzip.open(f,"rt",errors="replace") as fh:
            for i,line in enumerate(fh):
                r=i%4
                if r==1:
                    s=line.strip().upper(); n+=1; lens[len(s)]=lens.get(len(s),0)+1
                    for k,ps in SIG.items():
                        if any(x in s for x in ps): hits[k]+=1
                elif r==3: q.update(line.strip())
                if n>=400: break
    except Exception:
        rows.append(dict(project=proj,crop=meta.get(proj,{}).get("crop",""),verdict="READ_ERROR")); continue
    med=sorted(lens.items(), key=lambda x:-x[1])[0][0] if lens else 0
    its=100*hits["ITS"]/n; s16=100*hits["16S"]/n; ad=100*hits["adapter"]/n
    amo=100*hits["amoA"]/n; nif=100*hits["nifH"]/n
    # 判定
    if len(q)<=2: v="EXCLUDE_QUALITY"
    elif ad>30: v="EXCLUDE_ADAPTER"
    elif its>30 and its>s16: v="EXCLUDE_ITS"
    elif amo>20 or nif>20: v="EXCLUDE_FUNC"
    elif its>10 and its>=s16*0.5: v="CHECK_MIXED"
    elif med<100: v="EXCLUDE_SHORT"
    else: v="OK"
    rows.append(dict(project=proj, crop=meta.get(proj,{}).get("crop",""), mode=("paired" if paired else "single"),
                     qual=len(q), rlen=med, ITS=round(its), S16=round(s16), adapter=round(ad),
                     amoA=round(amo), nifH=round(nif),
                     samples=meta.get(proj,{}).get("ill_amp",""), verdict=v))
rows.sort(key=lambda r:(r["verdict"], r["crop"], r["project"]))
print("%-15s %-7s %-7s %4s %5s %5s %5s %5s  %s" % ("project","crop","mode","qual","rlen","ITS%","16S%","adpt%","判定"))
for r in rows:
    print("%-15s %-7s %-7s %4d %5d %5d %5d %5d  %s" % (r["project"],r["crop"],r["mode"],r["qual"],r["rlen"],r["ITS"],r["S16"],r["adapter"],r["verdict"]))
print()
from collections import Counter
c=Counter(r["verdict"] for r in rows)
for k,v in c.most_common(): print(f"  {k:<18} {v}")
with io.open(os.path.join(R,"out50k","batch3_screening.csv"),"w",encoding="utf-8-sig",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

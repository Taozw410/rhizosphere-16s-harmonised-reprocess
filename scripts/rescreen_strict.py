import csv, gzip, glob, os, io
R = r"D:\CodexWork\mic-re"
SIG_ITS = ["CTTGGTCATTTAGAGGAAGTAA","TCCGTAGGTGAACCTGCGG","GAATCATCGAATCTTTGAACGC",
           "GCGTTCAAAGATTCGATGATTC","GCTGCGTTCTTCATCGATGC","CGATGAAGAACGCAGCG","CCGCTTATTGATATGCTTAAG"]
SIG_16S = ["GTGCCAGCMGCCGCGGTAA","GTGYCAGCMGCCGCGGTAA","CCTACGGGNGGCWGCAG","ACTCCTACGGGAGGCAGCAG",
           "ATTAGAWACCCBDGTAGTCC","GACTACHVGGGTATCTAATCC","GGACTACNVGGGTWTCTAAT","TACGGGAGGCAGCAG"]
SIG_AD  = ["CTACACGACGCTCTTCCGATCT","AATGATACGGCGACCACCGAGATCTACAC","CAAGCAGAAGACGGCATACGAGAT"]
SIG_FN  = ["GGGGTTTCTACTGGTGGT","GCIWTYTAYGGIAARGGIGG"]
meta={}
for f in os.listdir(os.path.join(R,"meta")):
    if not f.endswith(".csv"): continue
    try:
        for r in csv.DictReader(io.open(os.path.join(R,"meta",f), encoding="utf-8-sig")):
            k=r.get("study") or r.get("project")
            if k and r.get("crop"): meta.setdefault(k,r)
    except Exception: pass
processed=set()
mi=os.path.join(R,"out50k","master_inventory.csv")
if os.path.exists(mi):
    for r in csv.DictReader(io.open(mi,encoding="utf-8-sig")): processed.add(r["project"])
kept=set(x.strip() for x in io.open(os.path.join(R,"meta","batch3_selected.txt"),encoding="utf-8-sig") if x.strip()) if os.path.exists(os.path.join(R,"meta","batch3_selected.txt")) else set()

rows=[]; seen=set()
for pbase in ["probe4","probe3","probe2","probe"]:
    base=os.path.join(R,"data",pbase)
    if not os.path.isdir(base): continue
    for proj in sorted(os.listdir(base)):
        if proj in seen: continue
        d=os.path.join(base,proj)
        if not os.path.isdir(d): continue
        fs=sorted(glob.glob(os.path.join(d,"*.fastq.gz")))
        if not fs: continue
        f=([x for x in fs if x.endswith("_1.fastq.gz")] or fs)[0]
        q=set(); lens={}; n=0; hits={"ITS":0,"16S":0,"ad":0,"fn":0}
        try:
            with gzip.open(f,"rt",errors="replace") as fh:
                for i,line in enumerate(fh):
                    r=i%4
                    if r==1:
                        s=line.strip().upper(); n+=1; lens[len(s)]=lens.get(len(s),0)+1
                        if any(x in s for x in SIG_ITS): hits["ITS"]+=1
                        if any(x in s for x in SIG_16S): hits["16S"]+=1
                        if any(x in s for x in SIG_AD): hits["ad"]+=1
                        if any(x in s for x in SIG_FN): hits["fn"]+=1
                    elif r==3: q.update(line.strip())
                    if n>=400: break
        except Exception: continue
        seen.add(proj)
        med=sorted(lens.items(), key=lambda x:-x[1])[0][0] if lens else 0
        its=100*hits["ITS"]/n; s16=100*hits["16S"]/n; ad=100*hits["ad"]/n; fn=100*hits["fn"]/n
        m=meta.get(proj,{}); crop=m.get("crop","?"); paired=any(x.endswith("_1.fastq.gz") for x in fs)
        if len(q)<6: v="EXCLUDE_QUALITY"
        elif ad>20: v="EXCLUDE_ADAPTER"
        elif its>20 and its>s16: v="EXCLUDE_ITS"
        elif fn>15: v="EXCLUDE_FUNC"
        elif med<100: v="EXCLUDE_SHORT"
        elif (not paired) and med>380: v="CHECK_LONGREAD"
        else: v="OK"
        rows.append(dict(project=proj, crop=crop, mode=("paired" if paired else "single"),
                         qual=len(q), rlen=med, ITS=round(its), S16=round(s16), adpt=round(ad),
                         func=round(fn), samples=m.get("ill_amp",""),
                         centre=(m.get("center","") or "")[:22],
                         used=("YES" if proj in processed or proj in kept else ""), verdict=v))
ok=[r for r in rows if r["verdict"]=="OK" and r["used"]!="YES"]
print("扫描候选总数 = %d   其中未被使用的 OK = %d" % (len(rows), len(ok)))
print()
print("%-15s %-7s %-7s %4s %5s %5s %5s  %s" % ("project","crop","mode","qual","rlen","ITS%","16S%","样本"))
for r in sorted(ok, key=lambda x:(x["crop"], -float(x["samples"] or 0))):
    print("%-15s %-7s %-7s %4d %5d %5d %5d  %s" % (r["project"],r["crop"],r["mode"],r["qual"],r["rlen"],r["ITS"],r["S16"],r["samples"]))
from collections import Counter
print()
print("OK 候选按作物：", dict(Counter(r["crop"] for r in ok)))
with io.open(os.path.join(R,"out50k","rescreen_strict.csv"),"w",encoding="utf-8-sig",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

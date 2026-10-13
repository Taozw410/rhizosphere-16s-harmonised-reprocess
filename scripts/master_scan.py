import re
# 全量扫描：所有已下载项目 + 已处理项目 + 标记基因判定
import csv, io, os, glob, gzip

R = r"D:\CodexWork\mic-re"
SIG_ITS = ["CTTGGTCATTTAGAGGAAGTAA","TCCGTAGGTGAACCTGCGG","GAATCATCGAATCTTTGAACGC",
           "GCGTTCAAAGATTCGATGATTC","GCTGCGTTCTTCATCGATGC","CGATGAAGAACGCAGCG",
           "GTTTCCGTAGGTGAACCTGCGG","CCGCTTATTGATATGCTTAAG"]
SIG_16S = ["GTGCCAGCAGCCGCGGTAA","GTGYCAGCMGCCGCGGTAA","CCTACGGGNGGCWGCAG","ACTCCTACGGGAGGCAGCAG",
           "ATTAGAWACCCBDGTAGTCC","GACTACHVGGGTATCTAATCC","GGACTACNVGGGTWTCTAAT","TACGGGAGGCAGCAG",
           "GTGCTCCCCCGCCAATTCCTT"]
SIG_OTH = {"adapter":["CTACACGACGCTCTTCCGATCT","AATGATACGGCGACCACCGAGATCTACAC","CAAGCAGAAGACGGCATACGAGAT"],
           "amoA":["GGGGTTTCTACTGGTGGT"], "nifH":["GCIWTYTAYGGIAARGGIGG"]}

def scan_asv(p):
    f = os.path.join(p, "asv_table.csv")
    if not os.path.exists(f): return None
    hits={"ITS":0,"16S":0}; oth={}; n=0
    with io.open(f, encoding="utf-8-sig") as fh:
        r=csv.reader(fh); next(r)
        for row in r:
            n+=1; s=row[0].upper()
            if any(x in s for x in SIG_ITS): hits["ITS"]+=1
            elif any(x in s for x in SIG_16S): hits["16S"]+=1
            for k,ps in SIG_OTH.items():
                if any(x in s for x in ps): oth[k]=oth.get(k,0)+1
            if n>=500: break
    gp=None
    tp=os.path.join(p,"asv_taxonomy.csv")
    if os.path.exists(tp):
        t=list(csv.DictReader(io.open(tp,encoding="utf-8-sig")))
        if t:
            g=sum(1 for x in t if (x.get("Genus") or "").strip() and (x.get("Genus") or "").strip().upper()!="NA")
            gp=100*g/len(t)
    return n,hits,oth,gp

def scan_reads(d):
    fs=sorted(glob.glob(os.path.join(d,"*.fastq.gz")))
    if not fs: return None
    f=([x for x in fs if x.endswith("_1.fastq.gz")] or fs)[0]
    hits={k:0 for k in ["ITS","16S","adapter","amoA","nifH"]}; n=0
    try:
        with gzip.open(f,"rt",errors="replace") as fh:
            for i,line in enumerate(fh):
                if i%4==1:
                    s=line.strip().upper(); n+=1
                    if any(x in s for x in SIG_ITS): hits["ITS"]+=1
                    if any(x in s for x in SIG_16S): hits["16S"]+=1
                    for k,ps in SIG_OTH.items():
                        if any(x in s for x in ps): hits[k]+=1
                if n>=400: break
    except Exception: return None
    return n,hits

# 收集所有项目
projs={}
for base,tag in [("out50k","W1"),("out_stage2","S2")]:
    p=os.path.join(R,base)
    if not os.path.isdir(p): continue
    for d in sorted(os.listdir(p)):
        if not re.match(r"^PRJ[A-Z]*\d+$", d): continue
        full=os.path.join(p,d)
        if not os.path.isdir(full): continue
        projs.setdefault(d,{})["out"]=os.path.join(base,d)
        projs.setdefault(d,{})["source"]=tag
for base in ["raw","raw50k","raw50k_longread","raw50k_w1done"]:
    p=os.path.join(R,"data",base)
    if not os.path.isdir(p): continue
    for d in sorted(os.listdir(p)):
        full=os.path.join(p,d)
        if os.path.isdir(full): projs.setdefault(d,{})[base]=True

rows=[]
for proj,info in sorted(projs.items()):
    outp = os.path.join(R, info["out"]) if "out" in info else None
    r = scan_asv(outp) if outp else None
    verdict="NOT_PROCESSED"; its=s16=gp=None; oth=""
    if r:
        n,hits,o,gp = r
        its=round(100*hits["ITS"]/n); s16=round(100*hits["16S"]/n)
        oth=",".join(f"{k}={v}" for k,v in o.items())
        if gp is not None and gp<20: verdict="EXCLUDE_MARKER"
        elif its>=50 and its>s16: verdict="EXCLUDE_MARKER"
        elif o and max(o.values())>=0.2*n: verdict="EXCLUDE_MARKER"
        elif (gp is not None and gp<60) or (its>=20 and (gp is None or gp<80)): verdict="CHECK_MIXED"
        else: verdict="OK_16S"
    else:
        # 未处理：用 reads 粗略判断
        rd = info.get("raw50k") and os.path.join(R,"data","raw50k",proj) or \
             (info.get("raw50k_longread") and os.path.join(R,"data","raw50k_longread",proj)) or \
             (info.get("raw") and os.path.join(R,"data","raw",proj))
        rr = scan_reads(rd) if rd else None
        if rr:
            n,h=rr
            its=round(100*h["ITS"]/n); s16=round(100*h["16S"]/n)
            ad=round(100*h["adapter"]/n)
            if its>30 and its>s16: verdict="NOT_PROCESSED_SUSPECT_ITS"
            elif ad>30: verdict="NOT_PROCESSED_SUSPECT_ADAPTER"
            else: verdict="NOT_PROCESSED"
    rows.append(dict(project=proj, source=info.get("source",""), out=info.get("out",""),
                     raw="Y" if info.get("raw") else "", raw50k="Y" if info.get("raw50k") else "",
                     longread="Y" if info.get("raw50k_longread") else "",
                     its=its if its is not None else "", s16=s16 if s16 is not None else "",
                     genus=(round(gp,1) if gp is not None else ""), other=oth, verdict=verdict))

with io.open(os.path.join(R,"out50k","master_inventory.csv"),"w",encoding="utf-8-sig",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

print("%-15s %-4s %5s %5s %7s %-22s %s" % ("project","src","ITS%","16S%","genus%","other","verdict"))
for r in rows:
    print("%-15s %-4s %5s %5s %7s %-22s %s" % (r["project"],r["source"],r["its"],r["s16"],r["genus"],r["other"],r["verdict"]))
print()
from collections import Counter
c=Counter(r["verdict"] for r in rows)
for k,v in c.most_common(): print(f"  {k:<26} {v}")
print("  总项目数:", len(rows))

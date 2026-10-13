import csv, gzip, glob, os, io
from collections import Counter
R = r"D:\CodexWork\mic-re"
meta={}
for f in os.listdir(os.path.join(R,"meta")):
    if not f.endswith(".csv"): continue
    try:
        for r in csv.DictReader(io.open(os.path.join(R,"meta",f), encoding="utf-8-sig")):
            k=r.get("study") or r.get("project")
            if k and r.get("crop"): meta.setdefault(k,r)
    except Exception: pass
rows=[]
seen=set()
for pb in ["probe","probe2","probe3","probe4"]:
    b=os.path.join(R,"data",pb)
    if not os.path.isdir(b): continue
    for proj in sorted(os.listdir(b)):
        if proj in seen: continue
        fs=sorted(glob.glob(os.path.join(b,proj,"*.fastq.gz")))
        if not fs: continue
        f=([x for x in fs if x.endswith("_1.fastq.gz")] or fs)[0]
        q=set()
        try:
            with gzip.open(f,"rt",errors="replace") as fh:
                n=0
                for i,line in enumerate(fh):
                    if i%4==3: q.update(line.strip()); n+=1
                    if n>=200: break
        except Exception: continue
        seen.add(proj)
        rows.append((proj, meta.get(proj,{}).get("crop","?"), len(q)))
from collections import defaultdict
bycrop=defaultdict(lambda: Counter())
for p,c,n in rows:
    bycrop[c][(">=6" if n>=6 else ("3-5" if n>=3 else "<=2"))] += 1
print("%-8s %6s %6s %6s %6s" % ("crop","total",">=6","3-5","<=2"))
for c in ("maize","rice","wheat","?"):
    if c not in bycrop: continue
    d=bycrop[c]; t=sum(d.values())
    print("%-8s %6d %6d %6d %6d" % (c,t,d[">=6"],d["3-5"],d["<=2"]))
print()
print("质量 >=6 的候选（按作物）：")
for c in ("maize","rice","wheat"):
    g=[p for p,cc,n in rows if cc==c and n>=6]
    print(f"  {c}: {len(g)} 个")
    if len(g)<=25: print("   ", ", ".join(g) if g else "(无)")

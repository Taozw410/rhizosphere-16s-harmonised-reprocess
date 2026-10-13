import gzip, glob, os, io, csv, sys
IUPAC = {"A":"A","C":"C","G":"G","T":"T","U":"T","R":"AG","Y":"CT","S":"GC","W":"AT","K":"GT","M":"AC","B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT"}
def match(seq, pat):
    if len(seq) < len(pat): return False
    return all(s in IUPAC.get(p,"") for s,p in zip(seq,pat))
PRIMERS = {"515F":"GTGYCAGCMGCCGCGGTAA","806R":"GGACTACNVGGGTWTCTAAT","341F":"CCTACGGGNGGCWGCAG",
           "805R":"GACTACHVGGGTATCTAATCC","338F":"ACTCCTACGGGAGGCAGCAG","1392R":"ACGGGCGGTGTGTRC",
           "27F":"AGAGTTTGATCMTGGCTCAG","ITS1F":"CTTGGTCATTTAGAGGAAGTAA"}
base = sys.argv[1]
meta = {r["study"]: r for r in csv.DictReader(io.open(r"D:\CodexWork\mic-re\meta\batch2_candidates.csv", encoding="utf-8-sig"))}
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
    hits=[n for n,p in PRIMERS.items() if sum(1 for s in seqs if match(s,p))>=150]
    med=sorted(lens.items(),key=lambda x:-x[1])[0][0] if lens else 0
    marker="16S"
    for s in seqs[:50]:
        if "GAATCATCGAATCTTTGAACGC" in s[:60] or "TCCGTAGGTGAACCTGCGG" in s[:60]: marker="ITS?"; break
    verdict = "EXCLUDE" if len(q)<=2 else ("CHECK" if len(q)<=5 else "OK")
    m = meta.get(proj, {})
    rows.append(dict(project=proj, crop=m.get("crop",""), centre=(m.get("center","") or "")[:26],
                     qual=len(q), rlen=med, marker=marker,
                     primers=",".join(hits) if hits else "-", verdict=verdict,
                     samples=m.get("ill_amp",""), marker_meta=m.get("marker","")))
rows.sort(key=lambda r: ({"OK":0,"CHECK":1,"EXCLUDE":2}[r["verdict"]], r["crop"], r["project"]))
print(f"{'project':<15}{'crop':<7}{'qual':>5}{'len':>6}  {'marker':<7}{'verdict':<9} primers")
for r in rows:
    print(f"{r['project']:<15}{r['crop']:<7}{r['qual']:>5}{r['rlen']:>6}  {r['marker']:<7}{r['verdict']:<9}{r['primers']}")
print()
for crop in ("maize","rice","wheat"):
    g=[r for r in rows if r["crop"]==crop]
    o=len([r for r in g if r["verdict"]=="OK"]); c=len([r for r in g if r["verdict"]=="CHECK"]); e=len([r for r in g if r["verdict"]=="EXCLUDE"])
    print(f"{crop:<8} OK={o:<3} CHECK={c:<3} EXCLUDE={e}")
with io.open(r"D:\CodexWork\mic-re\out50k\batch2_screening.csv","w",encoding="utf-8-sig",newline="") as f:
    w=csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

import re
# 标记基因可靠判据：属级注释率（主）+ ASV 特征序列（辅）
import csv, io, os
SIG_ITS = ["CTTGGTCATTTAGAGGAAGTAA","TCCGTAGGTGAACCTGCGG","GAATCATCGAATCTTTGAACGC",
           "GCGTTCAAAGATTCGATGATTC","GCTGCGTTCTTCATCGATGC","CGATGAAGAACGCAGCG",
           "GTTTCCGTAGGTGAACCTGCGG","CCGCTTATTGATATGCTTAAG"]
SIG_16S = ["GTGCCAGCAGCCGCGGTAA","GTGYCAGCMGCCGCGGTAA","CCTACGGGNGGCWGCAG","ACTCCTACGGGAGGCAGCAG",
           "ATTAGAWACCCBDGTAGTCC","GACTACHVGGGTATCTAATCC","GGACTACNVGGGTWTCTAAT","TACGGGAGGCAGCAG",
           "GTGCTCCCCCGCCAATTCCTT"]
SIG_OTH = {"adapter":["CTACACGACGCTCTTCCGATCT","AATGATACGGCGACCACCGAGATCTACAC","CAAGCAGAAGACGGCATACGAGAT"],
           "amoA":["GGGGTTTCTACTGGTGGT","GGACTTCACCCTGTATCC"], "nifH":["GCIWTYTAYGGIAARGGIGG"]}

def scan(d):
    f = os.path.join(d, "asv_table.csv")
    tax = os.path.join(d, "asv_taxonomy.csv")
    if not os.path.exists(f): return None
    hits = {"ITS":0,"16S":0}; oth = {}
    n = 0
    with io.open(f, encoding="utf-8-sig") as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            n += 1; s = row[0].upper()
            if any(p in s for p in SIG_ITS): hits["ITS"] += 1
            elif any(p in s for p in SIG_16S): hits["16S"] += 1
            for k, ps in SIG_OTH.items():
                if any(p in s for p in ps): oth[k] = oth.get(k,0)+1
            if n >= 500: break
    gp = None
    if os.path.exists(tax):
        t = list(csv.DictReader(io.open(tax, encoding="utf-8-sig")))
        if t:
            g = sum(1 for x in t if (x.get("Genus") or "").strip() and (x.get("Genus") or "").strip().upper()!="NA")
            gp = 100*g/len(t)
    return n, hits, oth, gp

rows = []
for base, src in [(r"D:\CodexWork\mic-re\out_stage2","S2"), (r"D:\CodexWork\mic-re\out50k","W1")]:
    if not os.path.isdir(base): continue
    for proj in sorted(os.listdir(base)):
        p = os.path.join(base, proj)
        if not os.path.isdir(p) or not re.match(r"^PRJ[A-Z]*\d+$", proj): continue
        r = scan(p)
        if not r: continue
        n, hits, oth, gp = r
        its_pct = 100*hits["ITS"]/n if n else 0
        s16_pct = 100*hits["16S"]/n if n else 0
        # 判定
        if gp is not None and gp < 20:
            verdict = "EXCLUDE_MARKER"
        elif its_pct >= 50 and its_pct > s16_pct:
            verdict = "EXCLUDE_MARKER"
        elif max(oth.values()) if oth else 0 >= 0.2*n:
            verdict = "EXCLUDE_MARKER"
        elif gp is not None and gp < 60:
            verdict = "CHECK_MIXED"
        elif its_pct >= 20 and gp is not None and gp < 80:
            verdict = "CHECK_MIXED"
        else:
            verdict = "OK_16S"
        rows.append(dict(project=proj, source=src, n=n, its=round(its_pct), s16=round(s16_pct),
                         other=",".join(f"{k}={v}" for k,v in oth.items()), genus=(round(gp,1) if gp is not None else ""),
                         verdict=verdict))
rows.sort(key=lambda r: (r["verdict"], r["project"]))
print("%-15s %-4s %5s %6s %6s %8s %-14s %s" % ("project","src","ASV","ITS%","16S%","genus%","other","判定"))
for r in rows:
    print("%-15s %-4s %5d %6d %6d %8s %-14s %s" % (r["project"],r["source"],r["n"],r["its"],r["s16"],r["genus"],r["other"],r["verdict"]))
print()
for v in ("OK_16S","CHECK_MIXED","EXCLUDE_MARKER"):
    g=[r for r in rows if r["verdict"]==v]
    print(f"  {v:<16} {len(g)}")
with io.open(r"D:\CodexWork\mic-re\out50k\marker_gene_qc.csv","w",encoding="utf-8-sig",newline="") as f:
    w=csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

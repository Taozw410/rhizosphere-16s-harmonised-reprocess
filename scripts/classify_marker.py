import re
import csv, io, os, sys
SIG = {
 "16S":  ["GTGCCAGCAGCCGCGGTAA","GTGYCAGCMGCCGCGGTAA","CCTACGGGNGGCWGCAG","ACTCCTACGGGAGGCAGCAG",
          "ATTAGAWACCCBDGTAGTCC","GACTACHVGGGTATCTAATCC","GGACTACNVGGGTWTCTAAT","TACGGGAGGCAGCAG",
          "GACTACCAGGGTATCTAATCC","GTGCTCCCCCGCCAATTCCTT"],
 "ITS":  ["CTTGGTCATTTAGAGGAAGTAA","TCCGTAGGTGAACCTGCGG","GAATCATCGAATCTTTGAACGC","GCGTTCAAAGATTCGATGATTC",
          "GCTGCGTTCTTCATCGATGC","CGATGAAGAACGCAGCG","GTTTCCGTAGGTGAACCTGCGG","CCGCTTATTGATATGCTTAAG"],
 "18S":  ["GTCTGGTGCCAGCAGCCGCGG"],
 "amoA": ["GGGGTTTCTACTGGTGGT","GGACTTCACCCTGTATCC"],
 "nifH": ["GCIWTYTAYGGIAARGGIGG","TACGGIAARGGIGG"],
 "adapter": ["CTACACGACGCTCTTCCGATCT","AATGATACGGCGACCACCGAGATCTACAC","CAAGCAGAAGACGGCATACGAGAT","AGATCGGAAGAGC"],
}
def classify(d):
    f = os.path.join(d, "asv_table.csv")
    if not os.path.exists(f): return None
    hits = {k:0 for k in SIG}; n = 0
    with io.open(f, encoding="utf-8-sig") as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            n += 1; s = row[0].upper()
            for k, pats in SIG.items():
                for p in pats:
                    if p in s: hits[k] += 1; break
            if n >= 500: break
    if n == 0: return None
    frac = {k: v/n for k,v in hits.items()}
    main = max(frac, key=frac.get)
    return n, frac, main

for base, label in [(r"D:\CodexWork\mic-re\out_stage2","阶段二"), (r"D:\CodexWork\mic-re\out50k","W1")]:
    if not os.path.isdir(base): continue
    print("="*70)
    print(label, base)
    rows = []
    for proj in sorted(os.listdir(base)):
        p = os.path.join(base, proj)
        if not os.path.isdir(p) or not re.match(r"^PRJ[A-Z]*\d+$", proj): continue
        r = classify(p)
        if not r: continue
        n, frac, main = r
        top = sorted([(k,v) for k,v in frac.items() if v>0.05], key=lambda x:-x[1])
        rows.append((proj, n, main, top, frac))
    print("%-15s %6s %-10s %s" % ("project","ASV","dominant","hits(>5%)"))
    for proj,n,main,top,frac in rows:
        print("%-15s %6d %-10s %s" % (proj,n,main, "  ".join("%s=%.0f%%"%(k,v*100) for k,v in top)))
    print()
    ok = [r for r in rows if r[2]=="16S" and r[4].get("ITS",0)<0.2]
    its = [r for r in rows if r[2]=="ITS"]
    other = [r for r in rows if r[2] not in ("16S","ITS")]
    mixed = [r for r in rows if r[2]=="16S" and r[4].get("ITS",0)>=0.2]
    print("  16S 为主: %d   含 ITS 混杂: %d   ITS 为主: %d   其他标记: %d" % (len(ok),len(mixed),len(its),len(other)))

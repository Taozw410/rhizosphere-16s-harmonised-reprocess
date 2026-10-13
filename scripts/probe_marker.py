import csv, io, os
OUT = r"D:\CodexWork\mic-re\out_stage2"
SIG = {
 "ITS2_fwd": "GAATCATCGAATCTTTGAACGC", "ITS2_rev": "GCGTTCAAAGATTCGATGATTC",
 "ITS1_fwd": "TCCGTAGGTGAACCTGCGG",     "ITS1_rev": "CCGCAGGTTCACCTACGGA",
 "ITS4_rc":  "GCTGCGTTCTTCATCGATGC",
 "16S_515":  "GTGCCAGCAGCCGCGGTAA",     "16S_806rc":"ATTAGAWACCCBDGTAGTCC",
 "18S_univ": "GTCTGGTGCCAGCAGCCGCGG",
 "nifH":     "GCIWTYTAYGGIAARGGIGG",    "amoA":     "GGGGTTTCTACTGGTGGT",
 "phoD":     "CAGTGGGACGACCACGAGGT",    "ITS_5.8S": "CGATGAAGAACGCAGCG",
}
projs = ["PRJNA506578","PRJNA549031","PRJNA595042","PRJNA639921","PRJNA888942",
         "PRJNA1159152","PRJNA770245","PRJNA576298","PRJNA1145186","PRJNA910480"]
print("%-15s %7s  %s" % ("project","ASVlen","首个 ASV 序列 60bp"))
for p in projs:
    f = os.path.join(OUT, p, "asv_table.csv")
    if not os.path.exists(f): continue
    with io.open(f, encoding="utf-8-sig") as fh:
        r = csv.reader(fh); next(r)
        row = next(r, None)
    if row:
        print("%-15s %7d  %s" % (p, len(row[0]), row[0][:60]))
print()
print("%-15s %5s  %s" % ("project","n","特征序列命中"))
for p in projs:
    f = os.path.join(OUT, p, "asv_table.csv")
    if not os.path.exists(f): continue
    hit = {k:0 for k in SIG}; n = 0
    with io.open(f, encoding="utf-8-sig") as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            n += 1; s = row[0].upper()
            for k,v in SIG.items():
                if v in s: hit[k] += 1
            if n >= 500: break
    top = sorted([(k,v) for k,v in hit.items() if v>0], key=lambda x:-x[1])[:4]
    print("%-15s %5d  %s" % (p, n, "  ".join("%s=%d" % kv for kv in top) if top else "(无命中)"))

# 标记基因全域检测（探针阶段用）：在 reads 的任意位置、双向搜索特征序列
import gzip, glob, os, sys
SIG = {
 "ITS": ["CTTGGTCATTTAGAGGAAGTAA","TCCGTAGGTGAACCTGCGG","GAATCATCGAATCTTTGAACGC",
         "GCGTTCAAAGATTCGATGATTC","GCTGCGTTCTTCATCGATGC","CGATGAAGAACGCAGCG","CCGCTTATTGATATGCTTAAG"],
 "16S": ["GTGCCAGCAGCCGCGGTAA","GTGYCAGCMGCCGCGGTAA","CCTACGGGNGGCWGCAG","ACTCCTACGGGAGGCAGCAG",
         "ATTAGAWACCCBDGTAGTCC","GACTACHVGGGTATCTAATCC","GGACTACNVGGGTWTCTAAT","TACGGGAGGCAGCAG"],
 "adapter": ["CTACACGACGCTCTTCCGATCT","AATGATACGGCGACCACCGAGATCTACAC","CAAGCAGAAGACGGCATACGAGAT"],
 "amoA": ["GGGGTTTCTACTGGTGGT"], "nifH": ["GCIWTYTAYGGIAARGGIGG"],
}
def check(f, nreads=400):
    hits = {k:0 for k in SIG}; n=0
    with gzip.open(f, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                s = line.strip().upper(); n += 1
                for k, ps in SIG.items():
                    if any(p in s for p in ps): hits[k] += 1
            if n >= nreads: break
    return n, hits
base = sys.argv[1]
print("%-15s %5s  %s" % ("project","n","特征命中率"))
for proj in sorted(os.listdir(base)):
    d = os.path.join(base, proj)
    if not os.path.isdir(d): continue
    fs = sorted(glob.glob(os.path.join(d, "*.fastq.gz")))
    if not fs: continue
    f = ([x for x in fs if x.endswith("_1.fastq.gz")] or fs)[0]
    n, hits = check(f)
    top = sorted([(k,v) for k,v in hits.items() if v>0], key=lambda x:-x[1])
    verdict = ""
    if hits.get("ITS",0)/n > 0.3 and hits.get("ITS",0) > hits.get("16S",0): verdict = "<-- 疑 ITS"
    elif hits.get("adapter",0)/n > 0.3: verdict = "<-- 接头未切"
    elif any(hits.get(k,0)/n>0.2 for k in ("amoA","nifH")): verdict = "<-- 疑功能基因"
    print("%-15s %5d  %s  %s" % (proj, n, "  ".join("%s=%.0f%%"%(k,100*v/n) for k,v in top[:4]), verdict))

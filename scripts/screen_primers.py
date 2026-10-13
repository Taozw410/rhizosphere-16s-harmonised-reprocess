# 半自动引物筛查：检查 reads 开头是否含常见 16S 引物
import gzip, glob, os

def ambig_match(seq, pat):
    """pat 允许 IUPAC 简并碱基；返回是否匹配前缀"""
    iupac = {
        "A":"A","C":"C","G":"G","T":"T","U":"T",
        "R":"AG","Y":"CT","S":"GC","W":"AT","K":"GT","M":"AC",
        "B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT",
    }
    if len(seq) < len(pat): return False
    for s, p in zip(seq, pat):
        if s not in iupac.get(p, ""): return False
    return True

PRIMERS = {
    "515F":  ("GTGYCAGCMGCCGCGGTAA", "A"),
    "515Fb": ("GTGCCAGCMGCCGCGGTAA", "A"),
    "806R":  ("GGACTACNVGGGTWTCTAAT", "A"),
    "806Rb": ("GGACTACHVGGGTWTCTAAT", "A"),
    "341F":  ("CCTACGGGNGGCWGCAG", "A"),
    "805R":  ("GACTACHVGGGTATCTAATCC", "A"),
    "926F":  ("AAACTYAAAKGAATTGACGG", "A"),
    "1392R": ("ACGGGCGGTGTGTRC", "A"),
    "27F":   ("AGAGTTTGATCMTGGCTCAG", "A"),
    "338F":  ("ACTCCTACGGGAGGCAGCAG", "A"),
    "ITS1F": ("CTTGGTCATTTAGAGGAAGTAA", "A"),
    "ITS2":  ("GCTGCGTTCTTCATCGATGC", "A"),
}

raw = r"D:\CodexWork\mic-re\data\raw"
print(f"{'project':<15}{'read':<7}{'n':>7}   best primer hits")
for proj in sorted(os.listdir(raw)):
    d = os.path.join(raw, proj)
    if not os.path.isdir(d): continue
    for f in sorted(glob.glob(os.path.join(d, "*.fastq.gz")))[:2]:
        tag = os.path.basename(f).replace(".fastq.gz", "")
        seqs = []
        with gzip.open(f, "rt", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    seqs.append(line.strip().upper())
                if len(seqs) >= 300: break
        hits = []
        for name, (p, _) in PRIMERS.items():
            c = sum(1 for s in seqs if ambig_match(s, p))
            if c >= 30:
                hits.append(f"{name}={c}/{len(seqs)}")
        print(f"{proj:<15}{tag[-4:]:<7}{len(seqs):>7}   " + (", ".join(hits) if hits else "(no known primer at read start)"))

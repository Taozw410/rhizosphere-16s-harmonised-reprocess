import gzip, glob, os

IUPAC = {"A":"A","C":"C","G":"G","T":"T","U":"T","R":"AG","Y":"CT","S":"GC","W":"AT",
         "K":"GT","M":"AC","B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT"}
def match(seq, pat):
    if len(seq) < len(pat): return False
    for s, p in zip(seq, pat):
        if s not in IUPAC.get(p, ""): return False
    return True

PRIMERS = {
    "515F":"GTGYCAGCMGCCGCGGTAA", "806R":"GGACTACNVGGGTWTCTAAT",
    "341F":"CCTACGGGNGGCWGCAG", "805R":"GACTACHVGGGTATCTAATCC",
    "338F":"ACTCCTACGGGAGGCAGCAG", "1392R":"ACGGGCGGTGTGTRC",
    "27F":"AGAGTTTGATCMTGGCTCAG", "ITS1F":"CTTGGTCATTTAGAGGAAGTAA",
}

base = r"D:\CodexWork\mic-re\data\probe"
print(f"{'project':<15}{'qual':>5}{'qmin':>5}{'qmax':>5}  {'len':>5}  {'status':<26} primers")
for proj in sorted(os.listdir(base)):
    d = os.path.join(base, proj)
    if not os.path.isdir(d): continue
    fs = sorted(glob.glob(os.path.join(d, "*.fastq.gz")))
    if not fs: continue
    f1 = [f for f in fs if f.endswith("_1.fastq.gz")] or [fs[0]]
    f = f1[0]
    q = set(); lens = {}; seqs = []
    with gzip.open(f, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            r = i % 4
            if r == 1:
                s = line.strip().upper(); lens[len(s)] = lens.get(len(s), 0) + 1
                if len(seqs) < 300: seqs.append(s)
            elif r == 3:
                q.update(line.strip())
            if i >= 4000: break
    lo, hi = (min(q), max(q)) if q else (0, 0)
    status = "PLACEHOLDER(reject)" if len(q) <= 2 else ("quantized(check)" if len(q) <= 5 else "ok")
    hits = [n for n, p in PRIMERS.items() if sum(1 for s in seqs if match(s, p)) >= 150]
    medlen = sorted(lens.items(), key=lambda x: -x[1])[0][0] if lens else 0
    paired = any(x.endswith("_2.fastq.gz") for x in fs)
    print(f"{proj:<15}{len(q):>5}{lo:>5}{hi:>5}  {medlen:>5}  {status:<26} {','.join(hits) if hits else '-'}  {'paired' if paired else 'single'}")

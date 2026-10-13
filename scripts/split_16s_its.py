import gzip, glob, os, shutil, json
d = r"D:\CodexWork\mic-re\data\raw50k\PRJNA1479906"
out = r"D:\CodexWork\mic-re\data\raw50k_excluded_its\PRJNA1479906"
os.makedirs(out, exist_ok=True)

ITS_SIG = ("GAATCATCGAATCTTTGAACGC", "TCCGTAGGTGAACCTGCGG", "GCTGCGTTCTTCATCGATGC",
           "GTGAATCATCGAATCTTTGAACGC", "AACTTTYRRCAAYGGATCWCT", "TCCTCCGCTTATTGATATGC")
def classify(f):
    with gzip.open(f, "rt", errors="replace") as fh:
        n = 0
        for i, l in enumerate(fh):
            if i % 4 == 1:
                s = l.strip().upper()
                if any(k in s[:60] for k in ITS_SIG):
                    return "ITS"
                if s.startswith("TTGCCAG") or "AGCCGCGGTAA" in s[:40]:
                    return "16S"
                n += 1
            if n >= 50: break
    return "unknown"

samples = {}
for f in sorted(glob.glob(os.path.join(d, "*.fastq.gz"))):
    s = os.path.basename(f)[:-len(".fastq.gz")]
    if s not in samples: samples[s] = classify(f)

its = [s for s, t in samples.items() if t == "ITS"]
six = [s for s, t in samples.items() if t == "16S"]
unk = [s for s, t in samples.items() if t == "unknown"]
print("16S 样本 :", len(six))
print("ITS 样本 :", len(its))
print("未知     :", len(unk), unk)

moved = 0
for s in its + unk:
    for f in glob.glob(os.path.join(d, s + "*.fastq.gz")):
        shutil.move(f, os.path.join(out, os.path.basename(f)))
        moved += 1
print("已移出文件数:", moved)
print("保留在 raw50k 的 16S 文件数:", len(glob.glob(os.path.join(d, "*.fastq.gz"))))

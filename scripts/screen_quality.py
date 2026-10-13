import gzip, glob, os, collections
raw = r"D:\CodexWork\mic-re\data\raw"
print(f"{'project':<15}{'file':<28}{'qual_uniq':>10}{'minq':>6}{'maxq':>6}   status")
for proj in sorted(os.listdir(raw)):
    d = os.path.join(raw, proj)
    if not os.path.isdir(d): continue
    fs = sorted(glob.glob(os.path.join(d, "*.fastq.gz")))
    if not fs: continue
    f = fs[0]
    vals = set()
    n = 0
    try:
        with gzip.open(f, "rb") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 3:
                    vals.update(line.strip())
                    n += 1
                if n >= 300: break
    except Exception as e:
        print(f"{proj:<15}{os.path.basename(f):<28}   read error: {e}")
        continue
    lo, hi = (min(vals), max(vals)) if vals else (0, 0)
    status = "PLACEHOLDER(常数质量)" if len(vals) <= 2 else "ok"
    print(f"{proj:<15}{os.path.basename(f):<28}{len(vals):>10}{lo:>6}{hi:>6}   {status}")

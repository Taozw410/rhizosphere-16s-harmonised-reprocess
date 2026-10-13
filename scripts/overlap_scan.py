import gzip
comp = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(comp)[::-1]

def head(f, n=150):
    out = []
    with gzip.open(f, "rt", errors="replace") as fh:
        for i, l in enumerate(fh):
            if i % 4 == 1: out.append(l.strip().upper())
            if len(out) >= n: break
    return out

def best_overlap(r1, r2c, min_ov=20, max_mm=3):
    # r2c 与 r1 在某个偏移处重叠：r2c 从 r1 的位置 k 开始匹配
    best = None
    for k in range(0, len(r1) - min_ov + 1):
        ov = min(len(r1) - k, len(r2c))
        a = r1[k:k+ov]; b = r2c[:ov]
        mm = sum(1 for x, y in zip(a, b) if x != y)
        if mm <= max_mm:
            best = (k, ov, mm); break
    return best

for proj, run in [("PRJNA430104","SRR6474320"), ("PRJNA401395","SRR6007430")]:
    d = r"D:\CodexWork\mic-re\data\raw50k\\" + proj
    r1 = head(d + "\\" + run + "_1.fastq.gz")
    r2 = head(d + "\\" + run + "_2.fastq.gz")
    hits = 0; offs = {}
    for a, b in zip(r1, r2):
        r = best_overlap(a, rc(b))
        if r:
            hits += 1
            offs[r[0]] = offs.get(r[0], 0) + 1
    print(f"{proj:<15} 可重叠配对 {hits}/{len(r1)}   重叠起点分布 {dict(sorted(offs.items())[:8])}")
    if offs:
        k = max(offs, key=offs.get)
        amp = len(r1[0]) + len(r2[0]) - (len(r1[0]) - k) + k*0
        print(f"{'':<15} 主要重叠起点 k={k}  推得扩增子长度 ≈ {len(r1[0]) + k} bp")

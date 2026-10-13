import gzip, collections
f = r"D:\CodexWork\mic-re\data\raw\PRJNA1049404\SRR27115409_1.fastq.gz"
q = collections.Counter(); n = 0; lens = collections.Counter(); ncount = 0
with gzip.open(f, "rt", errors="replace") as fh:
    for i, line in enumerate(fh):
        r = i % 4
        if r == 1:
            lens[len(line.strip())] += 1
            ncount += line.count("N")
        elif r == 3:
            q[line.strip()[:8]] += 1
            n += 1
        if n >= 2000: break
print("reads:", n)
print("length dist:", dict(lens))
print("total N bases in first 2000 reads:", ncount)
print("top quality prefixes:")
for k, v in q.most_common(6):
    print(f"   {k}  x{v}")

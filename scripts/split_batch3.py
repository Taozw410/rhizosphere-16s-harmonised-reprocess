import csv, io
sel = [x.strip() for x in io.open(r"D:\CodexWork\mic-re\meta\batch3_selected.txt", encoding="utf-8-sig").read().splitlines() if x.strip()]
info=[]
for p in sel:
    d = r"D:\CodexWork\mic-re\data\raw50k\\" + p
    import os, glob
    fs = glob.glob(os.path.join(d,"*.fastq.gz"))
    paired = any(x.endswith("_1.fastq.gz") for x in fs)
    n = len([x for x in fs if x.endswith("_1.fastq.gz")]) if paired else len(fs)
    info.append((p,n))
info.sort(key=lambda x:-x[1])
A=[];B=[];sA=sB=0
for p,n in info:
    if sA<=sB: A.append(p); sA+=n
    else: B.append(p); sB+=n
io.open(r"D:\CodexWork\mic-re\meta\batch3_streamA.txt","w",encoding="utf-8").write(",".join(A))
io.open(r"D:\CodexWork\mic-re\meta\batch3_streamB.txt","w",encoding="utf-8").write(",".join(B))
print("A:", len(A), "projects", sA, "samples")
print("B:", len(B), "projects", sB, "samples")
print("A =", ",".join(A))
print("B =", ",".join(B))

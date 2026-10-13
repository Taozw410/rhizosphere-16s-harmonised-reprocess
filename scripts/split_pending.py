import os, glob
R = r"D:\CodexWork\mic-re"
short = "PRJNA596848,PRJNA665860,PRJNA996548,PRJNA1282585,PRJNA1334676,PRJNA1476760,PRJNA601112,PRJNA753114,PRJNA816337".split(",")
info=[]
for p in short:
    d=os.path.join(R,"data","raw50k",p)
    fs=glob.glob(os.path.join(d,"*.fastq.gz"))
    paired=any(x.endswith("_1.fastq.gz") for x in fs)
    n=len([x for x in fs if x.endswith("_1.fastq.gz")]) if paired else len(fs)
    info.append((p,n))
info.sort(key=lambda x:-x[1])
A=[];B=[];sA=sB=0
for p,n in info:
    if sA<=sB: A.append(p); sA+=n
    else: B.append(p); sB+=n
open(os.path.join(R,"meta","pend_shortA.txt"),"w").write(",".join(A))
open(os.path.join(R,"meta","pend_shortB.txt"),"w").write(",".join(B))
print("A:", ",".join(A), "(%d samples)"%sA)
print("B:", ",".join(B), "(%d samples)"%sB)

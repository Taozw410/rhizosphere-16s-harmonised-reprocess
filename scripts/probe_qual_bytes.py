import gzip
f = r"D:\CodexWork\mic-re\data\raw\PRJNA1049404\SRR27115409_1.fastq.gz"
with gzip.open(f, "rb") as fh:
    for i in range(8):
        line = fh.readline()
        if i % 4 == 3:
            print("QUAL raw bytes (first 40):", line[:40])
            print("QUAL distinct values    :", sorted(set(line.strip())))
            print("ASCII range             :", min(line.strip()), max(line.strip()))
        elif i % 4 == 0:
            print("HEAD:", line[:70])
print()
print("--- another project: PRJNA401395 ---")
f2 = r"D:\CodexWork\mic-re\data\raw\PRJNA401395\SRR6007430_1.fastq.gz"
with gzip.open(f2, "rb") as fh:
    for i in range(8):
        line = fh.readline()
        if i % 4 == 3:
            print("QUAL raw bytes (first 40):", line[:40])
            print("QUAL distinct values    :", sorted(set(line.strip()))[:12], "...")
            print("ASCII range             :", min(line.strip()), max(line.strip()))
print()
print("--- single-end project: PRJNA608090 ---")
import glob
f3 = sorted(glob.glob(r"D:\CodexWork\mic-re\data\raw\PRJNA608090\*.fastq.gz"))[0]
with gzip.open(f3, "rb") as fh:
    for i in range(8):
        line = fh.readline()
        if i % 4 == 3:
            print("QUAL raw bytes (first 40):", line[:40])
            print("ASCII range             :", min(line.strip()), max(line.strip()))

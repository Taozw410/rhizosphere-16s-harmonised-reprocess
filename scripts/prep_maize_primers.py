# -*- coding: utf-8 -*-
"""为扩池的玉米候选检测 5' 引物，产出 pipeline 用的引物表（truncLen 由 plan_trunclen 后续填）。
usage: python prep_maize_primers.py
"""
import csv, glob, gzip, io, os, statistics

R = r"D:\CodexWork\mic-re"
PROJS = ["PRJNA633917", "PRJNA662201", "PRJNA967263"]
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
F = {"515F": "GTGYCAGCMGCCGCGGTAA", "515Fb": "GTGCCAGCMGCCGCGGTAA", "341F": "CCTACGGGNGGCWGCAG",
     "338F": "ACTCCTACGGGAGGCAGCAG", "27F": "AGAGTTTGATCMTGGCTCAG", "926F": "AAACTYAAAKGAATTGACGG"}
Rv = {"806R": "GGACTACNVGGGTWTCTAAT", "806Rb": "GGACTACHVGGGTWTCTAAT", "805R": "GACTACHVGGGTATCTAATCC",
      "1392R": "ACGGGCGGTGTGTRC", "534R": "ATTACCGCGGCTGCTGG", "907R": "CCGTCAATTCCTTTGAGTTT"}
ITS = {"ITS1F": "CTTGGTCATTTAGAGGAAGTAA", "ITS2": "GCTGCGTTCTTCATCGATGC"}


def prefix(s, p):
    return len(s) >= len(p) and all(c in IUPAC.get(q, "") for c, q in zip(s, p))


def carriage(path, table, n=200):
    hit = tot = 0
    buf = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            buf.append(line.rstrip("\n"))
            if len(buf) == 4:
                s = buf[1].upper(); tot += 1
                if any(prefix(s, p) for p in table.values()):
                    hit += 1
                buf = []
            if tot >= n:
                break
    return (100.0 * hit / tot) if tot else 0.0, tot


rows = []
for p in PROJS:
    files = sorted(glob.glob(os.path.join(R, "data", "raw", p, "*.fastq.gz")))
    f1 = ([x for x in files if x.endswith("_1.fastq.gz")] or files)[0]
    paired = any(x.endswith("_2.fastq.gz") for x in files)
    f2 = ([x for x in files if x.endswith("_2.fastq.gz")] or [None])[0]
    best_f, cf, n = None, 0.0, 0
    for name, seq in F.items():
        c, n = carriage(f1, {name: seq})
        if c > cf:
            best_f, cf = name, c
    ci, _ = carriage(f1, ITS)
    best_r, cr = None, 0.0
    if paired and f2:
        for name, seq in Rv.items():
            c, _ = carriage(f2, {name: seq})
            if c > cr:
                best_r, cr = name, c
    pf = F.get(best_f, "") if cf >= 50 else ""
    pr = Rv.get(best_r, "") if cr >= 50 else ""
    note = "maizeextra_%s%.0f%%_%s%.0f%%_ITS%.0f%%" % (best_f or "none", cf, best_r or "none", cr, ci)
    rows.append(dict(project=p, primer_f=pf, primer_r=pr, truncLen_f=0, truncLen_r=0, note=note))
    print("%-13s paired=%-5s R1引物=%s %.0f%%  R2引物=%s %.0f%%  ITS=%.0f%%" % (p, paired, best_f, cf, best_r, cr, ci))

out = os.path.join(R, "meta", "primer_table_maize_raw.csv")
with io.open(out, "w", encoding="ascii", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["project", "primer_f", "primer_r", "truncLen_f", "truncLen_r", "note"])
    w.writeheader(); w.writerows(rows)
print("wrote", out)

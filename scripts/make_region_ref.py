# -*- coding: utf-8 -*-
"""Build region-matched SILVA references by in-silico primer extraction.

For each reference sequence: find the forward primer (IUPAC-aware regex) near the 5' end, then the
reverse primer's reverse complement downstream, and keep the intervening amplicon (primers included,
matching what dada2 sees when primers are NOT trimmed off). Headers (taxonomy) are preserved verbatim.

Usage: python make_region_ref.py --ref <silva.fa.gz> --name V4_515F_806R --f GTGYCAGCMGCCGCGGTAA \
        --r GGACTACNVGGGTWTCTAAT --out <out.fa.gz> [--max-len 800] [--limit N]
"""
import gzip, io, sys, re, argparse

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMP = str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN")


def rc(s):
    return s.translate(COMP)[::-1]


def to_regex(p):
    return "".join("[" + IUPAC.get(ch, "N") + "]" for ch in p.upper())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True); ap.add_argument("--name", required=True)
    ap.add_argument("--f", required=True); ap.add_argument("--r", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--min-len", type=int, default=100)
    ap.add_argument("--max-len", type=int, default=800); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    reF = re.compile(to_regex(a.f))
    reR = re.compile(to_regex(rc(a.r)))
    n_in = n_out = 0
    with gzip.open(a.ref, "rt", errors="replace") as fin, gzip.open(a.out, "wt") as fout:
        hdr = None
        for line in fin:
            if line.startswith(">"):
                if hdr is not None and seq:
                    n_in += 1
                    mF = reF.search(seq)
                    if mF:
                        mR = reR.search(seq, mF.end())
                        if mR:
                            amp = seq[mF.start():mR.end()]
                            if a.min_len <= len(amp) <= a.max_len:
                                fout.write(hdr); fout.write(amp + "\n"); n_out += 1
                hdr = line; seq = ""; continue
            seq += line.strip()
            if a.limit and n_in >= a.limit:
                break
        if hdr is not None and seq:
            n_in += 1
            mF = reF.search(seq)
            if mF:
                mR = reR.search(seq, mF.end())
                if mR:
                    amp = seq[mF.start():mR.end()]
                    if a.min_len <= len(amp) <= a.max_len:
                        fout.write(hdr); fout.write(amp + "\n"); n_out += 1
    print("%s: scanned %d seqs, wrote %d region-matched seqs -> %s" % (a.name, n_in, n_out, a.out))


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Build a DECIPHER-format training FASTA for IDTAXA from the local SILVA 138.1 reference.

同一 SILVA 138.1 来源、按"唯一分类学字符串去冗余后等距抽样"，避免 DECIPHER 官方预训练模型不可达的问题。
usage: python make_idtaxa_train.py [--src <gzipped dada2 train fasta>] [--cap 6000] [--out <fa>]
"""
import argparse, gzip, io, os

R = r"D:\CodexWork\mic-re"
DEF_SRC = r"D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz"


def records(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        hdr = None
        seq = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq)
                hdr, seq = line[1:], []
            else:
                seq.append(line.strip())
        if hdr is not None:
            yield hdr, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=DEF_SRC)
    ap.add_argument("--cap", type=int, default=6000)
    ap.add_argument("--out", default=os.path.join(R, "ref_region", "idtaxa_train.fa"))
    a = ap.parse_args()

    seen = set()
    n_unique = 0
    for hdr, _ in records(a.src):
        t = hdr.strip().rstrip(";") + ";"
        if t not in seen:
            seen.add(t)
            n_unique += 1
    step = max(1, n_unique // a.cap)
    print("unique taxonomies: %d -> taking every %d-th (cap %d)" % (n_unique, step, a.cap))

    seen.clear()
    k = written = 0
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with io.open(a.out, "w", encoding="ascii", newline="\n") as out:
        for hdr, seq in records(a.src):
            t = hdr.strip().rstrip(";") + ";"
            if t in seen:
                continue
            seen.add(t)
            if k % step == 0 and written < a.cap:
                out.write(">seq%05d %s\n%s\n" % (written, t, seq.upper()))
                written += 1
            k += 1
            if written >= a.cap:
                break
    print("wrote %s  sequences=%d" % (a.out, written))


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Measure insert (amplicon) length from paired FASTQ by exact-offset overlap.

Convention: revcomp(R2) aligns to R1 such that rc(R2)[j] == R1[j + k].
  k >= 0  -> no read-through; overlap = len(R1) - k;  insert L = k + len(R2)
  k <  0  -> read-through (read longer than insert); overlap = len(R1); L = k + len(R2) < len(R1)
Reported columns: k_mode, ov_bp (= 2*R - L at that k), amp_len (= L), read_through flag.

Usage: python measure_amp.py <raw_root> <p1,p2,...> [--pairs 200] [--samples 2] [--out f.csv]
"""
import gzip, glob, os, sys, argparse, statistics, collections

COMP = str.maketrans("ACGTN", "TGCAN")
def rc(s):
    return s.translate(COMP)[::-1]

def read_pairs(f1, f2, n):
    out = []
    try:
        with gzip.open(f1, "rt", errors="replace") as a, gzip.open(f2, "rt", errors="replace") as b:
            i = 0
            while True:
                l1 = a.readline(); l2 = b.readline()
                if not l1 or not l2: break
                if i % 4 == 1:
                    out.append((l1.strip().upper(), l2.strip().upper()))
                    if len(out) >= n: break
                i += 1
    except Exception as e:
        sys.stderr.write("read error %s: %s\n" % (f1, e))
    return out

def read_single(f, n):
    out = []
    with gzip.open(f, "rt", errors="replace") as fh:
        i = 0
        while True:
            l = fh.readline()
            if not l: break
            if i % 4 == 1:
                out.append(l.strip().upper())
                if len(out) >= n: break
            i += 1
    return out

def best_offset(a, b, min_ov=15):
    """Return (k, ov, mm) for the best placement of b inside a; k may be negative."""
    La, Lb = len(a), len(b)
    best = None
    for k in range(-(Lb - min_ov), La - min_ov + 1):
        lo = max(0, -k); hi = min(Lb - 1, La - 1 - k)
        ov = hi - lo + 1
        if ov < min_ov: continue
        max_mm = max(1, ov // 50)
        mm = 0
        for j in range(lo, hi + 1):
            if a[j + k] != b[j]:
                mm += 1
                if mm > max_mm: break
        if mm <= max_mm:
            cand = (mm, -ov, k)
            if best is None or cand < best:
                best = cand
    if best is None: return None
    mm, negov, k = best
    return (k, -negov, mm)

def measure_project(raw_root, proj, pairs, n_samples):
    d = os.path.join(raw_root, proj)
    if not os.path.isdir(d):
        return dict(project=proj, status="NO_DIR", n_files=0, read_len=None, amp_len=None,
                    k_mode=None, ov_bp=None, read_through=None, amp_dist="", n_pairs=0, n_aligned=0)
    f1 = sorted(glob.glob(os.path.join(d, "*_1.fastq.gz")))
    f2 = sorted(glob.glob(os.path.join(d, "*_2.fastq.gz")))
    if f1 and f2:
        lens = []; rlens = []; ks = []; nal = 0; npr = 0
        for i in range(0, min(n_samples, len(f1))):
            pr = read_pairs(f1[i], f2[i], pairs)
            npr += len(pr)
            for a, b in pr:
                rlens.append(len(a))
                r = best_offset(a, rc(b), 15)
                if r:
                    k, ov, mm = r
                    nal += 1
                    ks.append(k)
                    lens.append(k + len(b))
        if not lens:
            return dict(project=proj, status="NO_ALIGN", n_files=len(f1), read_len=statistics.median(rlens),
                        amp_len=None, k_mode=None, ov_bp=None, read_through=None, amp_dist="",
                        n_pairs=npr, n_aligned=0)
        c = collections.Counter(lens)
        mode = c.most_common(1)[0][0]
        kmode = collections.Counter(ks).most_common(1)[0][0]
        R = statistics.median(rlens)
        return dict(project=proj, status="PAIRED_OK", n_files=len(f1), read_len=R,
                    amp_len=mode, k_mode=kmode, ov_bp=(min(mode, R) - max(kmode, 0)),
                    read_through=("YES" if kmode < 0 else "no"),
                    amp_dist="; ".join("%d(%d)" % (x, n) for x, n in c.most_common(5)),
                    n_pairs=npr, n_aligned=nal)
    fs = sorted(glob.glob(os.path.join(d, "*.fastq.gz")))
    if fs:
        seqs = read_single(fs[0], pairs)
        return dict(project=proj, status="SINGLE_END", n_files=len(fs),
                    read_len=(statistics.median([len(s) for s in seqs]) if seqs else None),
                    amp_len=None, k_mode=None, ov_bp=None, read_through=None, amp_dist="",
                    n_pairs=len(seqs), n_aligned=0)
    return dict(project=proj, status="NO_FASTQ", n_files=0, read_len=None, amp_len=None, k_mode=None,
                ov_bp=None, read_through=None, amp_dist="", n_pairs=0, n_aligned=0)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_root"); ap.add_argument("projects")
    ap.add_argument("--pairs", type=int, default=200)
    ap.add_argument("--samples", type=int, default=2)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    cols = ["project", "status", "n_files", "read_len", "amp_len", "k_mode", "ov_bp",
            "read_through", "amp_dist", "n_pairs", "n_aligned"]
    rows = []
    for p in [x.strip() for x in a.projects.split(",") if x.strip()]:
        r = measure_project(a.raw_root, p, a.pairs, a.samples)
        rows.append(r)
        print("%-14s %-10s R=%-5s L=%-5s k=%-5s ov=%-5s RT=%-4s aligned=%s/%s %s" % (
            r["project"], r["status"], r["read_len"], r["amp_len"], r["k_mode"], r["ov_bp"],
            r["read_through"], r["n_aligned"], r["n_pairs"], r["amp_dist"]))
    if a.out:
        with open(a.out, "w", encoding="utf-8", newline="") as fh:
            fh.write(",".join(cols) + "\n")
            for r in rows:
                fh.write(",".join('"%s"' % ("" if r.get(c) is None else r.get(c)) for c in cols) + "\n")
        print("wrote", a.out)

if __name__ == "__main__":
    main()


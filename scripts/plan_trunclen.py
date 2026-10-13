# -*- coding: utf-8 -*-
"""Plan per-project truncLen under strategy A (overlap must sit in good-quality cycles).

Robustness rules learned from batch7:
  * read length varies BETWEEN samples of one project -> cap truncLen at the MINIMUM read length.
  * the offset scan can return nonsense on noisy data -> require (a) enough aligned pairs,
    (b) insert length inside a plausible window; otherwise mark unreliable and keep truncLen 0.
  * if 2*min(R,L) - L < min_overlap the 50 bp target is unreachable -> keep reads untruncated
    (max overlap is what it is) and say so in the note.

Usage: python plan_trunclen.py <raw_root> <p1,p2,...> --out plan.csv [--pairs 120] [--qreads 25000]
"""
import gzip, glob, os, sys, argparse, statistics, collections, io, csv, math

COMP = str.maketrans("ACGTN", "TGCAN")
def rc(s):
    return s.translate(COMP)[::-1]

def read_pairs(f1, f2, n):
    out = []
    with gzip.open(f1, "rt", errors="replace") as a, gzip.open(f2, "rt", errors="replace") as b:
        i = 0
        while True:
            l1 = a.readline(); l2 = b.readline()
            if not l1 or not l2: break
            if i % 4 == 1:
                out.append((l1.strip().upper(), l2.strip().upper()))
                if len(out) >= n: break
            i += 1
    return out

def read_lengths(path, n=200):
    lens = []
    with gzip.open(path, "rt", errors="replace") as fh:
        i = 0
        while True:
            l = fh.readline()
            if not l: break
            if i % 4 == 1:
                lens.append(len(l.strip()))
                if len(lens) >= n: break
            i += 1
    return lens

def best_offset(a, b, min_ov=15, tol_div=12):
    La, Lb = len(a), len(b)
    best = None   # ((ov, -mm), k, ov, mm)
    for k in range(-(Lb - min_ov), La - min_ov + 1):
        lo = max(0, -k); hi = min(Lb - 1, La - 1 - k)
        ov = hi - lo + 1
        if ov < min_ov: continue
        max_mm = max(2, ov // tol_div)
        mm = 0
        for j in range(lo, hi + 1):
            if a[j + k] != b[j]:
                mm += 1
                if mm > max_mm: break
        if mm <= max_mm:
            # among valid placements prefer the LARGEST overlap (most constrained = the true one),
            # then the fewest mismatches. The old lexicographic (mm first) order could pick a
            # shifted pseudo-alignment when insert length ~= read length (batch7 PRJNA701534: L=35).
            cand = (ov, -mm)
            if best is None or cand > best[0]:
                best = (cand, k, ov, mm)
    if best is None: return None
    _, k, ov, mm = best
    return (k, ov, mm)

def qprofile(path, nmax):
    sums = None; n = 0
    with gzip.open(path, "rt", errors="replace") as fh:
        i = 0
        while True:
            l = fh.readline()
            if not l: break
            if i % 4 == 3:
                q = l.strip()
                if sums is None: sums = [0] * len(q)
                for j, ch in enumerate(q):
                    if j < len(sums): sums[j] += ord(ch) - 33
                n += 1
                if n >= nmax: break
            i += 1
    return ([s / n for s in sums] if n else []), n

def t25(prof, win=25, qthr=25.0):
    if not prof: return None
    for start in range(0, len(prof), win):
        if statistics.mean(prof[start:start + win]) < qthr:
            return start
    return len(prof)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_root"); ap.add_argument("projects")
    ap.add_argument("--pairs", type=int, default=120)
    ap.add_argument("--samples", type=int, default=3)
    ap.add_argument("--rlen-samples", type=int, default=30)  # 扫描全部样本，避免漏掉少数短读长样本组
    ap.add_argument("--qreads", type=int, default=25000)
    ap.add_argument("--qmin", type=float, default=25.0)
    ap.add_argument("--min-overlap", type=int, default=50)
    ap.add_argument("--overlap-target", type=int, default=60)
    ap.add_argument("--asv-root", default="")
    ap.add_argument("--primer-table", default="", help="adjust read/insert lengths for 5-primer trimming")
    ap.add_argument("--tag", default="plan")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    rows = []
    for p in [x.strip() for x in a.projects.split(",") if x.strip()]:
        d = os.path.join(a.raw_root, p)
        f1 = sorted(glob.glob(os.path.join(d, "*_1.fastq.gz")))
        f2 = sorted(glob.glob(os.path.join(d, "*_2.fastq.gz")))
        if not f1 or not f2:
            rows.append(dict(project=p, mode="single", read_len="", read_len_min="", insert_len="", t25_f="", t25_r="",
                             tf="", tr="", overlap="", n_aligned="", note="single-end: no merge, truncLen not planned"))
            print("%-14s single-end -> no truncLen planning" % p); continue
        # read lengths across several samples (project-internal spread matters!)
        med_per_sample = []
        rl_all = []
        for f in f1[:a.rlen_samples]:
            rl = read_lengths(f, 200)
            if rl:
                med_per_sample.append(statistics.median(rl)); rl_all += rl
        R_med = int(statistics.median(rl_all)) if rl_all else None
        pF = pR = 0
        if a.primer_table and os.path.exists(a.primer_table):
            for pr in list(csv.DictReader(io.open(a.primer_table, encoding="utf-8-sig"))):
                if pr.get("project") == p:
                    pF = len((pr.get("primer_f") or "").strip())
                    pR = len((pr.get("primer_r") or "").strip())
                    break
        # truncLen applies to every sample of the project, so cap at the smallest per-sample median;
        # isolated very short reads (< minLen=100) are dropped by the pipeline anyway.
        # truncLen applies to the whole project; use the 20th percentile of per-sample medians so that
        # one odd short-read sample (<=20% of samples) is dropped instead of crippling all samples.
        if med_per_sample:
            ms = sorted(med_per_sample)
            R_min = int(ms[max(0, int(0.20 * (len(ms) - 1)))]) if len(ms) >= 5 else int(ms[0])
        else:
            R_min = None
        # insert length from 2-3 samples
        lens = []; npr = 0
        for i in range(min(a.samples, len(f1))):
            for x, y in read_pairs(f1[i], f2[i], a.pairs):
                npr += 1
                r = best_offset(x, rc(y))
                if r:
                    k, ov, mm = r
                    lens.append(k + len(y))
        cnt = collections.Counter(lens)
        L = cnt.most_common(1)[0][0] if lens else None
        frac = (sum(cnt.values()) / npr) if npr else 0.0
        L_src = "raw_reads"
        # fallback / cross-check: ASV length mode from a finished run (AGENTS 10.4 step 2)
        if a.asv_root:
            t = os.path.join(a.asv_root, p, "asv_table.csv")
            if os.path.exists(t):
                c2 = collections.Counter()
                with io.open(t, encoding="utf-8-sig") as fh:
                    rr = csv.reader(fh); next(rr, None)
                    for row in rr:
                        if row: c2[len(row[0].strip().strip('"'))] += 1
                if c2:
                    L2 = c2.most_common(1)[0][0]
                    if (L is None) or (frac < 0.30) or (L < 100) or (L > 2 * (R_med or 0)):
                        L = L2; L_src = "asv_table"
                    elif abs(L2 - L) > 20:
                        # 已有跑完的 ASV 表时，行名长度主峰是权威插入长度（AGENTS §10.4）：
                        # 原始读段扫描在插入长度接近读长时会给出偏小的伪解，该偏小值会把真实重叠压成负数（合并全灭）。
                        L = L2; L_src = "asv_table(raw_reads=%d)" % L
        qf, _ = qprofile(f1[0], a.qreads)
        qr, _ = qprofile(f2[0], a.qreads)
        tF = t25(qf, qthr=a.qmin); tR = t25(qr, qthr=a.qmin)
        unreliable = (L is None) or (L < 100) or (L > 2 * (R_med or 0))
        if str(L_src).startswith("raw") and frac < 0.30:
            unreliable = True   # raw-read estimate without enough aligned pairs; ASV-based L is independent
        if unreliable:
            tf = tr = "0"; ov = ""; note = ("insert length unreliable (L=%s, aligned=%.0f%% of %d pairs)"
                                             % (L, 100 * frac, npr))
        else:
            # cutadapt runs BEFORE filterAndTrim (which drops reads shorter than truncLen), so plan in
            # post-trim coordinates; if the probe ends up NOT trimming, a smaller truncLen is still safe.
            R_eff_f, R_eff_r, L_eff = R_min - pF, R_min - pR, L - pF - pR
            cap = min(R_eff_f, R_eff_r, L_eff)
            # mergePairs(maxMismatch = 0) 的算术：重叠越长，出现至少一个错配的概率越大 -> 合并越少。
            # 因此目标不是"尽量保留长度"，而是"用尽量少的循环凑到刚够的重叠（默认 60bp），
            # 且不越过该读的质量拐点 T25"。错配预算比信息量更值钱（PRJDB8720 的 275bp 重叠合并率 0.1%）。
            qf = min(cap, (tF - pF)) if (tF and tF > pF) else cap
            qr = min(cap, (tR - pR)) if (tR and tR > pR) else cap
            L = L_eff
            ov_target = float(a.overlap_target)
            need = L + ov_target
            tf = int(min(qf, math.ceil(need / 2.0)))
            tr = int(min(qr, max(1, need - tf)))
            if tf + tr < need:
                tf = int(min(qf, need - tr))
            if tf + tr < need:
                tr = int(min(qr, need - tf))
            ov = tf + tr - L
            if abs(L - R_med) < 40 and str(L_src).startswith("raw"):
                # 插入长度≈读长时，偏移量搜索存在多解（伪解给出超大重叠），先保守不截断，交给 A/B 验证。
                tf = tr = 0
                ov = ""
                note = "insert≈read length (ambiguous alignment); keep untruncated, queued for A/B"
                rows.append(dict(project=p, mode="paired", read_len=R_med, read_len_min=R_min, insert_len=L,
                                 t25_f=tF, t25_r=tR, tf=tf, tr=tr, overlap=ov,
                                 n_aligned="%d/%d (%s)" % (len(lens), npr, L_src), note=note))
                print("%-14s R=%-4s L=%-5s primer=(%d,%d) -> keep untruncated (L≈R, ambiguous)" % (p, R_med, L, pF, pR))
                continue
            if 2 * cap - L < a.min_overlap:
                # truncating cannot reach the 50 bp target; keep reads as they are (0 = no truncation).
                # NOTE: cap itself would DISCARD reads (verified: filterAndTrim drops reads shorter than
                # truncLen), so 0 is the only safe fallback here.
                tf = tr = 0
                ov = ""
                note = "max overlap %d bp < %d bp (insert long); keep untruncated" % (2 * cap - L, a.min_overlap)
            else:
                note = "ok"
        rows.append(dict(project=p, mode="paired", read_len=R_med, read_len_min=R_min, insert_len=L,
                         t25_f=tF, t25_r=tR, tf=tf, tr=tr, overlap=ov,
                         n_aligned="%d/%d (%s)" % (len(lens), npr, L_src), note=note))
        print("%-14s R=%-4s Rmin=%-4s L=%-5s primer=(%d,%d) T25=(%s,%s) truncLen=(%s,%s) ov=%-4s %s" % (
            p, R_med, R_min, L, pF, pR, tF, tR, tf, tr, ov, note))
    if a.out:
        with io.open(a.out, "w", encoding="utf-8-sig", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["project", "mode", "read_len", "read_len_min", "insert_len",
                                               "t25_f", "t25_r", "tf", "tr", "overlap", "n_aligned", "note"])
            w.writeheader(); w.writerows(rows)
        print("wrote", a.out)

if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
# 对原始 FASTQ 做统一下限深度抽样（cap），双端保持配对
# 用法: subsample_fastq.py <src_dir> <dst_dir> <max_reads> [seed] [project1,project2,...]
# 加固（2026-10-08）：单个文件损坏（gzip CRC / 截断）不再中断整批——记录该样本并继续；
# 任何 FATAL 跳过都会写入 <dst_dir>/_subsample_report.csv 且以退出码 1 结束（链条据此判失败）。
# WARN（样本缺 _1/_2 一侧、非 *_1/*_2 命名的文件）只记录，不影响退出码。
import gzip, os, sys, glob, random, shutil, time, csv, datetime

src, dst, cap = sys.argv[1], sys.argv[2], int(sys.argv[3])
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 20261006
only = sys.argv[5].split(",") if len(sys.argv) > 5 and sys.argv[5].strip() else None

def count_reads(path):
    n = 0
    with gzip.open(path, "rb") as f:
        while True:
            buf = f.read(1 << 23)
            if not buf: break
            n += buf.count(b"\n")
    return n // 4

def write_subset(srcp, dstp, keep):
    with gzip.open(srcp, "rb") as fi, gzip.open(dstp, "wb", compresslevel=6) as fo:
        i = 0
        while True:
            l1 = fi.readline()
            if not l1: break
            l2 = fi.readline(); l3 = fi.readline(); l4 = fi.readline()
            if i in keep:
                fo.write(l1); fo.write(l2); fo.write(l3); fo.write(l4)
            i += 1
    return i

t0 = time.time()
tot_before = tot_after = 0
n_sampled = 0
skipped = []          # (kind, project, sample, reason)
for proj in sorted(os.listdir(src)):
    sd = os.path.join(src, proj)
    if not os.path.isdir(sd): continue
    if only is not None and proj not in only: continue
    dd = os.path.join(dst, proj); os.makedirs(dd, exist_ok=True)
    files = sorted(glob.glob(os.path.join(sd, "*.fastq.gz")))
    if not files: continue
    paired = any(f.endswith("_1.fastq.gz") for f in files)
    if paired:
        samples = {}
        for f in files:
            b = os.path.basename(f)
            if b.endswith("_1.fastq.gz"): samples.setdefault(b[:-len("_1.fastq.gz")], {})["1"] = f
            elif b.endswith("_2.fastq.gz"): samples.setdefault(b[:-len("_2.fastq.gz")], {})["2"] = f
    else:
        samples = {os.path.basename(f)[:-len(".fastq.gz")]: {"s": f} for f in files}

    for s, d in sorted(samples.items()):
        written = []
        try:
            rng = random.Random(seed)
            if paired:
                miss = [k for k in ("1", "2") if k not in d]
                if miss:
                    skipped.append(("WARN", proj, s, "missing side " + ",".join(miss)))
                    continue
                n = count_reads(d["1"]); tot_before += n
                dst1 = os.path.join(dd, s + "_1.fastq.gz"); dst2 = os.path.join(dd, s + "_2.fastq.gz")
                written = [dst1, dst2]
                if n <= cap:
                    shutil.copyfile(d["1"], dst1); shutil.copyfile(d["2"], dst2); m = n
                else:
                    keep = set(rng.sample(range(n), cap))
                    write_subset(d["1"], dst1, keep); write_subset(d["2"], dst2, keep); m = cap
                    n_sampled += 1
                tot_after += m
            else:
                n = count_reads(d["s"]); tot_before += n
                dst1 = os.path.join(dd, s + ".fastq.gz"); written = [dst1]
                if n <= cap:
                    shutil.copyfile(d["s"], dst1); m = n
                else:
                    keep = set(rng.sample(range(n), cap)); write_subset(d["s"], dst1, keep); m = cap
                    n_sampled += 1
                tot_after += m
        except Exception as e:
            skipped.append(("FATAL", proj, s, "%s: %s" % (type(e).__name__, e)))
            for f in written:
                try:
                    if os.path.exists(f): os.remove(f)
                except Exception:
                    pass

print("cap per sample : %s" % format(cap, ","))
print("samples sampled: %d" % n_sampled)
print("reads before   : %s" % format(tot_before, ","))
if tot_before:
    print("reads after    : %s  (%.1f%%)" % (format(tot_after, ","), 100.0 * tot_after / tot_before))
print("elapsed min    : %.1f" % ((time.time() - t0) / 60.0))
print("skipped        : %d (fatal %d, warn %d)" % (len(skipped),
      len([s for s in skipped if s[0] == "FATAL"]), len([s for s in skipped if s[0] == "WARN"])))
for s in skipped[:60]:
    print("   %s %-14s %-28s %s" % s)
if skipped:
    rep = os.path.join(dst, "_subsample_report.csv")
    with open(rep, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "project", "sample", "reason", "time"])
        for s in skipped:
            w.writerow([s[0], s[1], s[2], s[3], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")])
    print("wrote", rep)
if any(s[0] == "FATAL" for s in skipped):
    sys.exit(1)

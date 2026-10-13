import csv, os, gzip, glob, collections
raw = r"D:\CodexWork\mic-re\data\raw"
man = list(csv.DictReader(open(r"D:\CodexWork\mic-re\meta\w1_manifest_raw.csv", encoding="utf-8-sig")))

# expected bytes per run (manifest bytes = sum of all files of that run)
exp = {}
for r in man:
    if r["fastq"]:
        exp[r["run"]] = float(r["bytes"] or 0)

bad = []
ok = []
for proj in sorted(os.listdir(raw)):
    d = os.path.join(raw, proj)
    if not os.path.isdir(d): continue
    for f in sorted(glob.glob(os.path.join(d, "*.fastq.gz"))):
        base = os.path.basename(f)
        run = base.split("_")[0].replace(".fastq.gz", "")
        # gzip integrity: read to EOF
        integ = True
        try:
            with gzip.open(f, "rb") as fh:
                while fh.read(1 << 20):
                    pass
        except Exception as e:
            integ = False
        if not integ:
            bad.append((proj, base, "gzip-truncated", os.path.getsize(f)))
        else:
            ok.append((proj, base, os.path.getsize(f)))

# pair count check for paired runs
cnt = collections.defaultdict(list)
for proj, base, sz in ok:
    run = base.split("_")[0]
    cnt[(proj, run)].append(base)

print("OK files   :", len(ok))
print("BAD files  :", len(bad))
for b in bad:
    print("   ", b)

# per-run total vs expected
print("\n--- per-run size check ---")
mism = []
for (proj, run), files in sorted(cnt.items()):
    actual = sum(os.path.getsize(os.path.join(raw, proj, f)) for f in files)
    e = exp.get(run)
    if e and abs(actual - e) > 200:
        mism.append((proj, run, actual, e))
for m in mism[:20]:
    print("   MISMATCH {0}/{1}: actual={2:,} expected={3:,.0f}".format(*m))
print("run-size mismatches:", len(mism))

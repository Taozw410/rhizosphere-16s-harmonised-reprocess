import gzip, glob, os, collections
raw = r"D:\CodexWork\mic-re\data\raw"
for proj in sorted(os.listdir(raw)):
    d = os.path.join(raw, proj)
    if not os.path.isdir(d): continue
    fs = sorted(glob.glob(os.path.join(d, "*.fastq.gz")))
    n = len(fs); tot = sum(os.path.getsize(f) for f in fs)
    print(f"\n=== {proj}  files={n}  {tot/1e6:.1f} MB ===")
    for f in fs[:2]:
        with gzip.open(f, "rt", errors="replace") as fh:
            lens = []
            seqs = []
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    lens.append(len(line.strip()))
                    if len(seqs) < 3: seqs.append(line.strip())
                if i > 4000: break
        print(f"  {os.path.basename(f)}")
        print(f"    read length: min={min(lens)} max={max(lens)} mean={sum(lens)/len(lens):.0f}  n={len(lens)}")
        for s in seqs:
            print("    ", s[:60])
        # primer probe
        p515 = "GTGYCAGCMGCCGCGGTAA"; p515b = "GTGCCAGCMGCCGCGGTAA"
        p806 = "GGACTACNVGGGTWTCTAAT"; p806b = "GGACTACHVGGGTWTCTAAT"
        def match(seq, pat):
            if len(seq) < len(pat): return False
            s = seq.upper(); p = pat.upper()
            return all(pc in "N" or sc == pc for sc, pc in zip(s, p))
        with gzip.open(f, "rt", errors="replace") as fh:
            c515 = c806 = n0 = 0
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    s = line.strip().upper(); n0 += 1
                    if match(s, p515) or match(s, p515b): c515 += 1
                    if match(s, p806) or match(s, p806b): c806 += 1
                if i > 8000: break
        print(f"    primer at read start: 515F={c515}/{n0}  806R={c806}/{n0}")

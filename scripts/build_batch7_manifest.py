# -*- coding: utf-8 -*-
"""Build the combined download manifest for the new batch (Illumina only)."""
import csv, io, os
from collections import Counter
R = r"D:\CodexWork\mic-re"
out = os.path.join(R, "meta", "batch7_download_manifest.csv")
rows = []
for f in ("batch6_manifest.csv", "batch7_maize_manifest.csv"):
    p = os.path.join(R, "meta", f)
    if not os.path.exists(p):
        continue
    for r in csv.DictReader(io.open(p, encoding="utf-8-sig")):
        if not (r.get("fastq") or "").strip():
            continue
        if (r.get("platform") or "").strip().upper() != "ILLUMINA":
            continue
        rows.append(dict(project=r["project"].strip(), crop=(r.get("crop") or "").strip(),
                         run=r["run"].strip(), platform="ILLUMINA",
                         sample=(r.get("sample") or "").strip(), reads=(r.get("reads") or "").strip(),
                         fastq=r["fastq"].strip()))
with io.open(out, "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["project", "crop", "run", "platform", "sample", "reads", "fastq"])
    w.writeheader(); w.writerows(rows)
print("rows=%d  projects=%d  crops=%s" % (len(rows), len(set(r["project"] for r in rows)),
      dict(Counter(r["crop"] for r in rows))))
print("->", out)

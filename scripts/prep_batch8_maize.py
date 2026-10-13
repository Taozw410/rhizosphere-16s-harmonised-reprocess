# -*- coding: utf-8 -*-
"""Batch 8 prep: pick fresh maize candidates, fetch their ENA run manifests, for probe screening."""
import csv, io, os, re, time
import urllib.request, urllib.parse
from collections import Counter

R = r"D:\CodexWork\mic-re"
FIELDS = "run_accession,sample_accession,fastq_ftp,fastq_bytes,library_strategy,instrument_platform,read_count"
API = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def rows(path):
    try:
        return list(csv.DictReader(io.open(path, encoding="utf-8-sig")))
    except Exception:
        return []


used = set()
for r in rows(os.path.join(R, "out50k", "master_inventory.csv")):
    if r.get("project"):
        used.add(r["project"].strip())
for d in ["out50k", "out_stage2", "out"]:
    p = os.path.join(R, d)
    if os.path.isdir(p):
        for n in os.listdir(p):
            if re.match(r"^(PRJ|DRR|ERR|SRR)", n):
                used.add(n)
for d in ["raw", "raw50k", "raw50k_longread", "raw50k_w1done", "raw50k_excluded", "raw50k_excluded_its",
          "probe", "probe2", "probe3", "probe4", "probe5", "probe7"]:
    p = os.path.join(R, "data", d)
    if os.path.isdir(p):
        for n in os.listdir(p):
            if os.path.isdir(os.path.join(p, n)):
                used.add(n)
rec = os.path.join(R, "_recycle")
if os.path.isdir(rec):
    for n in os.listdir(rec):
        for tok in re.findall(r"PRJ[A-Z]*\d+", n):
            used.add(tok)
for f in ("rescreen_strict.csv", "marker_gene_qc.csv"):
    p = os.path.join(R, "out50k", f)
    for r in rows(p):
        if (r.get("verdict") or "").startswith("EXCLUDE"):
            used.add(r["project"])
for r in rows(os.path.join(R, "meta", "batch7_selected.txt")):
    pass
sel7 = [x.strip() for x in io.open(os.path.join(R, "meta", "batch7_selected.txt"), encoding="utf-8-sig").read().split() if x.strip()]
used |= set(sel7)

BAD = ("fungi", "fungal", "buckwheat", "metagenom", "metatranscriptom", "third-generation",
       "nanopore", "pacbio", "pepper", "soybean", "tomato", "potato", "peanut", "cotton",
       "barley", "banana", "sugarcane", "citrus", "apple", "tea ")


def score(r):
    s = 0
    n = int(r.get("ill_amp") or 0)
    mr = int(r.get("med_reads") or 0)
    if r.get("funcmarker") == "0":
        s += 40
    if r.get("nonplant") == "0":
        s += 20
    if 20 <= n <= 200:
        s += 30
    elif 200 < n <= 400:
        s += 22
    elif n > 400:
        s += 10
    elif n >= 15:
        s += 12
    if mr >= 50000:
        s += 20
    elif mr >= 20000:
        s += 12
    t = (r.get("title") or "").lower()
    if "rhizosphere" in t or "root" in t:
        s += 15
    if "16s" in t or "bacteri" in t or "microbi" in t:
        s += 15
    return s


elig = [r for r in rows(os.path.join(R, "meta", "eligible.csv")) if r.get("crop") == "maize"]
fresh = [r for r in elig if r["study"] not in used
         and not any(w in (r.get("title") or "").lower() for w in BAD)]
fresh.sort(key=lambda r: (-score(r), r["study"]))
take = fresh[:30]
print("fresh maize =", len(fresh), " probing =", len(take))


def fetch(study):
    q = urllib.parse.urlencode({"accession": study, "result": "read_run", "fields": FIELDS, "format": "tsv"})
    req = urllib.request.Request(API + "?" + q, headers={"User-Agent": "Mozilla/5.0 (research; maize-rhizosphere-16S)"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read().decode("utf-8", "replace")


out_rows = []
for r in take:
    study = r["study"]
    txt = None
    for _ in range(3):
        try:
            txt = fetch(study)
            break
        except Exception as e:
            print("  retry", study, e)
            time.sleep(3)
    if not txt:
        print("  FAIL", study)
        continue
    lines = [l for l in txt.strip().split("\n") if l.strip()]
    hdr = lines[0].split("\t")
    n = 0
    for l in lines[1:]:
        d = dict(zip(hdr, l.split("\t")))
        if (d.get("instrument_platform") or "").upper() != "ILLUMINA":
            continue
        if not (d.get("fastq_ftp") or "").strip():
            continue
        out_rows.append(dict(project=study, crop="maize", run=d["run_accession"], platform="ILLUMINA",
                             sample=d.get("sample_accession", ""), reads=d.get("read_count", ""),
                             fastq=d["fastq_ftp"].strip()))
        n += 1
    print("  %-14s runs=%-4d n=%-5s reads=%s" % (study, n, r.get("ill_amp"), r.get("med_reads")))

p = os.path.join(R, "meta", "batch8_maize_manifest.csv")
with io.open(p, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["project", "crop", "run", "platform", "sample", "reads", "fastq"])
    w.writeheader(); w.writerows(out_rows)
with io.open(os.path.join(R, "meta", "batch8_candidates.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["crop", "study", "ill_amp", "med_reads", "title"])
    for r in take:
        w.writerow(["maize", r["study"], r.get("ill_amp"), r.get("med_reads"), (r.get("title") or "")[:120]])
print("manifest rows =", len(out_rows), "projects =", len(set(x["project"] for x in out_rows)), "->", p)

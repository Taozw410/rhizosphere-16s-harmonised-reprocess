# -*- coding: utf-8 -*-
"""Fetch ENA run-level manifest (fastq URLs) for a list of study accessions."""
import csv, io, os, sys, time
import urllib.request
import urllib.parse

R = r"D:\CodexWork\mic-re"
FIELDS = "run_accession,sample_accession,fastq_ftp,fastq_bytes,library_strategy,instrument_platform,read_count"
API = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def fetch(study):
    q = urllib.parse.urlencode({"accession": study, "result": "read_run", "fields": FIELDS, "format": "tsv"})
    req = urllib.request.Request(API + "?" + q, headers={"User-Agent": "Mozilla/5.0 (research; rhizosphere-16S)"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read().decode("utf-8", "replace")


def main():
    sel = list(csv.DictReader(io.open(os.path.join(R, "meta", "batch7_selection.csv"), encoding="utf-8-sig")))
    maize = [r for r in sel if r["crop"] == "maize"]
    rows = []
    for r in maize:
        study = r["study"]
        for attempt in range(3):
            try:
                txt = fetch(study)
                break
            except Exception as e:
                print("  retry %s: %s" % (study, e))
                time.sleep(3)
        else:
            print("FAIL %s" % study)
            continue
        lines = [l for l in txt.strip().split("\n") if l.strip()]
        hdr = lines[0].split("\t")
        n = 0
        for l in lines[1:]:
            f = l.split("\t")
            d = dict(zip(hdr, f))
            ftp = (d.get("fastq_ftp") or "").strip()
            if not ftp:
                continue
            rows.append(dict(project=study, crop="maize", run=d["run_accession"],
                             platform=d.get("instrument_platform", ""), sample=d.get("sample_accession", ""),
                             reads=d.get("read_count", ""), strat=d.get("library_strategy", ""),
                             fastq=ftp))
            n += 1
        print("%-14s runs=%d" % (study, n))

    out = os.path.join(R, "meta", "batch7_maize_manifest.csv")
    with io.open(out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["project", "crop", "run", "platform", "sample", "reads", "strat", "fastq"])
        w.writeheader()
        w.writerows(rows)
    print("total runs =", len(rows), "->", out)


main()

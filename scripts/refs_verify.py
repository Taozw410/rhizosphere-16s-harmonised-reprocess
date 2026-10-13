"""Verify curated reference DOIs against CrossRef metadata.

No citation enters manuscript/refs/refs.md unless CrossRef returns metadata whose
title matches the expected keyword and whose year/authors are consistent.
Writes manuscript/refs/verified_refs.json + a markdown table for eyeballing.

Run: D:\\bioenv\\envs\\riz16s\\python.exe scripts\\refs_verify.py
"""

import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

MAILTO = "mic-re-refs@example.org"
OUTDIR = os.path.join("manuscript", "refs")

# (ref_id, doi, expected-title keyword(s) that must appear (case-insensitive, all of them))
CANDIDATES = [
    ("R01", "10.1073/pnas.1302837110", ["maize rhizosphere", "heritability"]),
    ("R02", "10.1038/s41396-021-00923-z", ["chronosequence", "breeding"]),
    ("R03", "10.1093/femsec/fiz040", ["rice roots", "rhizosphere"]),
    ("R04", "10.1111/nph.16890", ["host selection", "crop microbiome"]),
    ("R05", "10.1016/j.tplants.2020.04.015", ["macroecological", "rhizosphere"]),
    ("R06", "10.1111/geb.13266", ["study", "distance", "decay"]),
    ("R07", "10.1093/gigascience/gix029", ["quality", "sequence read archive"]),
    ("R08", "10.3389/fcimb.2021.720637", ["reprocessing", "16s"]),
    ("R09", "10.1016/j.tim.2018.11.003", ["contamination", "low microbial biomass"]),
    ("R10", "10.1186/s40168-018-0605-2", ["contaminant sequences", "marker-gene"]),
    ("R11", "10.1038/nmeth.3869", ["dada2", "amplicon"]),
    ("R12", "10.1186/s40793-025-00705-6", ["clustering and denoising"]),
    ("R13", "10.1111/1462-2920.16566", ["trimming", "read lengths"]),
    ("R14", "10.3389/fmicb.2020.550420", ["biased", "amplicon sequencing pipeline"]),
    ("R15", "10.3389/fmicb.2015.00771", ["primer and platform effects"]),
    ("R16", "10.1128/msphere.01202-20", ["primer", "pipelines", "parameters"]),
    ("R17", "10.1038/s41587-019-0209-9", ["qiime 2"]),
    ("R18", "10.1186/s40168-018-0521-5", ["idtaxa"]),
    ("R19", "10.1093/gbe/evae102", ["taxonomic classification"]),
    ("R20", "10.1016/s0168-6496(99)00083-5", ["acetate", "nitrate", "ferric iron", "sulfate"]),
    ("R21", "10.1038/ismej.2009.100", ["iron-reducing", "paddy"]),
    ("R22", "10.1128/aem.67.8.3586-3597.2001", ["methane oxidation", "rice rhizosphere"]),
    ("R23", "10.1016/j.soilbio.2018.05.025", ["methanotrophs", "methane oxidation", "taiwan"]),
    ("R24", "10.5194/bg-11-3353-2014", ["methane- and ammonia-oxidizing"]),
    ("R25", "10.1016/0168-6496(95)00066-6", ["sulfate", "paddy soil"]),
    ("R26", "10.3390/app11020688", ["faprotax", "soil bacteria"]),
    ("R27", "10.1186/s40793-020-00358-7", ["tax4fun2"]),
    ("R28", "10.1111/gcb.17581", ["alternate wetting and drying", "methane"]),
    ("R29", "10.1016/j.ejsobi.2012.10.008", ["upland", "irrigated rice"]),
    ("R30", "10.3390/microorganisms8091329", ["upland rice", "drought"]),
    ("R31", "10.1111/j.1466-8238.2011.00756.x", ["species replacement", "nestedness"]),
    ("R32", "10.1186/1471-2105-7-371", ["unifrac"]),
    ("R33", "10.1128/msystems.01066-22", ["core microbiome assignments"]),
    ("R34", "10.1186/s40168-017-0237-y", ["normalization", "differential abundance"]),
    ("R35", "10.1186/s13059-021-02306-1", ["siamcat"]),
    ("R36", "10.1128/msphere.00354-23", ["rarefaction"]),
    ("R37", "10.1186/s40168-018-0420-9", ["16s rrna gene copy number"]),
    ("R38", "10.1007/s11274-024-04000-3", ["methylotrophic bacteria", "paddy"]),
]


def norm(text: str) -> str:
    return re.sub(r"[^a-z0-9 ]+", " ", (text or "").lower())


def fetch(doi: str) -> dict | None:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi) + "?mailto=" + MAILTO
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "mic-re-refs/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode("utf-8")).get("message")
    except urllib.error.HTTPError as exc:
        print(f"  HTTP {exc.code} for {doi}", file=sys.stderr)
        return None
    except Exception as exc:  # noqa: BLE001
        print(f"  ERROR {exc!r} for {doi}", file=sys.stderr)
        return None


def summarize(msg: dict) -> dict:
    title = (msg.get("title") or [""])[0]
    journal = (msg.get("container-title") or [""])[0]
    year = None
    for key in ("published-print", "published-online", "issued", "created"):
        parts = (msg.get(key) or {}).get("date-parts") or []
        if parts and parts[0] and parts[0][0]:
            year = parts[0][0]
            break
    structured, authors = [], []
    for a in msg.get("author") or []:
        fam = a.get("family") or ""
        given = a.get("given") or ""
        structured.append({"family": fam, "given": given})
        if fam:
            authors.append(f"{fam} {given}".strip())
    return {
        "doi": msg.get("DOI", ""),
        "title": title,
        "journal": journal,
        "year": year,
        "volume": msg.get("volume", ""),
        "issue": msg.get("issue", ""),
        "page": msg.get("page", ""),
        "type": msg.get("type", ""),
        "n_authors": len(authors),
        "authors": authors,
        "authors_structured": structured,
        "first_author": authors[0] if authors else "",
    }


def main() -> int:
    os.makedirs(OUTDIR, exist_ok=True)
    out, bad = [], []
    for ref_id, doi, keywords in CANDIDATES:
        msg = fetch(doi)
        if not msg:
            bad.append((ref_id, doi, "no metadata returned"))
            print(f"MISS {ref_id} {doi}")
            continue
        info = summarize(msg)
        haystack = norm(info["title"]) + " " + norm(info["journal"])
        missing = [k for k in keywords if norm(k) not in haystack]
        info["ref_id"] = ref_id
        info["expected_keywords"] = keywords
        info["keyword_missing"] = missing
        out.append(info)
        if missing:
            bad.append((ref_id, doi, "missing keyword(s): " + "; ".join(missing)))
            print(f"CHECK {ref_id} {doi} :: {info['title']} :: missing {missing}")
        else:
            print(f"OK {ref_id} {doi} :: {info['title']} ({info['journal']}, {info['year']})")
    with open(os.path.join(OUTDIR, "verified_refs.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    lines = ["# CrossRef verification of curated references", "",
             "| ref | first author | year | journal | title | doi | status |",
             "| --- | --- | --- | --- | --- | --- | --- |"]
    for info in out:
        status = "ok" if not info["keyword_missing"] else "CHECK: " + "; ".join(info["keyword_missing"])
        lines.append(
            f"| {info['ref_id']} | {info['first_author']} | {info['year']} | {info['journal']} | "
            f"{info['title']} | {info['doi']} | {status} |")
    with open(os.path.join(OUTDIR, "_verified_table.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nverified {len(out)} / {len(CANDIDATES)}; problems: {len(bad)}")
    for ref_id, doi, why in bad:
        print(f"  {ref_id} {doi} :: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

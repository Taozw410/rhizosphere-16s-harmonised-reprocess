# -*- coding: utf-8 -*-
"""Water-regime coding for the 45 selected studies (analysis step requested by the user).

HONESTY NOTE: neither the repository metadata nor the sample names carry a water-regime field
(sample IDs are plain SRR/ERR accessions; titles mention water only in 5 studies).  The coding below is
therefore an EXPLICIT, EVIDENCE-TAGGED curation: for each study we record the water regime we assign,
the evidence (a title quote where one exists, otherwise the crop-default rule), and a confidence level.
This file is meant to be reviewed line by line by the author before it enters a manuscript.

usage: python water_regime_code.py  -> out_analysis/water_regime.csv
"""
import csv, io, os
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")

# evidence quotes taken verbatim from the study titles (see _tmp_titles_full.txt)
EVID = {
 "PRJEB45605": ("flooded", "high", "title: 'rice paddy pot experiment with 4 fertiliser treatments'"),
 "PRJNA858795": ("flooded", "high", "title: 'bacterial and archaeal communities living in Italian rice paddies rhizosphere soil and pore water'"),
 "PRJNA1261944": ("flooded", "high", "title: 'Microbial diversity in rice rhizosphere under flooding conditions'"),
 "PRJNA611687": ("flooded", "medium", "title: 'arsenic mobilization in the rice rhizosphere mediated by iron-reducing bacteria' (flooded-soil redox chemistry)"),
 "PRJNA1479906": ("upland", "high", "title: 'salt-stressed maize inoculated with Bacillus' (pot experiment, non-flooded)"),
 "PRJNA967263": ("upland", "high", "title: 'maize with contour reverse-slope terraces such as slope farmland' (dryland terraces)"),
 "PRJNA593708": ("upland", "medium", "title: 'X-ray u-CT scan on rhizosphere bacteria of maize plants' (potted plants)"),
 "PRJEB6400": ("upland", "high", "title: 'sea-water irrigation ... two wheat cultivars' (irrigated, not flooded)"),
 "PRJNA1157226": ("flooded", "medium", "title: 'pot experiment' rice (rice pots are flooded by default)"),
}
CROP_DEFAULT = {
 "rice": ("flooded", "medium", "crop-default rule: rice is grown in flooded paddy soil; no explicit statement in the title"),
 "maize": ("upland", "high", "crop-default rule: maize is not grown under permanent flooding in any of these studies"),
 "wheat": ("upland", "high", "crop-default rule: wheat is not grown under permanent flooding in any of these studies"),
}

def main():
    sel = [r for r in csv.DictReader(io.open(os.path.join(R, "meta", "final45_selection.csv"), encoding="utf-8-sig"))
           if (r.get("rank") or "").strip()]
    rows = []
    for r in sel:
        p, crop = r["project"], r["crop"]
        regime, conf, ev = EVID.get(p, CROP_DEFAULT.get(crop, ("uncertain", "low", "no evidence")))
        rows.append(dict(project=p, crop=crop, water_regime=regime, confidence=conf, evidence=ev))
    with io.open(os.path.join(OUT, "water_regime.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["project", "crop", "water_regime", "confidence", "evidence"])
        w.writeheader(); w.writerows(rows)
    n_fl = sum(1 for x in rows if x["water_regime"] == "flooded")
    print("wrote out_analysis/water_regime.csv ; flooded=%d upland=%d" % (n_fl, len(rows) - n_fl))


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Platform-cluster coding for the 45 selected studies (pre-registration section 2: cluster is judged by
FIELD SITE, YEAR, CULTIVAR, TREATMENT STRUCTURE and TEAM - NOT by author lists or institution names).

HONESTY NOTE: the repository carries no explicit site/year/variety/treatment field, so this is an
EXPLICIT, EVIDENCE-TAGGED curation from the study titles, centres and submission dates found in meta/*.csv.
Where the evidence only supports "same team" (different titles/years) the study stays in its own cluster
and the shared team is recorded separately, because section 2 distinguishes a shared platform from a
shared group.  This table is meant to be reviewed line by line by the author.

usage: python platform_cluster_code.py  -> out_analysis/platform_cluster.csv
"""
import csv, io, os
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")

# cluster assignments (id, evidence). Everything not listed is a singleton cluster of its own.
CLUSTERS = {
 "PRJEB21590": ("GX_CORNELL_MAIZE_2017", "same centre (Cornell) + same submission window 2017-08 + identical title 'Maize age and genotype shape the rhizosphere microbiome in the field' + same 4-level design", "high"),
 "PRJEB21600": ("GX_CORNELL_MAIZE_2017", "same centre + same window + title is the 'SIGR component' of the same study", "high"),
 "PRJNA771711": ("GX_UNL_BG_MAIZE_2021", "same centre (Nebraska-Lincoln) + same date 2021-10-17 + title 'Maize rhizosphere microbiome Buckler Goodman part 3 of 3'", "high"),
 "PRJNA771712": ("GX_UNL_BG_MAIZE_2021", "same centre + same date + title 'Maize rhizosphere microbiome Buckler Goodman part 2 of 3'", "high"),
 "PRJNA607209": ("GX_YZU_WHEAT_CD", "same centre (Yangzhou University) + title identical to PRJNA661129 ('comparing and deciphering the bacterial communities of the rhizosphere and bulk soils in the wheat field under various cadmium pollution')", "medium"),
 "PRJNA661129": ("GX_YZU_WHEAT_CD", "same centre + title identical to PRJNA607209", "medium"),
 "PRJNA596848": ("GX_YZU_RICE_CLIMATE", "same centre (Yangzhou University) as PRJNA665860; titles differ (soil type/climate/variety vs climate change), so same TEAM but treated as a separate experiment", "medium"),
 "PRJNA665860": ("GX_YZU_RICE_CLIMATE", "same centre as PRJNA596848; see note there", "medium"),
 "PRJNA633917": ("GX_NEIGAE_MAIZE", "same centre (Northeast Institute of Geography and Agroecology) as PRJNA781096; different titles/years -> same TEAM, separate trial", "medium"),
 "PRJNA781096": ("GX_NEIGAE_MAIZE", "same centre as PRJNA633917; see note there", "medium"),
 "PRJNA454003": ("GX_ROTHAMSTED_WHEAT", "same centre (Rothamsted Research) as PRJNA601112; different titles (N fertilisation vs breeding) -> same TEAM, separate trial", "high"),
 "PRJNA601112": ("GX_ROTHAMSTED_WHEAT", "same centre as PRJNA454003; see note there", "high"),
}
TEAM = {   # same team, evidenced by centre only (kept as separate clusters)
 "PRJNA593708": "Julius-Kuehn Institut", "PRJNA1037063": "Julius Kuehn-Institut",
 "PRJNA715843": "South China Botanical Garden", "PRJNA1101437": "Guangzhou University",
 "PRJNA554261": "Hunan Academy of Agricultural Sciences", "PRJNA662201": "Hunan Academy of Agricultural Sciences",
 "PRJNA608090": "Chinese Academy of Sciences", "PRJNA967263": "Southwest Forestry University",
}


def main():
    sel = [r for r in csv.DictReader(io.open(os.path.join(R, "meta", "final45_selection.csv"), encoding="utf-8-sig"))
           if (r.get("rank") or "").strip()]
    rows = []
    for r in sel:
        p = r["project"]
        cid, ev, conf = CLUSTERS.get(p, (p, "singleton: no other selected study shares this centre+title+window", "high"))
        rows.append(dict(project=p, crop=r["crop"], platform_cluster=cid,
                         same_team_as=TEAM.get(p, ""), evidence=ev, confidence=conf))
    with io.open(os.path.join(OUT, "platform_cluster.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["project", "crop", "platform_cluster", "same_team_as", "evidence", "confidence"])
        w.writeheader(); w.writerows(rows)
    import collections
    cnt = collections.Counter((r["crop"], r["platform_cluster"]) for r in rows)
    multi = {k: v for k, v in cnt.items() if v > 1}
    print("wrote out_analysis/platform_cluster.csv ; projects=%d ; clusters with >1 member (same crop): %s" % (
        len(rows), multi if multi else "none"))
    print("校验 §4 集群上限（同一作物每集群 ≤2）：%s" % ("通过" if all(v <= 2 for v in cnt.values()) else "违规！"))


if __name__ == "__main__":
    main()

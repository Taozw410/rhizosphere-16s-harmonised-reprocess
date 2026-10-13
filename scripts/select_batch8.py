# -*- coding: utf-8 -*-
"""Batch 8 maize: write selected list from screened probes."""
import csv, io, os
R = r"D:\CodexWork\mic-re"
rows = list(csv.DictReader(io.open(os.path.join(R, "out50k", "batch8_screening.csv"), encoding="utf-8-sig")))
ok = [r for r in rows if r["verdict"] == "OK"]
ok.sort(key=lambda r: r["project"])
with io.open(os.path.join(R, "meta", "batch8_selected.txt"), "w", encoding="utf-8-sig", newline="") as f:
    for r in ok:
        f.write(r["project"] + "\n")
print("batch8 selected =", len(ok))
for r in ok:
    print("  %-14s mode=%-7s rlen=%-4s qual=%-3s S16=%-3s ITS=%-3s ovl=%-5s n=%s" %
          (r["project"], r["mode"], r["rlen"], r["qual"], r["S16"], r["ITS"], r["ovl"], r["samples"]))

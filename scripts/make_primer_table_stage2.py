# 由探针筛查结果生成第二阶段引物表
import csv, io
scr = list(csv.DictReader(io.open(r"D:\CodexWork\mic-re\out50k\stage2_screening.csv", encoding="utf-8-sig")))
P341F, P805R = "CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"
P515F, P806R = "GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"
P338F, P806Rb = "ACTCCTACGGGAGGCAGCAG", "GGACTACHVGGGTWTCTAAT"
rows = ["project,primer_f,primer_r,truncLen_f,truncLen_r,note"]
for r in scr:
    if r["verdict"] != "OK":
        continue
    pr = r["primers"]
    if "341F" in pr or "805R" in pr:
        f, rv, note = P341F, P805R, "probe_341F_805R_V3V4"
    elif "338F" in pr:
        f, rv, note = P338F, P806Rb, "probe_338F_V3V4"
    elif "515F" in pr:
        f, rv, note = P515F, P806R, "probe_515F_806R_V4"
    elif "ITS1F" in pr or r["marker"] == "ITS?":
        f, rv, note = "CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC", "ITS_LIBRARY_check_marker"
    else:
        f, rv, note = P515F, P806R, "no_primer_detected_probe_decides"
    rows.append(f"{r['project']},{f},{rv},0,0,{note}")
io.open(r"D:\CodexWork\mic-re\meta\primer_table_stage2.csv", "w", encoding="utf-8-sig", newline="").write("\n".join(rows) + "\n")
print("written:", len(rows) - 1, "projects")
for x in rows[:6]: print(" ", x)

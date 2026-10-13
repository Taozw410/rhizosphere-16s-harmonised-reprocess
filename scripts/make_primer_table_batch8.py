import csv, io
sel = [x.strip() for x in io.open(r"D:\CodexWork\mic-re\meta\batch8_selected.txt", encoding="utf-8-sig").read().splitlines() if x.strip()]
P515F, P806R = "GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"
rows = ["project,primer_f,primer_r,truncLen_f,truncLen_r,note"]
for p in sel:
    rows.append(f"{p},{P515F},{P806R},0,0,batch8_no_primer_detected_probe_decides")
io.open(r"D:\CodexWork\mic-re\meta\primer_table_batch8.csv", "w", encoding="utf-8-sig", newline="").write("\n".join(rows) + "\n")
print("primer table written:", len(sel), "projects")

import re
# 对已完成项目做标记基因核查：属级注释率异常低 => 标记基因不符
import csv, io, os, sys
OUT = sys.argv[1] if len(sys.argv) > 1 else r"D:\CodexWork\mic-re\out_stage2"
THRESH = 20.0
rows = []
for d in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, d)
    if not os.path.isdir(p) or not re.match(r"^PRJ[A-Z]*\d+$", d): continue
    tp = os.path.join(p, "asv_taxonomy.csv")
    if not os.path.exists(tp): continue
    t = list(csv.DictReader(io.open(tp, encoding="utf-8-sig")))
    if not t: continue
    n = len(t)
    g = sum(1 for x in t if (x.get("Genus") or "").strip() and (x.get("Genus") or "").strip().upper() != "NA")
    gp = 100*g/n
    flag = "MARKER_MISMATCH" if gp < THRESH else ("LOW" if gp < 60 else "")
    rows.append((d, n, round(gp,1), flag))
rows.sort(key=lambda r: r[2])
print(f"{'project':<15}{'ASV':>7}{'genus%':>9}  flag")
for d,n,gp,f in rows:
    print(f"{d:<15}{n:>7}{gp:>9.1f}  {f}")
bad = [r for r in rows if r[3]=="MARKER_MISMATCH"]
print()
print("标记基因不符（属级注释率 <20%）：", len(bad))
for r in bad: print("  ", r[0], f"{r[2]}%")
low = [r for r in rows if r[3]=="LOW"]
print("偏低（20–60%）：", len(low))
for r in low: print("  ", r[0], f"{r[2]}%")
with io.open(os.path.join(OUT, "marker_check.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f); w.writerow(["project","n_asv","genus_pct","flag"]); w.writerows(rows)

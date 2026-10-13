# -*- coding: utf-8 -*-
"""把 6 个"连门都注释不上"的 artifact 型项目移入 _recycle（不真删），并写清单。"""
import os, io, csv, shutil, datetime
R = r"D:\CodexWork\mic-re"
TARGETS = ["PRJDB19884", "PRJEB65459", "PRJNA997935", "PRJEB104640", "PRJEB40690", "PRJNA394069"]
SLOTS = [("out_stage2", os.path.join(R, "out_stage2")),
         ("raw", os.path.join(R, "data", "raw")),
         ("raw50k", os.path.join(R, "data", "raw50k")),
         ("raw50k_longread", os.path.join(R, "data", "raw50k_longread"))]
dest_root = os.path.join(R, "_recycle", "unusable_20261010_artifact")
assert os.path.abspath(dest_root).startswith(os.path.abspath(os.path.join(R, "_recycle")))
rows = []
for p in TARGETS:
    moved = []
    for slot, base in SLOTS:
        src = os.path.join(base, p)
        if os.path.isdir(src):
            sz = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(src) for f in fs)
            dst = os.path.join(dest_root, p, slot)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.move(src, dst)
            moved.append((slot, sz))
    rows.append((p, moved))
    print("%-14s moved=%s  (%.1f MB)" % (p, ",".join(s for s, _ in moved), sum(s for _, s in moved) / 1e6))
with io.open(os.path.join(R, "out50k", "unusable_artifact_20261010.csv"), "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.writer(fh); w.writerow(["project", "moved_slots", "dest", "time"])
    for p, moved in rows:
        w.writerow([p, ";".join(s for s, _ in moved), os.path.join(dest_root, p), datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")])
print("manifest: out50k/unusable_artifact_20261010.csv")

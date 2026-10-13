# -*- coding: utf-8 -*-
"""Update truncLen rows in a DADA2 primer table (ASCII-only, BOM preserved).

Usage:
  python set_trunclen_rows.py <table.csv> "proj,tf,tr,note" ["proj2,tf2,tr2,note2" ...]

Rules: keep primer_f/primer_r untouched when the project already has a row; append a new row
with empty primers otherwise (empty primer = no cutadapt trimming, reads used as-is, which is
what the measured truncLen values assume). Refuses to write if any non-ASCII char is present.
"""
import sys, io, os

def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    path = sys.argv[1]
    entries = []
    for spec in sys.argv[2:]:
        parts = [p.strip() for p in spec.split(",")]
        if len(parts) < 3:
            sys.exit("bad entry: " + spec)
        proj, tf, tr = parts[0], parts[1], parts[2]
        note = parts[3] if len(parts) > 3 else "trunclen_set"
        entries.append((proj, tf, tr, note))

    with io.open(path, encoding="utf-8-sig") as fh:
        raw = fh.read()
    lines = [l for l in raw.replace("\r\n", "\n").split("\n")]
    header = lines[0].split(",")
    cols = {name: i for i, name in enumerate(header)}
    for need in ("project", "truncLen_f", "truncLen_r"):
        if need not in cols:
            sys.exit("table missing column: " + need + " (have " + ",".join(header) + ")")

    ncol = len(header)
    rows = [l.split(",") for l in lines[1:] if l.strip()]
    changed = []
    for proj, tf, tr, note in entries:
        hit = None
        for r in rows:
            if r[cols["project"]] == proj:
                hit = r
                break
        if hit is None:
            hit = [""] * ncol
            hit[cols["project"]] = proj
            rows.append(hit)
        old = (hit[cols["truncLen_f"]], hit[cols["truncLen_r"]])
        hit[cols["truncLen_f"]] = tf
        hit[cols["truncLen_r"]] = tr
        if "note" in cols:
            hit[cols["note"]] = note
        changed.append("%s %s,%s -> %s,%s" % (proj, old[0], old[1], tf, tr))

    out = ",".join(header) + "\r\n" + "\r\n".join(",".join(r) for r in rows) + "\r\n"
    bad = [(i, ch) for i, ch in enumerate(out) if ord(ch) > 127]
    if bad:
        sys.exit("refusing to write non-ASCII content at offset %d: %r" % bad[0])
    with io.open(path, "w", encoding="utf-8-sig", newline="") as fh:
        fh.write(out)
    for c in changed:
        print("updated", c)
    print("wrote", path, os.path.getsize(path), "bytes")

if __name__ == "__main__":
    main()

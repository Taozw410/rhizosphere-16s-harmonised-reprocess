"""Renumber [R##] citation tokens in the English manuscript by order of first appearance.

Reads  manuscript/07_manuscript_en.md           (source, tokens [R01]..[R38], internal annotation block)
Writes manuscript/07_manuscript_en_numbered.md  (submission candidate: numeric citations + reference list)
Writes manuscript/refs/refs_number_map.md       (R## -> number traceability, unreferenced refs flagged)

Run: D:\\bioenv\\envs\\riz16s\\python.exe scripts\\refs_number.py
"""

import io
import os
import re
import sys

SRC = os.path.join("manuscript", "07_manuscript_en.md")
DST = os.path.join("manuscript", "07_manuscript_en_numbered.md")
REFMD = os.path.join("manuscript", "refs", "refs.md")
MAPMD = os.path.join("manuscript", "refs", "refs_number_map.md")

TOKEN = re.compile(r"\[(R\d{2}(?:,R\d{2})*)\]")


def load_refs():
    """ref_id -> citation string, parsed from the CrossRef-verified table in refs.md."""
    out = {}
    with io.open(REFMD, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("| R"):
                continue
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) < 2:
                continue
            ref_id = cells[0].split()[0]
            if not re.fullmatch(r"R\d{2}", ref_id):
                continue  # skips the REF-## mapping table further down the file
            out[ref_id] = cells[1]
    return out


def main():
    refs = load_refs()
    with io.open(SRC, encoding="utf-8") as fh:
        text = fh.read()

    order = []
    for m in TOKEN.finditer(text):
        for rid in m.group(1).split(","):
            if rid not in order:
                order.append(rid)
    number = {rid: i + 1 for i, rid in enumerate(order)}

    def repl(m):
        nums = [number[rid] for rid in m.group(1).split(",")]
        return "[" + ",".join(str(n) for n in nums) + "]"

    body = TOKEN.sub(repl, text)
    body = "\n".join(l for l in body.split("\n") if not l.startswith("> "))
    body = body.rstrip("\n")

    lines = [body, "", "## References", ""]
    for rid in order:
        cite = refs.get(rid)
        if cite is None:
            print("FAIL: %s cited in text but absent from refs.md" % rid, file=sys.stderr)
            return 1
        lines.append("%d. %s" % (number[rid], cite))
    with io.open(DST, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    unused = sorted(set(refs) - set(order))
    map_lines = ["# Citation numbering map", "",
                 "| R## | number in manuscript |", "| --- | --- |"]
    for rid in order:
        map_lines.append("| %s | %d |" % (rid, number[rid]))
    map_lines += ["", "Not cited in the English manuscript: " + (", ".join(unused) if unused else "none")]
    with io.open(MAPMD, "w", encoding="utf-8") as fh:
        fh.write("\n".join(map_lines) + "\n")

    print("cited: %d refs; uncited: %d (%s)" % (len(order), len(unused), ", ".join(unused) or "-"))
    print("wrote %s and %s" % (DST, MAPMD))
    return 0


if __name__ == "__main__":
    sys.exit(main())

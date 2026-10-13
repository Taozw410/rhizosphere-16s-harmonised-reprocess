"""Stage a code-and-derived-data release for the manuscript (GitHub / Zenodo deposit).

Raw reads are NOT included (they are public in SRA/ENA; only accessions are listed).
Large intermediates, package caches and temporary files are excluded by rule.

    D:\\bioenv\\envs\\riz16s\\python.exe scripts\\build_release.py            # dry run: sizes only
    D:\\bioenv\\envs\\riz16s\\python.exe scripts\\build_release.py --build    # copy into release/<tag>/

Outputs release/<tag>/MANIFEST.md with per-file size and sha256 prefix.
"""

import hashlib
import io
import os
import shutil
import sys

TAG = "v1.0"
DEST_ROOT = os.path.join("release", TAG)

# (source path, destination subdirectory or "" for repository root)
INCLUDE_DIRS = [
    ("scripts", "scripts"),
    ("meta", "meta"),
    ("out_analysis", "out_analysis"),
    ("out50k", "selection"),
    ("supplementary", "supplementary"),
]
INCLUDE_FILES = [
    ("分析口径与实验预注册_20261008.md", "preregistration/分析口径与实验预注册_20261008.md"),
    ("out50k/选集报告.md", "selection/选集报告.md"),
    ("scripts/release_templates/README.md", "README.md"),
    ("scripts/release_templates/LICENSE", "LICENSE"),
    ("scripts/release_templates/LICENSE-DATA.md", "LICENSE-DATA.md"),
]

# Exclusions by name fragment or extension.
EXCLUDE_NAMES = ("__pycache__", "_bak", "pkgzips", "_tmp_", "_archive_", "_shared", "idtaxa_smoke",
                 "_dbg_", "_smoke_", "_help_", "_chk_", "_probe_", "release_templates")
EXCLUDE_EXT = (".rds", ".rda", ".npy", ".npz", ".zip", ".gz", ".fastq", ".fasta", ".fa", ".png.aux.xml")
# Figures are kept even though small binaries; PNG here means full-size figure exports, kept on purpose.
KEEP_EXT = (".png",)

# Only these out50k entries are staged (the directory holds gigabytes of intermediates).
SELECTION_ALLOW = ("选集报告.md", "power_simulation_report_final45.md", "idtaxa/idtaxa_vs_dada2.csv")


def excluded(rel):
    parts = rel.replace("\\", "/").split("/")
    if any(any(tok in p for tok in EXCLUDE_NAMES) for p in parts):
        return True
    ext = os.path.splitext(rel)[1].lower()
    if ext in KEEP_EXT:
        return False
    return ext in EXCLUDE_EXT


def collect():
    files = []
    for src, dst in INCLUDE_DIRS:
        if not os.path.isdir(src):
            continue
        for root, dirs, names in os.walk(src):
            dirs[:] = [d for d in dirs if not any(tok in d for tok in EXCLUDE_NAMES)]
            for n in names:
                full = os.path.join(root, n)
                rel = os.path.relpath(full, os.path.dirname(src) or ".")
                if src == "out50k":
                    rel_in_src = os.path.relpath(full, src).replace("\\", "/")
                    if rel_in_src not in SELECTION_ALLOW:
                        continue
                if excluded(rel):
                    continue
                files.append((full, os.path.join(dst, os.path.relpath(full, src))))
    for src, dst in INCLUDE_FILES:
        if os.path.isfile(src):
            files.append((src, dst))
    return files


def main():
    build = "--build" in sys.argv
    files = collect()
    total = 0
    per_dir = {}
    for full, _ in files:
        size = os.path.getsize(full)
        total += size
        key = full.replace("\\", "/").split("/")[0]
        per_dir[key] = per_dir.get(key, 0) + size
    print("release %s: %d files, %.1f MB" % (TAG, len(files), total / 1e6))
    for k, v in sorted(per_dir.items(), key=lambda kv: -kv[1]):
        print("   %-16s %8.1f MB" % (k, v / 1e6))
    if not build:
        print("\n(dry run; pass --build to stage into %s)" % DEST_ROOT)
        return 0

    if os.path.isdir(DEST_ROOT):
        print("refusing to overwrite existing %s" % DEST_ROOT)
        return 2
    rows = []
    for full, rel in files:
        dest = os.path.join(DEST_ROOT, rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(full, dest)
        with open(full, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()[:16]
        rows.append((rel.replace("\\", "/"), os.path.getsize(dest), digest))
    rows.sort()
    with io.open(os.path.join(DEST_ROOT, "MANIFEST.md"), "w", encoding="utf-8") as fh:
        fh.write("# Release manifest %s\n\n" % TAG)
        fh.write("Source: harmonised reprocessing of 45 public rhizosphere 16S projects.\n")
        fh.write("Raw reads are public in SRA/ENA and are not redistributed here.\n\n")
        fh.write("| file | bytes | sha256[:16] |\n| --- | --- | --- |\n")
        for rel, size, digest in rows:
            fh.write("| %s | %d | %s |\n" % (rel, size, digest))
    print("staged %d files into %s" % (len(rows), DEST_ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())

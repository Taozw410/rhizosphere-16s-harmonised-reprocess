"""Hygiene audit for manuscript text (AGENTS.md section 4.2 / 7.2).

Scans a markdown manuscript for internal information, placeholders and working notes that
must never reach a submission file, and writes a report. Docx/OOXML checks (line numbers,
PAGE field, comments, tracked changes, document properties) apply later to the docx build
and are NOT covered here.

Run: D:\\bioenv\\envs\\riz16s\\python.exe scripts\\audit_manuscript.py [path ...]
Default path: manuscript/07_manuscript_en_numbered.md
"""

import io
import os
import re
import sys

DEFAULT = os.path.join("manuscript", "07_manuscript_en_numbered.md")
REPORT = os.path.join("manuscript", "_audit", "hygiene_report.md")

# Tokens that are expected to still be present because the author has not supplied the value
# yet. They are reported separately and DO block submission, but they are not surprises.
KNOWN_PENDING = ["REPOSITORY-URL-AND-DOI"]

PATTERNS = [
    ("backtick", re.compile(r"`")),
    ("chinese", re.compile(r"[\u4e00-\u9fff]")),
    ("workspace_path", re.compile(r"(输出/|脚本/|公共/|广西版|全国版|AGENTS\.md|进度记录|_recycle|ref_region)")),
    ("internal_dir", re.compile(r"\b(out_analysis|out50k|out_stage2|out_abprime|raw50k_longread|manuscript/|meta/|scripts/|data/)\b")),
    ("internal_name", re.compile(r"(效应量数据提取表|analysis_[A-Za-z_]+\.json|\b[\w\-]+_筛_[\w\-]+\.csv)")),
    ("script_name", re.compile(r"\b(build|figures|analysis|export|probe|screen|run|refs|audit|manuscript|select|chain)_[A-Za-z_0-9]+\.(py|R|r|mjs|ps1)")),
    # Requires at least two characters before the dot, so author initials (D.R., J.R.) do not match.
    ("file_ext", re.compile(r"\b\w{2,}\.(?:csv|py|json|docx|md|tiff|eps|svg|R|r)\b")),
    ("old_journal", re.compile(r"Online Resource")),
    # Bracketed text, but numeric citations such as [1] or [16,17] are legitimate.
    ("placeholder", re.compile(r"(\[(?![0-9,;\s\u2013-]*\])[^\]]{0,90}\]|\bTBD\b|\bTODO\b|\bFIXME\b|to be completed|date of submission)", re.I)),
    ("work_note", re.compile(r"(should be confirmed|to be checked|author decision|this round|draft version|for internal use|AUTHOR_INPUT_NEEDED)", re.I)),
    ("identity", re.compile(r"([\w.+-]+@[\w-]+\.[\w.]+|\b\d{4}-\d{4}-\d{4}-[\dX]{4}\b|ORCID)", re.I)),
    ("stale_citation", re.compile(r"\[R\d{2}(?:,R\d{2})*\]")),
]


def audit(path):
    with io.open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().split("\n")
    hits, pending = [], []
    for i, line in enumerate(lines, 1):
        for name, rx in PATTERNS:
            for m in rx.finditer(line):
                snippet = m.group(0)
                rec = (os.path.basename(path), name, snippet, i, line.strip()[:160])
                if name == "placeholder" and any(k in snippet for k in KNOWN_PENDING):
                    pending.append(rec)
                else:
                    hits.append(rec)
    return hits, pending


def main():
    paths = sys.argv[1:] or [DEFAULT]
    all_hits, all_pending = [], []
    for p in paths:
        if not os.path.exists(p):
            print("MISSING FILE :: %s" % p)
            return 2
        h, pd = audit(p)
        all_hits += h
        all_pending += pd
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    out = ["# 卫生审核报告（markdown 层）", "",
           "扫描对象：" + "、".join(os.path.basename(p) for p in paths),
           "",
           "说明：docx/OOXML 层检查（行号、页脚 PAGE 域、批注、修订痕迹、文档属性）需在 docx 打样上另跑。",
           "",
           "## 命中（必须清零）", ""]
    if all_hits:
        out += ["| 文件 | 类别 | 片段 | 行 | 上下文 |", "| --- | --- | --- | --- | --- |"]
        out += ["| %s | %s | %s | %d | %s |" % r for r in all_hits]
    else:
        out += ["无命中。"]
    out += ["", "## 已知待作者填入（阻塞投稿，但非意外残留）", ""]
    if all_pending:
        out += ["| 文件 | 类别 | 片段 | 行 | 上下文 |", "| --- | --- | --- | --- | --- |"]
        out += ["| %s | %s | %s | %d | %s |" % r for r in all_pending]
    else:
        out += ["无。"]
    out += ["", "## 结论", "",
            "- 意外命中：%d（必须为 0）" % len(all_hits),
            "- 待作者填入：%d（投稿前必须替换为真实值，否则不得进入投稿包）" % len(all_pending)]
    with io.open(REPORT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    def safe(x):
        return str(x).encode("ascii", "replace").decode("ascii")

    for rec in all_hits:
        f, name, snippet, line_no, ctx = rec
        print("HIT  %s :: %s :: %s :: line %s :: %s" % (safe(f), safe(name), safe(snippet), line_no, safe(ctx)))
    for rec in all_pending:
        f, name, snippet, line_no, ctx = rec
        print("PEND %s :: %s :: %s :: line %s :: %s" % (safe(f), safe(name), safe(snippet), line_no, safe(ctx)))
    print("\nunexpected hits: %d; pending author input: %d; report: %s"
          % (len(all_hits), len(all_pending), REPORT))
    return 1 if (all_hits or all_pending) else 0


if __name__ == "__main__":
    sys.exit(main())

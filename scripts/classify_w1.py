import re
# 三分类清单：可用 / 需修正 / 排除
# 用法: classify_w1.py [输出目录]   默认 out50k
import csv, os, gzip, glob, io, datetime, sys

OUT   = sys.argv[1] if len(sys.argv) > 1 else r"D:\CodexWork\mic-re\out50k"
RAW   = r"D:\CodexWork\mic-re\data\raw50k"
EXC   = r"D:\CodexWork\mic-re\data\raw50k_excluded"
META  = r"D:\CodexWork\mic-re\meta"

def rd(p):
    if not os.path.exists(p): return None
    with io.open(p, encoding="utf-8-sig") as f: return list(csv.DictReader(f))

def rlen_of(proj):
    for base in (RAW, EXC):
        fs = sorted(glob.glob(os.path.join(base, proj, "*.fastq.gz")))
        if fs:
            with gzip.open(fs[0], "rt", errors="replace") as fh:
                for i, l in enumerate(fh):
                    if i % 4 == 1: return len(l.strip())
    return None

proj_meta = {}
for r in rd(os.path.join(META, "lightweight_60_final.csv")) or []:
    proj_meta[r["study"]] = r
for r in rd(os.path.join(META, "w1_replacement_candidates.csv")) or []:
    proj_meta.setdefault(r["study"], r)

primers = {r["project"]: r for r in (rd(os.path.join(META, "primer_table.csv")) or [])}
# ---- measured insert length -------------------------------------------------
# Deciding "is the read longer than the amplicon" must use a MEASURED insert length,
# not the proxy rule "read length >= 290". When the insert is longer than the read (L > R),
# a low filter-pass rate comes from 3' tail quality, and truncLen is only a trade-off.
# Source 1: meta/amp_len_*.csv (measured by measure_amp.py from raw reads; has read_through)
# Source 2: this project's own asv_table.csv ASV length mode (post-chimera; paired only)
insert_measured = {}
amp_status = {}
for _f in sorted(glob.glob(os.path.join(META, "amp_len_*.csv"))):
    for _r in (rd(_f) or []):
        _v = (_r.get("amp_len") or "").strip()
        amp_status[_r["project"]] = (_r.get("status") or "").strip()
        if _v and _v.upper() != "NA":
            try:
                insert_measured[_r["project"]] = (int(float(_v)), os.path.basename(_f))
            except ValueError:
                pass


def insert_len_of(proj, out, paired):
    """Return (measured insert length, source); (None, "") when it cannot be determined."""
    if not paired:
        return None, ""
    if proj in insert_measured:
        return insert_measured[proj]
    t = os.path.join(out, "asv_table.csv")
    if os.path.exists(t):
        c = {}
        try:
            with io.open(t, encoding="utf-8-sig") as fh:
                rr = csv.reader(fh)
                next(rr, None)
                for row in rr:
                    if not row:
                        continue
                    L = len(row[0].strip().strip('"'))
                    c[L] = c.get(L, 0) + 1
        except Exception:
            return None, ""
        if c:
            return max(c.items(), key=lambda kv: kv[1])[0], "asv_table"
    return None, ""


projects = []
for base, excluded_dir in ((RAW, False), (EXC, True)):
    if os.path.isdir(base):
        for d in sorted(os.listdir(base)):
            if os.path.isdir(os.path.join(base, d)) and re.match(r"^PRJ[A-Z]*\d+$", d):
                projects.append((d, excluded_dir))

rows = []
for proj, in_exc_dir in projects:
    out = os.path.join(OUT, proj)
    qc  = rd(os.path.join(out, "project_qc.csv"))
    qc  = qc[0] if qc else {}
    rt  = rd(os.path.join(out, "read_tracking.csv")) or []
    merge = None; keep = None; fpass = None
    if rt:
        tin   = sum(float(x["input"])    for x in rt)
        tfilt = sum(float(x["filtered"]) for x in rt)
        if tin:
            fpass = 100 * tfilt / tin
        # 合并率必须以"过滤后"为基数：merged/filtered（merged/input 会混入过滤损失）
        if tfilt and "merged" in rt[0]:
            merge = 100 * sum(float(x["merged"]) for x in rt) / tfilt
        if tin:
            keep = 100 * sum(float(x["nonchim"]) for x in rt) / tin
    qual = qc.get("quality_chars", "")
    # 过滤后剩余 reads（用于识别"全部被质量过滤剔除"的项目：这类项目连
    # read_tracking.csv 都写不出来，只能从 filter_stats.csv 判）
    fsum = None; fnsamp = 0
    _fsp = os.path.join(out, "filter_stats.csv")
    if os.path.exists(_fsp):
        try:
            _fr = list(csv.DictReader(io.open(_fsp, encoding="utf-8-sig")))
            fsum = sum(int(float(x.get("reads.out") or 0)) for x in _fr)
            fnsamp = len(_fr)
        except Exception:
            fsum = None; fnsamp = 0
    try: qual = int(qual)
    except Exception: qual = None
    asv = qc.get("n_asv", "")
    try: asv = int(asv)
    except Exception: asv = None
    mode = qc.get("mode", "")
    rlen = rlen_of(proj)
    paired = (mode == "paired") or (rlen is not None and mode == "")
    insert, insert_src = insert_len_of(proj, out, paired)

    reasons = []; tier = None
    if in_exc_dir:
        # 排除目录里可能混合多种排除原因，优先用引物表的 note 说明
        nt = (primers.get(proj, {}) or {}).get("note", "") or ""
        reason = "已排除"
        if "placeholder" in nt.lower() or nt == "":
            reason = "质量值占位符（常数），无法估计错误率"
        if "pairing" in nt.lower():
            reason = "配对数据损坏"
        if "merge" in nt.lower():
            reason = "合并率无法达标（已尽力调参）"
        if "EXCLUDED" in nt:
            reason = nt.replace("EXCLUDED_", "").replace("_", " ")
        tier = 3; reasons.append(reason)
    elif qual is not None and qual <= 2:
        tier = 3; reasons.append("质量值占位符（常数）")
    elif paired and rlen is not None and insert is not None and insert > 2 * rlen - 12 and merge is not None and merge < 40:
        # 插入片段长于两读有效覆盖 → 根本无可重叠（可变区选错），truncLen 无法修复
        tier = 3; reasons.append(f"双端无可检出重叠：实测插入长度 {insert}bp > 两读有效覆盖 {2*rlen}bp（读长 {rlen}bp），合并率 {merge:.1f}%，truncLen 无法修复（可变区选错）")
    elif fsum == 0 and fnsamp > 0:
        tier = 3; reasons.append(f"过滤后无 reads（{fnsamp} 个样本全部被质量过滤剔除），数据不可用")
    elif merge is not None and merge < 20 and (rlen is None or rlen < 290):
        # 读长未越界、且两读仍可重叠，却几乎合并不上 → 归档数据本身有问题
        tier = 3; reasons.append(f"配对异常：读长 {rlen}bp 正常但合并率仅 {merge:.1f}%，_1/_2 无法配对")
    elif paired and rlen is not None and insert is not None and rlen >= insert and (fpass is None or fpass < 60):
        # 读穿（读长 ≥ 实测插入长度）→ 截到 ≤ 插入长度即可消除读穿，必须重跑
        fp = f"{fpass:.0f}%" if fpass is not None else "未知"
        tier = 2; reasons.append(f"读长 {rlen}bp ≥ 实测插入长度 {insert}bp（读穿），过滤通过率 {fp}，需设 truncLen 后重跑")
    elif merge is not None and merge < 40:
        ov = (2 * rlen - insert) if (rlen is not None and insert is not None) else None
        if paired and rlen is not None and insert is not None and rlen < insert and ov is not None and ov >= 12:
            tier = 2; reasons.append(f"合并率偏低（{merge:.1f}%）：读长 {rlen}bp < 实测插入长度 {insert}bp（重叠仅 {ov}bp，落在 3' 尾段）；truncLen 属取舍，建议 A/B 实测")
        elif paired and amp_status.get(proj) == "NO_ALIGN" and merge < 10:
            tier = 3; reasons.append(f"双端无可检出重叠（{rlen}bp 配对读段 0 对齐），插入片段长于两读覆盖，合并率 {merge:.1f}%，truncLen 无法修复")
        else:
            tier = 2; reasons.append(f"合并率偏低（{merge:.1f}%），需复核引物/接头处理")
    elif merge is None and fpass is not None and fpass < 10 and qc:
        # 过滤通过率 <10% 且无合并数据：读段几乎全被质量过滤剔除，数据不可用
        tier = 3; reasons.append(f"过滤通过率仅 {fpass:.1f}%（读段几乎全部被质量过滤剔除），数据不可用")
    if tier is None and qual is not None and 3 <= qual <= 5:
        tier = 2; reasons.append(f"质量值仅 {qual} 档（被量化），建议做敏感性分析")
    if tier is None and asv is not None and asv > 10000:
        tier = 2; reasons.append(f"ASV 数 {asv} 偏高，需复核是否虚高")
    if tier is None and qc:
        tier = 1; reasons.append("质量值正常、合并率与 ASV 数在正常区间")
    if tier is None:
        tier = 0; reasons.append("尚未跑完")

    m = proj_meta.get(proj, {})
    rows.append(dict(project=proj, tier=tier, centre=(m.get("center","") or "")[:34], insert=(insert if insert is not None else ""),
                     year=(m.get("date","") or "")[:4], samples=qc.get("n_samples",""),
                     mode=mode or ("excluded" if in_exc_dir else ""), rlen=rlen or "",
                     qual=qual if qual is not None else "", asv=asv if asv is not None else "",
                     fpass=(f"{fpass:.1f}%" if fpass is not None else ""),
                     merge=(f"{merge:.1f}%" if merge is not None else ""),
                     keep=(f"{keep:.1f}%" if keep is not None else ""),
                     reason="；".join(reasons)))

order = {1:0, 2:1, 3:2, 0:3}
rows.sort(key=lambda r: (order[r["tier"]], r["project"]))
names = {1:"✅ 可用", 2:"⚠️ 需修正", 3:"❌ 排除", 0:"⏳ 待跑完"}

lines = ["# W1 数据可用性三分类清单\n",
         f"生成时间：{datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n",
         "**分类规则（按成因，不按症状）**\n",
         "- **指标口径**：过滤通过率 = filtered/input（诊断读长越界）；**合并率 = merged/filtered**（诊断配对/接头问题）。两者不可混为一谈。",
         "- **排除**：质量值占位符（≤2 档，不可修复）；合并率 <20% 且非读穿（归档数据本身有问题）；双端无可检出重叠 / 插入片段长于两读覆盖（可变区选错，truncLen 无法修复）",
         "- **需修正**：读长 ≥ 实测插入长度（读穿，需按项目设 truncLen）；合并率 <40%；质量值被量化（3–5 档）；ASV 数 >10000 需复核",
         "- **插入长度来源**：meta/amp_len_*.csv（measure_amp.py 从 raw reads 实测）优先，其次该项目 asv_table.csv 的 ASV 长度主峰（去嵌合后）；判据用实测值，不用读长 ≥290 这类代理指标。",
         "- **可用**：质量值正常、过滤通过率与合并率均在正常区间\n"]

for t in (1, 2, 3, 0):
    sub = [r for r in rows if r["tier"] == t]
    if not sub: continue
    lines.append(f"\n## {names[t]}（{len(sub)} 个）\n")
    lines.append("| 项目 | 机构 | 年 | 样本 | 类型 | 读长 | 质量档 | ASV | 过滤通过率 | 合并率 | 保留率 | 判定依据 |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sub:
        lines.append("| {project} | {centre} | {year} | {samples} | {mode} | {rlen} | {qual} | {asv} | {fpass} | {merge} | {keep} | {reason} |".format(**r))

tot = {t: len([r for r in rows if r["tier"] == t]) for t in (1,2,3,0)}
lines.append(f"\n---\n**总计：可用 {tot[1]} 个 / 需修正 {tot[2]} 个 / 排除 {tot[3]} 个 / 待跑完 {tot[0]} 个**\n")

with io.open(os.path.join(OUT, "W1_usability_tiers.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
with io.open(os.path.join(OUT, "W1_usability_tiers.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

print("可用", tot[1], "| 需修正", tot[2], "| 排除", tot[3], "| 待跑完", tot[0])
print("written: W1_usability_tiers.md")

import re
# 完整合并报告：合并 W1(out50k) + 第二阶段(out_stage2) + 全部筛查记录
import csv, os, io, sys, gzip, glob, datetime

R = r"D:\CodexWork\mic-re"
W1 = os.path.join(R, "out50k")
S2 = os.path.join(R, "out_stage2")
META = os.path.join(R, "meta")

def rd(p):
    if not os.path.exists(p): return []
    with io.open(p, encoding="utf-8-sig") as f: return list(csv.DictReader(f))

# ---- 项目级 QC ----
rows = []
def load(outdir, source):
    for d in sorted(os.listdir(outdir)):
        p = os.path.join(outdir, d)
        if not os.path.isdir(p) or not re.match(r"^PRJ[A-Z]*\d+$", d): continue
        qc = rd(os.path.join(p, "project_qc.csv"))
        if not qc: continue
        qc = qc[0]
        rt = rd(os.path.join(p, "read_tracking.csv"))
        fpass = merge = keep = None
        if rt:
            tin = sum(float(x["input"]) for x in rt)
            tf  = sum(float(x["filtered"]) for x in rt)
            if tin: fpass = 100*tf/tin
            if tf and "merged" in rt[0]:
                merge = 100*sum(float(x["merged"]) for x in rt)/tf
            if tin: keep = 100*sum(float(x["nonchim"]) for x in rt)/tin
        tax = rd(os.path.join(p, "asv_taxonomy.csv"))
        gp = None
        if tax:
            n=len(tax); g=sum(1 for t in tax if t.get("Genus","").strip() and t.get("Genus","").strip().upper()!="NA")
            gp = 100*g/n if n else None
        rows.append(dict(project=d, source=source, samples=qc.get("n_samples",""),
                         mode=qc.get("mode",""), qual=qc.get("quality_chars",""),
                         asv=qc.get("n_asv",""),
                         fpass=(f"{fpass:.1f}" if fpass is not None else ""),
                         merge=(f"{merge:.1f}" if merge is not None else ""),
                         keep=(f"{keep:.1f}" if keep is not None else ""),
                         genus=(f"{gp:.1f}" if gp is not None else "")))
load(W1, "W1")
load(S2, "S2")

# ---- 筛查记录（阶段2 + 批次2） ----
scr = []
pmeta = {}
for f in ("lightweight_60_final.csv","w1_replacement_candidates.csv","stage2_projects.csv","batch2_candidates.csv","batch2_selected.csv"):
    for r in rd(os.path.join(META, f)):
        k = r.get("study") or r.get("project")
        if k: pmeta.setdefault(k, r)

w1tier = {r["project"]: r for r in rd(os.path.join(W1,"W1_usability_tiers.csv"))}
s2scr  = {r["project"]: r for r in rd(os.path.join(W1,"stage2_screening.csv"))}
b2scr  = {r["project"]: r for r in rd(os.path.join(W1,"batch2_screening.csv"))}

out = []
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
out.append("# 甘蔗/作物根际微生物组公共数据重分析 —— 完整合并报告\n")
out.append(f"生成时间：{now}\n")
out.append("## 一、总体结论\n")
n_ok = sum(1 for r in rows if r["source"]=="S2") + sum(1 for r in rows if r["source"]=="W1")
out.append(f"- 已处理并产出 ASV 表的项目：**{len(rows)} 个**")
out.append(f"- 覆盖样本：**{sum(int(r['samples']) for r in rows if str(r['samples']).isdigit())} 个**")
out.append(f"- 累计 ASV：**{sum(int(r['asv']) for r in rows if str(r['asv']).isdigit())} 个**\n")

# 按作物统计
crop_of = {}
for r in rows:
    m = pmeta.get(r["project"], {})
    crop_of[r["project"]] = m.get("crop","")
out.append("### 按作物\n")
out.append("| 作物 | 项目数 | 样本 | ASV |")
out.append("|---|---|---|---|")
for c in ("maize","rice","wheat"):
    g=[r for r in rows if crop_of.get(r["project"])==c]
    if not g: continue
    s=sum(int(r["samples"]) for r in g if str(r["samples"]).isdigit())
    a=sum(int(r["asv"]) for r in g if str(r["asv"]).isdigit())
    out.append(f"| {c} | {len(g)} | {s} | {a} |")
out.append("")

# ---- 逐项目明细 ----
out.append("## 二、逐项目结果\n")
out.append("| 项目 | 来源 | 样本 | 类型 | 质量档 | ASV | 过滤通过率% | 合并率% | 保留率% | 属级注释率% |")
out.append("|---|---|---|---|---|---|---|---|---|---|")
for r in sorted(rows, key=lambda x:(crop_of.get(x["project"],""), x["project"])):
    out.append("| {project} | {source} | {samples} | {mode} | {qual} | {asv} | {fpass} | {merge} | {keep} | {genus} |".format(**r))
out.append("")

# ---- 三类缺陷命中 ----
out.append("## 三、公共数据缺陷统计\n")
ph = [p for p,r in s2scr.items() if r.get("verdict")=="EXCLUDE"] + [p for p,r in b2scr.items() if r.get("verdict")=="EXCLUDE"]
out.append("### 阶段二筛查（60 项目池）\n")
out.append("| 判定 | 项目数 |")
out.append("|---|---|")
for v,lab in [("OK","可直接使用"),("CHECK","需加敏感性分析"),("EXCLUDE","排除")]:
    n = len([1 for r in s2scr.values() if r.get("verdict")==v])
    out.append(f"| {v}（{lab}） | {n} |")
out.append("")
out.append("### 批次2替补候选筛查（44 个）\n")
out.append("| 判定 | 项目数 |")
out.append("|---|---|")
for v,lab in [("OK","可直接使用"),("CHECK","需加敏感性分析"),("EXCLUDE","排除")]:
    n = len([1 for r in b2scr.values() if r.get("verdict")==v])
    out.append(f"| {v}（{lab}） | {n} |")
out.append("")

# ---- 异常标记 ----
out.append("## 四、需要留意的异常项目\n")
out.append("| 项目 | 异常 | 说明 |")
out.append("|---|---|---|")
for r in rows:
    try: gp=float(r["genus"]) if r["genus"] else None
    except: gp=None
    try: a=int(r["asv"])
    except: a=None
    if gp is not None and gp < 60:
        out.append(f"| {r['project']} | 属级注释率仅 {gp:.0f}% | 远低于正常 85–98%，需核查引物/可变区 |")
    if a is not None and a > 20000:
        out.append(f"| {r['project']} | ASV 数 {a} 偏高 | 30 样本正常 5,000–15,000，疑引物未去或质量异常 |")
out.append("")

# ---- 已知缺陷类型 ----
out.append("## 五、本项目发现的 5 类公共数据缺陷\n")
out.append("| # | 缺陷类型 | 表现 | 检出方式 |")
out.append("|---|---|---|---|")
out.append("| 1 | 质量值占位符 | 质量值仅 1–2 档（常数） | `screen_quality.py`；DADA2 错误率学习直接失败 |")
out.append("| 2 | 配对数据损坏 | 读长正常但合并率 <10% | 逐样本合并率 + `overlap_scan.py` |")
out.append("| 3 | 读长超过扩增子 | 过滤通过率骤降；合并率受抑制 | ASV 长度分布反推扩增子长度 |")
out.append("| 4 | 质量值被量化 | 质量值 3–5 档；ASV 数虚高 | 质量值种类统计 |")
out.append("| 5 | 项目内标记基因混杂 | 同项目内一半 16S 一半 ITS | 逐样本过滤通过率双峰 + reads 开头序列比对 |")
out.append("")

with io.open(os.path.join(R, "out50k", "完整合并报告.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(out))
print("projects:", len(rows))
print("samples:", sum(int(r['samples']) for r in rows if str(r['samples']).isdigit()))
print("written: 完整合并报告.md")

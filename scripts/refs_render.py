"""Render manuscript/refs/refs.md from CrossRef-verified metadata + curated usage map.

Authors/volume/pages come only from CrossRef (verified_refs.json); the usage
column is the human curation. No bibliographic field is typed by hand.

Run: D:\\bioenv\\envs\\riz16s\\python.exe scripts\\refs_render.py
"""

import json
import os
import re

OUTDIR = os.path.join("manuscript", "refs")

# CrossRef records whose author fields are mis-split by the publisher deposit.
# Corrected ONLY against an independent source (OpenAlex record for the same DOI),
# never from memory; each correction is printed as a footnote in refs.md.
CORRECTIONS = {
    # CrossRef: family "O Bohn"/"D Kent" (initials leaked into the surname field);
    # OpenAlex same-DOI record: "Alonso Favela", "Martin Bohn", "Angela D. Kent".
    "R02": [{"family": "Favela", "given": "Alonso"},
            {"family": "Bohn", "given": "Martin"},
            {"family": "Kent", "given": "Angela D."}],
}
FOOTNOTES = {
    "R02": "CrossRef 元数据把中间名首字母并入姓氏（“O Bohn”“D Kent”）；此处按 OpenAlex 同 DOI 记录的作者字段校正为 Bohn, M. 与 Kent, A.D.。",
}

# ref_id -> (supports which part of the manuscript, one-line purpose)
USAGE = {
    "R01": ("P1; REF-13", "Field-scale evidence that maize rhizosphere composition is host-influenced and partly heritable"),
    "R02": ("P1; REF-13", "Crop breeding history shifts rhizosphere recruitment within one crop"),
    "R03": ("P1; REF-13", "Rice root/rhizosphere microbiome synthesis: crop-specific expectations"),
    "R04": ("P1; REF-18", "Host selection shapes crop microbiome assembly and network complexity"),
    "R05": ("P1; REF-14", "Macroecological framing of rhizosphere microbiomes; notes the missing cross-study scale"),
    "R06": ("P1; REF-03", "Study identity drives distance-decay differences: the cross-study comparability problem"),
    "R07": ("P2; REF-12", "Quantifying quality of public SRA data before reanalysis"),
    "R08": ("P2; REF-12", "Reprocessing public 16S studies: metadata gaps and robustness problems"),
    "R09": ("P2; REF-12", "Contamination issues in low-biomass studies (why negative controls matter)"),
    "R10": ("P2; 4.5", "decontam: statistical removal of contaminant sequences"),
    "R11": ("P3; 2.3; REF-11", "DADA2 ASV inference (the denoiser used here)"),
    "R12": ("P3; REF-11", "Benchmarking of clustering/denoising choices changes interpretation"),
    "R13": ("P3; REF-11", "Trimming/read-length choices change recovered mock-community composition"),
    "R14": ("P3; REF-10", "Pipeline choice biases environmental amplicon interpretation"),
    "R15": ("P3; REF-17", "Primer and platform effects on 16S tag sequencing"),
    "R16": ("P3; REF-10", "Primer/pipeline/parameter effects across 16S workflows"),
    "R17": ("P3; 2.3", "QIIME 2: reproducible, provenance-tracked amplicon analysis"),
    "R18": ("P3; REF-16", "IDTAXA: the second, independent taxonomy classifier used here"),
    "R19": ("P3; REF-16", "Comparison of machine-learning and database taxonomy classifiers"),
    "R20": ("P5; REF-02", "Electron-acceptor sequence (NO3- -> Fe(III) -> SO4(2-) -> CO2) in anoxic rice soil"),
    "R21": ("P5; REF-05", "Iron-reducing taxa identified directly in anoxic paddy soil"),
    "R22": ("P5; REF-04", "Methane oxidation and oxygen competition in the rice rhizosphere"),
    "R23": ("P5; REF-04", "Methanotroph community composition controls paddy methane oxidation"),
    "R24": ("P5; REF-04, REF-07", "Methanotroph-nitrifier competition couples C and N cycling in paddy soil"),
    "R25": ("P5; REF-05", "Sulfate turnover and sulfate-reducing bacteria in planted paddy soil"),
    "R26": ("P5; 3.10; REF-08", "How far FAPROTAX assignments can be trusted for soil bacteria"),
    "R27": ("3.10; REF-08", "Tax4Fun2: functional inference from 16S and its habitat dependence"),
    "R28": ("P5; REF-09", "Water management (AWD) changes methane/nitrous oxide emissions in rice fields"),
    "R29": ("P5; REF-09", "Upland versus irrigated rice rhizospheres: a design that separates water regime"),
    "R30": ("P5; REF-09", "Upland (aerobic) rice root microbiomes and drought tolerance"),
    "R31": ("3.4; methods", "Baselga turnover/nestedness partitioning used in the beta decomposition"),
    "R32": ("3.3; methods", "UniFrac distances used alongside Bray-Curtis and Aitchison"),
    "R33": ("4.1; 4.2; REF-01", "Core-microbiome assignments are definition- and method-dependent"),
    "R34": ("3.9; methods", "Normalization/differential-abundance choice depends on data characteristics"),
    "R35": ("3.8", "SIAMCAT: cross-study microbiome machine learning and transferability"),
    "R36": ("2.7; methods", "Rarefaction as the control for uneven sequencing effort"),
    "R37": ("3.10; caveat", "16S copy-number correction is unsolved: limits of read-based functional claims"),
    "R38": ("3.10; 4.3", "Methylotrophic bacteria are enriched in rice paddy bulk soil and rhizosphere"),
}


def initials(given: str) -> str:
    out = []
    for tok in re.split(r"[\s\-]+", (given or "").strip()):
        if tok:
            out.append(tok[0].upper() + ".")
    return "".join(out)


def authors_str(info: dict) -> str:
    structured = CORRECTIONS.get(info.get("ref_id")) or info.get("authors_structured") or []
    flipped = []
    for a in structured:
        fam = (a.get("family") or "").strip()
        if not fam:
            continue
        ini = initials(a.get("given") or "")
        flipped.append(f"{fam}, {ini}" if ini else fam)
    if not flipped:  # fall back to the flat strings if structured data is absent
        flipped = [a for a in (info.get("authors") or [])]
    if len(flipped) <= 3:
        return ", ".join(flipped)
    return ", ".join(flipped[:3]) + ", et al."


def main() -> int:
    with open(os.path.join(OUTDIR, "verified_refs.json"), encoding="utf-8") as fh:
        refs = json.load(fh)
    lines = [
        "# 文献表（CrossRef 核实，2026-10-13）",
        "",
        "> 生成方式：`scripts/refs_search.py`（OpenAlex 检索）→ 人工挑选 → `scripts/refs_verify.py`（逐条 CrossRef 核对题名/作者/年份）",
        "> → 本文件由 `scripts/refs_render.py` 从核实后的元数据渲染。**作者、卷期、页码全部来自 CrossRef，无手工键入**。",
        "> 未核实或题名不匹配者一律未收录；检索原始命中见 `raw/`、`_raw_index*.md`。",
        "> 正文引用一律用本表编号 `[R##]`，编号顺序即初次引用顺序（定稿时按目标刊要求调整为数字或作者-年份制）。",
        "",
        "| 编号 | 完整引文 | DOI | 支撑本文 | 用途 |",
        "| --- | --- | --- | --- | --- |",
    ]
    for info in refs:
        rid = info["ref_id"]
        usage = USAGE.get(rid, ("", ""))
        vol = info.get("volume") or ""
        issue = info.get("issue") or ""
        page = info.get("page") or ""
        loc = vol
        if issue:
            loc += f"({issue})"
        if page:
            loc += f", {page}"
        cite = (f"{authors_str(info)} ({info['year']}). {info['title'].rstrip('.')}. "
                f"*{info['journal']}*, {loc}. https://doi.org/{info['doi']}")
        flag = " †" if rid in CORRECTIONS else ""
        lines.append(f"| {rid}{flag} | {cite} | {info['doi']} | {usage[0]} | {usage[1]} |")
    if CORRECTIONS:
        lines += ["", "† 作者字段校正记录：", ""]
        for rid, note in FOOTNOTES.items():
            lines.append(f"- **{rid}**：{note}")
    lines += [
        "",
        "## 与讨论稿 `[REF-###]` 占位符的对应",
        "",
        "| 占位符 | 主题 | 本表对应编号 |",
        "| --- | --- | --- |",
        "| REF-01 | 核心微生物组的定义与尺度依赖 | R33（+R05 尺度框架） |",
        "| REF-02 | 稻田电子受体序列与产甲烷 | R20 |",
        "| REF-03 | 公共数据/多研究中研究层随机效应 | R06（+R35 跨研究整合） |",
        "| REF-04 | 稻田甲烷氧化与甲基营养 | R22、R23、R24（+R38 甲基营养） |",
        "| REF-05 | Geobacter/Desulfovibrio 与稻田铁硫还原 | R21、R25 |",
        "| REF-06 | 淹水对好氧异养的抑制 | R22、R29（间接；见备注） |",
        "| REF-07 | 稻田氮循环（硝酸呼吸） | R24（+R20 电子受体序列） |",
        "| REF-08 | 16S 功能预测的边界（潜力 ≠ 活性） | R26、R27、R37 |",
        "| REF-09 | 水旱轮作、落干、upland rice | R28、R29、R30 |",
        "| REF-10 | 扩增子区域对群落比较的影响 | R14、R15、R16 |",
        "| REF-11 | DADA2 参数（截断、合并率）影响 | R11、R12、R13 |",
        "| REF-12 | 公共数据质量缺陷与重分析陷阱 | R07、R08、R09、R10 |",
        "| REF-13 | 作物/基因型的根际效应 | R01、R02、R03 |",
        "| REF-14 | 跨地点/跨研究可复现性 | R05、R06 |",
        "| REF-15 | 合并率、读长与扩增子长度的关系 | R13（+R12） |",
        "| REF-16 | 分类学注释器偏差与交叉验证 | R18、R19 |",
        "| REF-17 | 扩增子区域/引物效应 | R15、R16（+R14） |",
        "",
        "## 备注（需要作者判断的两处）",
        "",
        "1. **REF-06（淹水抑制好氧异养）**：本次检索未找到直接以“淹水→好氧异养受抑”为结论的专文；",
        "   目前用 R22（根际氧竞争）与 R29（水田/旱作水稻根际对比）间接支撑。若作者要求更强支撑，",
        "   需补检 `(anoxia OR oxygen depletion) AND soil AND (aerobic heterotroph OR copiotroph)` 一类检索式。",
        "2. **REF-03（统计单元/研究层随机效应）**：微生物组领域缺一篇“必须把研究作为随机效应层”的方法学专文；",
        "   现用 R06（研究间差异驱动距离衰减）+ R35（跨研究整合）合并支撑，写作时应把这句写成“方法学共识/经验做法”，",
        "   不写成“某某证明”。",
    ]
    with open(os.path.join(OUTDIR, "refs.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote refs.md with {len(refs)} refs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

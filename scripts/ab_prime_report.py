# -*- coding: utf-8 -*-
"""A' 裁决对照：把 A' 臂（独立输出根）与 A 臂（现行结果）逐项目对比，输出四类读数。
usage: python ab_prime_report.py [--abprime out_abprime] [--out out50k/abprime_report.md]
"""
import argparse, csv, io, os
R = r"D:\CodexWork\mic-re"
P5 = ["PRJNA715843", "PRJNA593708", "PRJNA427813", "PRJNA596848", "PRJNA1310568"]
PILOT = os.path.join(R, "_recycle", "phase5b_rejected", "PRJNA1310568_185042")

def rd(p): return list(csv.DictReader(io.open(p, encoding="utf-8-sig"))) if os.path.exists(p) else []

def metrics(d):
    rt = rd(os.path.join(d, "read_tracking.csv"))
    m = dict(samples=0, merge=None, keep=None, n_asv=None, genus_r=None)
    if rt:
        tf = sum(float(x.get("filtered") or 0) for x in rt); tm = sum(float(x.get("merged") or 0) for x in rt)
        tin = sum(float(x.get("input") or 0) for x in rt); tn = sum(float(x.get("nonchim") or 0) for x in rt)
        m["samples"] = len(rt)
        if tf: m["merge"] = 100.0 * tm / tf
        if tin: m["keep"] = 100.0 * tn / tin
    ap = os.path.join(d, "asv_table.csv")
    if os.path.exists(ap):
        m["n_asv"] = sum(1 for _ in io.open(ap, encoding="utf-8-sig")) - 1
    tp = os.path.join(d, "asv_taxonomy.csv")
    if os.path.exists(tp):
        n = g = 0
        for row in csv.DictReader(io.open(tp, encoding="utf-8-sig")):
            n += 1; v = (row.get("Genus") or "").strip()
            if v and v.upper() != "NA": g += 1
        m["genus_r"] = 100.0 * g / n if n else None
    return m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abprime", default=os.path.join(R, "out_abprime"))
    ap.add_argument("--out", default=os.path.join(R, "out50k", "abprime_report.md"))
    a = ap.parse_args()
    rep = ["# A′ 裁决对照（预注册 §2）", "",
           "A 臂＝现行结果（out_stage2 / out50k）；A′ 臂＝truncLen 0,0 + maxEE_R2=5（out_abprime，含 2026-10-10 的 1310568 试点）。", "",
           "| 项目 | 臂 | 样本 | 合并率 | 最终保留率 | ASV 数 | 属级注释率 |",
           "| --- | --- | --- | --- | --- | --- | --- |"]
    verdicts = []
    for p in P5:
        for arm, d in (("A", os.path.join(R, "out_stage2", p)),
                       ("A'", PILOT if p == "PRJNA1310568" else os.path.join(a.abprime, p))):
            m = metrics(d)
            if not m["samples"] and not m["n_asv"]:
                rep.append("| %s | %s | - | - | - | - | - |" % (p, arm)); continue
            rep.append("| %s | %s | %d | %s | %s | %s | %s |" % (
                p, arm, m["samples"],
                ("%.1f%%" % m["merge"]) if m["merge"] is not None else "-",
                ("%.1f%%" % m["keep"]) if m["keep"] is not None else "-",
                m["n_asv"] if m["n_asv"] is not None else "-",
                ("%.1f%%" % m["genus_r"]) if m["genus_r"] is not None else "-"))
        ma = metrics(os.path.join(R, "out_stage2", p))
        mb = metrics(PILOT if p == "PRJNA1310568" else os.path.join(a.abprime, p))
        if ma["merge"] is not None and mb["merge"] is not None:
            d_m = mb["merge"] - ma["merge"]; d_k = (mb["keep"] or 0) - (ma["keep"] or 0)
            gok = (mb["genus_r"] is not None and ma["genus_r"] is not None and mb["genus_r"] >= ma["genus_r"] - 5 and mb["genus_r"] >= 60)
            if gok and d_m >= 5:
                v = "A′ 更好（合并率 +%.1fpp、保留率 +%.1fpp，属级不劣）→ 取 A′" % (d_m, d_k)
            elif abs(d_m) < 5:
                v = "两臂差异 <5pp → 维持 A（§2 规则）"
            else:
                v = "A′ 属级注释率不达标或下降 >5pp → 维持 A"
            verdicts.append("- %s：%s（Δ合并 %+.1fpp，Δ保留 %+.1fpp）" % (p, v, d_m, d_k))
    rep += ["", "## 逐项目判定（规则先写死，见 §2）", ""] + verdicts
    io.open(a.out, "w", encoding="utf-8").write("\n".join(rep) + "\n")
    print("wrote", a.out)
    for v in verdicts: print("  ", v)

if __name__ == "__main__":
    main()

import csv, io, os, glob
R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
def cols(p):
    return list(csv.reader(io.open(p, encoding="utf-8")))[0] if os.path.exists(p) else []
print("1) β：betadisper 有", "beta_all_report.md" if os.path.exists(os.path.join(OUT,"beta_all_report.md")) else "缺",
      "| 四距离:", all(os.path.exists(os.path.join(OUT,f)) for f in ("dist_unifrac_uw.csv","dist_unifrac_w.csv")))
print("2) α 列:", cols(os.path.join(OUT,"alpha_diversity.csv")))
print("   → 是否含 Chao1/ACE:", any("chao" in c.lower() or "ace" in c.lower() for c in cols(os.path.join(OUT,"alpha_diversity.csv"))))
ord_files = [os.path.basename(f) for f in glob.glob(os.path.join(OUT,"*")) if any(k in os.path.basename(f).lower() for k in ("pca","nmds","rda","cca","ordination","pcoa","hellinger"))]
print("3) 排序（PCA/NMDS/RDA/PCoA/Hellinger）产物:", ord_files if ord_files else "**无**")
bas = [os.path.basename(f) for f in glob.glob(os.path.join(OUT,"*")) if "baselga" in os.path.basename(f).lower() or "turnover" in os.path.basename(f).lower()]
print("4) Baselga/周转-嵌套分解产物:", bas if bas else "**无**")
sel = list(csv.reader(io.open(os.path.join(R,"meta","final45_selection.csv"), encoding="utf-8-sig")))[0]
print("5) 选集列（找 platform/cluster）:", [c for c in sel if c in ("cluster","platform_cluster","centre","title")])
ml = cols(os.path.join(OUT,"ml_report.md"))
print("6) ML 报告含校准/乐观度:", any(k in " ".join(ml).lower() for k in ("calibr","optimis","optimism")))
# 敏感性是否已传播到各模块
sens = {}
for f in ("beta_all_report.md","functional_faprotax.md","da_report.md","alpha_test_report.md","core_edge_report.md"):
    t = io.open(os.path.join(OUT,f), encoding="utf-8").read() if os.path.exists(os.path.join(OUT,f)) else ""
    sens[f] = ("1261944" in t, "一致属" in t or "agreed" in t or "consistent_genera" in t)
print("7) §18/§3 敏感性在各模块报告里出现情况:")
for k,v in sens.items(): print("   %-28s PRJNA1261944=%s  两法一致属=%s" % (k, v[0], v[1]))

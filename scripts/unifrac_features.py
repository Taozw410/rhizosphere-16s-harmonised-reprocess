# -*- coding: utf-8 -*-
"""UniFrac (unweighted + weighted) in the section-17 FEATURE space, vectorised with numpy.

The feature tree (NJ on DECIPHER distances) is unrooted: midpoint rooting is done here on the edge
table (BFS for the diameter, then root at its midpoint, splitting an edge if necessary), because
ape::midpoint() lives in phangorn, which is not installed.

Inputs : out_analysis/feature_matrix.csv, out_analysis/tree_edges.csv
Outputs: out_analysis/dist_unifrac_uw.csv, out_analysis/dist_unifrac_w.csv
usage: python unifrac_features.py
"""
import csv, io, os, sys, collections
import numpy as np

R = r"D:\CodexWork\mic-re"; OUT = os.path.join(R, "out_analysis")
META = ["project","base","crop","sample","all_reads","bacteria_reads","organelle_reads",
        "genus_reads","family_reads","unassigned_reads","n_asv"]
MINDEPTH = 1000


def load_matrix():
    rd = csv.reader(io.open(os.path.join(OUT, "feature_matrix.csv"), encoding="utf-8"))
    hdr = next(rd); feats = [c for c in hdr if c not in META]; gi = {c: i for i, c in enumerate(hdr)}
    meta, X = [], []
    for r in rd:
        if not r or int(float(r[gi["bacteria_reads"]] or 0)) < MINDEPTH:
            continue
        meta.append((r[gi["project"]], r[gi["crop"]], r[gi["sample"]]))
        X.append([float(r[gi[f]] or 0) for f in feats])
    return feats, meta, np.asarray(X, dtype=float)


def read_edges():
    adj, elen, tip = collections.defaultdict(list), {}, {}
    for r in csv.DictReader(io.open(os.path.join(OUT, "tree_edges.csv"), encoding="utf-8")):
        p, c, L = int(r["parent"]), int(r["child"]), float(r["length"])
        adj[p].append(c); adj[c].append(p)
        elen[(min(p, c), max(p, c))] = L
        if r["tip_label"]:
            tip[c] = r["tip_label"]
    return adj, elen, tip


def bfs(start, adj, elen):
    dist = {start: 0.0}; par = {start: None}; dq = collections.deque([start])
    while dq:
        u = dq.popleft()
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + elen[(min(u, v), max(u, v))]
                par[v] = u
                dq.append(v)
    return dist, par


def root_midpoint(adj, elen, tip):
    t0 = next(iter(tip))
    d1, _ = bfs(t0, adj, elen)
    a = max((n for n in d1 if n in tip), key=lambda n: d1[n])
    d2, par = bfs(a, adj, elen)
    b = max((n for n in d2 if n in tip), key=lambda n: d2[n])
    path = [b]
    while par[path[-1]] is not None:
        path.append(par[path[-1]])
    path.reverse()
    half = d2[b] / 2.0
    acc, root, split = 0.0, None, None
    for i in range(len(path) - 1):
        u, v = path[i], path[i + 1]; L = elen[(min(u, v), max(u, v))]
        if acc + L >= half:
            frac = (half - acc) / L if L else 0.0
            if frac <= 1e-9:
                root = u
            elif frac >= 1 - 1e-9:
                root = v
            else:
                split = (u, v, L, frac)
            break
        acc += L
    new_adj = collections.defaultdict(list)
    for k, vs in adj.items():
        new_adj[k] = list(vs)
    new_elen = dict(elen)
    if split:
        u, v, L, frac = split
        root = max(adj) + 1
        new_adj[u].remove(v); new_adj[v].remove(u)          # cut the edge that contains the midpoint
        new_adj[u].append(root); new_adj[root].append(u)
        new_adj[v].append(root); new_adj[root].append(v)
        new_elen[(min(u, root), max(u, root))] = L * frac
        new_elen[(min(v, root), max(v, root))] = L * (1 - frac)
        del new_elen[(min(u, v), max(u, v))]
    # DFS from the root -> directed edge list
    seen = {root}; stack = [root]; out = []
    while stack:
        u = stack.pop()
        for v in new_adj[u]:
            if v not in seen:
                seen.add(v)
                out.append((u, v, new_elen[(min(u, v), max(u, v))], tip.get(v, "")))
                stack.append(v)
    return out, len(seen), root


def descendants(edges):
    children = collections.defaultdict(list)
    for p, c, L, lab in edges:
        children[p].append((c, L, lab))
    desc = {}
    def walk(n):
        if n in desc:
            return desc[n]
        out = []
        for c, _, lab in children.get(n, []):
            if lab:
                out.append(lab)
            else:
                out.extend(walk(c))
        desc[n] = out
        return out
    sys.setrecursionlimit(100000)
    return [(L, walk(c)) for _, c, L, _ in edges]


def main():
    feats, meta, X = load_matrix()
    adj, elen, tip = read_edges()
    edges, n_nodes, root = root_midpoint(adj, elen, tip)
    E = descendants(edges)
    print("samples %d ; features %d ; tree tips %d ; nodes %d ; rooted at %s ; edges %d" %
          (len(meta), len(feats), len(tip), n_nodes, root, len(E)))
    fmap = {f: i for i, f in enumerate(feats)}
    # feature names with spaces were written to the FASTA with underscores (tax_bins/feature_reps.py),
    # so map the normalised form back to the original feature name.
    norm = {f.replace(" ", "_"): f for f in feats}
    L = np.array([e[0] for e in E])
    A = np.zeros((X.shape[0], len(E)))
    matched = 0
    for j, (_, cols) in enumerate(E):
        idx = [fmap[norm[t]] for t in cols if t in norm]
        if idx:
            A[:, j] = X[:, idx].sum(axis=1)
            matched += 1
    print("tree edges with at least one matched feature: %d / %d ; total matched tips: %d" %
          (matched, len(E), sum(1 for e in E for t in e[1] if t in norm)))
    if matched == 0:
        raise SystemExit("no tree tip matched a feature column - check name normalisation")
    B = (A > 0).astype(float)
    M = B * np.sqrt(L); G = M @ M.T; rs = (B * L).sum(axis=1)
    num = rs[:, None] + rs[None, :] - 2 * G
    den = rs[:, None] + rs[None, :] - G
    D_uw = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0); np.fill_diagonal(D_uw, 0.0)
    P = A / np.maximum(A.sum(axis=1, keepdims=True), 1)
    PL = P * L; rowsum = PL.sum(axis=1); D_w = np.zeros((X.shape[0],) * 2)
    for i0 in range(0, X.shape[0], 16):
        i1 = min(X.shape[0], i0 + 16)
        diff = np.abs(PL[i0:i1, None, :] - PL[None, :, :]).sum(axis=2)
        ssum = rowsum[i0:i1, None] + rowsum[None, :]
        D_w[i0:i1, :] = np.where(ssum > 0, diff / np.maximum(ssum, 1e-12), 0.0)
    D_w = (D_w + D_w.T) / 2.0; np.fill_diagonal(D_w, 0.0)
    print("unweighted UniFrac: min %.4f median %.4f max %.4f ; symmetric %s" %
          (D_uw[D_uw > 0].min(), np.median(D_uw[D_uw > 0]), D_uw.max(), np.allclose(D_uw, D_uw.T)))
    print("weighted   UniFrac: min %.4f median %.4f max %.4f" %
          (D_w[D_w > 0].min(), np.median(D_w[D_w > 0]), D_w.max()))
    for name, D in (("dist_unifrac_uw.csv", D_uw), ("dist_unifrac_w.csv", D_w)):
        with io.open(os.path.join(OUT, name), "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            names = ["%s|%s" % (m[0], m[2]) for m in meta]
            w.writerow([""] + names)
            for i, nm in enumerate(names):
                w.writerow([nm] + ["%.6f" % v for v in D[i]])
    print("wrote dist_unifrac_uw.csv / dist_unifrac_w.csv")


if __name__ == "__main__":
    main()

"""ContextDiffReg: condition-specific regulatory network inference and differential network analysis."""
from __future__ import annotations

from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.linear_model import ElasticNet, Ridge

METHODS = ("elasticnet", "ridge", "correlation")


def _to_dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x, dtype=float)


def infer_condition_network(expr, gene_names, prior_tf_to_targets, method="elasticnet", alpha=0.01, l1_ratio=0.5):
    """Infer a TF->target network for one condition.

    expr: Cell x gene matrix (log-normalised recommended).
    gene_names: gene symbols matching the columns of expr.
    prior_tf_to_targets: dict {TF: set(targets)}.
    method: elasticnet, ridge or correlation.
    """
    expr = _to_dense(expr)
    gene_names = list(gene_names)
    gene_idx = {g: i for i, g in enumerate(gene_names)}
    tfs = sorted(set(prior_tf_to_targets) & set(gene_idx))
    targets_union = {t for tg in prior_tf_to_targets.values() for t in tg} & set(gene_idx)
    rows = []

    for target in sorted(targets_union):
        ti = gene_idx[target]
        regs = [tf for tf in tfs if target in prior_tf_to_targets[tf]]
        if len(regs) == 0:
            continue
        x = expr[:, [gene_idx[r] for r in regs]]
        y = expr[:, ti]
        if x.shape[0] < 10 or np.std(y) < 1e-9:
            continue
        x = (x - x.mean(0)) / (x.std(0) + 1e-9)
        yc = (y - y.mean()) / (y.std() + 1e-9)

        if method == "correlation":
            coefs = np.corrcoef(x, yc, rowvar=False)[:-1, -1]
            coefs = np.nan_to_num(coefs)
        elif method == "ridge":
            model = Ridge(alpha=alpha)
            model.fit(x, yc)
            coefs = model.coef_
        else:
            model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=5000)
            model.fit(x, yc)
            coefs = model.coef_

        for reg, w in zip(regs, np.asarray(coefs).ravel()):
            if abs(w) > 1e-8 or method == "correlation":
                rows.append({"tf": reg, "target": target, "weight": float(w), "method": method})

    return pd.DataFrame(rows, columns=["tf", "target", "weight", "method"])


def differential_network(adata, prior_tf_to_targets, condition_key="condition", disease="disease", control="control",
                         method="elasticnet", n_bootstrap=100, seed=0):
    """Compare disease vs control networks and score differential edges.

    Returns edge table with tf, target, w_disease, w_control, delta, stability, significant.
    """
    rng = np.random.default_rng(seed)
    expr = _to_dense(adata.X)
    genes = list(adata.var_names)
    obs = adata.obs
    dis_mask = (obs[condition_key] == disease).values
    con_mask = (obs[condition_key] == control).values
    if dis_mask.sum() < 10 or con_mask.sum() < 10:
        raise ValueError("Both conditions need >=10 cells")

    w_dis = infer_condition_network(expr[dis_mask], genes, prior_tf_to_targets, method=method)
    w_con = infer_condition_network(expr[con_mask], genes, prior_tf_to_targets, method=method)

    base = w_dis.merge(w_con, on=["tf", "target"], how="outer", suffixes=("_dis", "_con"))
    base["w_disease"] = base["weight_dis"].fillna(0.0)
    base["w_control"] = base["weight_con"].fillna(0.0)
    base["delta"] = base["w_disease"] - base["w_control"]
    base = base[["tf", "target", "w_disease", "w_control", "delta"]].copy()

    boot = np.zeros(len(base))
    for _ in range(n_bootstrap):
        di = rng.integers(0, dis_mask.sum(), size=dis_mask.sum())
        ci = rng.integers(0, con_mask.sum(), size=con_mask.sum())
        bd = infer_condition_network(expr[dis_mask][di], genes, prior_tf_to_targets, method=method)
        bc = infer_condition_network(expr[con_mask][ci], genes, prior_tf_to_targets, method=method)
        bm = bd.merge(bc, on=["tf", "target"], how="outer", suffixes=("_d", "_c"))
        bm["d"] = bm["weight_d"].fillna(0.0) - bm["weight_c"].fillna(0.0)
        sign_map = dict(zip(zip(bm["tf"], bm["target"]), np.sign(bm["d"])))
        for i, (tf, tg) in enumerate(zip(base["tf"], base["target"])):
            s = sign_map.get((tf, tg), 0.0)
            if s != 0 and np.sign(s) == np.sign(base["delta"].iloc[i]):
                boot[i] += 1

    base["stability"] = boot / n_bootstrap
    base["significant"] = base["stability"] >= 0.9
    base = base.sort_values("delta", key=lambda s: s.abs(), ascending=False).reset_index(drop=True)
    return base


def rank_regulators(edge_table, prior_tf_to_targets=None, deg_genes=None, all_genes=None):
    """Rank differential regulators by network rewiring and downstream DEG enrichment."""
    sig = edge_table[edge_table["significant"]].copy() if "significant" in edge_table.columns else edge_table.copy()
    agg = sig.groupby("tf").agg(
        n_edges=("target", "size"),
        delta_sum=("delta", lambda s: s.abs().sum()),
        up_edges=("delta", lambda s: (s > 0).sum()),
        down_edges=("delta", lambda s: (s < 0).sum()),
    ).reset_index()

    if deg_genes is not None and all_genes is not None:
        from scipy.stats import fisher_exact
        deg_set = set(deg_genes)
        all_set = set(all_genes)
        rows = []
        for tf, g in sig.groupby("tf"):
            targets = set(g["target"])
            a = len(targets & deg_set)
            b = len(targets) - a
            c = len(deg_set - targets)
            d = len(all_set - targets - deg_set)
            odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({"tf": tf, "deg_p": p, "deg_odds": odds})
        if rows:
            agg = agg.merge(pd.DataFrame(rows), on="tf", how="left")

    agg["score"] = agg["delta_sum"] * agg["n_edges"]
    return agg.sort_values("score", ascending=False).reset_index(drop=True)
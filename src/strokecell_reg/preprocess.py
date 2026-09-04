"""Normalization, HVG, neighbors and pseudobulk aggregation."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc


def normalize_hvg(adata, target_sum: float = 1e4, n_top_genes: int = 2000,
                  flavor: str = "seurat_v3", batch_key=None):
    """Normalize, log-transform and select highly variable genes."""
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor, batch_key=batch_key)
    return adata


def build_neighbors(adata, n_pcs: int = 30, batch_key=None):
    """Scale HVG, run PCA and nearest neighbors (optional Harmony)."""
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=min(n_pcs, adata_hvg.n_vars - 1), svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    basis = "X_pca"
    if batch_key is not None:
        sc.external.pp.harmony_integrate(adata, key=batch_key, basis="X_pca", adjusted_basis="X_pca_harmony")
        basis = "X_pca_harmony"
    sc.pp.neighbors(adata, use_rep=basis, n_neighbors=15, n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]))
    return adata


def aggregate_pseudobulk(adata, sample_key: str = "sample", celltype_key: str = "cell_type",
                         condition_key: str = "condition"):
    """Build pseudobulk counts: rows = (sample, cell_type) combinations."""
    obs = adata.obs
    keys = [sample_key, celltype_key]
    df = pd.DataFrame(adata.X.toarray(), index=obs.index, columns=adata.var_names)
    df[sample_key] = obs[sample_key].values
    df[celltype_key] = obs[celltype_key].values
    pb = df.groupby(keys, observed=True).sum()
    meta = pb.index.to_frame(index=False)
    if condition_key in obs.columns:
        cond = obs[[sample_key, condition_key]].drop_duplicates().set_index(sample_key)[condition_key]
        meta[condition_key] = meta[sample_key].map(cond)
    return pb, meta

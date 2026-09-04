"""Quality control: metrics, filtering and doublet scores."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def compute_qc(adata, mito_prefix: str = "mt-") -> None:
    """Add QC metrics to adata.obs in place."""
    if sp.issparse(adata.X):
        adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
        adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    else:
        adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
        adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    mito = adata.var_names.astype(str).str.startswith(mito_prefix)
    if mito.sum() > 0:
        mito_counts = np.asarray(adata[:, mito].X.sum(axis=1)).ravel()
        adata.obs["pct_mito"] = mito_counts / np.maximum(adata.obs["n_counts"].values, 1) * 100
    else:
        adata.obs["pct_mito"] = 0.0
    adata.obs["log1p_n_counts"] = np.log1p(adata.obs["n_counts"].values)
    adata.obs["log1p_n_genes"] = np.log1p(adata.obs["n_genes"].values)


def filter_cells_genes(adata, min_genes: int = 200, max_genes: int = 8000,
                       min_counts: int = 500, max_pct_mito: float = 20.0, min_cells: int = 3):
    """Filter cells and genes. Returns a filtered copy."""
    compute_qc(adata)
    keep_cells = (
        (adata.obs["n_genes"] >= min_genes)
        & (adata.obs["n_genes"] <= max_genes)
        & (adata.obs["n_counts"] >= min_counts)
        & (adata.obs["pct_mito"] <= max_pct_mito)
    )
    keep_genes = np.asarray((adata.X > 0).sum(axis=0)).ravel() >= min_cells
    return adata[keep_cells.values, keep_genes].copy()


def doublet_scores(adata, n_pcs: int = 30):
    """Scrublet-style doublet score (PCA + kNN density)."""
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    x = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    x = np.log1p(x)
    pca = PCA(n_components=min(n_pcs, min(x.shape))).fit_transform(x)
    knn = NearestNeighbors(n_neighbors=min(50, max(x.shape[0] - 1, 2))).fit(pca)
    dist, _ = knn.kneighbors(pca)
    return np.median(dist, axis=1)

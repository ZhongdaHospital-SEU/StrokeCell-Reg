"""Data loading utilities for 10x Genomics count matrices and sample directories."""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad


def _open_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_10x(sample_dir: str, sample_name: str, prefix: str | None = None, use_var_names: str = "gene_symbols") -> ad.AnnData:
    """Load one 10x sample (barcodes.tsv.gz / genes-or-features.tsv.gz / matrix.mtx.gz)."""
    files = {f.name: f for f in Path(sample_dir).iterdir() if f.is_file()}
    matrix = None
    features = None
    barcodes = None
    if prefix:
        matrix = str(Path(sample_dir) / (prefix + "_matrix.mtx.gz"))
        if not Path(matrix).exists():
            matrix = str(Path(sample_dir) / (prefix + "_matrix.mtx"))
        barc = list(Path(sample_dir).glob(prefix + "_*barcodes.tsv*"))
        feat = list(Path(sample_dir).glob(prefix + "_*features.tsv*")) or list(Path(sample_dir).glob(prefix + "_*genes.tsv*"))
        if barc:
            barcodes = str(barc[0])
        if feat:
            features = str(feat[0])
    else:
        for name, f in files.items():
            base = name[:-3] if name.endswith(".gz") else name
            if base.endswith("matrix.mtx"):
                matrix = str(f)
            elif base.endswith(("genes.tsv", "features.tsv")):
                features = str(f)
            elif base.endswith("barcodes.tsv"):
                barcodes = str(f)
    if matrix is None or barcodes is None or features is None:
        raise FileNotFoundError("Could not find all 10x files in " + sample_dir)

    mat = scipy.io.mmread(matrix)
    if sp.isspmatrix(mat):
        mat = mat.tocsr()
    else:
        mat = sp.csr_matrix(mat)

    bc = pd.read_csv(_open_maybe_gz(barcodes), sep="\t", header=None)[0].astype(str)
    bc.name = None
    ft = pd.read_csv(_open_maybe_gz(features), sep="\t", header=None)
    var_ids = ft[0].astype(str).values
    var_symbols = ft[1].astype(str).values if ft.shape[1] >= 2 else var_ids.copy()
    var = pd.DataFrame(index=var_ids)
    var["gene_symbol"] = var_symbols

    if use_var_names == "gene_symbols":
        var.index = var["gene_symbol"]
        var.index.name = None
        keep = ~var.index.duplicated(keep="first")
        var = var[keep]
    else:
        keep = pd.Series(True, index=var.index)

    if mat.shape == (ft.shape[0], len(bc)):
        mat = mat.T
    if mat.shape != (len(bc), ft.shape[0]):
        raise ValueError("Shape mismatch after loading")
    keep_arr = np.asarray(keep.values) if hasattr(keep, "values") else np.asarray(keep, dtype=bool)
    mat = mat[:, keep_arr]

    adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=bc.astype(str)), var=var)
    adata.obs["sample"] = sample_name
    adata.obs_names_make_unique()
    return adata


def build_anndata_from_dirs(sample_spec: Dict[str, str], condition_map: Optional[Dict[str, str]] = None) -> ad.AnnData:
    """Concatenate multiple 10x sample directories into one AnnData."""
    adatas = [load_10x(d, s) for s, d in sample_spec.items()]
    adata = ad.concat(adatas, join="outer", fill_value=0, index_unique="-")
    if condition_map is not None:
        adata.obs["condition"] = adata.obs["sample"].map(condition_map)
    adata.var_names_make_unique()
    return adata


def write_h5ad(adata: ad.AnnData, path: str) -> None:
    adata.write_h5ad(path)


def read_h5ad(path: str) -> ad.AnnData:
    return ad.read_h5ad(path)

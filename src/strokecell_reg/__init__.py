"""StrokeCell-Reg: cell-type-specific regulatory reprogramming analysis for stroke scRNA-seq."""

__version__ = "0.1.0"

from .io import load_10x, build_anndata_from_dirs, read_h5ad, write_h5ad
from .qc import compute_qc, filter_cells_genes, doublet_scores
from .preprocess import normalize_hvg, build_neighbors, aggregate_pseudobulk
from .strokeregdb import build_strokeregdb, load_trrust, load_regnetwork, load_dorothea_tsv, prior_to_dict
from .contextdiffreg import infer_condition_network, differential_network, rank_regulators
from .enrich import enrichr, go_enrichment, kegg_enrichment, gsea_prerank
from .report import render_html, figure_block, table_block

__all__ = [
    "__version__",
    "load_10x", "build_anndata_from_dirs", "read_h5ad", "write_h5ad",
    "compute_qc", "filter_cells_genes", "doublet_scores",
    "normalize_hvg", "build_neighbors", "aggregate_pseudobulk",
    "build_strokeregdb", "load_trrust", "load_regnetwork", "load_dorothea_tsv", "prior_to_dict",
    "infer_condition_network", "differential_network", "rank_regulators",
    "enrichr", "go_enrichment", "kegg_enrichment", "gsea_prerank",
    "render_html", "figure_block", "table_block",
]
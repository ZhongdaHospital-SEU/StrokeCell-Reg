"""Command-line interface for StrokeCell-Reg."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def _load_sample_spec(path: str) -> dict:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if "samples" in data:
        return {s["name"]: s["dir"] for s in data["samples"]}, data.get("conditions", {})
    return data.get("sample_dirs", {}), data.get("conditions", {})


def cmd_preprocess(args):
    from .io import build_anndata_from_dirs
    from .qc import filter_cells_genes
    from .preprocess import normalize_hvg

    sample_spec, cond = _load_sample_spec(args.samples)
    adata = build_anndata_from_dirs(sample_spec, condition_map=cond)
    adata = filter_cells_genes(adata, min_genes=args.min_genes, max_genes=args.max_genes,
                               min_counts=args.min_counts, max_pct_mito=args.max_pct_mito)
    adata = normalize_hvg(adata, n_top_genes=args.n_top_genes)
    adata.write_h5ad(args.output)
    print(f"Wrote {args.output}: {adata.n_obs} cells x {adata.n_vars} genes")


def cmd_diffreg(args):
    import pandas as pd
    import anndata as ad
    from .contextdiffreg import differential_network, rank_regulators
    from .strokeregdb import prior_to_dict

    adata = ad.read_h5ad(args.input)
    if args.cell_type not in adata.obs.columns:
        raise ValueError("cell_type column not found")
    sub = adata[adata.obs[args.cell_type] == args.cell_type_label].copy()
    prior = pd.read_csv(args.prior, sep="\t", dtype=str)
    prior_tf = prior_to_dict(prior)
    edges = differential_network(sub, prior_tf, condition_key=args.condition_key,
                                 disease=args.disease, control=args.control,
                                 method=args.method, n_bootstrap=args.n_bootstrap)
    edges.to_csv(args.output, index=False)
    print(f"Wrote {args.output}: {len(edges)} edges, {int(edges['significant'].sum())} significant")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="strokecell-reg", description="StrokeCell-Reg analysis toolkit")
    sub = p.add_subparsers(dest="command", required=True)

    pp = sub.add_parser("preprocess", help="Build a filtered, normalized AnnData")
    pp.add_argument("--samples", required=True, help="JSON manifest of sample dirs and conditions")
    pp.add_argument("--output", required=True)
    pp.add_argument("--min-genes", type=int, default=200)
    pp.add_argument("--max-genes", type=int, default=8000)
    pp.add_argument("--min-counts", type=int, default=500)
    pp.add_argument("--max-pct-mito", type=float, default=20.0)
    pp.add_argument("--n-top-genes", type=int, default=2000)
    pp.set_defaults(func=cmd_preprocess)

    dr = sub.add_parser("diffreg", help="Run ContextDiffReg for one cell type")
    dr.add_argument("--input", required=True, help="AnnData with cell_type and condition obs")
    dr.add_argument("--prior", required=True, help="StrokeRegDB prior TSV")
    dr.add_argument("--cell-type", default="cell_type", help="obs column holding cell-type labels")
    dr.add_argument("--cell-type-label", required=True)
    dr.add_argument("--condition-key", default="condition")
    dr.add_argument("--disease", default="disease")
    dr.add_argument("--control", default="control")
    dr.add_argument("--method", default="elasticnet", choices=["elasticnet", "ridge", "correlation"])
    dr.add_argument("--n-bootstrap", type=int, default=100)
    dr.add_argument("--output", required=True)
    dr.set_defaults(func=cmd_diffreg)

    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
"""Functional enrichment via Enrichr (GO / KEGG / Reactome) and GSEA."""
from __future__ import annotations

import pandas as pd


def enrichr(gene_list, gene_set_library, top_n=20):
    """Run Enrichr enrichment for one library. Returns a tidy DataFrame."""
    import gseapy as gp
    res = gp.enrichr(gene_list=list(gene_list), gene_sets=gene_set_library, organism="mouse", outdir=None)
    df = res.results.copy()
    df = df.head(top_n)
    df["-log10_pval"] = -df["Adjusted P-value"].apply(lambda x: __import__("math").log10(max(float(x), 1e-300)))
    return df


def go_enrichment(gene_list, top_n=20):
    return enrichr(gene_list, "GO_Biological_Process_2023", top_n=top_n)


def kegg_enrichment(gene_list, top_n=20):
    return enrichr(gene_list, "KEGG_2021_Mouse", top_n=top_n)


def gsea_prerank(ranked_gene_scores, gene_set_library, outdir, seed=0, threads=4):
    """Run GSEA preranked on a Series/DataFrame of gene scores."""
    import gseapy as gp
    if isinstance(ranked_gene_scores, pd.Series):
        rnk = ranked_gene_scores.copy()
    else:
        rnk = ranked_gene_scores.iloc[:, 0]
    rnk = rnk[~rnk.index.duplicated(keep="first")].sort_values(ascending=False)
    res = gp.prerank(rnk=rnk, gene_sets=gene_set_library, outdir=outdir, seed=seed, threads=threads)
    return res.res2d
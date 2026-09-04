"""StrokeRegDB: build a stroke-weighted TF-target prior regulatory database.

Integrates DoRothEA, TRRUST and RegNetwork into a harmonised directed
TF -> target network with confidence levels and optional stroke relevance weights.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Optional

import pandas as pd

# Well-established stroke-relevant transcription factors (curated from literature).
STROKE_TFS = {
    "NFKB1", "RELA", "RELB", "STAT1", "STAT3", "STAT5A", "STAT5B",
    "HIF1A", "EPAS1", "SP1", "CEBPB", "CEBPA", "IRF1", "IRF3", "IRF7",
    "JUN", "FOS", "FOSL1", "FOSL2", "TP53", "PPARG", "NR3C1", "FOXO3",
    "CREB1", "ATF3", "ATF4", "EGR1", "EGR2", "KLF4", "KLF2", "MAF",
    "MAFB", "MAFF", "SRF", "MEF2C", "MEF2A", "SOX2", "SOX9", "OLIG2",
    "NFE2L2", "SMAD3", "SMAD4", "TCF4", "REST", "YAP1", "TEAD1", "MYC",
}

DO_LEVELS = {"A": 4, "B": 3, "C": 2, "D": 1, "E": 0}


def load_trrust(path: str, species: str = "mouse") -> pd.DataFrame:
    """Parse TRRUST raw TSV (headerless: TF, Target, Mode, PMID) into a long table."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
    out = pd.DataFrame({
        "tf": df[0].astype(str).str.upper(),
        "target": df[1].astype(str).str.upper(),
    })
    out["mode"] = df[2].astype(str) if df.shape[1] >= 3 else ""
    out["pmid"] = df[3].astype(str) if df.shape[1] >= 4 else ""
    out["confidence"] = 3
    out["source"] = "TRRUST"
    out["species"] = species
    return out


def load_regnetwork(path: str, species: str = "mouse") -> pd.DataFrame:
    """Parse a RegNetwork *.source file (cols: TF_symbol, TF_entrez, target_symbol, target_entrez)."""
    df = pd.read_csv(path, sep="\t", header=None, dtype=str, comment="#")
    ncol = df.shape[1]
    if ncol >= 3:
        out = pd.DataFrame({
            "tf": df[0].astype(str).str.upper(),
            "target": df[2].astype(str).str.upper(),
        })
        out["evidence"] = ""
    else:
        raise ValueError("RegNetwork file has fewer than 3 columns")
    out["mode"] = ""
    out["confidence"] = 2
    out["source"] = "RegNetwork"
    out["species"] = species
    return out


def load_dorothea_tsv(path: str, species: str = "mouse") -> pd.DataFrame:
    """Parse an exported DoRothEA TSV (columns tf, confidence, target, mor)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]
    out = pd.DataFrame({
        "tf": df["tf"].astype(str).str.upper(),
        "target": df["target"].astype(str).str.upper(),
    })
    out["confidence"] = df["confidence"].map(DO_LEVELS).fillna(2).astype(int)
    out["mode"] = df["mor"].astype(str) if "mor" in df.columns else ""
    out["source"] = "DoRothEA"
    out["species"] = species
    return out


def build_strokeregdb(
    sources: Iterable[pd.DataFrame],
    stroke_tfs: Optional[Iterable[str]] = None,
    min_confidence: int = 1,
) -> pd.DataFrame:
    """Merge curated sources into one deduplicated, weighted prior network.

    Returns a DataFrame with columns: tf, target, confidence, n_sources,
    stroke_tf, stroke_weight, mode.
    """
    frames = []
    for df in sources:
        need = ["tf", "target", "confidence", "source"]
        missing = [c for c in need if c not in df.columns]
        if missing:
            raise ValueError("Source missing columns: " + ", ".join(missing))
        frames.append(df[["tf", "target", "confidence", "source", "mode"]].copy())

    merged = pd.concat(frames, ignore_index=True)
    merged = merged[merged["confidence"] >= min_confidence]
    merged = merged.dropna(subset=["tf", "target"])

    agg = merged.groupby(["tf", "target"], as_index=False).agg(
        confidence=("confidence", "max"),
        n_sources=("source", "nunique"),
        modes=("mode", lambda s: ";".join(sorted({str(x) for x in s if x and str(x) != "nan"}))),
    )
    stroke_tfs = {str(x).upper() for x in (stroke_tfs or STROKE_TFS)}
    agg["stroke_tf"] = agg["tf"].isin(stroke_tfs)
    # stroke weight: TF in curated stroke set => x2, otherwise x1, scaled by confidence
    agg["stroke_weight"] = (agg["confidence"] + 1) * agg["stroke_tf"].map({True: 2, False: 1})
    agg = agg.sort_values(["stroke_weight", "n_sources", "confidence"], ascending=False).reset_index(drop=True)
    return agg


def prior_to_dict(prior: pd.DataFrame) -> Dict[str, set]:
    """Convert prior table to {tf: set(targets)}."""
    d: Dict[str, set] = {}
    for tf, tg in zip(prior["tf"], prior["target"]):
        d.setdefault(tf, set()).add(tg)
    return d

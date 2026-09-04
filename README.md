# StrokeCell-Reg

**StrokeCell-Reg** is a Python toolkit for decoding cell-type-specific gene
regulatory reprogramming in ischemic stroke from single-cell transcriptomics.
It couples a stroke-weighted transcription-factor target prior database,
**StrokeRegDB**, with **ContextDiffReg**, a module that infers disease-versus
control regulatory networks per cell type and scores differential edges by
bootstrap stability.

The toolkit was validated on synthetic single-cell data, on a mouse
middle-cerebral-artery-occlusion atlas, on six independent mouse stroke
cohorts, and on human intracerebral-hemorrhage immune cells.

The public repository is
[https://github.com/ZhongdaHospital-SEU/StrokeCell-Reg](https://github.com/ZhongdaHospital-SEU/StrokeCell-Reg).
The package source, the editable rendered figure set, and the StrokeRegDB
description are distributed here. Raw sequencing data are not included and
must be downloaded from GEO.

## Features

- Quality control, normalization, batch correction, and cell-type annotation
- Pseudobulk differential expression per cell type (DESeq2-style model)
- **StrokeRegDB**: a stroke-weighted TF-target prior integrating DoRothEA,
  TRRUST, and RegNetwork (mouse and human)
- **ContextDiffReg**: condition-specific regulatory-network inference with
  ridge / elastic-net / correlation variants and differential-network analysis
  with bootstrap stability scoring
- Functional enrichment, cell-cell communication, and a one-click HTML report
- Prior-free de novo regulon analysis, motif/ChIP anchoring, and perturbation
  validation

## Install (development)

`ash
pip install -e .
`

The package requires Python >= 3.10.

## Quick start (CLI)

`ash
# Preprocess raw matrices into an annotated AnnData object
strokecell-reg preprocess --samples samples.json --output qc.h5ad

# Run cell-type-specific differential regulatory analysis
strokecell-reg diffreg --input qc.h5ad --prior strokeregdb.tsv \
  --cell-type-label Microglia --output edges.tsv
`

## Figures

The editable, journal-compatible figure set used in the study is under
outputs/figures/ in SVG format.

## Data availability

All single-cell datasets are publicly available in the Gene Expression Omnibus:

- GSE174574, GSE225948, GSE142445, GSE197731, GSE245386, GSE319238 (mouse stroke)
- GSE166638 (human intracerebral hemorrhage)
- GSE162526 (PU.1 perturbation), GSE266422 (IRF8 knockout), GSE220041 (stroke microglia ATAC-seq)

Raw matrices are not distributed in this repository; download them from GEO.
**StrokeRegDB** is large and is distributed separately; contact the corresponding
author for access or a release link.

## Citation

If you use StrokeCell-Reg, please cite the method once the associated
publication is available.

## License

MIT. See LICENSE.

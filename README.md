# Soft Cell-Cell Communication via CellWalker-Weighted Ligand-Receptor Scoring

**Status:** work in progress

## Motivation

Ligand-receptor communication tools (CellChat, LIANA, CellPhoneDB) require every
cell to be assigned a single, discrete cell-type label before communication
scores are computed. In tumor microenvironments, differentiation trajectories,
and other transitional biology, that discreteness is a fiction — many cells sit
between identities (e.g. the partial-EMT malignant cells described in
[Puram et al. 2017](https://doi.org/10.1016/j.cell.2017.10.044)), and a hard
`argmax` label discards exactly the information that makes those cells
interesting.

[CellWalker / CellWalkR](https://github.com/PFPrzytycki/CellWalkR) already
computes what we need: a graph-diffusion-based **probability distribution**
over cell types for every cell, rather than a hard label. This project uses
that distribution directly as a weighting scheme for ligand-receptor scoring,
producing a **soft communication network** where a cell's contribution to each
cell type's signaling profile is proportional to how confidently it belongs to
that type — instead of all-or-nothing.

**Question:** do identity-ambiguous cells (e.g. p-EMT malignant cells) show
signaling behavior that is invisible to hard-clustered communication analysis?

## Approach

1. Run CellWalkR on a tumor scRNA-seq dataset using published marker gene sets
   to get, for every cell, a full probability vector across cell types (not
   just the top label).
2. Compute a per-cell **identity entropy** score from that distribution to
   flag ambiguous cells.
3. Implement a probability-weighted ligand-receptor scoring function
   (`R/soft_weighting.R`) that computes expected per-type ligand/receptor
   expression as a probability-weighted average across cells, instead of a
   mean within a hard cluster.
4. Build the equivalent **hard** network (standard argmax-label CellChat/LIANA
   style scoring) as a baseline.
5. Compare soft vs. hard networks: which edges appear only in the soft
   version, and do they concentrate around high-entropy (ambiguous) cells?

## Dataset

[GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322) —
Puram et al. 2017, head & neck squamous cell carcinoma scRNA-seq (~6,000
cells: malignant, fibroblast, T cell, B cell, myeloid, endothelial, mast).
Chosen because the original paper's p-EMT signature gives an independent way
to validate that CellWalker's ambiguity scores are picking up something real.

Backup dataset: GSE72056 (Tirosh et al. 2016, melanoma) — simpler logistics,
same cell-type diversity.

## Repo structure

```
R/                 Reusable functions (the actual novel contribution)
  soft_weighting.R   Probability-weighted L-R scoring + hard baseline + comparison
  utils.R            Entropy scoring, marker gene helpers, plotting helpers
scripts/            Ordered analysis pipeline
  01_download_data.R
  02_run_cellwalker.R
  03_soft_communication.R
  04_compare_networks.R
vignettes/          Worked walkthrough (getting-started.Rmd)
data/               Raw + processed data (gitignored except for small refs)
```

## Setup

```r
install.packages(c("Seurat", "Matrix", "igraph", "ggplot2", "GEOquery", "devtools"))
devtools::install_github("PFPrzytycki/CellWalkR", ref = "cellwalker2")
# LIANA gives a convenient curated ligand-receptor resource, even though we
# bypass its scoring engine:
devtools::install_github("saezlab/liana")
```

Then run `scripts/01_download_data.R` through `04_compare_networks.R` in order,
or step through `vignettes/getting-started.Rmd`.

## Status / TODO

- [ ] Confirm CellWalkR install path (CRAN version predates hierarchical
      features used here — need `cellwalker2` branch)
- [ ] Finalize marker gene sets per cell type (paper supplementary table vs.
      PanglaoDB cross-check)
- [ ] Implement permutation-based significance for soft LR scores
- [ ] p-EMT score correlation with entropy
- [ ] Figures for README (hard vs soft network diagrams)

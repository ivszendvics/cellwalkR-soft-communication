# 01_download_data.R
# Downloads GSE103322 (Puram et al. 2017, HNSCC scRNA-seq) from GEO.
# Backup dataset: GSE72056 (Tirosh et al. 2016, melanoma) -- swap the
# accession below if GSE103322 processed files turn out to be awkward
# to parse (GEO supplementary file formats vary a lot between series).

library(GEOquery)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

accession <- "GSE103322"

# GEOquery::getGEOSuppFiles() pulls the raw supplementary files (usually a
# processed counts/TPM matrix + a metadata/annotation file for this series).
# Inspect what comes down before assuming a format -- Puram et al. distribute
# a single tab-delimited matrix with cell type + tumor annotations baked into
# the first few rows, which needs manual splitting out (see comments below).
supp <- getGEOSuppFiles(accession, baseDir = "data/raw", makeDirectory = TRUE)
print(supp)

# --- Manual step you'll need after downloading ---
# Puram et al.'s matrix (GSE103322_HNSCC_all_data.txt.gz) has genes as rows,
# cells as columns, but the FIRST FEW ROWS are metadata (non-malignant cell
# type, malignant/non-malignant classification, patient ID, etc.) rather than
# gene expression. You'll need to:
#   1. Read the file with read.delim(..., check.names = FALSE)
#   2. Split off the top metadata rows into a separate cell-metadata data.frame
#   3. Coerce the remaining rows to a numeric expression matrix
#   4. Build a Seurat object from the matrix + metadata
#
# Example skeleton (uncomment and adjust once you've inspected the file):
#
# raw <- read.delim("data/raw/GSE103322/GSE103322_HNSCC_all_data.txt.gz",
#                    check.names = FALSE, row.names = 1)
# meta_rows <- c("non-malignant cell type", "malignant(1=no,2=yes)",
#                "tumor site", "patient")   # confirm exact row names in file
# cell_meta <- t(raw[rownames(raw) %in% meta_rows, ])
# expr <- as.matrix(raw[!rownames(raw) %in% meta_rows, ])
# storage.mode(expr) <- "numeric"
#
# saveRDS(list(expr = expr, meta = cell_meta), "data/gse103322_parsed.rds")

message("Downloaded to data/raw/", accession,
        ". Inspect the file structure before running 02_run_cellwalker.R -- ",
        "see comments in this script for the known quirks of this series.")

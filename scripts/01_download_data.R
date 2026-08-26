# 01_download_data.R
# Downloads GSE103322 (Puram et al. 2017, HNSCC scRNA-seq) from GEO, then
# parses the single matrix file GEO provides into an expression matrix +
# cell metadata data.frame.
# Backup dataset: GSE72056 (Tirosh et al. 2016, melanoma) -- swap the
# accession below if GSE103322 processed files turn out to be awkward
# to parse (GEO supplementary file formats vary a lot between series).

library(GEOquery)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

accession <- "GSE103322"
supp <- getGEOSuppFiles(accession, baseDir = "data/raw", makeDirectory = TRUE)
print(supp)

# --- Parse Puram et al.'s matrix ---
# GSE103322_HNSCC_all_data.txt.gz has genes as rows, cells as columns, but
# the first 6 lines are metadata, not expression (confirmed by inspecting
# the raw file):
#   1: header row of cell barcodes
#   2: "processed by Maxima enzyme" (always 1, not informative)
#   3: "Lymph node" (0/1: primary tumor vs. lymph node metastasis)
#   4: "classified  as cancer cell" (0/1 -- note the double space in the
#      literal label, a quirk of the file itself)
#   5: "classified as non-cancer cells" (0/1)
#   6: "non-cancer cell type" (string label for non-malignant cells, "0"
#      for malignant/unclassified cells)
# Rows 7+ are gene expression, already log2(TPM/10 + 1) normalized
# (Tirosh-lab convention for this era) -- NOT raw counts. Gene names carry
# literal single-quote characters in the file (e.g. 'RPS11') that need
# stripping. ~2,215 cells classified malignant and ~3,363 non-malignant
# out of 5,902 total match the paper's reported numbers exactly, which is
# a good sanity check that this parsing is correct.
#
# Patient ID is only present in the cell barcode itself (e.g.
# "HN28_P15_D06_S330_comb", "HNSCC26_P24_H05_S377_comb" -- both "HN" and
# "HNSCC" prefixes are used inconsistently for the same patients). A
# separate pool of 109 cells uses "HNSCC_combo1_HNSCC_ComboP1_*" barcodes
# with no extractable per-patient ID; these are kept with patient = "combo"
# since patient identity isn't used anywhere downstream of this script.

f <- file.path("data/raw", accession,
                paste0(accession, "_HNSCC_all_data.txt.gz"))

meta_lines <- readLines(gzfile(f), n = 6)
header <- strsplit(meta_lines[1], "\t")[[1]]
cell_names <- header[-1]

meta_fields <- lapply(meta_lines[2:6], function(l) trimws(strsplit(l, "\t")[[1]][-1]))

non_malignant_cell_type_raw <- meta_fields[[5]]
non_malignant_cell_type <- non_malignant_cell_type_raw
non_malignant_cell_type[non_malignant_cell_type == "-Fibroblast"] <- "Fibroblast"
non_malignant_cell_type[non_malignant_cell_type == "0"] <- NA

# Map the paper's non-malignant vocabulary onto this project's 7-type
# schema (matches marker_list in scripts/build_marker_genes.R). Macrophage
# and Dendritic are pooled into "Myeloid"; "myocyte" (19 cells) has no
# equivalent in the schema and is kept as its own label rather than
# silently dropped or miscoded.
type_map <- c(
  "Fibroblast" = "Fibroblast", "T cell" = "Tcell", "Endothelial" = "Endothelial",
  "B cell" = "Bcell", "Mast" = "Mast", "Macrophage" = "Myeloid",
  "Dendritic" = "Myeloid", "myocyte" = "Myocyte"
)
non_malignant_cell_type <- unname(type_map[non_malignant_cell_type])

classified_cancer <- meta_fields[[3]] == "1"
classified_noncancer <- meta_fields[[4]] == "1"

cell_type <- ifelse(classified_cancer, "Malignant", non_malignant_cell_type)

patient <- sub("^(HNSCC|HN)_?([0-9]+).*", "\\2", cell_names)
patient[!grepl("^(HNSCC|HN)_?[0-9]+", cell_names)] <- "combo"

cell_meta <- data.frame(
  row.names = cell_names,
  patient = patient,
  lymph_node = meta_fields[[2]] == "1",
  malignant = classified_cancer,
  non_malignant_cell_type = non_malignant_cell_type,
  cell_type = cell_type,
  stringsAsFactors = FALSE
)

message("Cell type breakdown:")
print(table(cell_meta$cell_type, useNA = "ifany"))

# --- Expression matrix (rows 7+) ---
# fread() is far faster than base read.delim() for a ~5,900 x 23,700 file.
expr_dt <- data.table::fread(f, skip = 6, header = FALSE, sep = "\t")
gene_names <- gsub("'", "", expr_dt[[1]])
stopifnot(!any(duplicated(gene_names)))

expr <- as.matrix(expr_dt[, -1])
storage.mode(expr) <- "numeric"
rownames(expr) <- gene_names
colnames(expr) <- cell_names

stopifnot(identical(colnames(expr), rownames(cell_meta)))

dir.create("data", showWarnings = FALSE)
saveRDS(list(expr = expr, meta = cell_meta), "data/gse103322_parsed.rds")

message("Saved data/gse103322_parsed.rds: ", nrow(expr), " genes x ", ncol(expr), " cells")

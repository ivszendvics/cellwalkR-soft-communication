# build_marker_genes.R
# Builds data/marker_genes.rds -- the curated marker_list that
# 02_run_cellwalker.R reads and converts (via markers_from_list()) into
# CellWalkR's expected markers table.
#
# --- Sourcing note ---
# The original plan was to transcribe cell-type marker genes directly from
# Puram et al. 2017's supplementary tables (Table S1/S5). In practice, both
# Cell.com (403) and PMC (reCAPTCHA) block automated/bot access to that PDF,
# so it couldn't be scraped programmatically -- someone with journal access
# should eventually cross-check this list against the real Table S1/S5.
#
# In the meantime, this list uses the well-established, broadly-canonical
# marker genes for each lineage (the same ones Puram et al. and essentially
# every HNSCC/PBMC scRNA-seq paper use for these broad cell types), cross-
# checked against PanglaoDB (https://panglaodb.se/markers.html) sensitivity/
# specificity tables for Fibroblasts, T cells, B cells, Macrophages,
# Endothelial cells, Mast cells, and Epithelial cells (malignant HNSCC cells
# are squamous-epithelial in origin). Genes were chosen for being
# lineage-specific (avoiding markers like PTPRC or CXCR4 that are shared
# across several immune types) and appearing consistently across sources.
#
# Cell types match the 7-type schema already used in scripts/02_run_cellwalker.R
# and the dataset description in README.md: Malignant, Fibroblast, Tcell,
# Bcell, Myeloid, Endothelial, Mast.

marker_list <- list(
  Malignant   = c("EPCAM", "KRT17", "KRT14", "KRT5", "SFN", "CDH1"),
  Fibroblast  = c("COL1A1", "COL1A2", "COL3A1", "DCN", "PDGFRA", "PDGFRB"),
  Tcell       = c("CD3D", "CD3E", "CD3G", "CD2", "TRAC", "IL7R"),
  Bcell       = c("CD79A", "CD79B", "MS4A1", "CD19", "IGHM", "BLNK"),
  Myeloid     = c("CD68", "LYZ", "AIF1", "CD14", "ITGAM", "CSF1R"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5", "ENG", "KDR"),
  Mast        = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC")
)

stopifnot(!any(duplicated(names(marker_list))))

dir.create("data", showWarnings = FALSE)
saveRDS(marker_list, "data/marker_genes.rds")

message("Saved marker_list for ", length(marker_list), " cell types to data/marker_genes.rds")
message("Genes per type: ", paste(names(marker_list), lengths(marker_list), sep = "=", collapse = ", "))

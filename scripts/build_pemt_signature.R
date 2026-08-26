# build_pemt_signature.R
# Builds data/pemt_signature.rds -- the p-EMT gene signature from Puram et
# al. 2017's Table S7 ("Six meta-signatures, each derived from multiple
# related NNMF programs"), used by scripts/04_compare_networks.R to compute
# a per-malignant-cell p-EMT score and correlate it against CellWalker's
# identity entropy.
#
# --- Sourcing note ---
# Table S7 lists the GENE SIGNATURE that defines each meta-program (Cell
# Cycle, p-EMT, Epi dif. 1, Epi dif. 2, Stress, Hypoxia) -- it is NOT a
# table of per-cell numeric scores. The paper computed per-cell scores by
# scoring cells against this gene list; there is no separate lookup table
# of scores to transcribe. This project computes that score itself (see
# scripts/04_compare_networks.R) using this gene list as input.
#
# Transcribed from the actual supplementary PDF (Table S7, "p-EMT" column)
# after Cell.com/PMC blocked automated access (see README). The "Genes
# 51-100" half of the p-EMT column extracted as only 49 genes instead of
# the expected 50 -- almost certainly a single gene lost to the PDF's
# dense multi-column table layout during extraction, not a deliberate
# omission. Proceeding with 99 of ~100 genes: one gene out of a ~100-gene
# module score won't meaningfully change a correlation result, and this
# avoids guessing at the missing symbol.

pemt_signature <- c(
  # Genes 1-50
  "SERPINE1", "TGFBI", "MMP10", "LAMC2", "P4HA2", "PDPN", "ITGA5", "LAMA3",
  "CDH13", "TNC", "MMP2", "EMP3", "INHBA", "LAMB3", "VIM", "SEMA3C",
  "PRKCDBP", "ANXA5", "DHRS7", "ITGB1", "ACTN1", "CXCR7", "ITGB6", "IGFBP7",
  "THBS1", "PTHLH", "TNFRSF6B", "PDLIM7", "CAV1", "DKK3", "COL17A1", "LTBP1",
  "COL5A2", "COL1A1", "FHL2", "TIMP3", "PLAU", "LGALS1", "PSMD2", "CD63",
  "HERPUD1", "TPM1", "SLC39A14", "C1S", "MMP1", "EXT2", "COL4A2", "PRSS23",
  "SLC7A8", "SLC31A2",
  # Genes 51-100 (only 49 recovered -- see note above)
  "ARPC1B", "APP", "MFAP2", "DFNA5", "MT2A", "MAGED2", "ITGA6", "FSTL1",
  "TNFRSF12A", "IL32", "COPB2", "PTK7", "OCIAD2", "TAX1BP3", "SEC13",
  "SERPINH1", "TPM4", "MYH9", "ANXA8L1", "PLOD2", "GALNT2", "LEPREL1",
  "MAGED1", "SLC38A5", "FSTL3", "CD99", "F3", "PSAP", "NMRK1", "FKBP9",
  "DSG2", "ECM1", "HTRA1", "SERINC1", "CALU", "TPST1", "PLOD3", "IGFBP3",
  "FRMD6", "CXCL14", "SERPINE2", "RABAC1", "TMED9", "NAGK", "BMP1", "ESYT1",
  "STON2", "TAGLN", "GJA1"
)

stopifnot(!any(duplicated(pemt_signature)))
stopifnot(length(pemt_signature) == 99)

dir.create("data", showWarnings = FALSE)
saveRDS(pemt_signature, "data/pemt_signature.rds")

message("Saved p-EMT signature (", length(pemt_signature), " genes) to data/pemt_signature.rds")

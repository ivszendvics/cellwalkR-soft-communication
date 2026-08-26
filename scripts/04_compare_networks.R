# 04_compare_networks.R
# Validation + figures:
#   1. Does CellWalker's entropy-based ambiguity score correlate with a
#      p-EMT score computed from Puram et al.'s published gene signature?
#      (independent sanity check that "ambiguous" cells aren't just noise)
#   2. Visualize soft vs. hard network differences.

library(ggplot2)

results <- readRDS("data/processed/communication_results.rds")
parsed  <- readRDS("data/gse103322_parsed.rds")

# --- p-EMT validation ---
# Puram et al.'s supplementary Table S7 does NOT publish a per-cell p-EMT
# score -- it publishes the p-EMT GENE SIGNATURE (the marker genes that
# define the program via NNMF), which is what's saved in
# data/pemt_signature.rds (see scripts/build_pemt_signature.R for sourcing
# detail, including a ~1-gene extraction gap: 99 of ~100 genes recovered
# from the actual PDF, since Cell.com/PMC block automated access to it).
# We score each malignant cell against that signature ourselves using
# Seurat::AddModuleScore() (average signature expression vs. a matched
# control gene set, binned by expression level -- the standard choice for
# this era of Tirosh-lab-style single-cell papers), then correlate that
# score against CellWalker's identity entropy.
dir.create("figures", showWarnings = FALSE)

pemt_signature <- readRDS("data/pemt_signature.rds")
cell_meta <- parsed$meta
malignant_cells <- rownames(cell_meta)[!is.na(cell_meta$cell_type) &
                                          cell_meta$cell_type == "Malignant"]
mal_expr <- parsed$expr[, malignant_cells]

mal_obj <- Seurat::CreateSeuratObject(counts = mal_expr)
mal_obj <- Seurat::SetAssayData(mal_obj, slot = "data", new.data = mal_expr)
mal_obj <- Seurat::AddModuleScore(mal_obj, features = list(pemt_signature), name = "pEMT")
pemt_score <- mal_obj$pEMT1
names(pemt_score) <- colnames(mal_obj)

entropy_mal <- results$entropy[malignant_cells]
stopifnot(identical(names(entropy_mal), names(pemt_score)))

pemt_cor <- cor.test(entropy_mal, pemt_score, method = "spearman")
message("Spearman correlation, identity entropy vs. p-EMT score (n = ",
        length(malignant_cells), " malignant cells): rho = ",
        round(pemt_cor$estimate, 3), ", p = ", format.pval(pemt_cor$p.value, digits = 3))

ggplot(data.frame(entropy = entropy_mal, pemt = pemt_score),
       aes(pemt, entropy)) +
  geom_point(alpha = 0.4) + geom_smooth(method = "lm") +
  labs(x = "p-EMT module score (this project's own scoring)",
       y = "CellWalker identity entropy",
       title = "Identity ambiguity vs. independently-derived p-EMT score",
       subtitle = paste0("Spearman rho = ", round(pemt_cor$estimate, 3),
                          ", p = ", format.pval(pemt_cor$p.value, digits = 3))) +
  theme_minimal()
ggsave("figures/pemt_entropy_correlation.png", width = 6, height = 4.5)
message("Saved figures/pemt_entropy_correlation.png")

# --- Soft vs hard network comparison plot ---
comparison <- results$comparison
comparison$edge <- paste(comparison$source, "->", comparison$target,
                          "(", comparison$ligand, "-", comparison$receptor, ")")

top_gained <- head(comparison[order(-comparison$diff), ], 15)

ggplot(top_gained, aes(x = reorder(edge, diff), y = diff)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = NULL, y = "Soft score - Hard score",
       title = "Interactions strengthened by probability-weighted scoring",
       subtitle = "Top edges gained when identity uncertainty is preserved") +
  theme_minimal()

ggsave("figures/soft_vs_hard_top_edges.png", width = 8, height = 5.5)

message("Saved figures/soft_vs_hard_top_edges.png")

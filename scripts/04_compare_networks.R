# 04_compare_networks.R
# Validation + figures:
#   1. Does CellWalker's entropy-based ambiguity score correlate with the
#      p-EMT score published in Puram et al.? (independent sanity check that
#      "ambiguous" cells aren't just noise)
#   2. Visualize soft vs. hard network differences.

library(ggplot2)

results <- readRDS("data/processed/communication_results.rds")
parsed  <- readRDS("data/gse103322_parsed.rds")

# --- p-EMT validation ---
# Puram et al. supply a per-malignant-cell p-EMT score (their Figure 4 /
# supplementary data). If you've pulled that column into cell_meta during
# 01_download_data.R parsing, correlate it against our entropy score for
# malignant cells only:
#
# pemt_score <- parsed$meta[names(results$entropy), "pEMT_score"]
# malignant_idx <- parsed$meta[names(results$entropy), "malignant"] == "yes"
# cor.test(results$entropy[malignant_idx], pemt_score[malignant_idx], method = "spearman")
#
# ggplot(data.frame(entropy = results$entropy[malignant_idx],
#                    pemt = pemt_score[malignant_idx]),
#        aes(pemt, entropy)) +
#   geom_point(alpha = 0.4) + geom_smooth(method = "lm") +
#   labs(x = "Published p-EMT score", y = "CellWalker identity entropy",
#        title = "Identity ambiguity tracks independently-derived p-EMT score")
# ggsave("figures/pemt_entropy_correlation.png", width = 6, height = 4.5)

# --- Soft vs hard network comparison plot ---
dir.create("figures", showWarnings = FALSE)

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
message("Uncomment the p-EMT block once pEMT_score is available in cell_meta.")

# 03_soft_communication.R
# Uses the CellWalker type-probability matrix from 02_run_cellwalker.R to
# compute both soft (probability-weighted) and hard (argmax-baseline)
# ligand-receptor communication networks.

devtools::load_all(".")  # loads R/soft_weighting.R and R/utils.R

parsed      <- readRDS("data/gse103322_parsed.rds")
type_probs  <- readRDS("data/processed/type_probs.rds")
expr_mat    <- parsed$expr  # genes x cells, normalized expression

stopifnot(identical(colnames(expr_mat), rownames(type_probs)))

# --- Identity ambiguity per cell ---
entropy   <- identity_entropy(type_probs)
ambiguous <- flag_ambiguous_cells(entropy, quantile_cutoff = 0.9)
message(sum(ambiguous), " cells flagged as high-ambiguity (top 10% entropy)")

# --- Ligand-receptor resource ---
# LIANA's resource name is "Consensus" (capitalized) -- select_resource() is
# case-sensitive and silently returns NULL on "consensus", not an error.
# The resource includes multi-subunit complexes as underscore-joined gene
# symbols (e.g. "ITGA4_ITGB7") in ~23% of rows; score_lr_pairs() expects a
# single gene symbol per side (matching a single row of the expression
# matrix), so complex entries are dropped here rather than expanded/summed.
# That leaves 3,607 simple pairs, 3,394 of which have both genes present in
# GSE103322's matrix.
consensus <- liana::select_resource("Consensus")[[1]]
is_complex <- grepl("_", consensus$source_genesymbol) | grepl("_", consensus$target_genesymbol)
lr_pairs <- consensus[!is_complex, c("source_genesymbol", "target_genesymbol")]
colnames(lr_pairs) <- c("ligand", "receptor")
lr_pairs <- unique(lr_pairs)

# --- Soft vs hard type-level expression ---
soft_expr <- weighted_type_expression(expr_mat, type_probs)
hard_expr <- hard_type_expression(expr_mat, type_probs)

# --- Score both networks ---
soft_scores <- score_lr_pairs(soft_expr, lr_pairs)
hard_scores <- score_lr_pairs(hard_expr, lr_pairs)

# --- Compare ---
comparison <- compare_soft_hard(soft_scores, hard_scores)

dir.create("data/processed", showWarnings = FALSE)
saveRDS(list(entropy = entropy, ambiguous = ambiguous,
             soft_scores = soft_scores, hard_scores = hard_scores,
             comparison = comparison),
        "data/processed/communication_results.rds")

message("Top interactions gained under soft weighting:")
print(head(comparison[, c("source", "target", "ligand", "receptor", "diff")], 10))

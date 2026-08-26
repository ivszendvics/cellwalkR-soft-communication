#' Soft cell-cell communication scoring
#'
#' Core idea: CellWalker/CellWalkR produces, for every cell, a probability
#' distribution over cell types (a diffusion score per cell x type), rather
#' than a single hard label. Standard ligand-receptor tools throw that away
#' by taking argmax before scoring. These functions keep the full
#' distribution and use it as a weighting scheme.
#'
#' Expected inputs throughout:
#'   expr_mat   : genes x cells expression matrix (normalized, e.g. Seurat's
#'                'data' slot), dgCMatrix or matrix
#'   type_probs : cells x cell-types probability matrix, rows summing to 1.
#'                This is what you get from CellWalkR's cell-to-label scores
#'                after row-normalizing.

#' Probability-weighted mean expression per cell type
#'
#' For each gene and each cell type, computes a probability-weighted average
#' expression across all cells, where a cell's weight for a given type is its
#' CellWalker probability of belonging to that type. This generalizes a
#' per-cluster mean (hard weighting: 1 if in cluster, 0 otherwise) to a soft
#' weighting.
#'
#' @param expr_mat genes x cells matrix
#' @param type_probs cells x types probability matrix (rownames = cell ids
#'   matching colnames(expr_mat), colnames = cell type names)
#' @return genes x types matrix of weighted mean expression
#' @export
weighted_type_expression <- function(expr_mat, type_probs) {
  stopifnot(all(colnames(expr_mat) == rownames(type_probs)))

  # Normalize weights per type so they act as a proper weighted mean
  # (columns of type_probs need not sum to 1 across cells, only rows do)
  col_sums <- colSums(type_probs)
  col_sums[col_sums == 0] <- 1  # avoid div by zero for unused types
  weights <- sweep(type_probs, 2, col_sums, "/")

  # genes x cells %*% cells x types -> genes x types
  as.matrix(expr_mat %*% weights)
}

#' Hard-label mean expression per cell type (baseline for comparison)
#'
#' Equivalent to what CellChat/LIANA compute internally: mean expression
#' within cells assigned to each type by argmax over type_probs.
#'
#' @inheritParams weighted_type_expression
#' @return genes x types matrix of mean expression per hard-assigned cluster
#' @export
hard_type_expression <- function(expr_mat, type_probs) {
  stopifnot(all(colnames(expr_mat) == rownames(type_probs)))
  hard_labels <- colnames(type_probs)[max.col(type_probs, ties.method = "first")]

  types <- colnames(type_probs)
  out <- matrix(0, nrow = nrow(expr_mat), ncol = length(types),
                 dimnames = list(rownames(expr_mat), types))
  for (ty in types) {
    idx <- which(hard_labels == ty)
    if (length(idx) == 0) next
    out[, ty] <- if (length(idx) == 1) expr_mat[, idx] else Matrix::rowMeans(expr_mat[, idx, drop = FALSE])
  }
  out
}

#' Score ligand-receptor pairs between cell types
#'
#' Given a genes x types expression matrix (from either
#' weighted_type_expression() or hard_type_expression()) and a ligand-receptor
#' pair table, compute a simple product score (ligand in source type x
#' receptor in target type) for every source-target-pair combination. This
#' mirrors the basic CellPhoneDB/CellChat scoring logic; swap in a fancier
#' score here later if desired.
#'
#' @param type_expr genes x types matrix
#' @param lr_pairs data.frame with columns `ligand`, `receptor` (gene symbols
#'   matching rownames(type_expr))
#' @return data.frame with columns source, target, ligand, receptor, score
#' @export
score_lr_pairs <- function(type_expr, lr_pairs) {
  types <- colnames(type_expr)
  lr_pairs <- lr_pairs[lr_pairs$ligand %in% rownames(type_expr) &
                          lr_pairs$receptor %in% rownames(type_expr), ]

  results <- vector("list", length(types)^2 * nrow(lr_pairs))
  i <- 1
  for (src in types) {
    for (tgt in types) {
      lig_expr <- type_expr[lr_pairs$ligand, src]
      rec_expr <- type_expr[lr_pairs$receptor, tgt]
      score <- lig_expr * rec_expr
      n <- nrow(lr_pairs)
      results[[i]] <- data.frame(
        source = rep(src, n), target = rep(tgt, n),
        ligand = lr_pairs$ligand, receptor = lr_pairs$receptor,
        score = score
      )
      i <- i + 1
    }
  }
  do.call(rbind, results)
}

#' Compare soft vs. hard communication networks
#'
#' Joins soft and hard LR score tables on (source, target, ligand, receptor)
#' and computes the difference, so you can rank interactions that only show
#' up (or show up much more strongly) once cell-identity uncertainty is
#' accounted for.
#'
#' @param soft_scores output of score_lr_pairs() run on weighted_type_expression()
#' @param hard_scores output of score_lr_pairs() run on hard_type_expression()
#' @return data.frame with soft_score, hard_score, diff, ranked by diff desc
#' @export
compare_soft_hard <- function(soft_scores, hard_scores) {
  key_cols <- c("source", "target", "ligand", "receptor")
  merged <- merge(soft_scores, hard_scores, by = key_cols,
                   suffixes = c("_soft", "_hard"))
  merged$diff <- merged$score_soft - merged$score_hard
  merged[order(-merged$diff), ]
}

#' Permutation-based significance for soft ligand-receptor scores
#'
#' Generalizes CellPhoneDB/CellChat-style permutation testing (shuffle which
#' cluster each cell belongs to, recompute, see how often chance beats the
#' observed score) to soft weights: shuffles which cell each row of
#' `type_probs` belongs to, breaking the link between a cell's expression and
#' its CellWalker identity distribution, and recomputes the soft LR score
#' under that null. The p-value for a given (source, target, ligand,
#' receptor) is the fraction of permutations whose null score is at least as
#' large as the observed one.
#'
#' @inheritParams weighted_type_expression
#' @param lr_pairs data.frame with columns `ligand`, `receptor`
#' @param n_perm number of permutations (default 1000)
#' @param p_adjust_method passed to `stats::p.adjust()` for the `q_value`
#'   column (default "BH"); use "none" to skip adjustment
#' @param seed optional seed for reproducibility
#' @return the observed `score_lr_pairs()` table with added `p_value` and
#'   `q_value` columns
#' @export
soft_lr_significance <- function(expr_mat, type_probs, lr_pairs, n_perm = 1000,
                                  p_adjust_method = "BH", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  observed_expr <- weighted_type_expression(expr_mat, type_probs)
  observed <- score_lr_pairs(observed_expr, lr_pairs)

  n_cells <- nrow(type_probs)
  ge_count <- rep(0L, nrow(observed))

  for (i in seq_len(n_perm)) {
    perm_probs <- type_probs[sample.int(n_cells), , drop = FALSE]
    rownames(perm_probs) <- rownames(type_probs)

    perm_expr <- weighted_type_expression(expr_mat, perm_probs)
    perm_scores <- score_lr_pairs(perm_expr, lr_pairs)

    # perm_scores has identical row order to observed: same types/lr_pairs
    # loop order in score_lr_pairs(), and expr_mat/lr_pairs are unchanged.
    ge_count <- ge_count + (perm_scores$score >= observed$score)
  }

  # Add-one smoothing avoids a p-value of exactly 0 from a finite sample.
  observed$p_value <- (ge_count + 1) / (n_perm + 1)
  observed$q_value <- stats::p.adjust(observed$p_value, method = p_adjust_method)
  observed
}

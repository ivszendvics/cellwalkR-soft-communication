#' Per-cell identity ambiguity (Shannon entropy over CellWalker type probs)
#'
#' Entropy 0 = cell is confidently one type. Higher entropy = the cell's
#' identity is split across multiple types (e.g. a p-EMT malignant cell
#' sitting between epithelial and mesenchymal signatures). Use this to flag
#' cells worth inspecting for identity-plasticity-driven signaling.
#'
#' @param type_probs cells x types probability matrix, rows summing to ~1
#' @return named numeric vector of entropy per cell, normalized to [0, 1]
#'   (1 = maximally ambiguous, i.e. uniform across all types)
#' @export
identity_entropy <- function(type_probs) {
  p <- type_probs
  p[p <= 0] <- NA  # avoid log(0); NA terms drop out of rowSums with na.rm
  raw_entropy <- -rowSums(p * log(p), na.rm = TRUE)
  max_entropy <- log(ncol(type_probs))
  result <- raw_entropy / max_entropy
  names(result) <- rownames(type_probs)  # rowSums should preserve these, but be explicit
  result
}

#' Flag high-ambiguity cells by entropy quantile
#'
#' @param entropy output of identity_entropy()
#' @param quantile_cutoff cells above this quantile are flagged (default: top 10%)
#' @return logical vector, same order/names as entropy
#' @export
flag_ambiguous_cells <- function(entropy, quantile_cutoff = 0.9) {
  entropy >= stats::quantile(entropy, quantile_cutoff, na.rm = TRUE)
}

#' Minimal curated ligand-receptor pair table (placeholder)
#'
#' A tiny starter set so the pipeline runs end-to-end before you swap in a
#' full curated resource (e.g. LIANA's consensus resource via
#' `liana::select_resource("consensus")`, or CellPhoneDB's list).
#'
#' @return data.frame with columns ligand, receptor
#' @export
example_lr_pairs <- function() {
  data.frame(
    ligand   = c("TGFB1", "VEGFA", "CXCL12", "JAG1", "TNF"),
    receptor = c("TGFBR2", "KDR", "CXCR4", "NOTCH1", "TNFRSF1A"),
    stringsAsFactors = FALSE
  )
}

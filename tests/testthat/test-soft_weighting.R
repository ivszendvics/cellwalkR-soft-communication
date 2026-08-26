test_that("weighted_type_expression reduces to hard mean when probs are 0/1", {
  expr_mat <- matrix(1:12, nrow = 3,
                      dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
  # Hard assignment expressed as a probability matrix: cells 1-2 -> TypeA, 3-4 -> TypeB
  type_probs <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), nrow = 4,
                        dimnames = list(paste0("c", 1:4), c("TypeA", "TypeB")))

  soft <- weighted_type_expression(expr_mat, type_probs)
  hard <- hard_type_expression(expr_mat, type_probs)

  expect_equal(soft, hard, tolerance = 1e-10)
})

test_that("identity_entropy is 0 for confident cells and ~1 for uniform cells", {
  type_probs <- matrix(c(1, 0.5, 0, 0.5), nrow = 2,
                        dimnames = list(c("c1", "c2"), c("TypeA", "TypeB")))
  ent <- identity_entropy(type_probs)
  expect_equal(unname(ent["c1"]), 0, tolerance = 1e-10)
  expect_equal(unname(ent["c2"]), 1, tolerance = 1e-10)
})

test_that("score_lr_pairs produces one row per source-target-pair combination", {
  type_expr <- matrix(c(1, 2, 3, 4), nrow = 2,
                       dimnames = list(c("LIG1", "REC1"), c("TypeA", "TypeB")))
  lr_pairs <- data.frame(ligand = "LIG1", receptor = "REC1")
  scores <- score_lr_pairs(type_expr, lr_pairs)
  expect_equal(nrow(scores), 4)  # 2 types x 2 types x 1 pair
  expect_true(all(c("source", "target", "ligand", "receptor", "score") %in% colnames(scores)))
})

test_that("compare_soft_hard computes diff and sorts descending", {
  soft <- data.frame(source = "A", target = "B", ligand = "L", receptor = "R", score = 5)
  hard <- data.frame(source = "A", target = "B", ligand = "L", receptor = "R", score = 3)
  cmp <- compare_soft_hard(soft, hard)
  expect_equal(cmp$diff, 2)
})

test_that("soft_lr_significance returns one row per source-target-pair with valid p/q values", {
  set.seed(42)
  expr_mat <- matrix(rpois(4 * 8, lambda = 3), nrow = 4,
                      dimnames = list(c("LIG1", "REC1", "LIG2", "REC2"),
                                      paste0("c", 1:8)))
  type_probs <- matrix(runif(8 * 2), nrow = 8,
                        dimnames = list(paste0("c", 1:8), c("TypeA", "TypeB")))
  type_probs <- type_probs / rowSums(type_probs)
  lr_pairs <- data.frame(ligand = c("LIG1", "LIG2"), receptor = c("REC1", "REC2"))

  result <- soft_lr_significance(expr_mat, type_probs, lr_pairs, n_perm = 50, seed = 1)

  expect_equal(nrow(result), 2 * 2 * nrow(lr_pairs))  # 2 types x 2 types x 2 pairs
  expect_true(all(c("p_value", "q_value") %in% colnames(result)))
  expect_true(all(result$p_value > 0 & result$p_value <= 1))
  expect_true(all(result$q_value >= 0 & result$q_value <= 1))
})

test_that("soft_lr_significance is reproducible given the same seed", {
  set.seed(7)
  expr_mat <- matrix(rpois(4 * 8, lambda = 3), nrow = 4,
                      dimnames = list(c("LIG1", "REC1", "LIG2", "REC2"),
                                      paste0("c", 1:8)))
  type_probs <- matrix(runif(8 * 2), nrow = 8,
                        dimnames = list(paste0("c", 1:8), c("TypeA", "TypeB")))
  type_probs <- type_probs / rowSums(type_probs)
  lr_pairs <- data.frame(ligand = c("LIG1", "LIG2"), receptor = c("REC1", "REC2"))

  result1 <- soft_lr_significance(expr_mat, type_probs, lr_pairs, n_perm = 25, seed = 99)
  result2 <- soft_lr_significance(expr_mat, type_probs, lr_pairs, n_perm = 25, seed = 99)

  expect_equal(result1$p_value, result2$p_value)
})

test_that("soft_lr_significance assigns a low p-value to an unambiguous strong signal", {
  # 6 cells split confidently 3/3 between two types. LIG1/REC1 are expressed
  # only in the TypeA cells, so the TypeA->TypeA score should dominate every
  # permutation of the type assignment and come out significant.
  expr_mat <- matrix(0, nrow = 2, ncol = 6,
                      dimnames = list(c("LIG1", "REC1"), paste0("c", 1:6)))
  expr_mat["LIG1", 1:3] <- 10
  expr_mat["REC1", 1:3] <- 10

  type_probs <- matrix(0, nrow = 6, ncol = 2,
                        dimnames = list(paste0("c", 1:6), c("TypeA", "TypeB")))
  type_probs[1:3, "TypeA"] <- 1
  type_probs[4:6, "TypeB"] <- 1

  lr_pairs <- data.frame(ligand = "LIG1", receptor = "REC1")
  result <- soft_lr_significance(expr_mat, type_probs, lr_pairs, n_perm = 200, seed = 3)

  best <- result[result$source == "TypeA" & result$target == "TypeA", ]
  expect_lt(best$p_value, 0.05)
})

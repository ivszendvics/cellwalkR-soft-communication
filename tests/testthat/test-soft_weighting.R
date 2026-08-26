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

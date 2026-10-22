test_that("getErrors returns MAE, RMSE and median segment variance", {
  data <- list(errors = c(0.1, -0.3, 0.2, 0, -0.1, 0.1),
               abs_seg = c(2, 2, 3, 3, 4, 4),
               abs_cn = c(1.8, 2.2, 2.5, 3.5, 3.9, 4.1),
               chromosomes = rep("1", 6))
  e <- getErrors(data)
  expect_named(e, c("clonality", "rmse", "MedianSegVar"))
  expect_equal(e[["clonality"]], mean(abs(data$errors)))
  expect_equal(e[["rmse"]], sqrt(mean(data$errors^2)))
  # per segment variances 0.08, 0.5 and 0.02
  expect_equal(e[["MedianSegVar"]], 0.08)
})

test_that("getErrors requires data", {
  expect_error(getErrors(NULL), "no data")
})

test_that("segment variance starts a new segment at the first changed bin", {
  # var(1, 1, 1) = 0 and var(5, 6, 7) = 1
  expect_equal(calculateSegmentVar(abs_seg = c(1, 1, 1, 5, 5, 5),
                                   abs_cn = c(1, 1, 1, 5, 6, 7)), 0.5)
})

test_that("segment variance splits segments at chromosome boundaries", {
  abs_seg <- c(2, 2, 2, 2)
  abs_cn <- c(1, 3, 2, 4)
  expect_equal(calculateSegmentVar(abs_seg, abs_cn), stats::var(abs_cn))
  expect_equal(calculateSegmentVar(abs_seg, abs_cn, chromosomes = c("1", "1", "2", "2")), 2)
})

test_that("segment variance ignores single-bin segments", {
  expect_equal(calculateSegmentVar(abs_seg = c(1, 2, 2), abs_cn = c(1, 2, 4)), 2)
})

test_that("segment variance accepts the matrices gridStats returns", {
  abs_seg <- matrix(c(1, 1, 1, 5, 5, 5), ncol = 1)
  abs_cn <- matrix(c(1, 1, 1, 5, 6, 7), ncol = 1)
  expect_equal(calculateSegmentVar(abs_seg, abs_cn), 0.5)
})

test_that("segment variance checks input lengths", {
  expect_error(calculateSegmentVar(abs_seg = c(1, 2), abs_cn = c(1, 2, 3)))
  expect_error(calculateSegmentVar(abs_seg = c(1, 2), abs_cn = c(1, 2), chromosomes = "1"))
})

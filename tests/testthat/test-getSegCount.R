test_that("getSegCount counts value changes and chromosome starts", {
  seg <- cbind(A = c(1, 1, 2, 2, 2, 2),
               B = c(1, 1, 1, 1, 1, 1))
  obj <- makeCN(seg, chromosome = c("1", "1", "1", "2", "2", "2"))
  # equal values either side of the chromosome boundary are separate segments
  expect_equal(getSegCount(obj), c(A = 3, B = 2))
})

test_that("getSegCount works on a single sample", {
  seg <- cbind(A = c(1, 1, 2, 2), B = c(1, 1, 1, 1))
  obj <- makeCN(seg, chromosome = rep("1", 4))
  expect_equal(getSegCount(obj[, "A"]), c(A = 2))
})

test_that("getSegCount ignores bins that are not in use", {
  seg <- matrix(c(1, 5, 1, 1), ncol = 1, dimnames = list(NULL, "A"))
  obj <- makeCN(seg, chromosome = rep("1", 4), use = c(TRUE, FALSE, TRUE, TRUE))
  expect_equal(getSegCount(obj), c(A = 1))
})

test_that("getSegCount matches the segments getSegTable returns", {
  seg <- cbind(A = c(1, 1, 2, 2, 3, 3, 3, 1),
               B = c(2, 2, 2, 2, 2, 2, 2, 2))
  obj <- makeCN(seg, chromosome = c(rep("1", 4), rep("2", 4)))
  segTable <- getSegTable(obj)
  expect_equal(as.vector(getSegCount(obj)),
               as.vector(table(segTable$sample)[colnames(seg)]))
})

test_that("getSegCount requires data", {
  expect_error(getSegCount(NULL), "no data")
})

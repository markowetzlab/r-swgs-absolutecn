test_that("getSegTable returns one row per segment with bin coordinates", {
  seg <- matrix(c(1, 1, 2, 2, 2, 2), ncol = 1, dimnames = list(NULL, "A"))
  obj <- makeCN(seg, chromosome = c("1", "1", "1", "2", "2", "2"))
  segTable <- getSegTable(obj)
  expect_equal(colnames(segTable), c("chromosome", "start", "end", "segVal", "sample"))
  expect_equal(segTable$chromosome, c("1", "1", "2"))
  expect_equal(segTable$start, c(1, 60001, 1))
  expect_equal(segTable$end, c(60000, 90000, 90000))
  expect_equal(segTable$segVal, c(1, 2, 2))
  expect_equal(segTable$sample, c("A", "A", "A"))
})

test_that("getSegTable rounds segment values and sets negatives to zero", {
  seg <- matrix(c(-0.2, -0.2, 2.12345, 2.12345), ncol = 1, dimnames = list(NULL, "A"))
  obj <- makeCN(seg, chromosome = rep("1", 4))
  expect_equal(getSegTable(obj)$segVal, c(0, 2.123))
})

test_that("getSegTable splits samples and skips unused bins", {
  seg <- cbind(A = c(1, 9, 1, 2), B = c(3, 9, 3, 3))
  obj <- makeCN(seg, chromosome = rep("1", 4), use = c(TRUE, FALSE, TRUE, TRUE))
  segTable <- getSegTable(obj)
  expect_equal(segTable$sample, c("A", "A", "B"))
  expect_equal(segTable$segVal, c(1, 2, 3))
  expect_equal(segTable$start, c(1, 90001, 1))
  expect_equal(segTable$end, c(90000, 120000, 120000))
})

test_that("getSegTable rejects other object types", {
  expect_error(getSegTable(data.frame(x = 1)), "segtable error")
})

test_that("collapseRDS returns a single object unchanged", {
  obj <- makeCN(cbind(A = c(1, 2)), chromosome = c("1", "1"))
  expect_identical(collapseRDS(list(obj)), obj)
})

test_that("collapseRDS combines single sample objects", {
  obj <- makeCN(cbind(A = c(1, 2, 2), B = c(3, 3, 1), C = c(2, 2, 2)), chromosome = rep("1", 3))
  combined <- collapseRDS(list(obj[, "A"], obj[, "B"], obj[, "C"]))
  expect_equal(Biobase::sampleNames(combined), c("A", "B", "C"))
  expect_equal(Biobase::assayDataElement(combined, "segmented"),
               Biobase::assayDataElement(obj, "segmented"))
})

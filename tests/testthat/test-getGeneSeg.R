# Five 30kb bins from 7,500,001 on chromosomes 1, 17 and X - the same
# coordinates on each chromosome so matches must also check the chromosome
# Bin indices are 1-5 (chr1), 6-10 (chr17) and 11-15 (chrX)
makeGeneObj <- function(use = NULL) {
  starts <- 7500001 + (0:4) * 30000
  seg <- matrix(2, nrow = 15, ncol = 1, dimnames = list(NULL, "A"))
  makeCN(seg, chromosome = rep(c("1", "17", "X"), each = 5),
         start = rep(starts, 3), use = use)
}

test_that("getGeneSeg defaults to the TP53 locus for hg19", {
  obj <- makeGeneObj()
  # TP53 hg19 17:7565097-7590863 spans the 3rd and 4th chr17 bins
  expect_equal(getGeneSeg(abs_data = obj, genome = "hg19"), 8:9)
})

test_that("getGeneSeg accepts targets with or without a chr prefix", {
  obj <- makeGeneObj()
  expect_equal(getGeneSeg(target = "17:7565097-7590863", abs_data = obj, genome = "hg19"), 8:9)
  expect_equal(getGeneSeg(target = "chr17:7565097-7590863", abs_data = obj, genome = "hg19"), 8:9)
})

test_that("getGeneSeg finds targets on chromosome X", {
  obj <- makeGeneObj()
  expect_equal(getGeneSeg(target = "X:7565097-7590863", abs_data = obj, genome = "hg19"), 13:14)
  expect_equal(getGeneSeg(target = "chrX:7565097-7590863", abs_data = obj, genome = "hg19"), 13:14)
})

test_that("getGeneSeg returns a single bin for an exact bin match", {
  obj <- makeGeneObj()
  expect_equal(getGeneSeg(target = "17:7560001-7589999", abs_data = obj, genome = "hg19"), 8)
})

test_that("getGeneSeg indexes bins after removing unused bins", {
  obj <- makeGeneObj(use = c(FALSE, rep(TRUE, 14)))
  expect_equal(getGeneSeg(abs_data = obj, genome = "hg19"), 7:8)
})

test_that("getGeneSeg rejects invalid targets", {
  obj <- makeGeneObj()
  expect_error(getGeneSeg(target = "17-7565097", abs_data = obj, genome = "hg19"),
               "not in format")
  expect_error(getGeneSeg(target = "18:1-100", abs_data = obj, genome = "hg19"),
               "no bins in use on chromosome 18")
  expect_error(getGeneSeg(target = "17:7590863-7565097", abs_data = obj, genome = "hg19"))
  expect_error(getGeneSeg(target = c("17:1-2", "17:3-4"), abs_data = obj, genome = "hg19"))
})

test_that("getGeneSeg validates data and genome", {
  obj <- makeGeneObj()
  expect_error(getGeneSeg(abs_data = NULL, genome = "hg19"), "no data")
  expect_error(getGeneSeg(abs_data = obj, genome = "hg20"))
})

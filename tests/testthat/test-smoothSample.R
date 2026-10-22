test_that("smoothSample validates its inputs", {
  obj <- makeCN(matrix(c(1, 2), ncol = 1), chromosome = c("1", "1"))
  expect_error(smoothSample(NULL), "no data")
  expect_error(smoothSample(obj, smooth = "TRUE"))
  expect_error(smoothSample(obj, smooth = TRUE, maxSegs = 22))
})

test_that("smoothSample returns the input when smoothing is off or not needed", {
  obj <- makeCN(matrix(c(1, 1, 2, 2), ncol = 1), chromosome = rep("1", 4))
  expect_identical(smoothSample(obj, smooth = FALSE), obj)
  expect_identical(smoothSample(obj, smooth = TRUE, maxSegs = 300), obj)
})

test_that("smoothSample resegments multi-sample objects over the segment limit", {
  set.seed(1)
  nbins <- 2000
  chromosome <- rep(c("1", "2"), each = nbins / 2)
  truth <- rep(rep(c(1, 1.5), length.out = 40), each = nbins / 40)
  cn <- cbind(A = truth + stats::rnorm(nbins, 0, 0.05),
              B = truth + stats::rnorm(nbins, 0, 0.05))
  # every bin starts a new segment before smoothing
  obj <- makeCN(cn, chromosome = chromosome, cn = cn)
  before <- getSegCount(obj)

  res <- suppressMessages(smoothSample(obj, smooth = TRUE, maxSegs = 30, seed = 1))
  after <- getSegCount(res)
  expect_s4_class(res, "QDNAseqCopyNumbers")
  expect_named(after, c("A", "B"))
  expect_true(all(after < before))
})

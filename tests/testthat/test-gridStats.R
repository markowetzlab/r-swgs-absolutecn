test_that("depthtocn converts relative depth to absolute copy number", {
  expect_equal(depthtocn(1, purity = 1, seqdepth = 0.5), 2)
  expect_equal(depthtocn(c(10, 20), purity = 1, seqdepth = 10), c(1, 2))
  # tumour CN 3 at purity 0.6 with single copy depth 10
  depth <- 10 * (0.6 * 3 + 2 * (1 - 0.6))
  expect_equal(depthtocn(depth, purity = 0.6, seqdepth = 10), 3)
})

test_that("gridStats recovers integer copy number at the true ploidy and purity", {
  purity <- 0.6
  seqdepth <- 50
  # mean CN of 2 so the true ploidy is 2
  absCN <- c(2, 2, 3, 3, 1, 1)
  rel <- matrix(seqdepth * (purity * absCN + 2 * (1 - purity)), ncol = 1)
  obj <- makeCN(rel, chromosome = rep(c("1", "2"), each = 3))

  g <- gridStats(obj, ploidy = 2, purity = purity)
  expect_named(g, c("cn", "seg", "seqdepth", "abs_seg", "abs_cn",
                    "integer_seg", "errors", "chromosomes"))
  expect_equal(g$seqdepth, seqdepth)
  expect_equal(as.vector(g$abs_seg), absCN)
  expect_equal(as.vector(g$integer_seg), absCN)
  expect_equal(as.vector(g$errors), rep(0, 6))
  expect_equal(g$chromosomes, rep(c("1", "2"), each = 3))
})

test_that("gridStats errors are the distance to the nearest integer state", {
  obj <- makeCN(matrix(c(1, 1, 2, 2, 2, 2), ncol = 1), chromosome = rep("1", 6))
  g <- gridStats(obj, ploidy = 2, purity = 1)
  # mean relative CN 10/6 at ploidy 2 scales values by 1.2
  expect_equal(as.vector(g$abs_seg), c(1.2, 1.2, 2.4, 2.4, 2.4, 2.4))
  expect_equal(as.vector(g$errors), c(0.2, 0.2, 0.4, 0.4, 0.4, 0.4))
})

test_that("gridStats validates its inputs", {
  obj <- makeCN(matrix(c(1, 2), ncol = 1), chromosome = c("1", "1"))
  expect_error(gridStats(NULL, ploidy = 2, purity = 1), "no data")
  expect_error(gridStats(obj, ploidy = "2", purity = 1))
  expect_error(gridStats(obj, ploidy = 2, purity = 0))
  expect_error(gridStats(obj, ploidy = 2, purity = 1.1))
  expect_error(gridStats(obj, ploidy = 1, purity = 1))
})

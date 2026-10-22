test_that("getDownsampleDepth at ploidy 2 and purity 1", {
  expect_equal(getDownsampleDepth(ploidy = 2, purity = 1, nbins_ref_genome = 1000),
               15 * 2 * 1000 * 1.098901)
})

test_that("getDownsampleDepth matches the original implementation", {
  ploidy <- 3.2
  purity <- 0.45
  nbins <- 90000
  original <- (((2 * (1 - purity) + purity * ploidy) / (ploidy * purity)) / purity) * 15 *
    (2 * (1 - purity) + purity * ploidy) * nbins * (1 / 0.91)
  expect_equal(getDownsampleDepth(ploidy = ploidy, purity = purity, nbins_ref_genome = nbins),
               original, tolerance = 1e-6)
})

test_that("getDownsampleDepth needs more reads at lower purity", {
  high <- getDownsampleDepth(ploidy = 2, purity = 0.9, nbins_ref_genome = 1000)
  low <- getDownsampleDepth(ploidy = 2, purity = 0.3, nbins_ref_genome = 1000)
  expect_gt(low, high)
})

test_that("getDownsampleDepth is vectorised over fits", {
  depths <- getDownsampleDepth(ploidy = c(2, 3), purity = c(1, 0.5), nbins_ref_genome = 1000)
  expect_equal(depths, c(getDownsampleDepth(ploidy = 2, purity = 1, nbins_ref_genome = 1000),
                         getDownsampleDepth(ploidy = 3, purity = 0.5, nbins_ref_genome = 1000)))
})

test_that("getDownsampleDepth validates its inputs", {
  expect_error(getDownsampleDepth(ploidy = 2, purity = 0, nbins_ref_genome = 1000))
  expect_error(getDownsampleDepth(ploidy = 2, purity = 1.5, nbins_ref_genome = 1000))
  expect_error(getDownsampleDepth(ploidy = "2", purity = 1, nbins_ref_genome = 1000))
  expect_error(getDownsampleDepth(ploidy = 2, purity = 1))
})

makePlotObj <- function() {
  seg <- matrix(rep(c(2, 3, 2, 1), each = 50), ncol = 1, dimnames = list(NULL, "A"))
  makeCN(seg, chromosome = rep(c("1", "2"), each = 100), cn = seg + 0.01)
}

test_that("plotProfile validates its inputs", {
  obj <- makePlotObj()
  expect_error(plotProfile(NULL), "no data")
  expect_error(plotProfile(obj, purity = 1, clonality = 0.1, rmse = 0.1))
  expect_error(plotProfile(obj, ploidy = 2, purity = 1.5, clonality = 0.1, rmse = 0.1))
  expect_error(plotProfile(obj, ploidy = 2, purity = 1, clonality = 0.1, rmse = 0.1, yrange = "10"),
               "yrange must be numeric")
})

test_that("plotProfile draws a profile", {
  withr::local_pdf(withr::local_tempfile(fileext = ".pdf"))
  expect_no_error(suppressMessages(
    plotProfile(makePlotObj(), ploidy = 2, purity = 1, clonality = 0.1, rmse = 0.1)
  ))
})

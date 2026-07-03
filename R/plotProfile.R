plotProfile <- function(relcn = NULL,
                        ploidy = NA,
                        purity = NA,
                        clonality = NA,
                        rmse = NA,
                        yrange = NULL) {
  if (is.null(relcn)) {
    stop("no data")
  }
  # Y axis range
  if (is.null(yrange)) {
    if (ploidy > 5) {
      yrange = 15
    } else {
      yrange = 10
    }
  } else{
    if (!is.numeric(yrange)) {
      stop("yrange must be a integer")
    }
  }
  # Plot abs fit
  sub <- paste0(
    "ploidy=",
    round(ploidy, 2),
    " | ",
    " purity=",
    round(purity, 2),
    " | ",
    " MAE=",
    round(clonality, 3),
    " | ",
    " RMSE=",
    round(rmse, 3)
  )

  QDNAseqmod::plot(
    relcn,
    doCalls = FALSE,
    showSD = TRUE,
    logTransform = FALSE,
    ylim = c(0, yrange),
    ylab = "Absolute tumour CN"
  )
  abline(h = 1:yrange - 1, col = "blue")
  mtext(sub, side = 1, line = 3.5)
}

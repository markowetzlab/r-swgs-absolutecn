getDownsampleDepth <- function(ploidy = NULL,
                               purity = NULL,
                               nbins_ref_genome = NULL,
                               rpc = 15,
                               ratio = 1.098901) {
  # original implementation
  # (((2*(1-purity)+purity*ploidy)/(ploidy*purity))/purity)*15
  # *(2*(1-purity)+purity*ploidy)*nbins_ref_genome*(1/0.91)
  if (any(is.null(ploidy),
          is.null(purity),
          is.null(nbins_ref_genome))) {
    stop("missing parameters")
  }

  cellploidy <- purity * ploidy + (2 * (1 - purity))
  relratio <- (cellploidy / (ploidy * purity)) / purity

  readRatio <- relratio * rpc * cellploidy * nbins_ref_genome * ratio
  return(readRatio)
}

#' getDownsampleDepth
#'
#' @param ploidy numeric ploidy value
#' @param purity numeric purity value
#' @param nbins_ref_genome number of bins in reference genome
#' @param rpc minimum read count per bin per copy
#' @param ratio availble/used read ratio
#'
#' @returns numeric
#' @export
getDownsampleDepth <- function(ploidy = NULL,
                               purity = NULL,
                               nbins_ref_genome = NULL,
                               rpc = 15,
                               ratio = 1.098901) {
  # original implementation
  # (((2*(1-purity)+purity*ploidy)/(ploidy*purity))/purity)*15
  # *(2*(1-purity)+purity*ploidy)*nbins_ref_genome*(1/0.91)
  stopifnot(is.numeric(ploidy),is.numeric(purity),
            is.numeric(nbins_ref_genome),
            is.numeric(rpc),is.numeric(ratio))

  stopifnot(purity > 0,purity <= 1)

  cellploidy <- purity * ploidy + (2 * (1 - purity))
  relratio <- (cellploidy / (ploidy * purity)) / purity

  readRatio <- relratio * rpc * cellploidy * nbins_ref_genome * ratio
  return(readRatio)
}

#' getSegCount
#'
#' Counts copy number segments per sample. A new segment starts where the
#' segmented value changes or a new chromosome begins, so equal values either
#' side of a chromosome boundary count as two segments (matching getSegTable
#' and the segment variance calculation)
#'
#' @param obj QDNASeq object of class "QDNAseqCopyNumbers"
#'
#' @returns named numeric vector of segment counts, one per sample
#' @export
getSegCount <- function(obj = NULL) {
  if (is.null(obj)) {
    stop("no data")
  }

  use <- Biobase::fData(obj)$use
  seg <- Biobase::assayDataElement(obj, "segmented")[use, , drop = FALSE]
  chromosomes <- Biobase::fData(obj)$chromosome[use]

  n <- nrow(seg)
  newChr <- chromosomes[-1] != chromosomes[-n]
  segCount <- apply(seg, 2, function(x) {
    newSeg <- c(TRUE, x[-1] != x[-n] | newChr)
    newSeg[is.na(newSeg)] <- TRUE
    sum(newSeg)
  })
  return(segCount)
}

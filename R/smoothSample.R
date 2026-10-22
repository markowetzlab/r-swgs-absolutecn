#' smoothSample
#'
#' Adjusted recursively to drop max segments below threshold using maximum StdDev
#' difference in means during segmentation splits in CBS
#'
#' @param relcn QDNASeq object
#' @param smooth Boolean to apply smoothing
#' @param maxSegs Maximum segments
#' @param seed Value for seed
#'
#' @returns QDNASeq object
#' @export
smoothSample <- function(relcn = NULL,
                          smooth = FALSE,
                          maxSegs = 300,
                          seed = NULL) {
  if (is.null(relcn)) {
    stop("no data")
  }

  stopifnot(is.logical(smooth))
  stopifnot(is.numeric(maxSegs),maxSegs > 22)

  # Check if smoothing needed
  if(smooth){
    currsamp <- relcn

    maxseg <- maxSegs
    sdadjust <- 1.5

    segnum <- getSegCount(currsamp)
    while(any(segnum > maxseg) && sdadjust < 5){
      currsamp <- QDNAseqmod::segmentBins(currsamp, transformFun = "sqrt",
                                          undo.SD = sdadjust, seeds = seed)
      segnum <- getSegCount(currsamp)
      sdadjust <- sdadjust + 0.5
    }
    relcn <- currsamp
  }
  return(relcn)
}

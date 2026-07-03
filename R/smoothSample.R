# segment smoothing function
# Adjusted recursively to drop max segments below threshold using maximum StdDev
# difference in means during segmentation splits in CBS
smoothSample <- function(relcn = NULL,
                          smooth = FALSE,
                          maxSegs = 300,
                          seed = NULL) {
  if (is.null(relcn)) {
    stop("segment smoothing provided with no data")
  }

  stopifnot(is.logical(smooth))
  stopifnot(is.numeric(maxSegs),maxSegs > 22)

  # Check if smoothing needed
  relative_tmp <- NULL
  if (smooth) {
    currsamp <- relcn

    maxseg <- maxSegs
    sdadjust <- 1.5

    condition <- Biobase::fData(currsamp)$use

    segments <- Biobase::assayDataElement(currsamp, "segmented")[condition, , drop =
                                                                   FALSE]
    segments <- apply(segments, 2, rle)
    segnum <- as.numeric(lapply(segments, function(x) {
      length(x$lengths)
    }))

    while (segnum > maxseg & sdadjust < 5) {
      currsamp <- QDNAseqmod::segmentBins(
        currsamp,
        transformFun = "sqrt",
        undo.SD = sdadjust,
        seeds = seed
      )

      segments <- Biobase::assayDataElement(currsamp, "segmented")[condition, , drop =
                                                                     FALSE]
      segments <- apply(segments, 2, rle)
      segnum <- as.numeric(lapply(segments, function(x) {
        length(x$lengths)
      }))

      sdadjust <- sdadjust + 0.5
    }
    #relative_tmp <- currsamp
    relcn <- currsamp
  }
  return(relcn)
}

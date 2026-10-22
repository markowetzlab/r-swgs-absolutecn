#' getSegTable
#'
#' @param x Object of class "QDNASeqCopyNumbers"
#'
#' @returns data.frame
#' @export
getSegTable <- function(x) {
  if (inherits(x, what = "QDNAseqCopyNumbers", which = F)) {
    sn <- Biobase::assayDataElement(x, "segmented")
    fd <- Biobase::fData(x)
    use <- fd$use
    fdfiltfull <- fd[use, ]
    sn <- sn[use, ]
    if (is.null(ncol(sn))) {
      sampleName <- Biobase::sampleNames(x)
      sn <- as.data.frame(sn)
      colnames(sn) <- sampleName
    }
    segTable <- c()
    for (smp in colnames(sn)) {
      for (chr in unique(fdfiltfull$chromosome)) {
        snfilt <- sn[fdfiltfull$chromosome == chr, colnames(sn) == smp]
        fdfilt <- fdfiltfull[fdfiltfull$chromosome == chr, ]
        sn.rle <- rle(snfilt)
        starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
        ends <- cumsum(sn.rle$lengths)

        segtmp <- do.call(rbind, lapply(1:length(sn.rle$lengths), function(i) {
          from <- fdfilt$start[starts[i]]
          to <- fdfilt$end[ends[i]]
          segValue <- sn.rle$values[i]
          c(fdfilt$chromosome[starts[i]], from, to, segValue)
        }))

        segTableRaw <- cbind(segtmp, sample = rep(smp, times = nrow(segtmp)))
        segTable <- rbind(segTable, segTableRaw)
      }
    }
    segTable <- as.data.frame(segTable)
    colnames(segTable) <- c("chromosome", "start", "end", "segVal", "sample")
    segTable$segVal <- round(as.numeric(segTable$segVal), 3)
    segTable$segVal[segTable$segVal < 0] <- 0
    segTable$start <- as.numeric(segTable$start)
    segTable$end <- as.numeric(segTable$end)
    return(segTable)
  } else {
    # NON QDNASEQ BINNED DATA FUNCTION
    stop("segtable error")
  }
}

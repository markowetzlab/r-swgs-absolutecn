#' getErrors
#'
#' @param data list of numeric vectors
#'
#' @returns named numeric vector
#' @export
getErrors <- function(data=NULL){
  if(is.null(data)){
    stop("no data")
  }
  # clonality is a legacy name for MAE
  clonality <- mean(abs(data$errors))
  # Root Mean Squared Error
  rmse <- sqrt(mean(data$errors^2))
  MedianSegVar <- calculateSegmentVar(abs_seg = data$abs_seg,
                                      abs_cn = data$abs_cn,
                                      chromosomes = data$chromosomes)
  return(c(clonality=clonality,
           rmse=rmse,
           MedianSegVar=MedianSegVar))
}

calculateSegmentVar <- function(abs_seg = NULL, abs_cn = NULL, chromosomes = NULL) {
  abs_seg <- as.numeric(abs_seg)
  abs_cn <- as.numeric(abs_cn)

  if (is.null(chromosomes)) {
    chromosomes <- rep(1, length(abs_seg))
  }
  stopifnot(length(chromosomes) == length(abs_seg),
            length(abs_cn) == length(abs_seg))

  # a new segment starts where segment value or chromosome changes
  n <- length(abs_seg)
  newSeg <- c(TRUE, abs_seg[-1] != abs_seg[-n] | chromosomes[-1] != chromosomes[-n])
  newSeg[is.na(newSeg)] <- TRUE
  segId <- cumsum(newSeg)

  segVar <- tapply(abs_cn, segId, stats::var)

  medianVar <- stats::median(segVar, na.rm = TRUE)
  return(medianVar)
}

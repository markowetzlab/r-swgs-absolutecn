#' gridStats
#'
#' @param obj QDNASeq object of class "QDNAseqCopyNumbers"
#' @param ploidy numeric ploidy value
#' @param purity numeric purity value
#'
#' @returns list
#' @export
gridStats <- function(obj=NULL,ploidy=NULL,purity=NULL){

  if(is.null(obj)){
    stop("no data")
  }

  stopifnot(is.numeric(ploidy),is.numeric(purity))
  stopifnot(purity > 0,purity <= 1)
  stopifnot(ploidy > 1)

  seg <- Biobase::assayDataElement(obj,"segmented")
  cn <- Biobase::assayDataElement(obj,"copynumber")
  chromosomes <- Biobase::fData(obj)$chromosome
  rel_ploidy <- mean(cn,na.rm=TRUE)
  cellploidy <- purity * ploidy + (2 * (1 - purity))
  seqdepth <- rel_ploidy / cellploidy

  abs_seg <- depthtocn(seg,purity,seqdepth)
  abs_cn <- depthtocn(cn,purity,seqdepth)
  integer_seg <- round(abs_seg,digits = 0)
  errors <- abs_seg - integer_seg
  return(list(cn=cn,
              seg=seg,
              seqdepth=seqdepth,
              abs_seg=abs_seg,
              abs_cn=abs_cn,
              integer_seg=integer_seg,
              errors=errors,
              chromosomes=chromosomes))
}

#' depthtocn
#'
#' converts read depth to copy number given purity and single copy depth
#'
#' @param x numeric vector
#' @param purity float
#' @param seqdepth integer
#'
#' @returns numeric vector
#' @export
depthtocn <- function(x, purity, seqdepth){
  (x / seqdepth - 2 * (1 - purity)) / purity
}

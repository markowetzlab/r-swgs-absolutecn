#' getGeneSeg
#'
#' @param target character string of genomic region in format "N:start-end"
#' @param abs_data QDNASeq object of class "QDNAseqCopyNumbers"
#' @param genome reference genome
#'
#' @returns numeric vector
#' @export
getGeneSeg <- function(target = NULL,
                         abs_data = NULL,
                         genome = NULL) {
  if (is.null(abs_data)) {
    stop("no data")
  }
  rlang::arg_match(genome,values = c("hg19","hg38"))

  if (is.null(target)) {
    if (genome == "hg19") {
      target <- c("17:7565097-7590863")
    } else if (genome == "hg38") {
      target <- c("17:7661779-7687538")
    }
  }

  to_use <- Biobase::fData(abs_data)$use
  cn_obj <- abs_data[to_use, ]
  bin_pos <- Biobase::fData(cn_obj)[, c("chromosome", "start", "end")]
  position <- as.numeric(stringr::str_split(
    string = target,
    pattern = ":|-",
    simplify = T
  ))

  mapply(assign,
         c("chr", "start", "end"),
         position,
         MoreArgs = list(envir = parent.frame()))

  gene_chr_pos <- bin_pos[bin_pos$chromosome == chr, ]

  minusStart <- gene_chr_pos$start - start
  min_start <- min(which(min(abs(minusStart)) == abs(minusStart)))

  minusEnd <- gene_chr_pos$end - end
  min_end <- max(which(min(abs(minusEnd)) == abs(minusEnd)))

  if (gene_chr_pos$start[min_start] > start & min_start != 1){
    min_start <- min_start - 1
  }
  if (gene_chr_pos$end[min_end] < end & min_end != length(gene_chr_pos$end)){
    min_end <- min_end + 1
  }

  index_min <- which(bin_pos$chromosome == chr & bin_pos$start == gene_chr_pos[min_start, 2])
  index_max <- which(bin_pos$chromosome == chr & bin_pos$end == gene_chr_pos[min_end, 3])
  gene_pos <- seq.int(index_min, index_max, 1)

  return(gene_pos)
}

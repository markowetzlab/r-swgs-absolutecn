#' getGeneSeg
#'
#' @param target character string of genomic region in format "chr:start-end"
#'  including chr prefix, autosomes, and X/Y
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

  if(is.null(target)) {
    if(genome == "hg19") {
      target <- c("17:7565097-7590863")
    } else if(genome == "hg38") {
      target <- c("17:7661779-7687538")
    }
  }

  to_use <- Biobase::fData(abs_data)$use
  cn_obj <- abs_data[to_use, ]
  bin_pos <- Biobase::fData(cn_obj)[, c("chromosome", "start", "end")]

  stopifnot(length(target) == 1)
  parts <- stringr::str_match(target, "^(?:chr)?([^:]+):([0-9]+)-([0-9]+)$")
  if(is.na(parts[1, 1])) {
    stop(paste0("target not in format 'chr:start-end' - ", target))
  }
  targetChr <- parts[1, 2]
  position <- as.numeric(parts[1, 3:4])
  names(position) <- c("start", "end")
  stopifnot(position["start"] <= position["end"])

  onChr <- as.character(bin_pos$chromosome) == targetChr
  gene_chr_pos <- bin_pos[onChr, ]
  if(nrow(gene_chr_pos) == 0) {
    stop(paste0("no bins in use on chromosome ", targetChr, " - ", target))
  }

  minusStart <- gene_chr_pos$start - position["start"]
  min_start <- min(which(min(abs(minusStart)) == abs(minusStart)))

  minusEnd <- gene_chr_pos$end - position["end"]
  min_end <- max(which(min(abs(minusEnd)) == abs(minusEnd)))

  if(gene_chr_pos$start[min_start] > position["start"] & min_start != 1){
    min_start <- min_start - 1
  }

  if(gene_chr_pos$end[min_end] < position["end"] & min_end != length(gene_chr_pos$end)){
    min_end <- min_end + 1
  }

  index_min <- which(onChr & bin_pos$start == gene_chr_pos[min_start, 2])
  index_max <- which(onChr & bin_pos$end == gene_chr_pos[min_end, 3])
  gene_pos <- seq.int(index_min, index_max, 1)

  return(gene_pos)
}

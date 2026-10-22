# Build a small QDNAseqCopyNumbers object for tests
# seg/cn are bins x samples matrices, bins are 30kb and numbered from 1 on each
# chromosome unless start positions are given
makeCN <- function(seg, chromosome, cn = seg, use = NULL, start = NULL) {
  seg <- as.matrix(seg)
  cn <- as.matrix(cn)
  if (is.null(colnames(seg))) {
    colnames(seg) <- paste0("S", seq_len(ncol(seg)))
  }
  colnames(cn) <- colnames(seg)
  if (is.null(use)) {
    use <- rep(TRUE, nrow(seg))
  }
  if (is.null(start)) {
    binIndex <- stats::ave(seq_along(chromosome), chromosome, FUN = seq_along)
    start <- (binIndex - 1) * 30000 + 1
  }
  end <- start + 30000 - 1
  ids <- paste0(chromosome, ":", start, "-", end)
  rownames(seg) <- ids
  rownames(cn) <- ids

  bins <- data.frame(chromosome = chromosome, start = start, end = end, use = use,
                     row.names = ids, stringsAsFactors = FALSE)
  # pData fields QDNAseqmod::plot reads
  pheno <- data.frame(name = colnames(seg), expected.variance = 0.01, used.reads = 1e6,
                      loess.span = 0.65, loess.family = "symmetric",
                      row.names = colnames(seg))

  classDef <- methods::getClass("QDNAseqCopyNumbers", where = asNamespace("QDNAseqmod"))
  obj <- methods::new(classDef, bins = bins, copynumber = cn, phenodata = pheno)
  Biobase::assayDataElement(obj, "segmented") <- seg
  obj
}

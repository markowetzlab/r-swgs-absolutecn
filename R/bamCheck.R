#' bamCheck
#'
#' @param x file path to BAM/CRAM file
#' @param filetype character string matching either "BAM" or "CRAM"
#' @param outname file path for output
#'
#' @returns none
#' @export
bamCheck <- function(x = NULL,filetype = NULL,outname = NULL) {
  if (is.null(x)) {
    stop("no BAMs provided")
  } else if(!file.exists(x)){
    stop(paste0("BAM file does not exist - ",x))
  }

  if (is.null(outname)) {
    stop("no outname provided")
  }

  # Remove existing bam.ok/bam.invalid files to prevent false pass
  if (grepl("_bam\\.ok$", outname)) {
    invalidOut <- sub("_bam\\.ok$", "_bam.invalid", outname)
  } else {
    invalidOut <- paste0(outname, ".invalid")
  }
  unlink(c(outname, invalidOut))

  rlang::arg_match(filetype,
                   values = c("BAM","CRAM"),
                   multiple = FALSE)

  ## Call samtools using cmdline - path shell-quoted as it could contain spaces
  xQ <- shQuote(x)
  cmd <- paste0("samtools quickcheck -q ", xQ)

  typecmd <- paste0("head -c 4 ", xQ)
  typecheck <- as.character(system(typecmd,intern = TRUE))

  if (system(cmd) == 0){
    if(typecheck == "CRAM" & filetype == "BAM"){
      log_vector <- paste0("Filetype specified as BAM, detected CRAM - ", x)
      check_vector <- TRUE
      # "\037\x8b\b\004" BAM file gzip magic number
    } else if(typecheck == "\037\x8b\b\004" & filetype == "CRAM"){
      log_vector <- paste0("Filetype specified as CRAM, detected BAM - ", x)
      check_vector <- TRUE
    } else {
      log_vector <- paste0("BAM/CRAM valid - ", x)
      check_vector <- FALSE
    }
  } else {
    log_vector <- paste0("BAM/CRAM invalid - ", x)
    check_vector <- TRUE
  }

  if(check_vector) {
    writeLines(text = as.character(log_vector), con = invalidOut)
  } else {
    writeLines(text = as.character(log_vector), con = outname)
  }
}

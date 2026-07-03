bamCheck <- function(x = NULL,filetype=NULL,outname = NULL) {
  if (is.null(x)) {
    stop("no BAMs provided")
  }
  if (is.null(filetype)) {
    stop("no filetype provided")
  }
  if (is.null(outname)) {
    stop("no outname provided")
  }

  ## Call samtools using cmdline
  cmd <- paste0("samtools quickcheck -q ", x)
  typecmd <- paste0("head -c 4 ", x)

  typecheck <- system(typecmd,intern = TRUE)

  if (system(cmd) == 0){
    if(typecheck == "CRAM" & filetype == "BAM"){
      log_vector <- paste0("Filetype specified as BAM, detected CRAM - ", x)
      check_vector <- TRUE
    } else if(typecheck == "" & filetype == "CRAM"){
      log_vector <- paste0("Filetype specified as CRAM, detected BAM - ", x)
      check_vector <- TRUE
    } else {
      log_vector <- paste0("BAM/CRAM valid - ", x)
      check_vector <- FALSE
    }

  } else {
    log_vector <- paste0("BAM/CRAM invalid or missing - ", x)
    check_vector <- TRUE
  }

  if(check_vector) {
    outname <- gsub(pattern = "ok", replacement = "invalid", outname)
    writeLines(text = as.character(log_vector), con = outname)
  } else {
    writeLines(text = as.character(log_vector), con = outname)
  }
}

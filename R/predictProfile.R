#' predictProfile
#'
#' @param qctable data.frame of absolute fits
#' @param model ML model object
#' @param method Method for autofitting
#' @param flagThreshold Threshold for flagging fits
#' @param errorMetric Error metric for fitting
#'
#' @returns Data.frame of fitted absolute fits
#' @export
predictProfile <- function(qctable = NULL,
                           model = NULL,
                           method = "randforest",
                           flagThreshold = 0.84,
                           errorMetric = "clonality"){

  if (is.null(model) & method == "randforest") {
    stop("no model")
  }

  if (!method %in% c("randforest", "errorOnly")) {
    stop(paste0("unknown method - use 'randforest' or 'errorOnly'"))
  }

  if (!errorMetric %in% c("clonality", "segvariance", "rmse")) {
    stop(paste0(
      "unknown error metric - use 'clonality', 'segvariance' or 'rmse'"
    ))
  }

  if (is.null(qctable) || nrow(qctable) == 0) {
    stop("no fits in qctable - check gridsearch filtering for this sample")
  }

  qctable$sample <- qctable$SAMPLE_ID

  switch(method, "randforest" = {
    newClass <- parsnip::predict.model_fit(model, qctable, type = "class")
    newProb <- round(parsnip::predict.model_fit(model, qctable, type = "prob"),
                     digits = 3)
    qctable <- cbind(qctable, newClass, newProb)

    qctable <- triageProfile(
      qctable = qctable,
      flagThreshold = flagThreshold,
      errorMetric = errorMetric
    )

    qctable <- qctable %>%
      dplyr::select(-sumt)
  }, "errorOnly" = {
    switch(
      errorMetric,
      "clonality" = {
        qctable <- qctable %>%
          dplyr::group_by(sample) %>%
          dplyr::arrange(sample, clonality) %>%
          dplyr::mutate(minVal = if(all(is.na(clonality))) NA_real_ else min(clonality, na.rm = TRUE),
                        use = !is.na(clonality) & clonality == minVal) %>%
          dplyr::select(-minVal)
      },
      "segvariance" = {
        qctable <- qctable %>%
          dplyr::group_by(sample) %>%
          dplyr::arrange(sample, segvariance) %>%
          dplyr::mutate(minVal = if(all(is.na(segvariance))) NA_real_ else min(segvariance, na.rm = TRUE),
                        use = !is.na(segvariance) & segvariance == minVal) %>%
          dplyr::select(-minVal)
      },
      "rmse" = {
        qctable <- qctable %>%
          dplyr::group_by(sample) %>%
          dplyr::arrange(sample, rmse) %>%
          dplyr::mutate(minVal = if(all(is.na(rmse))) NA_real_ else min(rmse, na.rm = TRUE),
                        use = !is.na(rmse) & rmse == minVal) %>%
          dplyr::select(-minVal)
      }
    )
    newCols <- c("pred_class",
                 "pred_FALSE",
                 "pred_TRUE",
                 "triageValue",
                 "flag")
    qctable <- qctable %>%
      dplyr::mutate(!!!stats::setNames(rep(NA, length(newCols)), newCols))
  })

  qctable <- qctable %>%
    dplyr::ungroup() %>%
    dplyr::select(-c("sample")) %>%
    dplyr::mutate(
      notes = paste0(
        "autofit|",
        "fitmethod=",
        method,
        "|flagThreshold=",
        flagThreshold,
        "|ErrorMetric=",
        errorMetric)) %>%
    dplyr::mutate(use = as.logical(use))

  # tie break for edge cases
  ## For cases with multiple where ploidy is n + 2CN and twice the purity
  ## have the same error metric values. Select the lowest ploidy of accepted fits.
  if(sum(qctable$use) > 1){
    # lowest ploidy among accepted fits only - rejected fits must not set the minimum
    minAcceptedPloidy <- qctable %>%
      dplyr::filter(use == TRUE) %>%
      dplyr::pull(ploidy) %>%
      min()
    qctable <- qctable %>%
      dplyr::mutate(use = ifelse(use == TRUE & ploidy == minAcceptedPloidy,TRUE,FALSE))
    if(sum(qctable$use) > 1){
      stop("multiple fits selected - unknown tiebreak failure")
    }
  }

  if(all(!as.logical(qctable$use))){
    # Default fit of ploidy 2 / purity 1. getDownsampleDepth is used to compute bin size
    nbins <- qctable$downsample_depth[1] /
      getDownsampleDepth(ploidy = qctable$ploidy[1], purity = qctable$purity[1],
                         nbins_ref_genome = 1)
    defaultDepth <- round(getDownsampleDepth(ploidy = 2, purity = 1,
                                             nbins_ref_genome = nbins))

    sampleCols <- intersect(c("TP53freq", "precPloidy", "precPurity"),
                            colnames(qctable))
    sampleVals <- lapply(qctable[sampleCols], function(x) unique(x))

    qctable <- qctable %>%
      dplyr::ungroup() %>%
      tibble::add_row(
        SAMPLE_ID = unique(.$SAMPLE_ID),
        PATIENT_ID = unique(.$PATIENT_ID),
        ploidy = 2.0,
        purity = 1,
        smooth = unique(.$smooth),
        segments = unique(.$segments),
        downsample_depth = defaultDepth,
        notes = unique(.$notes),
        !!!sampleVals,
        use = TRUE
      ) %>%
      dplyr::mutate(flag = ifelse(is.na(flag), NA, paste0("NOFIT|", flag))) %>%
      dplyr::mutate(notes = paste0("NO_FIT_FOUND|", notes))
  }
  return(qctable)
}

triageProfile <- function(qctable = NULL,
                          flagThreshold = 0.84,
                          errorMetric = "segvariance") {
  if(is.null(qctable)) {
    stop("no data")
  }

  qctable <- qctable %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(triageValue = abs(.pred_TRUE - .pred_FALSE)) %>%
    dplyr::mutate(sumt = sum(as.logical(.pred_class)))

  # Where multiple fits are predicted TRUE, select the lowest error metric
  switch(
    errorMetric,
    "clonality" = {
      qctable <- qctable %>%
        dplyr::group_by(sample) %>%
        dplyr::arrange(sample, clonality) %>%
        dplyr::mutate(predTrue = as.logical(.pred_class),
                      minTrue = if(any(predTrue & !is.na(clonality))) min(clonality[predTrue], na.rm = TRUE) else NA_real_,
                      use = ifelse(
                        sumt <= 1,
                        predTrue,
                        predTrue & !is.na(clonality) & clonality == minTrue
                      )) %>%
        dplyr::select(-predTrue,-minTrue)
    },
    "segvariance" = {
      qctable <- qctable %>%
        dplyr::group_by(sample) %>%
        dplyr::arrange(sample, segvariance) %>%
        dplyr::mutate(predTrue = as.logical(.pred_class),
                      minTrue = if(any(predTrue & !is.na(segvariance))) min(segvariance[predTrue], na.rm = TRUE) else NA_real_,
                      use = ifelse(
                        sumt <= 1,
                        predTrue,
                        predTrue & !is.na(segvariance) & segvariance == minTrue
                      )) %>%
        dplyr::select(-predTrue,-minTrue)
    },
    "rmse" = {
      qctable <- qctable %>%
        dplyr::group_by(sample) %>%
        dplyr::arrange(sample, rmse) %>%
        dplyr::mutate(predTrue = as.logical(.pred_class),
                      minTrue = if(any(predTrue & !is.na(rmse))) min(rmse[predTrue], na.rm = TRUE) else NA_real_,
                      use = ifelse(
                        sumt <= 1,
                        predTrue,
                        predTrue & !is.na(rmse) & rmse == minTrue
                      )) %>%
        dplyr::select(-predTrue,-minTrue)
    }
  )

  qctable <- qctable %>%
    dplyr::mutate(flag = ifelse(triageValue < flagThreshold, "LOWPROBFLAG", NA)) %>%
    dplyr::rename(pred_class = .pred_class,
                  pred_TRUE = .pred_TRUE,
                  pred_FALSE = .pred_FALSE)

  return(qctable)
}

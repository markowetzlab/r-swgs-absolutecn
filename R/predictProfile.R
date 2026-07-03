predictProfile <- function(qctable = NULL,
                           model = NULL,
                           method = "randforest",
                           flagThreshold = 0.74,
                           errorMetric = "clonality") {
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
          dplyr::mutate(use = ifelse(clonality == min(clonality), TRUE, FALSE))
      },
      "segvariance" = {
        qctable <- qctable %>%
          dplyr::group_by(sample) %>%
          dplyr::arrange(sample, segvariance) %>%
          dplyr::mutate(use = ifelse(segvariance == min(segvariance), TRUE, FALSE))
      },
      "rmse" = {
        qctable <- qctable %>%
          dplyr::group_by(sample) %>%
          dplyr::arrange(sample, rmse) %>%
          dplyr::mutate(use = ifelse(rmse == min(rmse), TRUE, FALSE))
      }
    )
    newCols <- c("pred_class",
                 "pred_FALSE",
                 "pred_TRUE",
                 "triageValue",
                 "flag")
    qctable <- qctable %>%
      dplyr::mutate(!!!setNames(rep(NA, length(newCols)), newCols))
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
        errorMetric
      )
    )

  # tie break for edge cases
  ## For cases with multiple where ploidy is n + 2CN and twice the purity
  ## have the same error metric values. Select the lowest ploidy of accepted fits.
  if(sum(qctable$use) > 1){
    qctable <- qctable %>%
      dplyr::mutate(use = ifelse(use == TRUE & ploidy == min(ploidy),TRUE,FALSE))
    if(sum(qctable$use) > 1){
      stop("multiple fits selected - unknown tiebreak failure")
    }
  }

  if (all(!as.logical(qctable$use))) {
    qctable <- qctable %>%
      dplyr::ungroup() %>%
      tibble::add_row(
        SAMPLE_ID = unique(.$SAMPLE_ID),
        PATIENT_ID = unique(.$PATIENT_ID),
        ploidy = 2.0,
        purity = 1,
        smooth = unique(.$smooth),
        segments = unique(.$segments),
        downsample_depth = 2756274,
        # default for pl=2,pu=1 @ 30kb
        use = "TRUE"
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

  switch(
    errorMetric,
    "clonality" = {
      qctable <- qctable %>%
        dplyr::group_by(sample) %>%
        dplyr::arrange(sample, clonality) %>%
        dplyr::mutate(use = ifelse(
          sumt <= 1,
          as.logical(.pred_class),
          ifelse(clonality == min(clonality), TRUE, FALSE)
        ))
    },
    "segvariance" = {
      qctable <- qctable %>%
        dplyr::arrange(sample, segvariance) %>%
        dplyr::mutate(use = ifelse(
          sumt <= 1,
          as.logical(.pred_class),
          ifelse(segvariance == min(segvariance), TRUE, FALSE)
        ))
    },
    "rmse" = {
      qctable <- qctable %>%
        dplyr::group_by(sample) %>%
        dplyr::arrange(sample, rmse) %>%
        dplyr::mutate(use = ifelse(
          sumt <= 1,
          as.logical(.pred_class),
          ifelse(rmse == min(rmse), TRUE, FALSE)
        ))
    }
  )

  qctable <- qctable %>%
    dplyr::mutate(flag = ifelse(triageValue < flagThreshold, "LOWPROBFLAG", NA)) %>%
    dplyr::rename(pred_class = .pred_class,
                  pred_TRUE = .pred_TRUE,
                  pred_FALSE = .pred_FALSE) %>%
    dplyr::mutate(use = pred_class) %>%
    dplyr::mutate(notes = ifelse(is.na(flag), NA, flag))

  return(qctable)
}

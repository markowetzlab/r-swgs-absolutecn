# Filtered fits for one sample - lowest clonality is 2.1, lowest rmse is 1.9
# and lowest segvariance is 3.0
makeQC <- function() {
  ploidy <- c(1.9, 2.1, 3.0, 3.9)
  purity <- c(0.5, 0.6, 0.8, 0.95)
  data.frame(SAMPLE_ID = "S1", PATIENT_ID = "P1",
             ploidy = ploidy,
             purity = purity,
             segments = 80,
             clonality = c(0.08, 0.05, 0.12, 0.20),
             rmse = c(0.07, 0.10, 0.15, 0.25),
             segvariance = c(0.03, 0.02, 0.01, 0.05),
             homozygousLoss = c(0, 0, 30000, 90000),
             downsample_depth = getDownsampleDepth(ploidy = ploidy, purity = purity,
                                                   nbins_ref_genome = 1000),
             smooth = FALSE, TP53freq = 0.5, precPloidy = NA, precPurity = NA)
}

selectedPloidy <- function(res) {
  res$ploidy[res$use]
}

test_that("errorOnly selects the fit with the lowest error metric", {
  qc <- makeQC()
  expect_equal(selectedPloidy(predictProfile(qc, method = "errorOnly", errorMetric = "clonality")), 2.1)
  expect_equal(selectedPloidy(predictProfile(qc, method = "errorOnly", errorMetric = "rmse")), 1.9)
  expect_equal(selectedPloidy(predictProfile(qc, method = "errorOnly", errorMetric = "segvariance")), 3.0)
})

test_that("errorOnly output has empty model columns and autofit notes", {
  res <- predictProfile(makeQC(), method = "errorOnly")
  expect_true(all(c("pred_class", "pred_FALSE", "pred_TRUE", "triageValue", "flag") %in% colnames(res)))
  expect_true(all(is.na(res$flag)))
  expect_false("sample" %in% colnames(res))
  expect_equal(unique(res$notes),
               "autofit|fitmethod=errorOnly|flagThreshold=0.84|ErrorMetric=clonality")
  expect_type(res$use, "logical")
})

test_that("equal error metrics are broken by the lowest ploidy", {
  qc <- makeQC()
  qc$clonality <- c(0.05, 0.08, 0.12, 0.05)
  res <- predictProfile(qc, method = "errorOnly")
  expect_equal(sum(res$use), 1)
  expect_equal(selectedPloidy(res), 1.9)
})

test_that("a sample with no usable fit gets a flagged ploidy 2 purity 1 default", {
  qc <- makeQC()
  qc$clonality <- NA
  res <- predictProfile(qc, method = "errorOnly")
  expect_equal(nrow(res), nrow(qc) + 1)
  expect_equal(sum(res$use), 1)

  nofit <- res[res$use, ]
  expect_equal(nofit$ploidy, 2)
  expect_equal(nofit$purity, 1)
  expect_equal(nofit$flag, "NOFIT")
  expect_true(startsWith(nofit$notes, "NO_FIT_FOUND|autofit|"))
  expect_equal(nofit$SAMPLE_ID, "S1")
  expect_equal(nofit$TP53freq, 0.5)
  expect_equal(nofit$downsample_depth,
               round(getDownsampleDepth(ploidy = 2, purity = 1, nbins_ref_genome = 1000)))
})

test_that("randforest adds model predictions and selects one fit", {
  res <- predictProfile(makeQC(), model = swgs_abs_rfranger_model, method = "randforest")
  expect_true(all(c("pred_class", "pred_FALSE", "pred_TRUE", "triageValue", "flag") %in% colnames(res)))
  expect_false("sumt" %in% colnames(res))
  expect_equal(res$pred_FALSE + res$pred_TRUE, rep(1, nrow(res)), tolerance = 0.002)
  expect_equal(res$triageValue, abs(res$pred_TRUE - res$pred_FALSE), tolerance = 0.002)
  expect_equal(sum(res$use), 1)

  predTrue <- as.logical(res$pred_class)
  if (any(predTrue)) {
    # lowest clonality among fits predicted TRUE
    expect_equal(res$clonality[res$use], min(res$clonality[predTrue]))
  } else {
    expect_equal(res$flag[res$use], "NOFIT")
  }
})

test_that("randforest flags fits below the probability threshold", {
  res <- predictProfile(makeQC(), model = swgs_abs_rfranger_model, method = "randforest",
                        flagThreshold = 0.84)
  fitted <- !grepl("NOFIT", res$flag)
  expect_equal(!is.na(res$flag[fitted]), res$triageValue[fitted] < 0.84)

  res <- predictProfile(makeQC(), model = swgs_abs_rfranger_model, method = "randforest",
                        flagThreshold = 0)
  expect_true(all(is.na(res$flag) | res$flag == "NOFIT"))
})

test_that("predictProfile validates its inputs", {
  qc <- makeQC()
  expect_error(predictProfile(qc, model = NULL, method = "randforest"), "no model")
  expect_error(predictProfile(qc, method = "svm"), "unknown method")
  expect_error(predictProfile(qc, method = "errorOnly", errorMetric = "mae"), "unknown error metric")
  expect_error(predictProfile(qc[0, ], method = "errorOnly"), "no fits in qctable")
  expect_error(predictProfile(NULL, method = "errorOnly"), "no fits in qctable")
})

# triageProfile is internal - the model output columns are set directly
makeTriage <- function(predClass = c("TRUE", "TRUE", "FALSE")) {
  data.frame(sample = "S1",
             ploidy = c(2, 3, 4),
             clonality = c(0.10, 0.05, 0.02),
             rmse = c(0.30, 0.20, 0.10),
             segvariance = c(0.10, 0.30, 0.05),
             .pred_class = factor(predClass, levels = c("FALSE", "TRUE")),
             .pred_FALSE = c(0.2, 0.45, 0.9),
             .pred_TRUE = c(0.8, 0.55, 0.1))
}

test_that("triageProfile picks the lowest error among fits predicted TRUE", {
  qc <- makeTriage()
  # ploidy 4 has the lowest errors but is predicted FALSE
  expect_equal(selectedPloidy(triageProfile(qc, errorMetric = "clonality")), 3)
  expect_equal(selectedPloidy(triageProfile(qc, errorMetric = "rmse")), 3)
  expect_equal(selectedPloidy(triageProfile(qc, errorMetric = "segvariance")), 2)
})

test_that("triageProfile defaults to clonality", {
  expect_equal(selectedPloidy(triageProfile(makeTriage())), 3)
})

test_that("triageProfile uses the single TRUE prediction or none", {
  expect_equal(selectedPloidy(triageProfile(makeTriage(c("TRUE", "FALSE", "FALSE")))), 2)
  expect_false(any(triageProfile(makeTriage(c("FALSE", "FALSE", "FALSE")))$use))
})

test_that("triageProfile flags low confidence predictions and renames columns", {
  res <- triageProfile(makeTriage(), flagThreshold = 0.5)
  # triage values are 0.6, 0.1 and 0.8
  expect_equal(res$flag[order(res$ploidy)], c(NA, "LOWPROBFLAG", NA))
  expect_true(all(c("pred_class", "pred_TRUE", "pred_FALSE") %in% colnames(res)))
  expect_error(triageProfile(NULL), "no data")
})

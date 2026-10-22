test_that("column name data are character vectors without duplicates", {
  for (cols in list(fittingColumnNames, dsFittingColumnNames, rel2absColumnNames)) {
    expect_type(cols, "character")
    expect_false(anyDuplicated(cols) > 0)
    expect_true("SAMPLE_ID" %in% cols)
  }
})

test_that("the autofit model uses features produced by gridsearch", {
  predictors <- swgs_abs_rfranger_model$fit$forest$independent.variable.names
  expect_s3_class(swgs_abs_rfranger_model, "model_fit")
  expect_setequal(predictors, c("segments", "clonality", "homozygousLoss",
                                "segvariance", "ploidy", "purity"))
  expect_true(all(predictors %in% fittingColumnNames))
  expect_equal(swgs_abs_rfranger_model$lvl, c("FALSE", "TRUE"))
})

test_that("the autofit model does not capture the training environment", {
  spec <- swgs_abs_rfranger_model$spec
  quos <- Filter(rlang::is_quosure, c(spec$args, spec$eng_args, spec$method$fit$args))
  for (q in quos) {
    expect_identical(rlang::quo_get_env(q), rlang::empty_env())
  }
  # the captured environment held the training data and tuning results (3.8 MB)
  expect_lt(length(serialize(swgs_abs_rfranger_model, NULL)), 1e6)
})

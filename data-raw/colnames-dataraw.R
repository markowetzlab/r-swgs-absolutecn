## colnames datasets
fittingColumnNames <- c(
  "SAMPLE_ID",
  "ploidy",
  "purity",
  "segments",
  "clonality",
  "rmse",
  "segvariance",
  "downsample_depth",
  "powered",
  "TP53cn",
  "expected_TP53_AF",
  "homozygousLoss"
)

usethis::use_data(fittingColumnNames, overwrite = TRUE)
rm(fittingColumnNames)

rel2absColumnNames <- c(
  "SAMPLE_ID",
  "PATIENT_ID",
  "ploidy",
  "purity",
  "segments",
  "TP53cn",
  "expected_TP53_AF",
  "TP53freq",
  "clonality",
  "rmse",
  "segvariance"
)
usethis::use_data(rel2absColumnNames, overwrite = TRUE)
rm(rel2absColumnNames)

dsFittingColumnNames <- c(
  "SAMPLE_ID",
  "PATIENT_ID",
  "ploidy",
  "purity",
  "precPloidy",
  "precPurity",
  "TP53freq",
  "segments.pre",
  "segments.post",
  "clonality.pre",
  "clonality.post",
  "rmse.pre",
  "rmse.post",
  "segvariance.pre",
  "segvariance.post",
  "paired.ends",
  "total.reads",
  "used.reads",
  "expected.variance",
  "loess.span",
  "loess.family",
  "downsample_depth",
  "powered",
  "TP53cn.pre",
  "TP53cn.post",
  "expected_TP53_AF.pre",
  "expected_TP53_AF.post",
  "smooth",
  "homozygousLoss",
  "rank_clonality",
  "pl_diff",
  "new_state_n",
  "new_state",
  "use",
  "notes",
  "pred_class",
  "pred_FALSE",
  "pred_TRUE",
  "triageValue",
  "flag"
)
usethis::use_data(dsFittingColumnNames, overwrite = TRUE)
rm(dsFittingColumnNames)


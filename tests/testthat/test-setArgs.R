# Argument lists match those passed by workflow/argo-auto-workflow.yaml in
# swgs-absolutecn so a removed or renamed option fails here first
parseArgs <- function(script, args) {
  optparse::parse_args(setArgs(script = script), args = args)
}

test_that("setArgs requires a script name", {
  expect_error(setArgs(), "no script name provided")
})

test_that("bamCheck arguments parse", {
  opts <- parseArgs("bamCheck", c("--bam", "S1.bam", "--fileType", "BAM", "--okOut", "S1_bam.ok"))
  expect_equal(opts$bam, "S1.bam")
  expect_equal(opts$fileType, "BAM")
  expect_equal(opts$okOut, "S1_bam.ok")
})

test_that("qdnaseqMod arguments parse", {
  opts <- parseArgs("qdnaseqMod", c("--bin", "30", "--sampleName", "007", "--meta", "meta.tsv",
                                    "--bam", "S1.bam", "--rdsOut", "S1.rds", "--genome", "hg19",
                                    "--autoSmooth", "FALSE", "--smoothThreshold", "300",
                                    "--use_seed", "TRUE", "--seed_val", "9999",
                                    "--pairedEnd", "FALSE"))
  expect_equal(opts$bin, 30)
  # sample IDs stay character so leading zeros are kept
  expect_equal(opts$sampleName, "007")
  expect_false(opts$autoSmooth)
  expect_true(opts$use_seed)
  expect_equal(opts$seed_val, "9999")
})

test_that("gridsearch arguments parse", {
  opts <- parseArgs("gridsearch", c("--rds", "S1.rds", "--meta", "meta.tsv", "--pdfOut", "S1.pdf",
                                    "--tsvOut", "S1.tsv", "--plotOut", "S1.png",
                                    "--filtOut", "S1_filt.tsv", "--fitOut", "S1_fit.tsv",
                                    "--bin", "30", "--project", "PRJ", "--genome", "hg38",
                                    "--pl_min", "1.6", "--pl_max", "8", "--pu_min", "0.15",
                                    "--pu_max", "1", "--sampleName", "S1", "--af_cutoff", "0.15",
                                    "--filter_underpowered", "FALSE", "--filter_homozygous", "TRUE",
                                    "--homozygous_prop", "5e7", "--homozygous_threshold", "0.4"))
  expect_equal(opts$pl_min, 1.6)
  expect_equal(opts$af_cutoff, 0.15)
  expect_false(opts$filter_underpowered)
  expect_true(opts$filter_homozygous)
  expect_equal(opts$homozygous_prop, 5e7)
  expect_equal(opts$homozygous_threshold, 0.4)
  expect_equal(opts$project, "PRJ")
})

test_that("autofit arguments parse", {
  opts <- parseArgs("autofit", c("--tsv", "S1.tsv", "--tsvOut", "S1_fitted.tsv",
                                 "--flagThreshold", "0.84", "--fitMethod", "randforest",
                                 "--errorMetric", "clonality"))
  expect_equal(opts$flagThreshold, 0.84)
  expect_equal(opts$fitMethod, "randforest")
})

test_that("downsampleBams arguments parse", {
  opts <- parseArgs("downsampleBams", c("--bam", "S1.bam", "--meta", "S1.tsv", "--rds", "S1.rds",
                                        "--bamOut", "S1_ds.bam", "--sampleName", "S1",
                                        "--prplpu", "FALSE", "--fileType", "CRAM",
                                        "--reference", "ref.fa"))
  expect_false(opts$prplpu)
  expect_equal(opts$reference, "ref.fa")
})

test_that("qdnaseqRelToAbs arguments parse", {
  opts <- parseArgs("qdnaseqRelToAbs", c("--rds", "S1.rds", "--meta", "S1.tsv", "--genome", "hg19",
                                         "--sampleName", "S1", "--tsvOut", "S1.tsv",
                                         "--rdsOut", "S1_abs.rds", "--tabOut", "S1_seg.tsv",
                                         "--plotOut", "S1.png"))
  expect_equal(opts$tabOut, "S1_seg.tsv")
  # the script does not use bin or project
  expect_error(parseArgs("qdnaseqRelToAbs", c("--bin", "30")), "bin")
  expect_error(parseArgs("qdnaseqRelToAbs", c("--project", "PRJ")), "project")
})

test_that("splitInput arguments parse", {
  opts <- parseArgs("splitInput", c("--tsv", "all.tsv", "--tsvOut", "S1.tsv", "--sampleName", "S1"))
  expect_equal(opts$tsvOut, "S1.tsv")
})

test_that("combineResults takes comma-separated paths as one string", {
  opts <- parseArgs("combineResults", c("--tsv", "a.tsv,b.tsv", "--rds", "a.rds,b.rds",
                                        "--tab", "a_seg.tsv,b_seg.tsv", "--tsvOut", "all.tsv",
                                        "--rdsOut", "all.rds", "--tabOut", "all_seg.tsv"))
  expect_equal(opts$tsv, "a.tsv,b.tsv")
  expect_equal(strsplit(opts$rds, ",")[[1]], c("a.rds", "b.rds"))
})

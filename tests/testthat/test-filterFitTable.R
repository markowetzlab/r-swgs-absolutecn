# Gridsearch fits for one sample - the best fit per ploidy is
# 2.0/0.6 (0.05), 2.1/0.6 (0.07), 3.0/0.8 (0.08) and 4.0/0.9 (0.20)
makeFits <- function(sample = "S1") {
  data.frame(SAMPLE_ID = sample,
             ploidy = c(2.0, 2.0, 2.1, 3.0, 4.0),
             purity = c(0.5, 0.6, 0.6, 0.8, 0.9),
             clonality = c(0.10, 0.05, 0.07, 0.08, 0.20),
             powered = c(1, 1, 1, 0, 1),
             homozygousLoss = c(0, 0, 0, 30000, 90000),
             expected_TP53_AF = c(0.50, 0.55, 0.52, 0.70, NA))
}

makeMeta <- function(sample = "S1", TP53freq = 0.5) {
  data.frame(SAMPLE_ID = sample, PATIENT_ID = "P1", TP53freq = TP53freq, smooth = FALSE,
             precPloidy = NA, precPurity = NA, file = paste0(sample, ".bam"))
}

runFilter <- function(table = makeFits(), metadata = makeMeta(), filter_underpowered = FALSE,
                      filter_homozygous = FALSE, af_cutoff = 0.15, homozygous_prop = 1e8,
                      ranks = 10) {
  filterFitTable(table = table, metadata = metadata,
                 filter_underpowered = filter_underpowered,
                 filter_homozygous = filter_homozygous, af_cutoff = af_cutoff,
                 homozygous_prop = homozygous_prop, ranks = ranks)
}

test_that("filterFitTable keeps the best fit per ploidy within the TP53 cutoff", {
  res <- runFilter()
  # 3.0 is 0.20 from TP53freq and 4.0 has no expected AF
  expect_equal(sort(res$filtered$ploidy), c(2.0, 2.1))
  expect_equal(res$filtered$purity[res$filtered$ploidy == 2.0], 0.6)
  expect_false("file" %in% colnames(res$filtered))
})

test_that("filterFitTable prunes ploidies within 0.3 of each other", {
  res <- runFilter()
  expect_equal(res$pruned$ploidy, 2.0)
  expect_equal(res$pruned$purity, 0.6)
  expect_true(all(is.na(res$pruned$use)))
  expect_true(all(is.na(res$pruned$notes)))
})

test_that("af_cutoff = 1 disables the TP53 filter", {
  res <- runFilter(af_cutoff = 1)
  expect_equal(sort(res$filtered$ploidy), c(2.0, 2.1, 3.0, 4.0))
  expect_equal(res$pruned$ploidy, c(2.0, 3.0, 4.0))
})

test_that("samples without a TP53 VAF skip the TP53 filter", {
  res <- runFilter(metadata = makeMeta(TP53freq = NA))
  expect_equal(sort(res$filtered$ploidy), c(2.0, 2.1, 3.0, 4.0))
})

test_that("ranks limits the number of ploidy states and keeps ties", {
  res <- runFilter(af_cutoff = 1, ranks = 2)
  expect_equal(sort(res$filtered$ploidy), c(2.0, 2.1))

  tied <- makeFits()
  tied$clonality[tied$ploidy == 3.0] <- 0.07
  res <- runFilter(table = tied, af_cutoff = 1, ranks = 2)
  expect_equal(sort(res$filtered$ploidy), c(2.0, 2.1, 3.0))
})

test_that("ranks is applied within each sample", {
  other <- makeFits("S2")
  other$clonality <- other$clonality / 10
  res <- runFilter(table = rbind(makeFits(), other),
                   metadata = rbind(makeMeta(), makeMeta("S2")),
                   af_cutoff = 1, ranks = 2)
  expect_equal(as.vector(table(res$filtered$SAMPLE_ID)), c(2, 2))
})

test_that("filter_underpowered removes unpowered fits", {
  res <- runFilter(af_cutoff = 1, filter_underpowered = TRUE)
  expect_equal(sort(res$filtered$ploidy), c(2.0, 2.1, 4.0))

  unpowered <- makeFits()
  unpowered$powered <- 0
  expect_error(runFilter(table = unpowered, filter_underpowered = TRUE),
               "filter_underpowered")
})

test_that("filter_homozygous removes fits above the homozygous loss limit", {
  res <- runFilter(af_cutoff = 1, filter_homozygous = TRUE, homozygous_prop = 30000)
  expect_equal(sort(res$filtered$ploidy), c(2.0, 2.1, 3.0))

  expect_error(runFilter(filter_homozygous = TRUE, homozygous_prop = -1),
               "filter_homozygous")
})

test_that("filterFitTable validates its inputs", {
  expect_error(runFilter(table = NULL), "no data")
  expect_error(runFilter(metadata = NULL), "no metadata")
  expect_error(runFilter(af_cutoff = 0))
  expect_error(runFilter(af_cutoff = 1.5))
  expect_error(runFilter(ranks = 1))
  expect_error(runFilter(ranks = 2.5))
  expect_error(runFilter(filter_underpowered = "FALSE"))
  expect_error(runFilter(homozygous_prop = "1"))
})

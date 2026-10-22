test_that("bamCheck stops when the BAM or output name is missing", {
  expect_error(bamCheck(x = NULL, filetype = "BAM", outname = "S1_bam.ok"), "no BAMs provided")

  missing <- file.path(withr::local_tempdir(), "missing.bam")
  expect_error(bamCheck(x = missing, filetype = "BAM", outname = "S1_bam.ok"),
               missing, fixed = TRUE)

  bam <- withr::local_tempfile(fileext = ".bam", lines = "not a bam")
  expect_error(bamCheck(x = bam, filetype = "BAM", outname = NULL), "no outname provided")
})

test_that("bamCheck removes stale outputs before checking the file type", {
  bam <- withr::local_tempfile(fileext = ".bam", lines = "not a bam")
  dir <- withr::local_tempdir()
  okOut <- file.path(dir, "S1_bam.ok")
  invalidOut <- file.path(dir, "S1_bam.invalid")
  writeLines("stale", okOut)
  writeLines("stale", invalidOut)

  expect_error(bamCheck(x = bam, filetype = "SAM", outname = okOut))
  expect_false(file.exists(okOut))
  expect_false(file.exists(invalidOut))
})

test_that("bamCheck writes an invalid file for a corrupt BAM", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("samtools")), "samtools not available")

  bam <- withr::local_tempfile(fileext = ".bam", lines = "not a bam")
  dir <- withr::local_tempdir()
  okOut <- file.path(dir, "S1_bam.ok")
  invalidOut <- file.path(dir, "S1_bam.invalid")

  bamCheck(x = bam, filetype = "BAM", outname = okOut)
  expect_false(file.exists(okOut))
  expect_true(file.exists(invalidOut))
  expect_match(readLines(invalidOut), "BAM/CRAM invalid", fixed = TRUE)
})

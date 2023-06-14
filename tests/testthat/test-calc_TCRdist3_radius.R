test_that("multiplication works", {
  library(reticulate)
  file_path <- testthat::test_path("testdata", "clonosets_vdjtools_format.tsv")
  TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7)
  clonoset(TCRgrObject) <- clonoset(TCRgrObject)[1:100,]
  TCRgrObject <- calc_TCRdist3_radius(TCRgrObject)
  expect_equal('tcrdist3.radius' %in% colnames(clonoset(TCRgrObject)), TRUE)
})

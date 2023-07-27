test_that("TCRgrapherCounts works", {
  file_path <- testthat::test_path("testdata", "clonosets_vdjtools_format.tsv")
  TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7)
  clonoset(TCRgrObject) <- clonoset(TCRgrObject)[1:100,]
  TCRgrapherCountsObject <- TCRgrapherCounts(TCRgrObject, cluster_id = TRUE)
  expect_equal(nrow(TCRgrapherCountsObject@edges), 29)
})

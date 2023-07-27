test_that("find_TCR_components_by_bfs works", {
  file_path <- testthat::test_path("testdata", "clonosets_vdjtools_format.tsv")
  TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7)
  clonoset(TCRgrObject) <- clonoset(TCRgrObject)[1:100,]
  TCRgrObject <- find_TCR_components_by_bfs(TCRgrObject)
  expect_equal(length(unique((TCRgrObject@clonoset$cluster_id))), 71)
  expect_equal(nrow(TCRgrObject@edges), 29)
})

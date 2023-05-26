test_that("reading a single file works", {
  file_path <- testthat::test_path("testdata", "clonosets_vdjtools_format.tsv")
  TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7)
  expect_equal(attr(TCRgrObject@metadata, "class"), c("data.table", "data.frame"))
  expect_equal(attr(TCRgrObject@clonoset, "class"), c("data.table", "data.frame"))
  expect_equal(nrow(TCRgrObject@metadata), 1)
  expect_equal(ncol(TCRgrObject@metadata), 2)
  expect_equal(colnames(TCRgrObject@metadata), c('file', 'sample_id'))
  expect_equal(TCRgrObject@metadata$file, file_path)
  expect_equal(TCRgrObject@metadata$sample_id, file_path)
  clonoset_names <- c('count', 'cdr3nt', 'cdr3aa', 'bestVGene', 'bestJGene')
  expect_equal(all(clonoset_names %in% colnames(TCRgrObject@clonoset)), TRUE)
})

test_that("reading files from the directory without metadata works", {
  file_path <- testthat::test_path("testdata/test_dir_with_clonosets", "")
  TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7)
  expect_equal(attr(TCRgrObject@metadata, "class"), c("data.table", "data.frame"))
  expect_equal(attr(TCRgrObject@clonoset, "class"), c("data.table", "data.frame"))
  expect_equal(ncol(TCRgrObject@metadata), 2)
  expect_equal(nrow(TCRgrObject@metadata), length(list.files(file_path)))
  expect_equal(colnames(TCRgrObject@metadata), c('file', 'sample_id'))
  expect_equal(TCRgrObject@metadata$file, list.files(file_path))
  clonoset_names <- c('count', 'cdr3nt', 'cdr3aa', 'bestVGene', 'bestJGene')
  expect_equal(all(clonoset_names %in% colnames(TCRgrObject@clonoset)), TRUE)
})

test_that("reading files from the directory with metadata works", {
  file_path <- testthat::test_path("testdata/test_dir_with_clonosets", "")
  metadata_path <- testthat::test_path("testdata/metadata.tsv")
  TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7, metadata_path, 1, 2)
  expect_equal(attr(TCRgrObject@metadata, "class"), c("data.table", "data.frame"))
  expect_equal(attr(TCRgrObject@clonoset, "class"), c("data.table", "data.frame"))
  #expect_equal(ncol(TCRgrObject@metadata), 3)
  expect_equal(nrow(TCRgrObject@metadata), length(list.files(file_path)))
  expect_equal(all(c('file', 'sample_id') %in% colnames(TCRgrObject@metadata)), TRUE)
  expect_equal(TCRgrObject@metadata$file, list.files(file_path))
  clonoset_names <- c('count', 'cdr3nt', 'cdr3aa', 'bestVGene', 'bestJGene')
  expect_equal(all(clonoset_names %in% colnames(TCRgrObject@clonoset)), TRUE)
})

test_that("getter for a clonoset works", {
  file_path <- testthat::test_path("testdata", "clonosets_vdjtools_format.tsv")
  TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7)
  test_clonoset <- clonoset(TCRgrObject)
  test_row <- data.table(count = 3, freq = 0.00016481705307109108,
                         cdr3nt = 'TGTGCCAGCTCACAAGACAGACTAAACACAGAAGTCTTCTTT',
                         cdr3aa = 'CASSQDRLNTEVFF', bestVGene = 'TRBV12-2',
                         d = '', bestJGene = 'TRBJ1-1',
                         sample_id = 'testdata/clonosets_vdjtools_format.tsv',
                         clone_id = 1)
  expect_equal(test_clonoset[1,], test_row)
})

test_that("setter for a clonoset works", {
  file_path <- testthat::test_path("testdata", "clonosets_vdjtools_format.tsv")
  TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7)
  clonoset(TCRgrObject) <- clonoset(TCRgrObject)[,test_col := 'test',]
  test_row <- data.table(count = 3, freq = 0.00016481705307109108,
                         cdr3nt = 'TGTGCCAGCTCACAAGACAGACTAAACACAGAAGTCTTCTTT',
                         cdr3aa = 'CASSQDRLNTEVFF', bestVGene = 'TRBV12-2',
                         d = '', bestJGene = 'TRBJ1-1',
                         sample_id = 'testdata/clonosets_vdjtools_format.tsv',
                         clone_id = 1, test_col = 'test')
  expect_equal(clonoset(TCRgrObject)[1,], test_row)
})

test_that("validity method works", {
  file_path <- testthat::test_path("testdata", "clonosets_vdjtools_format.tsv")
  TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7)
  expect_error(clonoset(TCRgrObject)[1,'count'] <- -1)
})

test_that("subsetting method works", {
  file_path <- testthat::test_path("testdata/test_dir_with_clonosets", "")
  metadata_path <- testthat::test_path("testdata/metadata.tsv")
  TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7, metadata_path, 1, 2)
  TCRgrObject <- subset(TCRgrObject, c('Sample_18', 'Sample_19'))
  expect_equal(length(unique(TCRgrObject@clonoset$sample_id)), 2)
  expect_equal(nrow(TCRgrObject@metadata), 2)
})

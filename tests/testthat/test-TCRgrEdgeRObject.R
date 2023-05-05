test_that("basics", {
  s1 <- data.table(count = c(10, 5, 1), cdr3nt = c('NT1', 'NT2', 'NT3'),
                   cdr3aa = c('AA1', 'AA2', 'AA3'), bestVGene = c('V1', 'V2', 'V3'),
                   bestJGene = c('J1', 'J2', 'J3'), sample_id = c('s1', 's1', 's1'))
  s2 <- data.table(count = c(10, 5, 1), cdr3nt = c('NT1', 'NT2', 'NT3'),
                   cdr3aa = c('AA1', 'AA2', 'AA3'), bestVGene = c('V1', 'V2', 'V3'),
                   bestJGene = c('J1', 'J2', 'J3'), sample_id = c('s2', 's2', 's2'))
  s3 <- data.table(count = c(10, 8, 5, 1), cdr3nt = c('NT1', 'NT4', 'NT2', 'NT3'),
                   cdr3aa = c('AA1', 'AA1', 'AA2', 'AA3'), bestVGene = c('V4', 'V4', 'V5', 'V6'),
                   bestJGene = c('J1', 'J4', 'J2', 'J3'), sample_id = c('s3', 's3', 's3', 's3'))
  TCRgrObject <- new('TCRgrapher',
                     metadata  = data.table(file = c('file1.txt', 'file2.txt', 'file3.txt'),
                                            sample_id = c('s1', 's2', 's3')),
                     clonoset = rbindlist(list(s1, s2, s3))
  )

  count_table <- list()
  feature_info <- list()
  metadata <- TCRgrObject@metadata
  clonoset <- TCRgrObject@clonoset
  v_gene = TRUE

  for(sample_t in metadata$sample_id){
    sample <- clonoset[sample_id == sample_t]
    if(!v_gene){
      sample <- sample[, .(count = sum(count)), by = .(cdr3aa)]
    } else {
      sample <- sample[, .(count = sum(count)), by = .(cdr3aa, bestVGene)]
    }
    setnames(sample, 'count', sample_t)
    count_table <- append(count_table, list(sample))
  }
  count_table <- Reduce(function(x, y) merge.data.table(x, y, by=c('cdr3aa', 'bestVGene'), all = TRUE),
                        count_table)
  # TODO
})

library(tcrgrapher)
library(data.table)
library(edgeR)
# for finding connectivity components of the graph
library(igraph)

# Let's find expanded clonotypes after vaccination. There are 8 mice and 5 time
# points: 2 points before vaccination and 3 after vaccination

# Data loading. See ?TCRgrapher
# One of the ways to read the data
file_path <- testthat::test_path("testdata/vaccinated_mice", "")
metadata_path <- testthat::test_path("testdata/metadata_vaccination.tsv")
TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7, metadata_path, 1, 2)

# see general information about the object
TCRgrObject
# see metadata
metadata(TCRgrObject)
# add column to metadata
metadata(TCRgrObject)$vaccination <- 'before'
# change some values
metadata(TCRgrObject)[metadata(TCRgrObject)$time > 2, 'vaccination'] <- 'after'

# create a count table with aggregation by V segments
TCRgrCounts <- TCRgrapherCounts(TCRgrObject)
# run EdgeR pipeline
# There are two comparison levels: before and after vaccination
edgeR_res <- edgeR_pipeline(TCRgrCounts, 'vaccination')
unique(edgeR_res$comparison)
head(edgeR_res)

# example to show more than two comparison levels
edgeR_res_time <- edgeR_pipeline(TCRgrCounts, 'time')
unique(edgeR_res_time$comparison)

# To take a subset we should specify samples that we want to take
samples <-  metadata(TCRgrCounts)[time <= 3, 'sample_id']
TCRgrCounts_3 <- subset(TCRgrCounts, unlist(samples))
# check count table columns and metadata sample_id values
colnames(count_table(TCRgrCounts_3))
metadata(TCRgrCounts_3)$sample_id
# run edgeR pipeline
edgeR_res_3 <- edgeR_pipeline(TCRgrCounts_3, 'vaccination')
head(edgeR_res_3)

# create a count table with aggregation by V and J segments
TCRgrCounts_VJ <- TCRgrapherCounts(TCRgrObject, j_gene = TRUE)
# take the subset again
samples <-  metadata(TCRgrCounts_VJ)[time <= 3, 'sample_id']
TCRgrCounts_VJ_3 <- subset(TCRgrCounts_VJ, unlist(samples))
# run edgeR pipeline
edgeR_res_VJ_3 <- edgeR_pipeline(TCRgrCounts_VJ_3, 'vaccination')
head(edgeR_res_VJ_3)

# create a count table with aggregation by one amino acid sequence
TCRgrCounts_aa <- TCRgrapherCounts(TCRgrObject, v_gene = FALSE, j_gene = FALSE)
# take the subset again
samples <-  metadata(TCRgrCounts_aa)[time <= 3, 'sample_id']
TCRgrCounts_aa_3 <- subset(TCRgrCounts_aa, unlist(samples))
# run edgeR pipeline
edgeR_res_aa_3 <- edgeR_pipeline(TCRgrCounts_aa_3, 'vaccination')
head(edgeR_res_aa_3)

# create a count table by clusters
# the example will run for a veeery long time
samples <-  metadata(TCRgrObject)[time <= 1, 'sample_id']
TCRgrObject_3 <- subset(TCRgrObject, unlist(samples))
TCRgrCounts_clusters <- TCRgrapherCounts(TCRgrObject_3, cluster_id = TRUE)
# if you want to run edgeR pipeline with aggregation by clusters, you can do the following
edgeR_res_clusters <- edgeR_pipeline(TCRgrCounts_clusters, 'vaccination')
head(edgeR_res_clusters)


library(tcrgrapher)
library(data.table)
library(reticulate)

TCRgrObject <- TCRgrapher('Tub_mice_merged_data.tsv', 2, 24, 33, 6, 8)
clonoset(TCRgrObject)$bestVGene <- sapply(strsplit(clonoset(TCRgrObject)$bestVGene, '\\*'), function(x) x[[1]][1])
clonoset(TCRgrObject)$bestJGene <- sapply(strsplit(clonoset(TCRgrObject)$bestJGene, '\\*'), function(x) x[[1]][1])

TCRgrObject <- ALICE_pipeline(TCRgrObject, cores = cores, thres_counts = 1,
                              N_neighbors_thres = 0)

TCRgrObject <- tcrgrapher::run_TCRNET(TCRgrObject,
                                      'TCRnet_control/C57BL6_all_genes_with_counts.pool.aaVJ.table.txt',
                                      command = 'java -jar vdjtools-1.2.1.jar')

TCRgrObject <- run_GLIPH2(TCRgrObject, 'irtools.centos')

TCRgrObject <- calc_TCRdist3_radius(TCRgrObject, cores = cores)

df <- apply(clonoset(TCRgrObject), 2, as.character)

write.table(df, quote = F, sep = '\t', row.names = F, file ='result.tsv')

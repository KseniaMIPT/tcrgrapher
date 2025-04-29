library(tcrgrapher)
library(data.table)
library(reticulate)
library(stringr)

# You may need to specify path to python with installed OLGA to run ALICE
Sys.setenv(PATH="/home/klupyr/.conda/envs/statbiophys/bin/")

# specify your number of cores
cores = 2

# data loading
TCRgrObject <- TCRgrapher('Tub_mice_data.tsv', 2, 24, 33, 6, 8)
clonoset(TCRgrObject)$bestVGene <- sapply(strsplit(clonoset(TCRgrObject)$bestVGene, '\\*'), function(x) x[[1]][1])
clonoset(TCRgrObject)$bestJGene <- sapply(strsplit(clonoset(TCRgrObject)$bestJGene, '\\*'), function(x) x[[1]][1])

# run ALICE
TCRgrObject <- ALICE_pipeline(TCRgrObject, cores = cores, thres_counts = 1, N_neighbors_thres = 0)

# run TCRNET
# you need to write your path instead of
# 'TCRnet_control/olga_generated/C57BL6_all_genes_with_counts.pool.aaVJ.table.txt'
# and your own path to vdjtools
TCRgrObject <- tcrgrapher::run_TCRNET(TCRgrObject,
                                        'TCRnet_control/olga_generated/C57BL6_all_genes_with_counts.pool.aaVJ.table.txt',
                                        command = 'java -jar vdjtools-1.2.1.jar')

# run TCRdist3
TCRgrObject <- calc_TCRdist3_radius(TCRgrObject, cores = cores)

# run GLIPH2
# specify your own path to GLIPH2 downloded from http://50.255.35.37:8080/tools
TCRgrObject <- run_GLIPH2(TCRgrObject, '~/soft/irtools.centos')

write.table(clonoset(TCRgrObject), quote = F, sep = '\t', row.names = F,
            file = paste0('result.tsv'))

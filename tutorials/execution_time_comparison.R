# Write your path to the directory with ALICE https://github.com/pogorely/ALICE/tree/master
setwd('/home/klupyr/soft/ALICE/')
source('/home/klupyr/soft/ALICE/ALICE.R')

# Write your path to the directory with OLGA
Sys.setenv(PATH="/home/klupyr/.conda/envs/statbiophys/bin/")

S1d15<-fread("sample/S1_d15_V9_J2_7.tsv")
# S1d0<-fread("sample/S1_d0_V9_J2_7.tsv")
S1<-list(d15=S1d15)

# ALICE 16 cores
start.time <- Sys.time()
S1_alice<-ALICE_pipeline_OLGA(DTlist=S1,cores=16)
end.time <- Sys.time()
time.taken_16_ALICE_no_OLGA <- end.time - start.time
print(time.taken_16_ALICE_no_OLGA)

# ALCIE 1 core
start.time <- Sys.time()
S1_alice<-ALICE_pipeline_OLGA(DTlist=S1,cores=1)
end.time <- Sys.time()
time.taken_1_ALICE_no_OLGA <- end.time - start.time
print(time.taken_1_ALICE_no_OLGA)

# ALICE 16 cores with OLGA
start.time <- Sys.time()
S1_alice<-ALICE_pipeline_OLGA(DTlist=S1,cores=16)
end.time <- Sys.time()
time.taken_16_ALICE <- end.time - start.time
print(time.taken_16_ALICE)

# ALCIE 1 core with OLGA
start.time <- Sys.time()
S1_alice<-ALICE_pipeline_OLGA(DTlist=S1,cores=1)
end.time <- Sys.time()
time.taken_1_ALICE <- end.time - start.time
print(time.taken_1_ALICE)

# libraries
library(tcrgrapher)
library(data.table)
library(reticulate)
library(stringr)

# data loading
TCRgrObject <- TCRgrapher("sample/S1_d15_V9_J2_7.tsv", 2, 4, 5, 6, 7)
nrow(clonoset(TCRgrObject))

TCRgrObject_2 <- TCRgrapher("sample/S1_d15_V9_J2_7.tsv", 2, 4, 5, 6, 7)
clonoset(TCRgrObject_2) <- do.call("rbind", replicate(10, clonoset(TCRgrObject_2), simplify = FALSE))
nrow(clonoset(TCRgrObject_2))

TCRgrObject_3 <- TCRgrapher("sample/S1_d15_V9_J2_7.tsv", 2, 4, 5, 6, 7)
clonoset(TCRgrObject_3) <- do.call("rbind", replicate(100, clonoset(TCRgrObject_3), simplify = FALSE))
nrow(clonoset(TCRgrObject_3))

samples <- list(TCRgrObject, TCRgrObject_2, TCRgrObject_3)

dt <- c()

for(nb_of_cores in c(1, 16)){
  for(size_of_sample in 1:3){
    if((nb_of_cores == 1 & size_of_sample == 1) | (nb_of_cores == 16)){
      sample_t <- samples[[size_of_sample]]

      # ALICE
      Sys.setenv(PATH="/home/klupyr/.conda/envs/statbiophys/bin/")
      start.time <- Sys.time()
      TCRgrObject <- tcrgrapher::ALICE_pipeline(sample_t, cores = nb_of_cores, N_neighbors_thres = 1, chain = 'humanTRB')
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      dt <- rbind(dt, c(nb_of_cores, size_of_sample, 'ALICE', time.taken))

      # TCRdist3
      start.time <- Sys.time()
      TCRgrObject <- calc_TCRdist3_radius(sample_t, cores = nb_of_cores)
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      print(paste(c(nb_of_cores, size_of_sample, 'TCRdist3', time.taken)))
      # dt <- rbind(dt, c(nb_of_cores, size_of_sample, 'TCRdist3', time.taken))

      # TCRNET
      Sys.setenv(PATH="/home/klupyr/.conda/envs/vdjtools/bin")
      start.time <- Sys.time()
      cmd = paste0('java -Xmx128G -XX:ActiveProcessorCount=', nb_of_cores, ' -jar /home/klupyr/soft/vdjtools-1.2.1.jar')
      TCRgrObject <- tcrgrapher::run_TCRNET(sample_t,
                                            '/home/klupyr/TCRnet_control/olga_generated/C57BL6_all_genes_with_counts.pool.aaVJ.table.txt',
                                            command = cmd)
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      dt <- rbind(dt, c(nb_of_cores, size_of_sample, 'TCRNET', time.taken))

      # GLIPH2
      TCRgrObject <- TCRgrapher("/home/klupyr/soft/ALICE/sample/S1_d15_V9_J2_7.tsv", 2, 4, 5, 6, 7)
      start.time <- Sys.time()
      TCRgrObject <- run_GLIPH2(TCRgrObject, '~/soft/irtools.centos')
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      dt <- rbind(dt, c(nb_of_cores, size_of_sample, 'GLIPH2', time.taken))
    }
  }
}

write.table(dt, file = '~/time.tsv', sep = '\t', quote = F, row.names = F)

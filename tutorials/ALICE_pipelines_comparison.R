# Write your path to the directory with ALICE https://github.com/pogorely/ALICE/tree/master
#setwd('/home/username/ALICE/')
#source('/home/username/ALICE/ALICE.R')

# Write your path to the directory with OLGA
#Sys.setenv(PATH="/home/username/.conda/envs/env_with_olga/bin/")

S1d15<-fread("sample/S1_d15_V9_J2_7.tsv")
S1d0<-fread("sample/S1_d0_V9_J2_7.tsv")
#S1<-list(d0=S1d0,d15=S1d15)
S1<-list(d15=S1d15)

# ALICE 16 cores
start.time <- Sys.time()
S1_alice<-ALICE_pipeline_OLGA(DTlist=S1,cores=16)
end.time <- Sys.time()
time.taken_16_ALICE_no_OLGA <- end.time - start.time
time.taken_16_ALICE_no_OLGA

# ALCIE 1 core
start.time <- Sys.time()
S1_alice<-ALICE_pipeline_OLGA(DTlist=S1,cores=16)
end.time <- Sys.time()
time.taken_1_ALICE_no_OLGA <- end.time - start.time
time.taken_1_ALICE_no_OLGA

# ALICE 16 cores with OLGA
start.time <- Sys.time()
S1_alice<-ALICE_pipeline_OLGA(DTlist=S1,cores=16)
end.time <- Sys.time()
time.taken_16_ALICE <- end.time - start.time
time.taken_16_ALICE

# ALCIE 1 core with OLGA
start.time <- Sys.time()
S1_alice<-ALICE_pipeline_OLGA(DTlist=S1,cores=16)
end.time <- Sys.time()
time.taken_1_ALICE <- end.time - start.time
time.taken_1_ALICE

# TCRgrapher 16 cores
library(tcrgrapher)
library(data.table)

TCRgrObject <- TCRgrapher("sample/S1_d15_V9_J2_7.tsv", 2, 4, 5, 6, 7)
start.time <- Sys.time()
TCRgrObject <- tcrgrapher::ALICE_pipeline(TCRgrObject, cores = 16, N_neighbors_thres = 1, chain = 'humanTRB')
end.time <- Sys.time()
time.taken_16_tcgrapher <- end.time - start.time
time.taken_16_tcgrapher

# TCRgrapher 1 cores

TCRgrObject <- TCRgrapher("sample/S1_d15_V9_J2_7.tsv", 2, 4, 5, 6, 7)
start.time <- Sys.time()
TCRgrObject <- tcrgrapher::ALICE_pipeline(TCRgrObject, cores = 1, thres_counts = 1,
                                          N_neighbors_thres = 1, p_adjust_method = "BH",
                                          chain = 'humanTRB', stats = 'OLGA', model= '-')
end.time <- Sys.time()
time.taken_1_tcgrapher <- end.time - start.time
time.taken_1_tcgrapher

print(time.taken_16_ALICE_no_OLGA)
print(time.taken_1_ALICE_no_OLGA)
print(time.taken_1_ALICE)
print(time.taken_16_ALICE)
print(time.taken_1_tcgrapher)
print(time.taken_16_tcgrapher)

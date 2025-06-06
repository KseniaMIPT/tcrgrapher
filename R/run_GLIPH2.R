#     TCRgrapher: R package for identifying condition associated T cell clonotypes
#     Copyright (C) 2024 Kseniia Lupyr
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @importFrom data.table merge.data.table
NULL
#' take_min_pval_gliph2
#'
#' The function takes clusters.csv from gliph2 output and assigns the minimum p-value
#' for every clonotype from the TCRgrapher clonotype table.
#'
#' @param gliph2_clusters clusters.csv from gliph2 output
#' @param TCRgrObject TCRgrapher object
#' @return The function returns a vector of minimum p-values
#' @export
take_min_pval_gliph2 <- function(gliph2_clusters, TCRgrObject){
  gliph2_clusters <- gliph2_clusters[,min(Fisher_score),by=.(TcRa, V)]
  gliph2_clusters <- merge.data.table(clonoset(TCRgrObject), gliph2_clusters,
                                      by.x = c('cdr3aa', 'bestVGene'),
                                      by.y = c('TcRa', 'V'), all.x = TRUE, sort=FALSE)
  gliph2_clusters <- gliph2_clusters$V1
  gliph2_clusters[is.na(gliph2_clusters)] <- 1
  return(gliph2_clusters)
}

#' run_GLIPH2
#'
#' The function takes the TCRgrapher object as an input and runs GLIPH2 analysis.
#' In GLIPH2 output, one sequence can be included in different clusters and
#' have different p-values. At the last step, the function assigns the minimum p-value
#' for every clonotype from the TCRgrapher clonotype table. If there are different types
#' of patterns in the analysis, the output will have a separate column for each
#' type of pattern. The GLIPH2 executable file must be installed.
#' See details at http://50.255.35.37:8080/tools
#'
#' @param TCRgrObject TCRgrapher object
#' @param gliph2_path path to the GLIPH2 executable file
#' @param control The possible values are 'TCRgrapher' and 'GLIPH2'. The default
#' value is 'TCRgrapher'. If you choose 'GLIPH2' the standard GLIPH2 control from the
#' server will be used. If you would like to use the control generated in our lab,
#' choose 'TCRgrapher'
#' @param kmer_min_depth One of the GLIPH2 parameters. The default value is 100000000.
#' We recommend not using clusters based on k-mers. The default value is chosen
#' to avoid its computation. If you want to use them, choose the small value
#' (kmer_min_depth = 3, for example). Keep in mind that it will also increase
#' computational time.
#' @return The function returns the same TCRgrapher object with additional columns
#' in the clonotype table: 'min_gliph2_fisher_score_global' and 'min_gliph2_fisher_score_(k)mer'
#' if k-mers were present.
#' @export
run_GLIPH2 <- function(TCRgrObject, gliph2_path, control = 'TCRgrapher',
                       kmer_min_depth = '100000000'){
  # convert to gliph2 format
  data <- clonoset(TCRgrObject)
  data$cond <- 'Subject:condition'
  data$freq <- data$count / sum(data$count)
  data <- data[,c('cdr3nt', 'bestVGene', 'bestJGene', 'cdr3aa', 'cond', 'freq')]

  path <- getwd()
  input_file <- paste0(path, '/input.tsv')
  write.table(data, file = input_file, sep='\t', quote = F, row.names = F, col.names = F)

  # control made from https://zenodo.org/record/6339774
  if(control == 'TCRgrapher'){
    control <- system.file('extdata/GLIPH2_control/', package = "tcrgrapher", mustWork = TRUE)
  } else if(control == 'GLIPH2'){
    control <- system.file('extdata/GLIPH2_control_server/', package = "tcrgrapher", mustWork = TRUE)
  }
  text = paste0('out_prefix=output\n',
                '\ncdr3_file=', input_file, '\n',
                'refer_file=', control, 'ref_CD48_ms.txt\n',
                'v_usage_freq_file=', control, 'ref_V_CD48_ms.txt\n',
                'cdr3_length_freq_file=', control, 'ref_L_CD48_ms.txt\n',
                'local_min_pvalue=0.001\n',
                'p_depth = 1000\n',
                'global_convergence_cutoff = 1\n',
                'simulation_depth=1000\n',
                'kmer_min_depth=', kmer_min_depth, '\n',
                'local_min_OVE=10\n',
                'algorithm=GLIPH2\n',
                'all_aa_interchangeable=1\n')
  script_path <- paste0(path, '/script_gliph2.txt')
  write(text, file = script_path)

  # rm old files if they exist
  old_files <- c('output_parameter.txt', 'output_cluster.csv', 'output_HLA.csv',
                 'output_score.txt', 'output_cluster.txt', 'script_gliph2.txt',
                 'input.tsv', 'output_kmer.log', 'output_kmer.txt')
  path_files <- list.files(path)
  for(file in old_files){
    if(file %in% path_files){
      system(paste0('rm ', path, '/', file))
    }
  }
  # run gliph2
  system(paste0(gliph2_path, ' -c ', script_path))
  # read output
  output_dt <- read.delim(paste0(path, '/output_cluster.csv'), sep=',')
  output_dt <- setDT(output_dt)
  # column with global scores
  output_dt_global <- output_dt[str_detect(type, 'global')]
  output_dt_kmers <- output_dt[!str_detect(type, 'global')]
  clonoset(TCRgrObject)$min_gliph2_fisher_score_global <-
    take_min_pval_gliph2(output_dt_global, TCRgrObject)
  # columns with kmers scores
  if(nrow(output_dt_kmers) != 0){
    for(k in unique(nchar(output_dt_kmers$pattern))){
      output_dt_kmer <- output_dt_kmers[nchar(pattern) == k]
      clonoset(TCRgrObject)[, paste0('min_gliph2_fisher_score_', k, 'mer')] <-
        take_min_pval_gliph2(output_dt_kmer, TCRgrObject)
    }
  }
  return(TCRgrObject)
}

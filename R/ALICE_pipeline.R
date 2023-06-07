#' @importFrom stats p.adjust
#' @importFrom stats ppois
#' @importFrom stats approxfun
#' @importFrom utils write.table
#' @importFrom stringdist stringdistmatrix
#' @importFrom data.table :=
#' @importFrom data.table rbindlist
#' @importFrom data.table fread
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
NULL
"OLGAVJ_MOUSE_TRB"
"OLGAVJ_HUMAN_TRB"
"OLGAVJ_HUMAN_TRA"
"C57BL6_MOUSE_TRB"

# Secondary functions -----------------------------------------------------

calculate_nb_of_neighbors_one_side <- function(DT) {
  # Every sequence with one mismatch is a neighbor
  tmp <- stringdistmatrix(DT$cdr3aa, DT$cdr3aa, method = "hamming")
  apply(tmp,
        MARGIN = 1,
        function(x) { sum(x <= 1, na.rm = T) })
}

calculate_nb_of_neighbors <- function(DT) {
  # give unique number 'leftgr' to the rows with the same VJ combination and left
  # CDR3aa part
  DT[, leftgr := .GRP, .(substr(cdr3aa, 1, floor(nchar(cdr3aa) / 2)), bestVGene, bestJGene)]
  # give unique number 'rightgr' to the rows with the same VJ combination and right
  # CDR3aa part
  DT[, rightgr := .GRP, .(substr(cdr3aa, floor(nchar(cdr3aa) / 2) + 1, nchar(cdr3aa)), bestVGene, bestJGene)]

  DT[, D_left := calculate_nb_of_neighbors_one_side(.SD), .(leftgr)]
  DT[, D_right := calculate_nb_of_neighbors_one_side(.SD), .(rightgr)]
  DT[, D_id := .N, .(cdr3aa)]
  DT[, D := (D_left + D_right - D_id - 1), ]
  DT <- subset(DT, select = -c(rightgr, leftgr, D_left, D_right, D_id))
  DT[D < 0, D := 0, ]
  DT
}

all_other_variants_one_mismatch_regexp <- function(str) {
  # All one mismatch variants with X (regexp!)
  unique(as.vector(sapply(
    2:(nchar(str) - 1),
    function(x) {
      tmp <- str
      substr(tmp, x, x) <- "X"
      tmp
    }
  )))
}

parallel_wrapper_beta <- function(DT, cores = 1, chain = "mouseTRB",
                                  stats = 'OLGA', model='-') {
  # Calculate generation probability with OLGA or SONIA

  # add ind column for sequence combining
  if (!("ind" %in% colnames(DT))) DT[, ind := 1:.N, ]

  tmp_names <- paste0("tmp", 1:cores, ".tsv")
  tmp_names_out <- paste0("tmp_out", 1:cores, ".tsv")
  path <- tempdir()
  path <- paste0(path, '/')

  for (f in c(paste0(path, tmp_names), paste0(path,tmp_names_out))) if (file.exists(f)) file.remove(f)

  DTt <- split(DT, sort((1:nrow(DT) - 1) %% cores + 1))

  for (i in 1:length(DTt)) {
    write.table(as.data.frame(DTt[[i]][, .(cdr3aa, bestVGene, bestJGene, ind), ]),
                quote = F, row.names = F, sep = "\t", file = paste0(path, tmp_names[i]), col.names = F
    )
  }

  if(model != '-'){
    chain=model
  }

  if(stats == 'OLGA'){
    commands <- paste0(
      "olga-compute_pgen --", chain,
      " --display_off --time_updates_off --seq_in 0 --v_in 1 --j_in 2 -d 'tab' -i ",
      path, "tmp", 1:cores, ".tsv -o ", path, "tmp_out", 1:cores, ".tsv"
    )
  } else if (stats == 'SONIA'){
    commands <- paste0('sonia-evaluate --', chain,
                       " --seq_in 0 --v_in 1 --j_in 2 -d 'tab' --ppost -i ", path, "tmp",
                       1:cores, ".tsv -o ", path, "tmp_out", 1:cores, ".tsv")
  }

  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  foreach(i=1:cores) %dopar% {
    system(commands[i], intern=TRUE)
  }
  parallel::stopCluster(cl)

  stats_output <- rbindlist(lapply(paste0(path, tmp_names_out), fread))

  if(stats == 'OLGA'){
    DT$Pgen <- stats_output$V2
  }else if (stats == 'SONIA'){
    DT$Pgen <- stats_output$Pgen
    DT$Q <- stats_output$Q
    DT$Ppost <- stats_output$Ppost
  }
  DT
}

# Main function -----------------------------------------------------------
#' ALICE_pipeline
#'
#' The function takes a TCRgrapher object as an input and performs neighborhood
#' enrichment analysis using the ALICE algorithm. See ?TCRgrapher.
#'
#' @param TCRgrObject TCRgrapher object that contains clonoset table
#' @param Q_val selection factor. 1/Q sequences pass selection in a thymus. The
#' default value for mouses 6.27. If a human model is taken and Q is not changed
#' manually Q = 27 is used
#' @param cores number of used cores, 1 by default
#' @param thres_counts Only sequences with number of counts above this threshold
#' are taken into account
#' @param N_neighbors_thres Only sequences with number of neighbors above
#' threshold are used to calculate generation probability
#' @param p_adjust_method One of the methods from p.adjust from the stats package
#' possible options: "bonferroni", "holm", "hochberg", "hommel", "BH" or "fdr",
#' "BY", "none". "BH" is a default method.
#' @param chain Statistical model selection. Possible options: "mouseTRB",
#' "humanTRB", "humanTRA".
#' @param stats Tool that will be used for generation probability calculation.
#' Possible options: "OLGA", "SONIA". "SONIA" also calculates Q for every sequence.
#' @param model Standard OLGA generation probability model is used by default.
#' To set your one generation probability model, write "set_custom_model_VDJ
#' <path_to_folder_with_model>". A generation probability model is usually IGOR output.
#' A folder should contain the following files: V_gene_CDR3_anchors.csv,
#' J_gene_CDR3_anchors.csv, model_marginals.txt, model_params.txt. Some models
#' can be found in the folder "model"
#' @return Function returns TCRgrapher object filtered by number
#' of counts and number of neighbors with additional columns. Additional columns
#' are the following
#' \itemize{
#' \item{"ALICE.D"}{"Number of neighbors in clonoset. Neighbor is a similar sequence
#' with one mismatch"}
#' \item{"ALICE.VJ_n_total"}{"Number of unique sequences with given VJ combination"}
#' \item{"ALICE.Pgen"}{"Probability to be generated computed by OLGA"}
#' \item{"ALICE.p_value"}{"p value under null hypothesis that number of sequence's
#' neighbors is not more than in the random model"}
#' \item{"ALICE.p_adjust"}{"p value with multiple testing correction"}
#' \item{"ALICE.log_p_value"}{"log of p_value. p values are usually very small
#' and can be rounded to zero"}
#' }
#' @export
ALICE_pipeline <- function(TCRgrObject, Q_val = 6.27, cores = 1, thres_counts = 1,
                           N_neighbors_thres = 0, p_adjust_method = "BH",
                           chain = 'mouseTRB', stats = 'OLGA', model = '-'){
  if(!("TCRgrapher" %in% attr(TCRgrObject, 'class'))){
    stop("The function takes TCRgrapher object as an input. See ?TCRgrapher",
         call. = FALSE)
  }
  # TODO other requirements

  DT <- clonoset(TCRgrObject)
  message("checking for unproductive sequences if it hasn't been made earlier")
  DT <- DT[!grepl(cdr3aa, pattern = "*", fixed = T) & ((nchar(cdr3nt) %% 3) == 0)]

  model_marginals <- list('mouseTRB'=OLGAVJ_MOUSE_TRB,
                          'humanTRB'=OLGAVJ_HUMAN_TRB,
                          'humanTRA'=OLGAVJ_HUMAN_TRA)

  model_Q_val <- list('mouseTRB'=6.27,
                      'humanTRB'=27,
                      'humanTRA'=27)

  message('model selection')
  if(model == '-'){
    OLGAVJ <- model_marginals[chain][[1]]
    Q_val <- model_Q_val[chain][[1]]
  } else {
    path_to_model <- unlist(base::strsplit(model, ' '))[2]
    params <- read.table(paste0(path_to_model, 'model_params.txt'))
    V_names <- sapply(strsplit(params$V1[grep('TRBV', params$V1)], ";"), `[`, 1)
    V_names <- sub('%TRBV', 'TRBV', V_names)
    J_names <- sapply(strsplit(params$V1[grep('TRBJ', params$V1)], ";"), `[`, 1)
    J_names <- sub('%TRBJ', 'TRBJ', J_names)
    marginals <- read.table(paste0(path_to_model, 'model_marginals.txt'))
    V_prob <- as.numeric(unlist(strsplit(substring(marginals$V1[3], 2, nchar(marginals$V1[3])), ',')))
    J_prob <- as.numeric(unlist(strsplit(substring(marginals$V1[6], 2, nchar(marginals$V1[6])), ',')))
    OLGAVJ <- V_prob %*% t(J_prob)
    rownames(OLGAVJ) <- V_names
    colnames(OLGAVJ) <- J_names
  }

  message('filtering sequences by number of counts')
  DT <- DT[count >= thres_counts,]
  stopifnot(nrow(DT) != 0)
  message('filtering V and J for present in model')
  DT <- DT[bestVGene %in% rownames(OLGAVJ) & bestJGene %in% colnames(OLGAVJ)]
  stopifnot(nrow(DT) != 0)
  message('calculating number of neighbors for every sequence')
  DT <- calculate_nb_of_neighbors(DT)
  DT[, VJ_n_total := .N, .(bestVGene, bestJGene)]
  message('filtering sequences by number of neighbors')
  DT <- DT[D >= N_neighbors_thres][, ind := 1:.N, ]
  stopifnot(nrow(DT) != 0)
  message('generating all possible sequences with one mismatch')
  DT_with_mismatch <- DT[, .(bestVGene, bestJGene,
                             cdr3aa = all_other_variants_one_mismatch_regexp(cdr3aa)
  ), ind]
  message('generation probability calculation')
  DT <- parallel_wrapper_beta(DT = DT, cores = cores, chain = chain,
                              stats = stats, model=model)
  DT_with_mismatch <- parallel_wrapper_beta(DT = DT_with_mismatch, cores = cores,
                                            chain = chain,  stats = stats, model=model)
  if(stats == 'OLGA'){
    # Pgen - probability to be generated computed by OLGA
    # Pgen_sum - sum of Pgen of all sequences similar to the given with one mismatch
    DT$Pgen_sum <- DT_with_mismatch[, sum(Pgen), ind]$V1
    # Pgen_sum_corr - Pgen_sum without probabilities of the main sequence
    # I removed the correction y VJ combination from Pogorelyy's script because
    # olga canculates conditional probability by V and J segments.
    DT[, Pgen_sum_corr := Pgen_sum - Pgen * (nchar(cdr3aa) - 2), ]
    DT[, p_val := ppois(D-1,
                        # normalize for conditioning on observing a variant
                        lambda = (Q_val * VJ_n_total * Pgen_sum_corr) / (1 - exp(-(Q_val * VJ_n_total * Pgen_sum_corr))),
                        lower.tail = F)]
    # p values are toooo small
    DT[, log_p_val := ppois(D-1,
                        # normalize for conditioning on observing a variant
                        lambda = (Q_val * VJ_n_total * Pgen_sum_corr) / (1 - exp(-(Q_val * VJ_n_total * Pgen_sum_corr))),
                        lower.tail = F,
                        log.p = TRUE)]
  } else if (stats == 'SONIA'){
    DT$Ppost_sum <- DT_with_mismatch[, sum(Ppost), ind]$V1
    DT[, Ppost_sum_corr := Ppost_sum - Ppost * (nchar(cdr3aa) - 2), ]
    DT[, p_val := ppois(D-1,
                        # normalize for conditioning on observing a variant
                        lambda =  (VJ_n_total * Ppost_sum_corr) / (1 - exp(-(VJ_n_total * Ppost_sum_corr))),
                        lower.tail = F)]
    # p values are toooo small
    DT[, log_p_val := ppois(D-1,
                        # normalize for conditioning on observing a variant
                        lambda =  (VJ_n_total * Ppost_sum_corr) / (1 - exp(-(VJ_n_total * Ppost_sum_corr))),
                        lower.tail = F,
                        log.p = TRUE)]
  }

  DT[is.na(ALICE.p_value), 'ALICE.p_value'] <- 1
  DT[is.na(ALICE.log_p_value), 'ALICE.log_p_value'] <- 0
  DT[, p_adjust := p.adjust(p_val, method = p_adjust_method)]

  # deletion of unnecessary columns
  DT <- subset(DT, select = -c(ind, Pgen_sum, Pgen_sum_corr))

  setnames(DT, c('D', 'VJ_n_total', 'Pgen', 'p_val', 'p_adjust', 'log_p_val'),
           c('ALICE.D', 'ALICE.VJ_n_total', 'ALICE.Pgen', 'ALICE.p_value',
             'ALICE.p_adjust', 'ALICE.log_p_value'))

  clonoset(TCRgrObject) <- DT
  return(TCRgrObject)
}

# Additional functions ---------------------------------------------------------

#' pval_with_abundance
#'
#' The function calculates the p-value, taking into account the abundance of every
#' clonotype. The method is described in Pogorelyy et al. 2019. In the case of
#' very large datasets with high numbers of neighbors, it may not be possible to
#' calculate the corrected p-value, especially pval_with_abundance_counts.
#'
#' @param clonoset clonoset after the analysis with the ALICE pipeline.
#' To get it, use clonoset(TCRgrapher object)
#' @return Function returns the same clonoset with additional columns:
#' \itemize{
#' \item{pval_with_abundance_log2_counts}{recalculated p-value considering count
#'  number of every clonotype. Log2 is used for count normalization}
#' \item{pval_with_abundance_counts}{recalculated p-value considering count number
#'  of every clonotype. There is no count normalization}
#'  \item{log_pval_with_abundance_log2_counts}{log of pval_with_abundance_log2_counts}
#' \item{log_pval_with_abundance_counts}{log of pval_with_abundance_counts}
#'  }
#' @export
pval_with_abundance <- function(clonoset) {
  # check that ALICE analysis was performed
  if(!(all(c('count', 'ALICE.D', 'ALICE.p_value') %in% colnames(clonoset)))){
    stop("A clonoset must include ALICE.D and ALICE.p_value columns",
         call. = FALSE)
  }
  counts <- clonoset[,count]
  log_counts <- log2(counts)
  clonoset[, log2_counts := log2(count),]
  all_numbers_of_neighbors <- unique(clonoset[,ALICE.D])
  all_numbers_of_neighbors <- all_numbers_of_neighbors[all_numbers_of_neighbors != 0]

  clonoset[ALICE.D == 0, ALICE.pval_with_abundance_log2_counts := 1]
  clonoset[ALICE.D == 0, ALICE.pval_with_abundance_counts := 1]

  for(nb in sort(all_numbers_of_neighbors)){
    # log2
    PDF_f <- approxfun(density(log_counts^nb))
    try(clonoset[ALICE.D == nb,
             ALICE.pval_with_abundance_log2_counts := PDF_f(clonoset[ALICE.D == nb,
                                                                     log2_counts]) * clonoset[ALICE.D == nb,
                                                                                              ALICE.p_value]])
    try(clonoset[ALICE.D == nb,
            ALICE.log_pval_with_abundance_log2_counts := log(PDF_f(clonoset[ALICE.D == nb,
                                                                    log2_counts])) + clonoset[ALICE.D == nb,
                                                                                             ALICE.log_p_value]])
    # just counts
    PDF_f <- approxfun(density(counts^nb))
    try(clonoset[ALICE.D == nb,
             ALICE.pval_with_abundance_counts := PDF_f(clonoset[ALICE.D == nb,
                                                                count]) * clonoset[ALICE.D == nb,
                                                                                   ALICE.p_value]])
    try(clonoset[ALICE.D == nb,
             ALICE.log_pval_with_abundance_counts := log(PDF_f(clonoset[ALICE.D == nb,
                                                                count])) + clonoset[ALICE.D == nb,
                                                                                   ALICE.log_p_value]])
  }
  return(clonoset)
}

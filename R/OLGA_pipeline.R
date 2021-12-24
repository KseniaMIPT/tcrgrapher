#' @importFrom stats p.adjust
#' @importFrom stats ppois
#' @importFrom utils write.table
#' @importFrom stringdist stringdistmatrix
#' @importFrom data.table :=
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
NULL
"OLGAVJ_MOUSE_TRB"
"OLGAVJ_HUMAN_TRB"
"OLGAVJ_HUMAN_TRA"

# Secondary functions -----------------------------------------------------

calculate_nb_of_neighbors_one_side <- function(df) {
  # Calculate nb of neighbors above threshold
  # Every sequence with one mismatch is a neighbor. Only sequences with number
  # of counts above threshold are taken into account
  tmp <- stringdistmatrix(df$cdr3aa, df$cdr3aa, method = "hamming")
  apply(tmp,
    MARGIN = 1,
    function(x) { sum(x <= 1, na.rm = T) })
}

calculate_nb_of_neighbors <- function(df, stats) {
  # TODO: another stats within all sequences not only within VJ combination
  if(stats == 'ALICE'){
    # Filter data by nb of neighbors
    # # Calculate nb of neighbors using computational trick = usage of sequences
    # # with identical left or right parts
    #
    # df[, leftgr := .GRP, .(substr(cdr3aa, 1, floor(nchar(cdr3aa) / 2)),
    #                        bestVGene, bestJGene)]
    # df[, rightgr := .GRP, .(substr(cdr3aa, floor(nchar(cdr3aa) / 2) + 1, nchar(cdr3aa)),
    #                         bestVGene, bestJGene)]

    # df[, D_left := calculate_nb_of_neighbors_one_side(.SD), .(leftgr)]
    # df[, D_right := calculate_nb_of_neighbors_one_side(.SD),.(rightgr)]
    df[ ,ind := 1:.N, ]
    tmp <- stringdistmatrix(df$cdr3aa, df$cdr3aa, method = "hamming")
    df$D <- apply(tmp,
                  MARGIN = 1,
                  function(x) { sum(x <= 1, na.rm = T) }) - 1

    df$D_counts <- apply(tmp,
                         MARGIN = 1,
                         function(x) {sum(df$Read.count[x <= 1])})
    df$D_counts <- df$D_counts - df$Read.count

    df$D_log2_counts <- apply(tmp,
                         MARGIN = 1,
                         function(x) {sum(log2(df$Read.count[x <= 1]))})
    df$D_log2_counts <- df$D_log2_counts - df$Read.count

    # df[, D_id := .N, .(cdr3aa)]
    # df[, D := (D_left + D_right - D_id - 1), ]
    # df <- subset(df, select = -c(rightgr, leftgr, D_left, D_right, D_id))
    df[D < 0, D := 0, ]
    df
  }
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

parallel_wrapper_beta <- function(df, cores = 1, chain = "mouseTRB",
                                  stats = 'OLGA') {
  # Calculate generation probability with OLGA or SONIA

  # add ind column for sequence combining
  if (!("ind" %in% colnames(df))) df[, ind := 1:.N, ]

  fn <- paste0("tmp", 1:cores, ".tsv")
  fn2 <- paste0("tmp_out", 1:cores, ".tsv")
  path <- tempdir()
  path <- paste0(path, '/')

  for (f in c(paste0(path, fn), paste0(path,fn2))) if (file.exists(f)) file.remove(f)

  dfl <- split(df, sort((1:nrow(df) - 1) %% cores + 1))

  for (i in 1:length(dfl)) {
    write.table(as.data.frame(dfl[[i]][, .(cdr3aa, bestVGene, bestJGene, ind), ]),
      quote = F, row.names = F, sep = "\t", file = paste0(path, fn[i],  col.names = F)
    )
  }

  if(stats == 'OLFA'){
    commands <- paste0(
      "olga-compute_pgen --", chain,
      " --display_off --time_updates_off --seq_in 0 --v_in 1 --j_in 2 -d 'tab' -i ",
      path, "tmp", 1:cores, ".tsv -o ", path, "tmp_out", 1:cores, ".tsv"
    )
  } else if (stats == 'SONIA'){
    commands <- paste0('sonia-evaluate --', chain, " --ppost -i ", path, "tmp",
                       1:cores, ".tsv -o ", path, "tmp_out", 1:cores, ".tsv")
  }

  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  foreach(i=1:cores) %dopar% {
    system(commands[i])
  }
  parallel::stopCluster(cl)

  fnt <- c()
  for(file in fn2){
    fnt <- rbind(fnt, fread(paste0(path, file)))
  }

  if(stats == 'OLGA'){
    df$Pgen <- fnt$V2
  }else if (stats == 'SONIA'){
    df$Pgen <- fnt$Pgen
    df$Q <- fnt$Q
    df$Ppost <- fnt$Ppost
  }
  df
}

#' @export
pval_with_abundance <- function(df) {
  counts <- df[,1]
  log_counts <- log2(counts)
  neighbors <- df[,'D']
  df$pval_with_abundance <- -1
  for(d in unique(neighbors)){
    PDF_f <- approxfun(density((log_counts)^d))
    df[df$D == d, 'pval_with_abundance_log2_counts'] <- PDF_f(df[df$D == d, 'D_log2_counts']) * df[df$D == d, 'p_val']
    PDF_f <- approxfun(density((counts)^d))
    df[df$D == d, 'pval_with_abundance_counts'] <- PDF_f(df[df$D == d, 'D_counts'])* df[df$D == d, 'p_val']
  }
  df
}

# Main function -----------------------------------------------------------
#' pipeline_OLGA
#'
#' Main function that takes a table with CDR3 sequences as an input. The table
#' should have the following columns (Order of the columns are not important but
#' the following names are necessary.)
#' \itemize{
#' \item{"Read.count"}{"Number of unique reads per cdr3 sequence"}
#' \item{"cdr3nt"}{"CDR3 nucleotide sequence"}
#' \item{"cdr3aa"}{"CDR3 aminoacid sequence"}
#' \item{"bestVGene"}{"TRBV segment"}
#' \item{"bestJGene"}{"TRBJ segment"}
#' }
#'
#' @param df data.table
#' @param Q_val selection factor. 1/Q sequences pass selection in the thymus. The
#' default value for mouses 6.27. If a human model is taken and Q is not changed
#' manually Q = 27 is used
#' @param cores number of used cores, 1 by default
#' @param thres_counts Only sequences with number of counts above this threshold
#' are taken into account
#' @param N_neighbors_thres Only sequences with number of neighbors above
#' threshold are used to calculate generation probability
#' @param p_adjust_method One of the method from p.adjust from stats package
#' possible options: "bonferroni", "holm", "hochberg", "hommel", "BH" or "fdr",
#' "BY", "none". "BH" is a default method.
#' @param chain Statistical model selection. Possible options: "mouseTRB",
#' "humanTRB", "humanTRA".
#' @param stats OLGA or SONIA
#' @return Function returns the same table that was in input filtered by number
#' of counts and number of neighbors with additional columns. Additional columns
#' are the following
#' \itemize{
#' \item{"D"}{"Number of neighbors in clonoset. Neighbor is a similar sequence
#' with one mismatch"}
#' \item{"VJ_n_total"}{"Number of unique sequences with given VJ combination"}
#' \item{"Pgen"}{"Probability to be generated computed by OLGA"}
#' \item{"Pgen_sum_corr"}{"Sum of Pgen of all sequences similar to the given
#' with one mismatch"}
#' \item{"Pgen_by_VJ"}{"Conditional probability. Sum of probabilities to be
#' generated with given VJ combination"}
#' \item{"p_val"}{"p value under null hypothesis that number of sequence's
#' neighbors is not more than in the random model"}
#' \item{"p_adjust"}{"p value with multiple testing correction"}
#' }
#' @export
pipeline_OLGA <- function(df, Q_val = 6.27, cores = 1, thres_counts = 1,
                          N_neighbors_thres = 1, p_adjust_method = "BH",
                          chain = 'mouseTRB', stats = 'OLGA') {

  # TODO: check table format
  # TODO: check arguments' values
  # I can use stopifnot() to check what should be true

  # checking for unproductive sequences if it hasn't been made earlier
  df <- df[!grepl(cdr3aa, pattern = "*", fixed = T) & ((nchar(cdr3nt) %% 3) == 0)]

  if (chain == 'mouseTRB'){
    OLGAVJ = OLGAVJ_MOUSE_TRB
  } else if (chain == 'humanTRB'){
    OLGAVJ = OLGAVJ_HUMAN_TRB
    if(Q_val == 6.27){
      Q_val = 27
    }
  } else if (chain == 'humanTRA'){
    OLGAVJ = OLGAVJ_HUMAN_TRA
    if(Q_val == 6.27){
      Q_val = 27
    }
  } else {
    stop('There is no such model')
  }

  # filter sequences by number of counts
  df <- df[Read.count > thres_counts,]
  stopifnot(nrow(df) != 0)
  # filter V and J for present in model
  df <- df[bestVGene %in% row.names(OLGAVJ) & bestJGene %in% colnames(OLGAVJ)]
  stopifnot(nrow(df) != 0)

  df <- calculate_nb_of_neighbors(df, stats = stats)

  df[, VJ_n_total := .N, .(bestVGene, bestJGene)]
  df <- df[D >= N_neighbors_thres][, ind := 1:.N, ]
  stopifnot(nrow(df) != 0)

  df_with_mismatch <- df[, .(bestVGene, bestJGene,
    cdr3aa = all_other_variants_one_mismatch_regexp(cdr3aa)
  ), ind]

  df <- parallel_wrapper_beta(df = df, cores = cores, chain = chain)
  df_with_mismatch <- parallel_wrapper_beta(df = df_with_mismatch,
                                                 cores = cores, chain = chain)
  if(stats == 'OLGA'){
    # Pgen - probability to be generated computed by OLGA
    # Pgen_sum - sum of Pgen of all sequences similar to the given with one mismatch
    df$Pgen_sum <- df_with_mismatch[, sum(Pgen), ind]$V1
    # Pgen_sum_corr - Pgen_sum without probabilities of the main sequence
    df[, Pgen_sum_corr := Pgen_sum - Pgen * (nchar(cdr3aa) - 2), ]
    # Bayes' rule
    df[, Pgen_by_VJ := 1 * Pgen_sum_corr / OLGAVJ[cbind(bestVGene, bestJGene)], ]
    df[, p_val := ppois(D-1, lambda = Q_val * VJ_n_total * Pgen_by_VJ, lower.tail = F)]
  } else if (stats == 'SONIA'){
    df$Ppost_sum <- df_with_mismatch[, sum(Ppost), ind]$V1
    df[, Ppost_sum_corr := Ppost_sum - Ppost * (nchar(cdr3aa) - 2), ]
    f[, Ppost_by_VJ := 1 * Ppost_sum_corr / OLGAVJ[cbind(bestVGene, bestJGene)], ]
    df[, p_val := ppois(D-1, lambda =  VJ_n_total * Ppost_by_VJ, lower.tail = F)]
  }

  df[, p_adjust := p.adjust(p_val, method = p_adjust_method)]
  # add cluster IDs
  # df <- find_cluster(df)

  # deletion of unnecessary columns
  df <- subset(df, select = -c(ind))
  return(df)
}

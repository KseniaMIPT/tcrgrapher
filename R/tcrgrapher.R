#' @importFrom stats p.adjust
#' @importFrom stats ppois
#' @importFrom stats approxfun
#' @importFrom utils write.table
#' @importFrom stringdist stringdistmatrix
#' @importFrom data.table :=
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph components
NULL
"OLGAVJ_MOUSE_TRB"
"OLGAVJ_HUMAN_TRB"
"OLGAVJ_HUMAN_TRA"
"C57BL6_MOUSE_TRB"

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

calculate_nb_of_neighbors <- function(df) {
  # give unique number 'leftgr' to the rows with the same VJ combination and left
  # CDR3aa part
  df[, leftgr := .GRP, .(substr(cdr3aa, 1, floor(nchar(cdr3aa) / 2)),
                         bestVGene, bestJGene)]
  # give unique number 'rightgr' to the rows with the same VJ combination and right
  # CDR3aa part
  df[, rightgr := .GRP, .(substr(cdr3aa, floor(nchar(cdr3aa) / 2) + 1, nchar(cdr3aa)),
                          bestVGene, bestJGene)]

  df[, D_left := calculate_nb_of_neighbors_one_side(.SD), .(leftgr)]
  df[, D_right := calculate_nb_of_neighbors_one_side(.SD), .(rightgr)]
  df[, D_id := .N, .(cdr3aa)]
  df[, D := (D_left + D_right - D_id - 1), ]
  df <- subset(df, select = -c(rightgr, leftgr, D_left, D_right, D_id))
  df[D < 0, D := 0, ]
  df
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
                                  stats = 'OLGA', model='-') {
  # Calculate generation probability with OLGA or SONIA

  # add ind column for sequence combining
  if (!("ind" %in% colnames(df))) df[, ind := 1:.N, ]

  tmp_names <- paste0("tmp", 1:cores, ".tsv")
  tmp_names_out <- paste0("tmp_out", 1:cores, ".tsv")
  path <- tempdir()
  path <- paste0(path, '/')

  for (f in c(paste0(path, tmp_names), paste0(path,tmp_names_out))) if (file.exists(f)) file.remove(f)

  dft <- split(df, sort((1:nrow(df) - 1) %% cores + 1))

  for (i in 1:length(dft)) {
    write.table(as.data.frame(dft[[i]][, .(cdr3aa, bestVGene, bestJGene, ind), ]),
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
    commands <- paste0('sonia-evaluate --', chain, " --ppost -i ", path, "tmp",
                       1:cores, ".tsv -o ", path, "tmp_out", 1:cores, ".tsv")
  }

  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  foreach(i=1:cores) %dopar% {
    system(commands[i], intern=TRUE)
  }
  parallel::stopCluster(cl)

  stats_output <- c()
  for(file in tmp_names_out){
    stats_output <- rbind(stats_output, fread(paste0(path, file)))
  }

  if(stats == 'OLGA'){
    df$Pgen <- stats_output$V2
  }else if (stats == 'SONIA'){
    df$Pgen <- stats_output$Pgen
    df$Q <- stats_output$Q
    df$Ppost <- stats_output$Ppost
  }
  df
}

# Main function -----------------------------------------------------------
#' tcrgrapher
#'
#' Main function that takes a table with CDR3 sequences as an input. The table
#' should have the following columns (Order of the columns are not important but
#' the following names are necessary)
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
#' @param stats Tool that will be used for generation probability calculation.
#' Possible options: "OLGA", "SONIA". "SONIA" also calculate Q for every sequence.
#' @param model Standard OLGA generation probability model is used by default.
#' To set your one generation probability model write "set_custom_model_VDJ
#' <path_to_folder_with_model>". Generation probability model is usually IGOR output.
#' Folder should contain the following files: V_gene_CDR3_anchors.csv,
#' J_gene_CDR3_anchors.csv, model_marginals.txt, model_params.txt. Some models
#' one can find in the folder "model"
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
tcrgrapher <- function(df, Q_val = 6.27, cores = 1, thres_counts = 1,
                          N_neighbors_thres = 1, p_adjust_method = "BH",
                          chain = 'mouseTRB', stats = 'OLGA', model = '-') {

  message("checking for unproductive sequences if it hasn't been made earlier")
  df <- df[!grepl(cdr3aa, pattern = "*", fixed = T) & ((nchar(cdr3nt) %% 3) == 0)]

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
  df <- df[Read.count >= thres_counts,]
  stopifnot(nrow(df) != 0)
  message('filtering V and J for present in model')
  df <- df[bestVGene %in% rownames(OLGAVJ) & bestJGene %in% colnames(OLGAVJ)]
  stopifnot(nrow(df) != 0)
  message('calculating number of neighbors for every sequence')
  df <- calculate_nb_of_neighbors(df)
  df[, VJ_n_total := .N, .(bestVGene, bestJGene)]
  message('filtering sequences by number of neighbors')
  df <- df[D >= N_neighbors_thres][, ind := 1:.N, ]
  stopifnot(nrow(df) != 0)
  message('generating all possible sequences with one mismatch')
  df_with_mismatch <- df[, .(bestVGene, bestJGene,
    cdr3aa = all_other_variants_one_mismatch_regexp(cdr3aa)
  ), ind]
  message('generation probability calculation')
  df <- parallel_wrapper_beta(df = df, cores = cores, chain = chain,
                              stats = stats, model=model)
  df_with_mismatch <- parallel_wrapper_beta(df = df_with_mismatch, cores = cores,
                                            chain = chain,  stats = stats, model=model)
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
    df[, Ppost_by_VJ := 1 * Ppost_sum_corr / OLGAVJ[cbind(bestVGene, bestJGene)], ]
    df[, p_val := ppois(D-1, lambda =  VJ_n_total * Ppost_by_VJ, lower.tail = F)]
  }

  df[, p_adjust := p.adjust(p_val, method = p_adjust_method)]

  # deletion of unnecessary columns
  df <- subset(df, select = -c(ind))
  return(df)
}

# Additional functions ---------------------------------------------------------

#' pval_with_abundance
#'
#' Function calculates p-value taking into account abundance of every clonotype
#'
#' @param df output of tcrgrapher function
#' @return Function returns the same data.table with additional columns:
#' \itemize{
#' \item{"pval_with_abundance_log2_counts"}{"Recalculated p-value considering
#' count number of every clonotype. Log2 is used for count normalization"}
#' #' \item{"pval_with_abundance_counts"}{"Recalculated p-value considering
#' count number of every clonotype. There is no count normalization"}
#' @export
pval_with_abundance <- function(df) {
  counts <- df[,1]
  log_counts <- log2(counts)
  neighbors <- df[,'D']
  df$pval_with_abundance <- -1
  for(d in unique(neighbors)){
    PDF_f <- approxfun(density((log_counts)^d))
    df[df$D == d,
       'pval_with_abundance_log2_counts'] <- PDF_f(df[df$D == d,
                                                      'D_log2_counts']) * df[df$D == d, 'p_val']
    PDF_f <- approxfun(density((counts)^d))
    df[df$D == d,
       'pval_with_abundance_counts'] <- PDF_f(df[df$D == d,
                                                 'D_counts'])* df[df$D == d, 'p_val']
  }
  df
}

#' make_graph
#'
#' Function makes graph from tcrgrapher output with igraph package. Every node
#' of the graph is an unique clonotype from the table (one line). Edges
#' connects clonotypes with one amino acid mismatch or identical clonotypes if
#' they were in separate lines.
#'
#' @param df output of tcrgrapher function
#' @return Function returns an igraph graph object
#' @export
make_graph <- function(df){
  adj_matrix <- stringdistmatrix(df$cdr3aa, df$cdr3aa, method = "hamming")
  adj_matrix <- 1*(adj_matrix <= 1)
  rownames(adj_matrix) <- df$cdr3aa
  colnames(adj_matrix) <- df$cdr3aa
  diag(adj_matrix) <- 0
  g <- graph_from_adjacency_matrix(
    adj_matrix,
    mode = "undirected",
    weighted = NULL,
    diag = TRUE,
    add.colnames = NULL,
    add.rownames = NA
  )
  g
}

#' find_cluster
#'
#' Function takes tcrgrapher output and returns the same table with additional
#' column "cluster_id". All clusters of neighbours with one mismatch have unique
#' id. Function uses "components" function from igraph package.
#'
#' @param df output of tcrgrapher function
#' @param g make_graph output
#' @return Function returns the same data.table with additional column
#' "cluster_id"
#' @export
find_cluster <- function(df, g){
  components <- components(g)
  df$cluster_id <- components$membership
  df
}

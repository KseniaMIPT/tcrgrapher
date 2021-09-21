#' @importFrom stats p.adjust
#' @importFrom stats ppois
#' @importFrom utils write.table
#' @importFrom stringdist stringdistmatrix
#' @importFrom data.table :=

# Secondary functions -----------------------------------------------------

calculate_nb_of_neighbors <- function(df, N_neighbors_thres = 1) {
  # Calculate nb of neighbors above threshold
  # Every count is a one neighbor
  if (nrow(df) > 1) {
    tmp <- stringdistmatrix(df$cdr3aa, df$cdr3aa, method = "hamming")
    apply(tmp,
      MARGIN = 1,
      function(x) {
        sum(x[df$Read.count > N_neighbors_thres] <= 1, na.rm = T)
      }
    ) - 0
  } else {
    1
  }
}

filter_data_by_nb_of_neighbors <- function(df, N_neighbors_thres = 1) {
  # Filter data by nb of neighbors
  # Calculate nb of neighbors using computational trick = usage of sequences
  # with identical left or right parts

  df[, leftgr := .GRP, .(substr(cdr3aa, 1, floor(nchar(cdr3aa) / 2)),
    bestVGene, bestJGene)]
  df[, rightgr := .GRP, .(substr(cdr3aa, floor(nchar(cdr3aa) / 2) + 1, nchar(cdr3aa)),
    bestVGene, bestJGene)]

  # TODO: тут так-то отдельный трешхолд, надо разобраться с ними

  df[
    , D_left := calculate_nb_of_neighbors(.SD, N_neighbors_thres = N_neighbors_thres),
    .(leftgr)
  ]
  df[
    , D_right := calculate_nb_of_neighbors(.SD, N_neighbors_thres = N_neighbors_thres),
    .(rightgr)
  ]
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

olga_parallel_wrapper_beta <- function(df, cores = 1, chain = "mouseTRB",
                                       withoutVJ = F) {
  # Calculate generation probability with OLGA

  # add ind column for sequence combining
  if (!("ind" %in% colnames(df))) df[, ind := 1:.N, ]

  fn <- paste0("tmp", 1:cores, ".tsv")
  fn2 <- paste0("tmp_out", 1:cores, ".tsv")
  path <- tempdir()

  for (f in c(paste0(path, fn), paste0(path,fn2))) if (file.exists(f)) file.remove(f)

  dfl <- split(df, sort((1:nrow(df) - 1) %% cores + 1))

  for (i in 1:length(dfl)) {
    write.table(as.data.frame(dfl[[i]][, .(cdr3aa, bestVGene, bestJGene, ind), ]),
      quote = F, row.names = F, sep = "\t", file = paste0(path, fn[i])
    )
  }

  olga_commands <- paste0(
    "olga-compute_pgen --", chain,
    " --display_off --time_updates_off	--seq_in 0 --v_in 1 --j_in 2 --lines_to_skip 1 -d 'tab' -i ",
    path, "tmp", 1:cores, ".tsv -o ", path, "tmp_out", 1:cores, ".tsv"
  )
  if (withoutVJ) {
    olga_commands <- paste0(
      "olga-compute_pgen --", chain,
      " --display_off --time_updates_off --seq_in 0 --lines_to_skip 1 -d 'tab' -i ",
      path, "tmp", 1:cores, ".tsv -o ", path, "tmp_out", 1:cores, ".tsv"
    )
  }

  system(olga_commands, wait = T)
  #system(paste0(olga_commands, collapse = " & "), wait = T)
  system("echo done", wait = T)

  # TODO вот тут беда была
  fnt <- fread(paste0(path, fn2)) #do.call(rbind, lapply(paste0(path, fn2), fread))
  df$Pgen <- fnt$V2
  df
}

# TODO дописать return

# Main function -----------------------------------------------------------
#' pipeline_OLGA
#'
#' Main fucntion that takes table with cdr3 sequences as an input. Table should
#' have the following columns (names of the colums are not important but the
#' followig order is necessary)
#' \itemize{
#' \item{"Read.count"}{"Number of unique reads per cdr3 sequence"}
#' \item{"freq"}{"Clonotype frequency in the clonoset"}
#' \item{"cdr3nt"}{"CDR3 nucleotide sequence"}
#' \item{"bestVGene"}{"TRBV segment"}
#' \item{"bestVGene"}{"TRBD segment"}
#' \item{"bestJGene"}{"TRBJ segment"}
#' \item{"VEnd"}{"Position of the end of V segment in CDR3 sequence"}
#' \item{"DStart""}{"Position of the start of D segment in CDR3 sequence"}
#' \item{"DEnd"}{"Position of the end of D segment in CDR3 sequence"}
#' \item{"JStart"}{"Position of the start of J segment in CDR3 sequence"}
#' }
#'
#' @param df data.table
#' @param Q selection factor. 1/Q sequences pass selection in thymus. Default
#' value for mouses 6.27
#' @param cores number of used cores, 1 by default
#' @param Read_thres threshold 1
#' @param Read_thres2 threshold 2
#' @param N_neighbors_thres threshold 3
#' @param p_adjust_method one of the method from p.adjust from stats package
#' possible options: "bonferroni", "holm", "hochberg", "hommel", "BH" or "fdr",
#' "BY", "none". "BH" is a default method.
#' @return Function returns tha same table that was in input with additional
#' columns
#' \itemize{
#' \item{"D"}{"Number of neighbors in clonoset. Neighbor is a simular sequence
#' with one mismatch"}
#' }
#' @export
pipeline_OLGA <- function(df, Q = 6.27, cores = 1, Read_thres = 0,
                          Read_thres2 = 1, N_neighbors_thres = 1,
                          p_adjust_method = "BH") {
  colnames(df) <- c(
    "Read.count", "freq", "cdr3nt", "cdr3aa", "bestVGene", "bestDGene",
    "bestJGene", "VEnd", "DStart", "DEnd", "JStart"
  )

  # TODO: check table format
  # TODO: check arguments' values
  # I can use stopifnot() to check what should be true

  # checking for unproductive sequences if it hasn't been made earlier
  df <- df[!grepl(cdr3aa, pattern = "*", fixed = T) & ((nchar(cdr3nt) %% 3) == 0)]
  # filter V and J for present in model
  df <- df[Read.count > Read_thres,][bestVGene %in% row.names(OLGAVJ) &
    bestJGene %in% colnames(OLGAVJ)]

  df <- filter_data_by_nb_of_neighbors(df, N_neighbors_thres = 1)

  df[Read.count > Read_thres2, n_total := .N, .(bestVGene, bestJGene)]
  df <- df[D >= N_neighbors_thres][, ind := 1:.N, ]

  df_with_mismatch <- df[, .(bestVGene, bestJGene,
    cdr3aa = all_other_variants_one_mismatch_regexp(cdr3aa)
  ), ind]

  df <- olga_parallel_wrapper_beta(df = df, cores = cores)
  df_with_mismatch <- olga_parallel_wrapper_beta(df = df_with_mismatch,
    cores = cores)

  df$Pgen1 <- df_with_mismatch[, sum(Pgen), ind]$V1
  df[, Pgen3 := Pgen1 - Pgen * (nchar(cdr3aa) - 2), ]
  df[, space := Pgen3 / (OLGAVJ[cbind(bestVGene, bestJGene)]), ]
  df[, space_n := Pgen3 / (OLGAVJ[cbind(bestVGene, bestJGene)]), ]
  df[, p_val := ppois(D, lambda = 3 * Q * n_total * Pgen3 /
    (OLGAVJ[cbind(bestVGene, bestJGene)]), lower.tail = F), ]

  df[, p_adjust := p.adjust(p_val, method = p_adjust_method)]
  return(df)
}

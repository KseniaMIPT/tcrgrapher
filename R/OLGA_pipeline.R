#' @importFrom stats p.adjust
#' @importFrom stats ppois
#' @importFrom utils write.table
#' @importFrom stringdist stringdistmatrix
#' @importFrom data.table :=

"VDJT"
"OLGA_V_J_mouse_beta.rda"

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
    ) - 0 # я не поняла зачем это действие нужно, но оставила на всякий случай
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
                                       withoutVJ = F, prompt = T) {
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

  # Pause before read
  Sys.sleep(60)

  if (prompt) {
    readline(prompt = "Press [enter] to continue")
  }
  # TODO вот тут беда была
  fnt <- fread(paste0(path, fn2)) #do.call(rbind, lapply(paste0(path, fn2), fread))
  df$Pgen <- fnt$V2
  df
}

# Main function -----------------------------------------------------------
#' Main function
#'
#' @param df data.table
#' @param Q selection factor
#' @param cores number of cores
#' @param prompt smth
#' @param Read_thres threshold 1
#' @param Read_thres2 threshold 2
#' @param N_neighbors_thres threshold 3
#' @export
pipeline_OLGA <- function(df, Q = 6.27, cores = 1, prompt = F, Read_thres = 0,
                          Read_thres2 = 1, N_neighbors_thres = 1) {
  colnames(df) <- c(
    "Read.count", "freq", "cdr3nt", "cdr3aa", "bestVGene", "d",
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

  df <- olga_parallel_wrapper_beta(df = df, cores = cores, prompt = prompt)
  df_with_mismatch <- olga_parallel_wrapper_beta(
    df = df_with_mismatch,
    cores = cores, prompt = prompt
  )

  df$Pgen1 <- df_with_mismatch[, sum(Pgen), ind]$V1
  df[, Pgen3 := Pgen1 - Pgen * (nchar(cdr3aa) - 2), ]
  df[, space := Pgen3 / (OLGAVJ[cbind(bestVGene, bestJGene)]), ]
  df[, space_n := Pgen3 / (OLGAVJ[cbind(bestVGene, bestJGene)]), ]
  df[, p_val := ppois(D, lambda = 3 * Q * n_total * Pgen3 /
    (OLGAVJ[cbind(bestVGene, bestJGene)]), lower.tail = F), ]

  df[, p_adjust := p.adjust(p_val, method = "BH")]
  return(df)
}

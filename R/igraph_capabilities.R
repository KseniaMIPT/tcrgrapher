#' make_TCR_graph
#'
#' Function makes graph from tcrgrapher output with igraph package. Every node
#' of the graph is an unique clonotype from the table (one line). Edges
#' connects clonotypes with one amino acid mismatch or identical clonotypes if
#' they were in separate lines.
#'
#' @param df output of tcrgrapher function
#' @return Function returns an igraph graph object
#' @export
make_TCR_graph <- function(df, v_gene = TRUE){
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \"igraph\" must be installed to use this function.",
      call. = FALSE
    )
  }
  adj_matrix <- stringdistmatrix(df$cdr3aa, df$cdr3aa, method = "hamming")
  adj_matrix <- 1*(adj_matrix <= 1)
  rownames(adj_matrix) <- df$cdr3aa
  colnames(adj_matrix) <- df$cdr3aa
  diag(adj_matrix) <- 0
  if(v_gene){
    adj_matrix[outer(df$v_segment, df$v_segment, FUN = '!=')] <- 0
  }
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
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \"igraph\" must be installed to use this function.",
      call. = FALSE
    )
  }
  components <- components(g)
  df$cluster_id <- components$membership
  df
}

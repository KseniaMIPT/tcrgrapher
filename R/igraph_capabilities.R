#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
NULL

#' make_TCR_graph
#'
#' Function makes a graph from the clonoset table with the igraph package. Every node
#' of the graph is an unique clonotype from the table (one line). Edges
#' connects clonotypes with one amino acid mismatch or identical clonotypes if
#' they are in separate lines.
#'
#' @param clonoset clonoset clonoset table. To get from TCRgrapher object use clonoset(<TCRgrapher object>)
#' @return Function returns an igraph graph object
#' @export
make_TCR_graph <- function(clonoset, v_gene = TRUE, j_gene = FALSE){
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \"igraph\" must be installed to use this function.",
      call. = FALSE
    )
  }
  # graph from adjacency matrix
  # adj_matrix <- stringdistmatrix(clonoset$cdr3aa, clonoset$cdr3aa, method = "hamming")
  # adj_matrix <- 1*(adj_matrix <= 1)
  # rownames(adj_matrix) <- clonoset$clone_id
  # colnames(adj_matrix) <- clonoset$clone_id
  # diag(adj_matrix) <- 0
  # if(v_gene){
  #   adj_matrix[outer(clonoset$bestVGene, clonoset$bestVGene, FUN = '!=')] <- 0
  # }
  # if(j_gene){
  #   adj_matrix[outer(clonoset$bestJGene, clonoset$bestJGene, FUN = '!=')] <- 0
  # }
  # g <- graph_from_adjacency_matrix(
  #   adj_matrix,
  #   mode = "undirected"
  # )
  g <- make_empty_graph(directed = FALSE)
  n <- nrow(clonoset)
  g <- add_vertices(g, n)
  to_check <- rep(TRUE, n)
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  for(v in 1:n){
    if(to_check[v]){
      to_check[v] <- FALSE
      if(v_gene){
        V_current <- clonoset$bestVGene[v]
        candidates <- clonoset[to_check & bestVGene == V_current, .(cdr3aa, clone_id)]
        if(j_gene){
          J_current <- clonoset$bestJGene[v]
          candidates <- candidates[bestJGene == J_current]
        }
      } else {
        candidates <- clonoset[to_check, .(cdr3aa, clone_id)]
      }
      if(nrow(candidates) != 0){
        neighbors <- stringdistmatrix(clonoset$cdr3aa[v], candidates$cdr3aa, method = "hamming") <= 1
        for(x in candidates[neighbors[1,], clone_id]){
          g <- add_edges(g, c(v, x))
        }
      }
    }
    setTxtProgressBar(pb, v/n*100)
  }
  g
}

#' find_TCR_components
#'
#' Function takes clonoset table and returns the same table with additional
#' column "cluster_id". All clusters of neighbors with one mismatch have unique
#' id. Function uses "components" function from igraph package.
#'
#' @param clonoset clonoset table. To get from TCRgrapher object use clonoset(<TCRgrapher object>)
#' @param g make_TCR_graph output
#' @return Function returns the same clonoset table with additional column "cluster_id"
#' @export
find_TCR_components <- function(clonoset, g){
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \"igraph\" must be installed to use this function.",
      call. = FALSE
    )
  }
  components <- components(g)
  clonoset$cluster_id <- components$membership
  clonoset
}

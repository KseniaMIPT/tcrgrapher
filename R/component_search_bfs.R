#' @importFrom stringdist stringdistmatrix
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
NULL
#' bfs_for_TCRs
#'
#' The function searches for all neighbors of the source (src) clonotype. A neighbor
#' is a clonotype that differs by one amino acid mismatch.
#'
#' @param clonoset To get clonoset use clonoset(your_TCR_grapher_object)
#' @param src clone_id of the source clonotype
#' @param comp_id cluster_id for the visited vertexes. All other vertexes will have
#' -1 cluster_id if there are none. If some cluster_ids have already been set,
#' they will be kept.
#' @return The function returns the clonoset where all neighrors of the source
#' clonotype have the given id.
#' @export
bfs_for_TCRs <- function(clonoset, src, comp_id, v_gene = TRUE, j_gene = FALSE){
  if(!('cluster_id' %in% colnames(clonoset))){
    clonoset[, cluster_id := -1]
  }
  queue <- c(src)
  clonoset[src, cluster_id := comp_id]
  V_genes <- if (v_gene) clonoset$bestVGene[src] else unique(clonoset$bestVGene)
  J_genes <- if (j_gene) clonoset$bestJGene[src] else unique(clonoset$bestJGene)
  while(length(queue) != 0){
    cur <- queue[1]
    if(length(queue) > 1){
      queue <- queue[2:length(queue)]
    } else {
      queue <- c()
    }
    candidates <- clonoset[clonoset$cluster_id == -1 &
                             bestVGene %in% V_genes &
                             bestJGene %in% J_genes,
                           .(cdr3aa, clone_id)]
    if(nrow(candidates) != 0){
      neighbors <- stringdistmatrix(clonoset[,cdr3aa][cur],
                                    candidates[,cdr3aa], method = "hamming") <= 1
      for(n in candidates[neighbors[1,], clone_id]){
        clonoset[n, cluster_id := comp_id]
        clonoset[cur, cluster_id := comp_id]
        queue <- c(queue, n)
      }
    }
  }
  clonoset
}

#' find_TCR_components_by_bfs
#'
#' The function searches for connectivity components of the graph, where every node is a
#' clonotype with a unique clone_id and edges connect clonotypes that differ by one
#' amino acid mismatch.
#'
#' @param clonoset To get clonoset use clonoset(your_TCR_grapher_object)
#' @param v_gene Boolean value. If 'v_gene' is 'TRUE', only clonotypes with the
#' same V genes will be part of one connectivity component. Default value is 'TRUE'.
#' @param j_gene Boolean value. If 'j_gene' is 'TRUE', only clonotypes with the
#' same J genes will be part of one connectivity component. Default value is 'FALSE'.
#' @return The function returns the clonoset with the new column "cluster_id", where
#' every connectivity component has a unique value.
#' @export
find_TCR_components_by_bfs <- function(clonoset, v_gene = TRUE, j_gene = FALSE){
  src <- 1
  comp_id <- 1
  clonoset[, cluster_id := -1]
  pb <- txtProgressBar(min = 0, max = log(nrow(clonoset)), style = 3)
  while(sum(clonoset[,cluster_id] == -1) != 0){
    clonoset <- bfs_for_TCRs(clonoset, src, comp_id, v_gene, j_gene)
    src <- clonoset[cluster_id == -1, clone_id][1]
    comp_id <- comp_id + 1
    setTxtProgressBar(pb, log(src))
  }
  clonoset
}

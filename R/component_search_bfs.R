#' @importFrom stringdist stringdistmatrix
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
NULL

#' @param clonoset clonoset
#' @param src clone_id of the source clonotype
#'
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

find_TCR_components_by_bfs <- function(TCRgrObject, v_gene = TRUE, j_gene = FALSE){
  src <- 1
  comp_id <- 1
  clonoset <- clonoset(TCRgrObject)
  clonoset[, cluster_id := -1]
  pb <- txtProgressBar(min = 0, max = nrow(clonoset), style = 3)
  while(sum(clonoset[,cluster_id] == -1) != 0){
    clonoset <- bfs_for_TCRs(clonoset, src, comp_id, v_gene, j_gene)
    src <- clonoset[cluster_id == -1, clone_id][1]
    comp_id <- comp_id + 1
    setTxtProgressBar(pb, src)
  }
  clonoset(TCRgrObject) <- clonoset
  TCRgrObject
}

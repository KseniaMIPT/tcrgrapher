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
#' @param edges matrix to which new edges will be added
#' @return The function returns a list with a clonoset where all neighbors of the source
#' clonotype have the given id. The second item of the list is the matrix with edges.
#' @export
bfs_for_TCRs <- function(clonoset, src, comp_id, v_gene = TRUE, j_gene = FALSE, edges = c()){
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
        edges <- rbind(edges, c(cur, n))
      }
    }
  }
  return(list(clonoset, edges))
}

#' find_TCR_components_by_bfs
#'
#' The function searches for connectivity components of the graph, where every node is a
#' clonotype with a unique clone_id and edges connect clonotypes that differ by one
#' amino acid mismatch. It returns the same TCRgrapher object with list of edges and
#' with additional column cluster_id in the clonoset table.
#'
#' @param TCRgrObject input object
#' @param v_gene Boolean value. If 'v_gene' is 'TRUE', only clonotypes with the
#' same V genes will be part of one connectivity component. Default value is 'TRUE'.
#' @param j_gene Boolean value. If 'j_gene' is 'TRUE', only clonotypes with the
#' same J genes will be part of one connectivity component. Default value is 'FALSE'.
#' @return The function returns TCRgrapher object. The clonoset has the new column "cluster_id", where
#' every connectivity component has a unique value. TCRgrapher  object contains
#' found edges. To see them use edges(<your object>)
#' @export
find_TCR_components_by_bfs <- function(TCRgrObject, v_gene = TRUE, j_gene = FALSE){
  clonoset <- TCRgrObject@clonoset
  src <- 1
  comp_id <- 1
  clonoset[, cluster_id := -1]
  edges <- c()
  pb <- txtProgressBar(min = 0, max = log(nrow(clonoset)), style = 3)
  while(sum(clonoset[,cluster_id] == -1) != 0){
    output_bfs <- bfs_for_TCRs(clonoset, src, comp_id, v_gene, j_gene, edges)
    clonoset <- output_bfs[[1]]
    edges <- output_bfs[[2]]
    src <- clonoset[cluster_id == -1, clone_id][1]
    comp_id <- comp_id + 1
    setTxtProgressBar(pb, log(src))
  }

  TCRgrObject@clonoset <- clonoset
  TCRgrObject@edges <- as.matrix(edges)
  TCRgrObject
  #return(list(clonoset, edges))
}

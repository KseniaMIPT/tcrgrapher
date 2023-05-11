#' @import methods
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom data.table :=

# S4 class
setClass('TCRgrapherCounts',
         contains = 'TCRgrapher',
         slots = c(
           count_table = 'data.frame',
           feature_info = 'data.table'
         ))

#' TCRgrapherCounts object
#'
#' Function takes TCRgrapher object as an input and creates a count table
#' where each row corresponds to unique feature (unique amino acid sequence or
#' cluster of sequences) and each column corresponds to the sample.
#'
#' @param
#' @return Function returns TCRgrapherCounts object
#' @export
TCRgrapherCounts <- function(TCRgrObject, v_gene = TRUE, j_gene = FALSE){
  # requirements
  if(!("TCRgrapher" %in% attr(TCRgrObject, 'class'))){
    stop("The function takes TCRgrapher object as an input. See ?TCRgrapher",
         call. = FALSE)
  }
  if(!is.logical(v_gene)){
    stop("'v_gene' must be TRUE or FALSE",
         call. = FALSE)
  }
  if(!is.logical(j_gene)){
    stop("'v_gene' must be TRUE or FALSE",
         call. = FALSE)
  }
  if(!is.logical(by_clusters)){
    stop("'by_clusters' must be TRUE or FALSE",
         call. = FALSE)
  }

  metadata <- TCRgrObject@metadata
  clonoset <- TCRgrObject@clonoset
  grouping <- c('cdr3aa', 'bestVGene', 'bestJGene')[c(TRUE, v_gene, j_gene)]

  formula <- as.formula(paste0(paste(grouping, collapse = ' + '), ' ~ sample_id'))
  count_table <- dcast(clonoset, formula, value.var = 'count', fun.aggregate = sum)
  count_table <- as.data.frame(count_table)

  rownames(count_table) <- Reduce(paste, count_table[,grouping])
  count_table <- count_table[,-(1:length(grouping))]

  feature_info <- cbind(clonoset, Reduce(paste, clonoset[,grouping]))

  new('TCRgrapherCounts', clonoset = TCRgrObject@clonoset,
      metadata = TCRgrObject@metadata, count_table = count_table,
      feature_info = feature_info)
}

#' @import methods
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom data.table :=
NULL

# -----------------------------------------------------------------------------

#' TCRgrapherCounts object
#'
#' Function takes a TCRgrapher object as an input and creates a count table
#' where each row corresponds to an unique amino acid sequence and each column
#' corresponds to the sample.
#'
#' @param TCRgrObject TCRgrapher object. It can be constructed by calling the
#' TCRgrapher( ) function.
#' @param v_gene Boolean value. If 'v_gene' is 'TRUE', clonotypes with the same amino
#' acid sequences but different V genes are presented in different rows.
#' Default value is 'TRUE'.
#' @param j_gene Boolean value. If 'j_gene' is 'TRUE', clonotypes with the same amino
#' acid sequences but different J genes are presented in different rows.
#' Default value is 'FALSE'.
#' @return Function returns TCRgrapherCounts object that contains clonoset and
#' metadata from TCRgrapher object, count table and feature info table
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
  metadata <- TCRgrObject@metadata
  clonoset <- TCRgrObject@clonoset
  grouping <- c('cdr3aa', 'bestVGene', 'bestJGene')[c(TRUE, v_gene, j_gene)]

  formula <- as.formula(paste0(paste(grouping, collapse = ' + '), ' ~ sample_id'))
  count_table <- dcast(clonoset, formula, value.var = 'count', fun.aggregate = sum)
  count_table <- as.data.frame(count_table)

  rownames(count_table) <- Reduce(if(length(grouping) == 1) c else paste,
                                  count_table[,grouping])

  count_table <- count_table[,-(1:length(grouping))]

  feature_info <- cbind(clonoset, Reduce(if(length(grouping) == 1) c else paste,
                                         clonoset[,grouping, with=FALSE]))
  setnames(feature_info, 'V2', 'feature')
  feature_info$grouping <- paste(grouping, collapse = ' ')

  new('TCRgrapherCounts', clonoset = TCRgrObject@clonoset,
      metadata = TCRgrObject@metadata, count_table = count_table,
      feature_info = feature_info)
}

#' @export
setGeneric("count_table", function(x) standardGeneric("count_table"))

#' @export
setMethod("count_table", "TCRgrapher", function(x) x@count_table)

#' @export
setGeneric("count_table<-", function(x, value) standardGeneric("count_table<-"))

#' @export
setMethod("count_table<-", "TCRgrapher", function(x, value) {
  x@count_table <- value
  validObject(x)
  x
})

#' @export
setGeneric("feature_info", function(x) standardGeneric("feature_info"))

#' @export
setMethod("feature_info", "TCRgrapher", function(x) x@feature_info)

#' @export
setGeneric("feature_info<-", function(x, value) standardGeneric("feature_info<-"))

#' @export
setMethod("feature_info<-", "TCRgrapher", function(x, value) {
  x@feature_info <- value
  validObject(x)
  x
})

#' @export
TCRgrCountsValidity <- function(object){
  if(!all(apply(object@count_table, 2, is.numeric))){
    return("All count_table values must be numeric")
  }
  if(!all(apply(TCRgrCounts@count_table, 2, function(x) sum(x < 0) == 0))){
    return("All count_table values must be positive")
  }
  if(!all(rownames(object@count_table) %in% object@feature_info$feature)){
    return("count_table rownames must match feature column from feature_info")
  }
  if(!all(unique(object@metadata$sample_id) %in% object@feature_info$sample_id)){
    return("'sample_id' metadata column and 'sample_id' feature_info column must have the same values")
  }
  if(!all(unique(object@metadata$sample_id) %in% colnames(object@count_table))){
    return("'sample_id' metadata column and count_table column names must have the same values")
  }
  return(TRUE)
}

#' @export
setValidity("TCRgrapherCounts", TCRgrCountsValidity)

#' @export
setMethod("show", "TCRgrapherCounts", function(object) {
  cat(is(object)[[1]], "\n\n",
      "  metadata: ", nrow(object@metadata), " rows and ",
      ncol(object@metadata), " columns\n",
      "  clonoset: ", nrow(object@clonoset), " rows and ",
      ncol(object@clonoset), " columns\n",
      "  count_table: ", nrow(object@count_table), " unique features and ",
      ncol(object@count_table), " samples\n\n",
      "use clonoset(<your_object>) to see clonoset table \n",
      "use metadata(<your_object>) to see metadata table \n",
      "use count_table(<your_object>) to see count_table \n",
      "use feature_info(<your_object>) to see feature_info table \n",
      sep = ""
  )
})

#' take_subset_from_count_table
#'
#' If you would like to take part of the samples into analysis, it is important
#' to take a subset from a count table correctly, because number of rows influence
#' the result.
#'
#' @param count_table a data frame where each row corresponds to unique
#' amino acid clonotype and each column corresponds to the sample.
#' @param samples a vector of names that matches the names of the columns to be
#' taken into analysis
#' @return count_table with specified columns and without rows with zero counts
#' @examples
#' # metadata_CD8 <- metadata[metadata$Cell_Population == 'CD8',]
#' # count_table_CD8 <- take_subset_from_count_table(count_table, metadata_CD8$samples)
#' @export
take_subset_from_count_table <- function(count_table, samples){
  count_table <- count_table[,samples]
  count_table <- count_table[apply(count_table, 1, sum) != 0,]
  count_table
}

#' @export
setMethod("subset", "TCRgrapherCounts", function(x, samples) {
  x@metadata <- x@metadata[sample_id %in% samples]
  x@clonoset <- x@clonoset[sample_id %in% samples]
  x@count_table <- take_subset_from_count_table(x@count_table, samples)
  x@feature_info <- x@feature_info[sample_id %in% samples]
  x
})

# transform_from_clonotypes_to_clusters <- function(TCRgrCounts){
#   # temporary columns
#   count_table <- TCRgrCounts@count_table
#
#   count_table$cdr3aa <- sapply(strsplit(rownames(count_table), ' '), function(x) x[1])
#   if(v_gene){
#     count_table$bestVGene <- sapply(strsplit(rownames(count_table), ' '), function(x) x[2])
#   } else {
#     count_table$bestVGene <- ''
#   }
#   # functions from igraph_capabilities.R
#   g <- make_TCR_graph(count_table)
#   count_table <- find_cluster(count_table, g)
#   # using aggregation from data.table
#   count_table <- setDT(count_table)
#   components_table <- count_table[,.(cdr3aa, bestVGene, cluster_id)]
#   count_table <- count_table[,-c('cdr3aa', 'bestVGene')]
#   count_table <- count_table[,lapply(.SD, sum), by='cluster_id']
#   count_table <- as.data.frame(count_table)
#   # now we have clusters_id instead of clonotypes!
#   rownames(count_table) <- count_table$cluster_id
#   count_table$cluster_id <- NULL
#   # return TODO
#   list('count_table'=count_table, 'components_table'=components_table)
# }

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
#' @param cluster_id Boolean value. Default value is 'FALSE'. If 'cluster_id' is 'TRUE',
#' clonotypes with the same cluster_id will be grouped and will be presented in
#' one row. In feature_info table feature' column will contain cluster_id and
#' genes if corresponding parameters are chosen. If 'cluster_id' is already in the
#' clonotype table it will be taken from there. Otherwise, it will be defined using
#' 'make_TCR_graph' and 'find_TCR_components' functions. See documentation.
#' @return Function returns TCRgrapherCounts object that contains clonoset and
#' metadata from TCRgrapher object, count table and feature info table
#' @export
TCRgrapherCounts <- function(TCRgrObject, v_gene = TRUE, j_gene = FALSE, cluster_id = FALSE){
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
  grouping <- c('cdr3aa', 'cluster_id', 'bestVGene', 'bestJGene')[c(!cluster_id, cluster_id, v_gene, j_gene)]
  # find components if it hasn't been made earlier
  if(!('cluster_id' %in% colnames(clonoset)) & cluster_id){
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop(
        "Package \"igraph\" must be installed to find cluster_ids.",
        call. = FALSE
      )
    }
    #g <- make_TCR_graph(clonoset, v_gene = v_gene, j_gene = j_gene)
    clonoset <- find_TCR_components_by_bfs(clonoset)
  }

  formula <- as.formula(paste0(paste(grouping, collapse = ' + '), ' ~ sample_id'))
  count_table <- dcast(clonoset, formula, value.var = 'count', fun.aggregate = sum)
  count_table <- as.data.frame(count_table)

  if(length(grouping) == 1){
    rownames(count_table) <- count_table[,grouping]
    feature_info <- cbind(clonoset, feature = clonoset[,grouping, with=FALSE])
  } else {
    rownames(count_table) <- Reduce(paste, count_table[,grouping])
    feature_info <- cbind(clonoset, feature = Reduce(paste, clonoset[,grouping, with=FALSE]))
  }

  count_table <- subset(count_table, select = (length(grouping) + 1):(ncol(count_table)))
  feature_info$grouping <- paste(grouping, collapse = ' ')

  new('TCRgrapherCounts', clonoset = clonoset,
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
  if(!all(apply(object@count_table, 2, function(x) sum(x < 0) == 0))){
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

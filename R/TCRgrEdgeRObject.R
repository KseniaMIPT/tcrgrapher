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

# constructor function
# TODO
# create_count_table <- function(TCRgrObject, v_gene = TRUE, by_clusters = FALSE,
#                                custom_groupping = FALSE){
#   # requirements
#   if(!("TCRgrapher" %in% attr(TCRgrObject, 'class'))){
#     stop("The function takes TCRgrapher object as an input. See ?TCRgrapher",
#          call. = FALSE)
#   }
#
#   count_table <- list()
#   feature_info <- list()
#   metadata <- TCRgrObject@metadata
#   clonoset <- TCRgrObject@clonoset
#
#   for(sample_t in metadata$sample_id){
#     sample <- clonoset[sample_id == sample_t]
#     if(!v_gene){
#       sample <- sample[, .(count = sum(count)), by = .(cdr3aa)]
#     } else {
#       sample <- sample[, .(count = sum(count)), by = .(cdr3aa, bestVGene)]
#     }
#
#     setnames(sample, "count", "sample_id")
#     count_table <- merge.data.table(count_table, sample,
#                                     by=c('cdr3aa', 'bestVGene'), all=TRUE)
#   }
#
#   count_table <- Reduce(count_table, function(x, y) merge.data.table(x, y, by=c('cdr3aa', 'bestVGene'), all = TRUE))
#
#   count_table <- as.data.frame(count_table)
#   if(v_gene){
#     rownames(count_table) <- paste(count_table$cdr3aa, count_table$bestVGene)
#   } else {
#     rownames(count_table) <- count_table$cdr3aa
#   }
#   count_table <- count_table[,-(1:2)]
#
#
#   new('TCRgrapherCounts', clonoset = TCRgrObject@clonoset,
#       metadata = TCRgrObject@metadata, count_table = count_table,
#       feature_info = feature_info)
# }

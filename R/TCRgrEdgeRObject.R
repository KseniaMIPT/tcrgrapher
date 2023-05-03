#' @import methods
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom data.table :=

# S4 class
setClass('TCRgrapherCounts',
         contains = 'TCRgrapher',
         slots = c(
           count_table = 'data.table',
           feature_info = 'data.table'
         ))

# constructor function
# TODO
# create_count_table <- function(..., TCRgrapherObject = NA, v_gene = TRUE){
#   if(is.na(TCRgrapherObject)){
#     TCRgrapher(...)
#   } else {
#     TCRgrapherObject
#   }
#   # count_table <- data.table(cdr3aa='0', bestVGene='0')
#   # metadata <- metadata(TCRgrObject)
#   # for(i in 1:nrow(metadata)){
#   #   sample_id <- metadata$sample_id[i]
#   #   sample <- clonoset(TCRgrObject)[sample == sample_id]
#   #   if(!v_gene){
#   #     sample$bestVGene <- ''
#   #   }
#   #   sample <- sample[, .(count = sum(count)), by = .(cdr3aa, bestVGene)]
#   #   setnames(sample, "count", sample_id)
#   #   count_table <- merge.data.table(count_table, sample,
#   #                                   by=c('cdr3aa', 'bestVGene'), all=TRUE)
#   # }
#   # count_table <- count_table[cdr3aa != '0']
#   # count_table[is.na(count_table), ] <- 0
#   # count_table <- as.data.frame(count_table)
#   # if(v_gene){
#   #   rownames(count_table) <- paste(count_table$cdr3aa, count_table$bestVGene)
#   # } else {
#   #   rownames(count_table) <- count_table$cdr3aa
#   # }
#   # count_table <- count_table[,-(1:2)]
#   # count_table
# }

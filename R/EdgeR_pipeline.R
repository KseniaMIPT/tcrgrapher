#' edgeR_pipeline
#'
#' The function performs statistical analysis by edgeR to identify significantly
#' expanded clonotypes. First, it filters the data using relatively mild conditions
#' that better suit TCR repertoire data. Second, it normalize counts and estimate
#' dispersion by standard edgeR methods. Finally, it performs all pairwise comparisons
#' between groups and compares each group vs all others using glm and QL F-test.
#'
#' @param TCRgrCounts a TCRgrapherCounts object. For more details see ?TCRgrapherCounts
#' @param comparison a name of the column that specifies levels of comparison in the metadata
#' @param min.count parameter from edgeR::filterByExpr function. Minimum count
#' required for at least some samples. Default value is 1.
#' @param min.total.count parameter from edgeR::filterByExpr function.
#' Minimum total count required. Default value is 5.
#' @param large.n parameter from edgeR::filterByExpr function. Number of samples
#'  per group that is considered to be “large”. Default value is 1.
#' @param min.prop parameter from edgeR::filterByExpr function. Minimum
#' proportion of samples in the smallest group that express the gene. Default
#' value is 0.5.
#' @param normalization_method method to be used to edgeR::calcNormFactors function
#' @return data.table with statistics for every clonotype performed for each pair
#' of groups and for every sample vs all others. 'comparison' column shows
#' comparisons in the format 'Group1 vs Group2'.
#' In this case, LFC values are positive if the clonotype is expanded in 'Group1' and
#' negative if it is expanded in 'Group2'.
#' @examples
#' # find significantly expanded clonotypes after vaccination
#' # my_metadata table contains column 'vaccination_status' with groups of comparison
#' # edgeR_pipeline(my_count_table, my_metadata, 'vaccination_status')
#' @export
edgeR_pipeline <- function(TCRgrCounts, comparison, min.count = 1,
                           min.total.count = 5, large.n = 1, min.prop = 0.5,
                           normalization_method ="TMM"){
  # Important conditions
  if(!requireNamespace("edgeR", quietly = TRUE)){
    stop("Package \"edgeR\" must be installed and loaded to use this function.",
         call. = FALSE)
  }
  if(!("TCRgrapherCounts" %in% attr(TCRgrCounts, 'class'))){
    stop("The function takes TCRgrapherCounts object as an input. See ?TCRgrapherCounts and ?TCRgrapher",
         call. = FALSE)
  }

  setDTthreads(threads = 0)
  count_table <- TCRgrCounts@count_table
  metadata <- TCRgrCounts@metadata

  if(!('sample_id' %in% colnames(metadata))){
    stop("There is no column 'sample_id' in metadata.", call. = FALSE)
  }
  if(!all(colnames(count_table) %in% metadata$sample_id)){
    stop("The names of the count table columns and the 'sample_id' metadata column do not match.",
         call. = FALSE)
  }
  if(!is.character(comparison) | !(comparison %in% colnames(metadata))){
    stop("Specify the factor by which the comparison will be made.
    It should be character corresponding to particular metadata column name.",
         call. = FALSE)
  }
  # Data preparation
  metadata <- as.data.frame(metadata)
  metadata[,comparison] <- as.factor(metadata[,comparison])
  rownames(metadata) <- metadata$sample_id
  design <- model.matrix(as.formula(paste('~ 0 +', comparison)), data=metadata)
  sample <- DGEList(counts=count_table)
  # Data filtration
  keep <- filterByExpr(sample, design,
                       min.count = min.count, min.total.count = min.total.count,
                       large.n = large.n, min.prop = min.prop)
  sample <- sample[keep, , keep.lib.sizes=FALSE]
  # Normalization and dispersion estimation
  sample <- calcNormFactors(sample, method=normalization_method)
  sample <- estimateDisp(sample, design)
  # Fit a quasi-likelihood negative binomial generalized log-linear model to count data
  fit <- glmQLFit(sample, design)
  # All pairwise comparisons by QL F-test
  sign_result <- c()
  comparison_levels <- colnames(design)
  nb_of_comparison_levels <- length(comparison_levels)
  # pairwise comparisons
  for(i in 1:(nb_of_comparison_levels-1)){
    nb_of_rows = nb_of_comparison_levels - i
    comparison_matrix_pairwise <- cbind(matrix(0, nrow = nb_of_rows, ncol = i-1),
                                        matrix(1, nrow = nb_of_rows, ncol = 1),
                                        -1*diag(nb_of_rows))
    for(j in 1:nb_of_rows){
      qlf <- glmQLFTest(fit, contrast = comparison_matrix_pairwise[j,])
      topTags <- topTags(qlf, n = nrow(count_table), p.value = 1)$table
      topTags$comparison <- paste(comparison_levels[i], 'vs', comparison_levels[i+j])
      topTags$feature <- rownames(topTags)
      sign_result <- rbind(sign_result, topTags)
    }
  }
  # each group vs all others
  if(nb_of_comparison_levels > 2){
    comparison_matrix_vs_all <- matrix(-1/(nb_of_comparison_levels-1),
                                       nb_of_comparison_levels,
                                       nb_of_comparison_levels)
    comparison_matrix_vs_all[row(comparison_matrix_vs_all) == col(comparison_matrix_vs_all)] <- 1
    for(i in 1:nb_of_comparison_levels){
      qlf <- glmQLFTest(fit,contrast = comparison_matrix_vs_all[i,])
      topTags <- topTags(qlf, n = nrow(count_table), p.value = 1)$table
      topTags$comparison <- paste(comparison_levels[i], 'vs all')
      topTags$feature <- rownames(topTags)
      sign_result <- rbind(sign_result, topTags)
    }
  }
  sign_result <- setDT(sign_result)
  setcolorder(sign_result, c("feature", setdiff(names(sign_result), "feature")))
  sign_result
}

#' filter_edgeR_res
#'
#' The function takes "vs all" comparisons and checks if all pairwise comparisons
#' are consistent with the given "vs all" comparison. 'the_worst_pairwise_p' column
#' shows the worst p value in all pairwise comparisons.
#'
#' @param res_dt the output of edgeR_pipeline function
#' @return subset of res_dt with only 'vs all' comparisons and additional columns
#' @export
filter_edgeR_res <- function(res_dt){
  res_dt_to_check <- res_dt[str_detect(res_dt$comparison, 'vs all'),]
  res_dt_to_check$consistent <- FALSE
  res_dt_to_check$the_worst_pairwise_p <- 1
  for(i in 1:nrow(res_dt_to_check)){
    feature_t <- res_dt_to_check[i,feature]
    level <- unlist(str_split(res_dt_to_check[i, comparison], ' vs all'))[1]
    dt_t <- res_dt[feature == feature_t]
    dt_t <- dt_t[str_detect(dt_t$comparison, level),]
    up_level <- sapply(str_split(dt_t$comparison, ' vs '), function(x) x[[1]])
    res_dt_to_check[i, 'consistent'] <- nrow(dt_t[up_level == level & logFC > 0]) + nrow(dt_t[up_level != level & logFC < 0]) == nrow(dt_t)
    dt_t <- dt_t[!str_detect(dt_t$comparison, 'vs all')]
    res_dt_to_check[i, 'the_worst_pairwise_p'] <- max(dt_t$PValue)
  }
  return(res_dt_to_check)
}

#' wilcox_pipeline
#'
#' The function performs all pairwise comparisons between groups and compares
#' each group vs all others using wilcox test to identify significantly
#' expanded clonotypes of clusters of clonotypes.
#'
#' @param TCRgrCounts a TCRgrapherCounts object. For more details see ?TCRgrapherCounts
#' @param comparison a name of the column that specifies levels of comparison in the metadata
#' @return data.table with statistics for every clonotype or cluster of clonotypes
#' performed for each pair of groups and for every group vs. all others.
#' 'comparison' column shows comparisons in the format 'Group1 vs Group2'.
#' P-values correspond to the alternative hypothesis that the median of the first
#' group is higher than the median of the second group.
#' @export
wilcox_pipeline <- function(TCRgrObject, comparison){
  if(!requireNamespace("stats", quietly = TRUE)){
    stop("Package \"stats\" must be installed and loaded to use this function.",
         call. = FALSE)
  }
  metadata <- as.data.frame(metadata(TCRgrObject))
  counts <- count_table(TCRgrObject)
  sign_result <- c()
  # pairwise comparisons
  pairs <- combn(unique(metadata[,comparison]), 2)
  for(i in 1:ncol(pairs)){
    pair <- pairs[,i]
    curr_comparison <- (metadata[,comparison] == pair[1])*1 + (metadata[,comparison] == pair[2])*(-1)
    for(feature in rownames(counts)){
      p_val_t <- wilcox.test(as.numeric(counts[feature, curr_comparison == 1]),
                             as.numeric(counts[feature, curr_comparison == -1]),
                             alternative = 'greater')$p.value
      comparison_t <- paste(pair[1], 'vs', pair[2])
      sign_result <- rbind(sign_result, cbind(feature, comparison_t, p_val_t))
    }
  }
  # each group vs all others
  if(length(unique(metadata[,comparison])) > 2){
    comp_levels <- unique(metadata[,comparison])
    for(i in 1:length(comp_levels)){
      curr_comparison <- (metadata[,comparison] == comp_levels[i])*1 + (metadata[,comparison] != comp_levels[i])*(-1)
      for(feature in rownames(counts)){
        p_val_t <- wilcox.test(as.numeric(counts[feature, curr_comparison == 1]),
                               as.numeric(counts[feature, curr_comparison == -1]),
                               alternative = 'greater')$p.value
        comparison_t <- paste(comp_levels[i], 'vs all')
        sign_result <- rbind(sign_result, cbind(feature, comparison_t, p_val_t))
      }
    }
  }
  colnames(sign_result) <- c('feature', 'comparison', 'p_value_greater')
  sign_result
}

#' filter_wilcox_res
#'
#' The function takes "vs all" comparisons and checks if all pairwise comparisons
#' are consistent with the given "vs all" comparison.
#'
#' @param res_dt the output of wilcox_pipeline function
#' @return subset of res_dt with only 'vs all' comparisons and additional column
#' 'consistent'
#' @export
filter_wilcox_res <- function(res_dt){
  res_dt_to_check <- res_dt[str_detect(res_dt$comparison, 'vs all'),]
  res_dt_to_check$consistent <- FALSE
  for(i in 1:nrow(res_dt_to_check)){
    feature_t <- res_dt_to_check[i,feature]
    level <- unlist(str_split(res_dt_to_check[i, comparison], ' vs all'))[1]
    dt_t <- res_dt[feature == feature_t]
    dt_t <- dt_t[str_detect(dt_t$comparison, level),]
    up_level <- sapply(str_split(dt_t$comparison, ' vs '), function(x) x[[1]])
    res_dt_to_check[i, 'consistent'] <- nrow(dt_t[up_level == level & p_value_greater < 0.5]) + nrow(dt_t[up_level != level & p_value_greater > 0.5]) == nrow(dt_t)
    dt_t <- dt_t[!str_detect(dt_t$comparison, 'vs all')]
  }
  return(res_dt_to_check)
}

#' heatmap_expanded
#'
#' The function takes an output of edgeR_pipeline or wilcox_pipeline functions
#' and creates a heatmap using the ComplexHeatmap library. It is recommended
#' to use filter_edgeR_res or filter_wilcox_res previously and filter results by
#' some threshold.
#' @param TCRgrCounts TCRgrapherCounts object. For more details see ?TCRgrapherCounts
#' @param expanded_test_res table with features (clonotypes or clusters of clonotypes)
#' @param comparison a name of the column that specifies levels of comparison in the metadata
#' to draw. It could be an output of edgeR_pipeline or wilcox_pipeline.
#' @export
heatmap_expanded <- function(TCRgrCounts, expanded_test_res, comparison){
  if(!requireNamespace("ComplexHeatmap", quietly = TRUE)){
    stop("Package \"ComplexHeatmap\" must be installed and loaded to use this function.",
         call. = FALSE)
  }
  count_table <- count_table(TCRgrCounts)
  count_table <- count_table[unique(expanded_test_res$feature),]
  data <- log(count_table)
  colnames(data) <- colnames(data)
  data[data == -Inf] <- -1
  data <- as.matrix(data)
  pht <- ComplexHeatmap::pheatmap(data, cluster_cols = FALSE,
                                  column_split = metadata(TCRgrCounts)[,comparison])
  pht
}

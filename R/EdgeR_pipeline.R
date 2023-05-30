#' edgeR_pipeline
#'
#' Function performs statistical analysis by edgeR to identify significantly
#' expanded clonotypes. First, it filters the data using relatively mild conditions
#' that better suit TCR repertoire data. Second, it normalize counts and estimate
#' dispersion by standard edgeR methods. Finally, it performs all pairwise comparisons
#' between groups and compares each group vs all others using glm and QL F-test.
#'
#' @param count_table a data frame where each row corresponds to unique
#' amino acid clonotype and each column corresponds to the sample.
#' @param comparison name of the column that specifies levels of comparison in the metadata
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
#' In this case, LFC values are negative if the clonotype is expanded in 'Group1' and
#' positive if it is expanded in 'Group2'.
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
                                        matrix(-1, nrow = nb_of_rows, ncol = 1),
                                        diag(nb_of_rows))
    for(j in 1:nb_of_rows){
      qlf <- glmQLFTest(fit,contrast = comparison_matrix_pairwise[j,])
      topTags <- topTags(qlf, n = nrow(count_table), p.value = 1)$table
      topTags$comparison <- paste(comparison_levels[i], 'vs', comparison_levels[i+j])
      topTags$feature <- rownames(topTags)
      sign_result <- rbind(sign_result, topTags)
    }
  }
  # each group vs all others
  if(nb_of_comparison_levels > 2){
    comparison_matrix_vs_all <- matrix(-1/(nb_of_comparison_levels-1), nb_of_comparison_levels, nb_of_comparison_levels)
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

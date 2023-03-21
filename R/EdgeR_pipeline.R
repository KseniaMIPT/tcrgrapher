#' create_count_table
#'
#' Function takes clonotype tables from specified directory and creates a table
#' where each row corresponds to unique amino acid clonotype and each column
#' corresponds to the sample.
#'
#' @param dir path to the directory where the clonotype tables are located
#' @param metadata table should contain 'sample' column with names of the files
#' in directory specified earlier. Samples should contain 'cdr3aa' and 'count'
#' columns.
#' @return Function returns count table with unique amino acid sequences as
#' row names and sample names as column names
#' @export
create_count_table <- function(dir, metadata){
  count_table <- data.table(cdr3aa='0')
  for(sample_name in metadata$sample){
    sample <- fread(paste0(dir, sample_name))
    sample <- sample[, .(count = sum(count)), by = 'cdr3aa']
    setnames(sample, "count", sample_name)
    count_table <- merge.data.table(count_table, sample, by='cdr3aa', all=TRUE)
  }
  count_table <- count_table[cdr3aa != '0']
  count_table[is.na(count_table), ] <- 0
  count_table <- as.data.frame(count_table)
  rownames(count_table) <- count_table$cdr3aa
  count_table <- count_table[,-1]
  count_table
}

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

#' edgeR_pipeline
#'
#' Function performs statistical analysis by edgeR to identify significantly
#' expanded clonotypes. First, it filters the data using relatively mild conditions
#' that better suit TCR repertoire data. Second, it normalize counts and estimate
#' dispersion by standard edgeR methods. Finally, it performs all comparisons
#' between groups using glm and QL F-test.
#'
#' @param count_table a data frame where each row corresponds to unique
#' amino acid clonotype and each column corresponds to the sample.
#' @param metadata A data frame should contain the column 'sample'. The names of
#' the count table columns and the 'sample' metadata column must be the same,
#' and they must be in the same order. Also, metadata must contain a column that
#' specifies levels of comparison. The name of the column must be the same as
#' in 'comparison' variable.
#' @param min.count parameter from edgeR::filterByExpr function. Minimum count
#' required for at least some samples. Default value is 1.
#' @param min.total.count parameter from edgeR::filterByExpr function.
#' Minimum total count required. Default value is 5.
#' @param large.n parameter from edgeR::filterByExpr function. Number of samples
#'  per group that is considered to be “large”. Default value is 1.
#' @param min.prop parameter from edgeR::filterByExpr function. Minimum
#' proportion of samples in the smallest group that express the gene. Default
#' value is 0.5.
#' @param alpha cutoff value for adjusted p-values. Only clonotypes with adjusted
#' p-values equal or lower than specified are returned
#' @param normalization_method method to be used to edgeR::calcNormFactors function
#' @return data.frame with statistics for every clonotype performed for each pair
#' of groups and for every sample vs all others. 'comparison' column shows
#' comparisons in the format 'Group1 vs Group2'.
#' In this case, LFC values are negative if the clonotype is expanded in 'Group1' and
#' positive if it is expanded in 'Group2'.
#' @examples
#' # find significantly expanded clonotypes after vaccination
#' # my_metadata table contains column 'vaccination_status' with groups of comparison
#' # to get information for every clonotype alpha = 1 is specified
#' # edgeR_pipeline(my_count_table, my_metadata, 'vaccination_status', alpha = 1)
#' @export
edgeR_pipeline <- function(count_table, metadata, comparison, min.count = 1,
                           min.total.count = 5, large.n = 1, min.prop = 0.5,
                           alpha = 0.05, normalization_method ="TMM"){
  # Important conditions
  if(!requireNamespace("edgeR", quietly = TRUE)){
    stop("Package \"edgeR\" must be installed to use this function.",
      call. = FALSE)
  }
  if(!('sample' %in% colnames(metadata))){
    stop("There is no column 'sample' in metadata.", call. = FALSE)
  }
  if(colnames(count_table) != metadata$sample){
    stop("The names of the count table columns and the 'sample' metadata column do not match.",
         call. = FALSE)
  }
  if(!is.character(comparison) | !(comparison %in% colnames(metadata))){
    stop("Specify the factor by which the comparison will be made.
    It should be character corresponding to particular metadata column name.",
         call. = FALSE)
  }
  # Data preparation
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata$sample
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
  nb_of_comparisons <- length(comparison_levels)
  # pairwise comparisons
  for(i in 1:(nb_of_comparisons-1)){
    nb_of_rows = nb_of_comparisons - i
    comparison_matrix_pairwise <- cbind(matrix(0, nrow = nb_of_rows, ncol = i-1),
                               matrix(-1, nrow = nb_of_rows, ncol = 1),
                               diag(nb_of_rows))
    for(j in 1:nb_of_rows){
      qlf <- glmQLFTest(fit,contrast = comparison_matrix_pairwise[j,])
      topTags <- topTags(qlf, n = nrow(count_table), p.value = alpha)$table
      topTags$comparison <- paste(comparison_levels[i], 'vs', comparison_levels[i+j])
      topTags$feature <- rownames(topTags)
      sign_result <- rbind(sign_result, topTags)
    }
  }
  # each group vs all others
  comparison_matrix_vs_all <- matrix(-1/(nb_of_comparisons-1), nb_of_comparisons, nb_of_comparisons)
  comparison_matrix_vs_all[row(comparison_matrix_vs_all) == col(comparison_matrix_vs_all)] <- 1
  for(i in 1:nb_of_comparisons){
    qlf <- glmQLFTest(fit,contrast = comparison_matrix_vs_all[i,])
    topTags <- topTags(qlf, n = nrow(count_table), p.value = alpha)$table
    topTags$comparison <- paste(comparison_levels[i], 'vs all')
    topTags$cluster <- rownames(topTags)
    sign_result <- rbind(sign_result, topTags)
  }
  sign_result
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ALICE_pipeline.R
\name{pval_with_abundance}
\alias{pval_with_abundance}
\title{pval_with_abundance}
\usage{
pval_with_abundance(clonoset, by_counts = TRUE, by_log2_counts = TRUE)
}
\arguments{
\item{clonoset}{clonoset after the analysis with the ALICE pipeline.
To get it, use clonoset(TCRgrapher object)}

\item{by_counts}{A logical value indicates whether to calculate corrected
p-values using counts or not}

\item{by_log2_counts}{A logical value indicates whether to calculate corrected
p-values using log2 counts or not}
}
\value{
Function returns the same clonoset with additional columns:
\itemize{
\item{pval_with_abundance_log2_counts}{recalculated p-value considering count
number of every clonotype. Log2 is used for count normalization}
\item{pval_with_abundance_counts}{recalculated p-value considering count number
of every clonotype. There is no count normalization}
\item{log_pval_with_abundance_log2_counts}{log of pval_with_abundance_log2_counts}
\item{log_pval_with_abundance_counts}{log of pval_with_abundance_counts}
}
}
\description{
The function calculates the p-value, taking into account the abundance of every
clonotype. The method is described in Pogorelyy et al. 2019. In the case of
very large datasets with high numbers of neighbors, it may not be possible to
calculate the corrected p-value, especially pval_with_abundance_counts.
}

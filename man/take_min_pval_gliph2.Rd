% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_GLIPH2.R
\name{take_min_pval_gliph2}
\alias{take_min_pval_gliph2}
\title{take_min_pval_gliph2}
\usage{
take_min_pval_gliph2(gliph2_clusters, TCRgrObject)
}
\arguments{
\item{gliph2_clusters}{clusters.csv from gliph2 output}

\item{TCRgrObject}{TCRgrapher object}
}
\value{
The function returns a vector of minimum p-values
}
\description{
The function takes clusters.csv from gliph2 output and assigns the minimum p-value
for every clonotype from the TCRgrapher clonotype table.
}

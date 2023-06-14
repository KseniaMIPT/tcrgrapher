#' calc_TCRdist3_radius
#'
#' The function takes TCRgrapherObject as an input and finds the optimal radius
#' with tcrdist3 python package.
#'
#' @param TCRgrObject See ?TCRgrapher
#' @param cores the number of cores to use
#' @export
calc_TCRdist3_radius <- function(TCRgrObject, cores, organism = 'mouse',
                                 chain = 'beta'){
  if(!requireNamespace("reticulate", quietly = TRUE)){
    stop("Package \"reticulate\" must be installed and loaded to use this function.",
         call. = FALSE)
  }
  clonoset <- clonoset(TCRgrObject)
  if(!('freq' %in% colnames(clonoset))){
    clonoset$freq <- clonoset$count / sum(clonoset$count)
    message('Frequency values were calculated from count column')
  }
  source_python(system.file('TCRdist3.py', package = 'tcrgrapher'))
  res <- tcrdist_radii(clonoset, cores, organism, chain)
  clonoset(TCRgrObject)$tcrdist3.radius <- res$radius
  TCRgrObject
}

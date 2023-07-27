#' calc_TCRdist3_radius
#'
#' The function takes TCRgrapher object as an input and finds the optimal radius
#' with tcrdist3 python package.
#'
#' @param TCRgrObject See ?TCRgrapher
#' @param cores the number of cores to use. The default value is 1.
#' @param organism Possible options: mouse, human. The default value is "mouse"
#' @param chain Possible options: alpha, beta. The default value is "beta"
#' @param max_radius Only distances that are less than or equal to max_radius are stored
#' @return TCRgrapher object. A clonoset contains additional column "tcrdist3.radius"
#' @export
calc_TCRdist3_radius <- function(TCRgrObject, cores = 1, organism = 'mouse',
                                 chain = 'beta', max_radius = 50){
  message('Python package tcrdist3 must be installed. Installation: "pip install
          git+https://github.com/kmayerb/tcrdist3.git@0.2.2". See documentation
          https://tcrdist3.readthedocs.io/en/latest/')
  if(!requireNamespace("reticulate", quietly = TRUE)){
    stop("Package \"reticulate\" must be installed and loaded to use this function.",
         call. = FALSE)
  }
  clonoset <- clonoset(TCRgrObject)
  if(!('freq' %in% colnames(clonoset))){
    clonoset$freq <- clonoset$count / sum(clonoset$count)
    message('Frequency values were calculated from count column')
  }
  source_python(system.file('Python/TCRdist3.py', package = 'tcrgrapher'))
  res <- tcrdist_radii(clonoset, cores, organism, chain, max_radius)
  clonoset(TCRgrObject)$tcrdist3.radius <- as.numeric(res)
  return(TCRgrObject)
}

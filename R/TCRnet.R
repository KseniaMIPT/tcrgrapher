#' run_TCRnet
#'
#' The function runs TCRnet (CalcDegreeStats command from vdjtools). VDJtools must
#' be installed. See documentation at https://vdjtools-doc.readthedocs.io/en/master/
#'
#' @param command A command or a path to run on the command line. The default
#' value is 'vdjtools'. For example, the user can specify 'java -jar /your_path/vdjtools-1.2.1.jar'
#' @param background A path to the background sample
#' @param search_scope A vector with three numbers: allowed number of substitutions (s),
#' indels (id) and total number of mismatches (t). The default value is c(1,0,1)
#' @param grouping The default value is 'vj'. Possible values: 'vj', 'vjl', 'v', 'dummy'
#' @export
run_TCRnet <- function(TCRgrObject, background_path, command = 'vdjtools',
                       search_scope = c(1,0,1), grouping = 'vj'){
  background <- fread(background_path)
  search_scope <- paste(search_scope, collapse = ',')
  dt <- clonoset(TCRgrObject)
  if(!('freq' %in% colnames(dt))){
    dt$freq <- dt$count / sum(dt$count)
  }
  dt$d <- '-'
  dt <- dt[,c('count', 'freq', 'cdr3nt', 'cdr3aa', 'bestVGene', 'd', 'bestJGene')]
  setnames(dt, c('bestVGene', 'bestJGene'), c('v', 'j'))
  path <- tempdir()
  input_path <- paste0(path, '/input.tsv')
  write.table(dt, file = input_path, sep = '\t', quote = F, row.names = F)
  output_path <- paste0(path, '/output.tsv')
  cmd <- paste(command, 'CalcDegreeStats -b', background, '-o', search_scope,
                '-g', grouping, '-g2', grouping, input_path, output_path)
  system(cmd)
  res <- read.delim(output_path)
  res <- res[, 8:ncol(res)]
  colnames(res) <- paste0('TCRNET.', colnames(res))
  clonoset(TCRgrObject) <- cbind(clonoset(TCRgrObject), res)
  return(TCRgrObject)
}

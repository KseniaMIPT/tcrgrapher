#     TCRgrapher: R package for identifying condition associated T cell clonotypes
#     Copyright (C) 2024 Kseniia Lupyr
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' run_TCRNET
#'
#' The function runs TCRNET (CalcDegreeStats command from vdjtools). VDJtools must
#' be installed. See documentation at https://vdjtools-doc.readthedocs.io/en/master/
#'
#' @param command A command or a path to run on the command line. The default
#' value is 'vdjtools'. For example, the user can specify 'java -jar /your_path/vdjtools-1.2.1.jar'
#' @param background_path A path to the background sample
#' @param search_scope A vector with three numbers: allowed number of substitutions (s),
#' indels (id) and total number of mismatches (t). The default value is c(1,0,1)
#' @param grouping The default value is 'vj'. Possible values: 'vj', 'vjl', 'v', 'dummy'
#' @export
run_TCRNET <- function(TCRgrObject, background_path, command = 'vdjtools',
                       search_scope = c(1,0,1), grouping = 'vj'){
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
  cmd <- paste(command, 'CalcDegreeStats -b', background_path, '-o', search_scope,
                '-g', grouping, '-g2', grouping, input_path, path)
  system(cmd)
  res <- fread(paste0(path, '/input.txt'))
  res <- res[, 12:ncol(res)]
  colnames(res) <- paste0('TCRNET.', colnames(res))
  clonoset(TCRgrObject) <- cbind(clonoset(TCRgrObject), res)
  return(TCRgrObject)
}

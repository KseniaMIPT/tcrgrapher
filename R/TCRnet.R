#' run_TCRnet
#'
#' The function runs TCRnet (CalcDegreeStats command from vdjtools). VDJtools must
#' be installed. See documentation at https://vdjtools-doc.readthedocs.io/en/master/
#'
#' @param command A command or a path to run on the command line. The default
#' value is 'vdjtools'. For example, the user can specify 'java -jar /your_path/vdjtools-1.2.1.jar'
#' @param background A user may specify a path to the background sample.
#' The package includes several background samples. The default value is
#' 'olga_gen_mouseTRB' which corresponds to the olga generated sample of mouse TRB chains.
#' The result was pooled by amino acid sequences, V and J segments.
#' @param search_scope A vector with three numbers: allowed number of substitutions (s),
#' indels (id) and total number of mismatches (t). The default value is c(1,0,1)
#' @param grouping The default value is 'vj'. Possible values: 'vj', 'vjl', 'v', 'dummy'
#' @export
run_TCRnet <- function(input_path, output_path, command = 'vdjtools', background = 'olga_gen_mouseTRB',
                       search_scope = c(1,0,1), grouping = 'vj'){
  if(background == 'olga_gen'){
    background <- system.file('extdata/background_sample/C57BL6_all_genes_with_counts.pool.aaVJ.table.txt',
                              package = 'tcrgrapher')
  }
  search_scope <- paste(search_scope, collapse = ',')
  cmd <- paste(command, 'CalcDegreeStats -b', background, '-o', search_scope,
                '-g', grouping, '-g2', grouping, input_path, output_path)
  system(cmd)
  # TODO work with TCRgrapher objects
}

#' make_gen_model
#'
#'The function uses the IGoR software to infer a new generation model. IGoR should
#' be installed. The recommended version is v1.3.0. The function takes the path to
#' the file with non-functional nucleotide CDR3 sequences and writes a generation
#' model based on these sequences to the working directory. To use more features,
#' one can use IGoR directly. See the documentation at https://qmarcou.github.io/IGoR/#version.
#'
#' @param nonfunc_path A path to the file with nonfunctional nucleotide sequences.
#' The file should contain only sequences (one sequence per line).
#' @param command a command or a path to run on the command line. The default
#' value is 'igor'. For example, the user can specify '/user/.local/bin/igor'.
#' @param wd a name of the working directory. The default value is 'igor_output'.
#' @param batch_name The default value is 'new_inference'
#' @param species Possible values: mouse, human. If "human" is set, the default
#' IGoR genomic template is used. If "mouse" is set, the default TCRgrapher
#' genomic template is used. See Details.
#' @param chain Possible values for human model: alpha, beta. Possible values for
#' mouse model: beta.If "human" is set, the default IGoR genomic template is used.
#' If "mouse" is set, the default TCRgrapher genomic template is used. See Details.
#' @details
#' Additional details...How tcgrapher genomic template was made...# TODO
#'
#' @export
make_gen_model <- function(nonfunc_path = 'test', command = 'igor', wd = 'igor_output',
                           batch_name = 'new_inference', species = 'mouse',
                           chain = 'beta'){
  cmd_read <- paste0(command, ' -set_wd ', wd, ' -batch ', batch_name,
                     ' -read_seqs ', nonfunc_path)
  read_human <- paste0(' -species human -chain ', chain)
  cmd_align <- paste0(command, ' -set_wd ', wd, ' -batch ', batch_name, ' -align --all --ntCDR3')
  set_anchors_and_genomic <- paste0(' -set_CDR3_anchors --V ',
                            system.file('extdata/mouse_TRB_reference/V_gene_CDR3_anchors.csv', package = 'tcrgrapher'),
                            ' --J ', system.file('extdata/mouse_TRB_reference/J_gene_CDR3_anchors.csv', package = 'tcrgrapher'),
                            ' -set_genomic --V ', system.file('extdata/mouse_TRB_reference/genomicVs.fasta', package = 'tcrgrapher'),
                            ' --D ', system.file('extdata/mouse_TRB_reference/genomicDs.fasta', package = 'tcrgrapher'),
                            ' --J ', system.file('extdata/mouse_TRB_reference/genomicJs.fasta', package = 'tcrgrapher'))
  cmd_infer <- paste0(command, ' -set_wd ', wd, ' -batch ',  batch_name,
                                   ' -infer --N_iter 10')
  infer_mouse_TRB <- paste0(' -set_custom_model ',
                            system.file('extdata/mouse_TRB_reference/model_params.txt', package = 'tcrgrapher'),
                            system.file('extdata/mouse_TRB_reference/model_marginals.txt', package = 'tcrgrapher'))
  if(species == 'mouse' & chain == 'beta'){
    system(cmd_read)
    system(paste0(cmd_align, set_anchors_and_genomic))
    system(paste0(cmd_infer, set_anchors_and_genomic, infer_mouse_TRB))
  } else {
    system(paste0(cmd_read, read_human))
    system(cmd_align)
    system(cmd_infer)
  }
}

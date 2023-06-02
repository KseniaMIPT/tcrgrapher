#' make_gen_model
#'
#' The function uses the IGoR software to infer a new generation model.
#'
#' @param nonfunc_path A path to the file with nonfunctional nucleotide sequences.
#' The file should contain only sequences (one sequence per line).
#' @param command a command or a path to run on the command line. The default
#' value is 'igor'. For example, the user can specify '/user/.local/bin/igor'.
#' @param wd a name of the working directory. The default value is 'igor_output'.
#' @param batch_name The default value is 'new_inference'
#' @export
make_gen_model <- function(nonfunc_path = 'test', command = 'igor', wd = 'igor_output',
                           batch_name = 'new_inference'){
  # cmd_read <- paste0(command, ' -set_wd ', wd, ' -batch ', batch_name,
  #                    ' -read_seqs ', nonfunc_path)
  # cmd_align <- paste0(command, ' -set_wd ', wd, ' -batch ', batch_name,
  #                     ' -align --all --ntCDR3 -set_CDR3_anchors --V ', my_reference/V_gene_CDR3_anchors_01.csv --J my_reference/J_gene_CDR3_anchors_01.csv -set_genomic --V my_reference/genomicVs_01.fasta --D my_reference/genomicDs_01.fasta --J my_reference/genomicJs_01.fasta)
  # cmd_align_my_mouse_ref <- paste0(command, ' -set_wd igor_my_ref_Th_01 -batch new -infer --N_iter 10 -set_CDR3_anchors --V my_reference/V_gene_CDR3_anchors_01.csv --J my_reference/J_gene_CDR3_anchors_01.csv -set_genomic --V my_reference/genomicVs_01.fasta --D my_reference/genomicDs_01.fasta --J my_reference/genomicJs_01.fasta -set_custom_model /home/ksenia/sc_LCMV/olga/model_params.txt /home/ksenia/sc_LCMV/olga/model_marginals.txt
}

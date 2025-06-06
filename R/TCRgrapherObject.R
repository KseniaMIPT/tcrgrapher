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

#' @import methods
#' @import data.table

# secondary functions
check_path <- function(files_path){
  if(substring(files_path, nchar(files_path)) != '/'){
    files_path <- paste0(files_path, '/')
  }
  files_path
}

#' TCRgrapher object
#'
#' TCRgrapher is a S4 class that can be constructed by calling the TCRgrapher( )
#' function. See "Details". TCRgrapher contains clonoset and metadata as
#' data.tables. You can get them by calling metadata(x) or clonoset(x).
#'
#' There are three ways to initialize an object.
#'
#' (1) Specify the path to the one file without metadata
#'
#' If the path for only one file is specified, metadata will be produced automatically
#' and will contain one row and two columns: "file" and "sample_id".
#' Clonoset will have an additional column "sample_id" with one unique value and a
#' column "clone_id" with an unique id for every row.
#' Positions of "count", "cdr3nt", "cdr3aa", "V gene", "J gene" clonoset's columns must be specified.
#'
#' (2) Specify the path to the directory with files without metadata
#'
#' If the path for the directory with files is specified, and metadata is not specified
#' metadata will be produced automatically and will have a number of rows equal
#' to the number of samples in the directory. All files will be merged into one
#' data.table with additional columns "sample_id" and "clone_id" with
#' unique ids for every row. Positions of "count", "cdr3nt", "cdr3aa",
#' "V gene", "J gene" clonoset's columns must be specified.
#'
#' (3) Specify the path to the directory with files and metadata
#'
#' If the path for the directory with files and the path to the metadata are specified,
#' only files from the metadata column "file" will be taken into consideration.
#' Clonosets will be merged into one data.table with additional columns "sample_id"
#'  and "clone_id" with unique ids for every row.
#' Values in the "sample_id" column will be the same as in the "sample_id"
#' metadata column. Positions of "count", "cdr3nt", "cdr3aa", "V gene", "J gene" clonoset's
#' columns and "file", "sampl_id" metadata's columns must be specified.
#'
#' TCRgrapherCounts is a subclass of TCRgrapher that store additional data:
#' count_table and feature_info table
#'
#' use clonoset(<your_object>) to see clonoset table
#' use metadata(<your_object>) to see metadata table
#' use edges(<your_object>) to see edges if they are present
#' use subset(<your_object>, <your sample_ids>) to take a subset
#'
#'
#' @examples
#' # # (1) Path to the one file without metadata
#' # file_path <- paste0(find.package('tcrgrapher'), "tests/testthat/testdata/clonosets_vdjtools_format.tsv")
#' # TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7) # positions of clonoset's columns
#' # # metadata has the same values in both columns. Let's change one of them
#' # metadata(TCRgrObject)[1,'sample_id'] <- 'sample_1'
#'
#' # # (2) Path to the directory with files without metadata
#' # dir_path <- paste0(find.package('tcrgrapher'), "tests/testthat/testdata/test_dir_with_clonosets")
#' # TCRgrObject <- TCRgrapher(dir_path, 1, 3, 4, 5, 7)
#'
#' # # (3) Path to the directory with files and metadata
#' # dir_path <- paste0(find.package('tcrgrapher'), "tests/testthat/testdata/test_dir_with_clonosets")
#' # metadata_path <- testthat::test_path("testdata/metadata.tsv")
#' # TCRgrObject <- TCRgrapher(dir_path, 1, 3, 4, 5, 7, # positions of clonoset's columns
#'                         #  metadata_path, 1, 2)      # positions of metadtata's columns
#'
#'
#' @param files_path Path to the file with a clonoset or path to the directory
#' with files. The file should contain a header and the following values: count,
#' cdr3nt, cdr3aa, V gene and J gene. The position of each column should be specified further.
#' Indexing is 1-based. Additional columns will be kept.
#' @param count_column Integer value that specifies the position of the count column in the clonoset.
#' @param cdr3nt_column Integer value that specifies the position of the cdr3nt column in the clonoset.
#' @param cdr3aa_column Integer value that specifies the position of the cdr3aa column in the clonoset.
#' @param v_gene_column Integer value that specifies the position of the V gene column in the clonoset.
#' @param j_gene_column Integer value that specifies the position of the J gene column in the clonoset.
#' @param metadata_path Path to the metadata file. The file must contain a header,
#' column with file names in the directory specified earlier, and a column with an unique
#' name for every file. The default value is NA. If the metadata path is not specified,
#' metadata will be generated.
#' @param files_column Integer value that specifies the position of the column with
#' file names in metadata
#' @param ids_column Integer value that specifies the position of the column with
#' unique names for every file in metadata
#' @seealso #TODO
#' @export
TCRgrapher <- function(files_path, count_column, cdr3nt_column,
                       cdr3aa_column, v_gene_column, j_gene_column,
                       metadata_path = NA, files_column = 1, ids_column = 2){

  if(is.na(metadata_path)){
    # without metadata
    files <- list.files(files_path)
    if(length(files) == 0){
      # only one file
      clonoset <- fread(files_path, header = TRUE)
      clonoset$sample_id <- files_path
      metadata <- data.table(file = files_path, sample_id = files_path)
    } else {
      # directory with files
      files_path <- check_path(files_path)
      clonoset <- rbindlist(
        lapply(
          files,
          function(file) fread(paste0(files_path, file), header = TRUE)[, sample_id := file, ]
        )
      )
      metadata <- data.table(file = files, sample_id = files)
    }
  } else {
    # with metadata
    metadata <- fread(metadata_path, header = TRUE)
    colnames(metadata)[c(files_column, ids_column)] <- c('file', 'sample_id')
    files_path <- check_path(files_path)
    metadata_files <- paste0(files_path, metadata$file)
    clonoset <- rbindlist(
      lapply(
        1:nrow(metadata),
        function(i) fread(metadata_files[i], header = TRUE)[,sample_id := metadata$sample_id[i],]
      )
    )
  }
  clonoset_indexes <-c(
    count_column, cdr3nt_column, cdr3aa_column, v_gene_column, j_gene_column
  )
  clonoset_names <- c('count', 'cdr3nt', 'cdr3aa', 'bestVGene', 'bestJGene')
  colnames(clonoset)[clonoset_indexes] <- clonoset_names
  clonoset[, clone_id := 1:.N, ]
  new('TCRgrapher', clonoset = clonoset, metadata = metadata)
}

#' @export
setGeneric("clonoset", function(x) standardGeneric("clonoset"))

#' @export
setMethod("clonoset", "TCRgrapher", function(x) x@clonoset)

#' @export
setGeneric("clonoset<-", function(x, value) standardGeneric("clonoset<-"))

#' @export
setMethod("clonoset<-", "TCRgrapher", function(x, value) {
  x@clonoset <- value
  validObject(x)
  x
})

#' @export
setGeneric("metadata", function(x) standardGeneric("metadata"))

#' @export
setMethod("metadata", "TCRgrapher", function(x) x@metadata)

#' @export
setGeneric("metadata<-", function(x, value) standardGeneric("metadata<-"))

#' @export
setMethod("metadata<-", "TCRgrapher", function(x, value) {
  x@metadata <- value
  validObject(x)
  x
})

#' @export
setGeneric("edges", function(x) standardGeneric("edges"))

#' @export
setMethod("edges", "TCRgrapher", function(x) x@edges)

#' @export
setGeneric("edges<-", function(x, value) standardGeneric("edges<-"))

#' @export
setMethod("edges<-", "TCRgrapher", function(x, value) {
  x@edges <- value
  validObject(x)
  x
})

#' @export
TCRgrValidity <- function(object){
  if(!is.data.table(object@clonoset)){
    return("Clonoset must be a data.table")
  }
  if(!is.data.table(object@metadata)){
    return("Metadata must be a data.table")
  }
  if(!all(c('count', 'cdr3nt', 'cdr3aa', 'bestVGene', 'bestJGene') %in% colnames(object@clonoset))){
    return("Clonoset must contain the following columns: 'count', 'cdr3nt', 'cdr3aa', 'bestVGene', 'bestJGene', 'sample_id'")
  }
  if(!all(c('file', 'sample_id') %in% colnames(object@metadata))){
    return("Metadata must contain the following columns: 'file', 'sample_id'")
  }
  if(anyDuplicated(object@metadata) != 0){
    return("Metadata must contain unique rows")
  }
  if(!all(unique(object@clonoset$sample_id) %in% object@metadata$sample_id)){
    return("'sample_id' clonoset column and 'sample_id' metadata column must have the same values")
  }
  counts <- object@clonoset$count
  if(!all(is.numeric(counts), counts > 0, counts %% 1 == 0)){
    return("Counts should be positive integers")
  }
  if(!is.character(object@clonoset$cdr3nt)){
    return("cdr3nt column must contain character values")
  }
  if(!is.character(object@clonoset$cdr3aa)){
    return("cdr3aa column must contain character values")
  }
  if(!is.character(object@clonoset$bestVGene)){
    return("bestVGene column must contain character values")
  }
  if(!is.character(object@clonoset$bestJGene)){
    return("bestJGene column must contain character values")
  }
  if(ncol(object@edges) == 2){
    if(!all(c(object@edges[,1], object@edges[,2]) %in% object@clonoset$clone_id)){
      return("Edges must connect nodes that are present in the clone_id column")
    }
  }
  return(TRUE)
}

#' @export
setValidity("TCRgrapher", TCRgrValidity)

#' @export
setMethod("show", "TCRgrapher", function(object) {
  if(ncol(object@edges) != 2){
    text_for_edges <- "  edges are not defined\n"
  } else {
    text_for_edges <- paste(" ", nrow(object@edges) ,"edges\n\n")
  }
  cat(is(object)[[1]], "\n\n",
      "  metadata: ", nrow(object@metadata), " rows and ",
      ncol(object@metadata), " columns\n",
      "  clonoset: ", nrow(object@clonoset), " rows and ",
      ncol(object@clonoset), " columns\n",
      text_for_edges,
      "use clonoset(<your_object>) to see clonoset table \n",
      "use metadata(<your_object>) to see metadata table \n",
      "use edges(<your_object>) to see edges if they are present \n",
      "use subset(<your_object>, <your sample_ids>) to take a subset\n",
      sep = ""
  )
})

#' #' @export
#' setGeneric("subset", function(x, samples) standardGeneric("subset"))

#' Subseting for TCRgrapher objects
#'
#' The function takes subset from clonotype table, metadata and edges if they are
#' present. Samples to keep should be specified. 'clone_id' column will be updated.
#'
#' @param x TCRgrapher object
#' @param subset Vector with sample ids that should be kept
#'
#' @export
setMethod("subset", "TCRgrapher", function(x, samples) {
  x@metadata <- x@metadata[sample_id %in% samples]
  x@clonoset <- x@clonoset[sample_id %in% samples]
  x@clonoset[, clone_id := 1:.N, ]
  if(ncol(x@edges == 2)){
    x@edges <- x@edges[x@edges[,1] %in% x@clonoset$clone_id & x@edges[,2] %in% x@clonoset$clone_id,]
  }
  x
})

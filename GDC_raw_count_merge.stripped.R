GDC_raw_count_merge <- function( id_list="my_id_list", my_rot="no", pseudo_fudge=NA, order_rows=TRUE,  order_columns=TRUE, debug=FALSE, verbose=FALSE, remove_tag=".htseq.counts")

{
  ### MAIN ###
  ###### load the neccessary packages
  if ( is.element("matlab", installed.packages()[,1]) == FALSE ){ install.packages("matlab") }
  library(matlab)

  if(debug==TRUE){print("Made it here 1")}

  my_ids <- flatten_list(as.list(scan(file=id_list, what="character")))

  if(debug==TRUE){print("Made it here 2")}

  # assuming all files are the same format, take a look at the number of files
  # and the first file to create a large enough data matrix
  compiled_matrix = matrix(,nrow = length(readLines(my_ids[1])),
                           ncol = (length(my_ids)))
  row_names <- unlist(as.list(scan(file=my_ids[1], what="character")))
  row_names <- row_names[is.na(suppressWarnings(as.numeric(row_names)))]

  # add column and row names (columns are files, rows are genes)
  rownames(compiled_matrix) <- row_names
  # remove tag from colnames
  col_names <- gsub(remove_tag, "", my_ids)
  colnames(compiled_matrix) <- col_names

  # assuming that every file has the same genes, reads through the files
  # and builds up the data matrix in "cimpiled_matrix"
  for (i in 1:length(my_ids))
  {
    if(verbose==TRUE){print(paste("Processing sample (", i, ")"))}
    read_lines <- unlist(as.list(scan(file=my_ids[i], what="character")))
    counts <- read_lines[!is.na(suppressWarnings(as.numeric(read_lines)))]
    compiled_matrix[,i] = counts
  }

  if(debug==TRUE){print("Made it here 3")}

  # rotate the matrix if that option is selected
  if( identical(my_rot, "yes")==TRUE ){
    compiled_matrix <- rot90(rot90(rot90(compiled_matrix)))
  }

  if(debug==TRUE){print("Made it here 4")}

  # output the matrix as a flat file
  fileout_name <- paste0(id_list, ".merged_data.txt")
  export_data(compiled_matrix, fileout_name)
}



### SUBS ###
export_data <- function(data_object, file_name){
  write.table(data_object, file=file_name, sep="\t", col.names = NA, row.names = TRUE, quote = FALSE, eol="\n")
}



flatten_list <- function(some_list){
  flat_list <- unlist(some_list)
  flat_list <- gsub("\r","",flat_list)
  flat_list <- gsub("\n","",flat_list)
  flat_list <- gsub("\t","",flat_list)
}

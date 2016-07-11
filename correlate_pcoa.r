source("~/git/Ronald-and-Mert/render_calculated_pcoa.r")

correlate_pcoa <- function(metadata_file,pcoa_file){
  #in: metadata file, pcoa file
  #out: correlation matrix file
  import_metadata_v2 <- function(group_table){
    #import metadata and changes every field to numeric values
    meta <- read.table(
      file=group_table,row.names=1,header=TRUE,sep="\t",
      comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE
    )
    col_names <- names(meta)
    meta[,col_names] <- lapply(meta[,col_names],as.factor)
    meta[,col_names] <- lapply(meta[,col_names],as.numeric)
    print(ncol(meta))
    print(nrow(meta))
    return(meta)
  }
  
  my_metadata <- import_metadata_v2(metadata_file)
  my_pcoa <- load_pcoa_data(pcoa_file)
  #print(my_pcoa$eigen_vectors)
  #print(nrow(my_pcoa$eigen_values))
  cor_summary_matrix <- matrix(NA,ncol(my_metadata),ncol(my_pcoa$eigen_vectors))
  for(i in 1:ncol(my_metadata)){
    for(j in 1:ncol(my_pcoa$eigen_vectors)){
      cor_summary_matrix[i,j] <- cor(my_metadata[,i],my_pcoa$eigen_vectors[,j],use="pairwise.complete.obs")
    }
  }
  write.table(cor_summary_matrix,file = "cor_summary.txt", sep = "\t",quote = F,row.names = F,col.names = F)
}
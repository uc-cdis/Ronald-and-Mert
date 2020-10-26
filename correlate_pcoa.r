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
    #fix column types
    meta[,col_names] <- lapply(meta[,col_names],as.factor)
    meta[,col_names] <- lapply(meta[,col_names],as.numeric)
    return(meta)
  }

  my_metadata <- import_metadata_v2(metadata_file)
  my_pcoa <- load_pcoa_data(pcoa_file)
  #create empty matrices with correct sizes
  cor_summary_matrix_pearson <- matrix(NA,ncol(my_metadata),ncol(my_pcoa$eigen_vectors))
  cor_summary_matrix_spearman <- matrix(NA,ncol(my_metadata),ncol(my_pcoa$eigen_vectors))
  cor_summary_matrix_pearson_test<- matrix(NA,ncol(my_metadata),ncol(my_pcoa$eigen_vectors))
  cor_summary_matrix_spearman_test<- matrix(NA,ncol(my_metadata),ncol(my_pcoa$eigen_vectors))

  for(i in 1:ncol(my_metadata)){
    for(j in 1:ncol(my_pcoa$eigen_vectors)){
      #fill the cor matricies
      cor_summary_matrix_pearson[i,j] <- cor(my_metadata[,i],my_pcoa$eigen_vectors[,j],use="pairwise.complete.obs",method = "pearson")
      cor_summary_matrix_spearman[i,j] <- cor(my_metadata[,i],my_pcoa$eigen_vectors[,j],use="pairwise.complete.obs",method = "spearman")
      if(!is.na(my_metadata[,i])){
        #fill the p-value matricies
        cor_summary_matrix_spearman_test[i,j] <- (cor.test(my_metadata[,i],my_pcoa$eigen_vectors[,j],method = "spearman")$p.value)
        cor_summary_matrix_pearson_test[i,j] <- (cor.test(my_metadata[,i],my_pcoa$eigen_vectors[,j],method = "pearson")$p.value)
        }
    }
  }
  #write all of the matricies
  write.table(cor_summary_matrix_pearson,file = "cor_summary_pearson.txt", sep = "\t",quote = F,row.names = F,col.names = F)
  write.table(cor_summary_matrix_spearman,file = "cor_summary_spearman.txt", sep = "\t",quote = F,row.names = F,col.names = F)
  write.table(cor_summary_matrix_spearman_test,file = "cor_summary_spearman_test.txt", sep = "\t",quote = F,row.names = F,col.names = F)
  write.table(cor_summary_matrix_pearson_test,file = "cor_summary_pearson_test.txt", sep = "\t",quote = F,row.names = F,col.names = F)
}

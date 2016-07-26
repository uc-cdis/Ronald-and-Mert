library(biomaRt)
library(dplyr)

get_y_freq_from_gene_df <- function(gene_df){
  mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  gene_df <- read.table(file = gene_df,sep = "\t",stringsAsFactors = F,fill = T)
  result <- NULL
  colnames(gene_df) <- gene_df[1,]
  y_list <- c()
  for (col in gene_df){
    ensembl_genes <- gsub( "\\..*", "", col[-1])

    print(col[1])
    result <- getBM(
      filters= "ensembl_gene_id", 
      attributes= c("ensembl_gene_id","chromosome_name"),
      values= ensembl_genes,
      mart= mart)
    y_genes <- result[which(result$chromosome_name=="Y"),]

    if (length(y_genes$ensembl_gene_id)!=0){
    y_list <- c(y_list,paste(as.character(unlist(y_genes[1])),col[1],sep = " "))
    }
    check <- gsub( "\\..*", "", col[-1]) %in% y_genes[,1]
  }
  y_gene_mat <- unlist(y_list)
  
  gene <- sapply(y_gene_mat,function (l) strsplit(l,split = " ")[[1]][1])
  project_name <- sapply(y_gene_mat,function (l) strsplit(l,split = " ")[[1:2]][1])
  y_gene_df <- cbind(gene,project_name)

  final <- aggregate(project_name~gene,data = y_gene_df,FUN = paste,collapse = ",")
  write.table(final,file = "y_gene_frequencies.txt",sep = "\t",row.names = F)
}
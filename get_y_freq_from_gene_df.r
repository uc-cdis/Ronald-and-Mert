library(biomaRt)
library(dplyr)

#creates a frequency table of y genes by project name
#input: data frame with project names as column names and gene ID's under corresponding project names
#output: frequency table showing which y gene occured in which projects

get_y_freq_from_gene_df <- function(gene_df){
  mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  gene_df <- read.table(file = gene_df,
                        sep = "\t",
                        stringsAsFactors = F,
                        fill = T)
  result <- NULL
  colnames(gene_df) <- gene_df[1,]
  y_list <- c()
  for (col in gene_df){
    #do it for every project
    ensembl_genes <- gsub( "\\..*", "", col[-1])

    print(col[1])
    result <- getBM(
      filters= "ensembl_gene_id", 
      attributes= c("ensembl_gene_id","chromosome_name"),
      values= ensembl_genes,
      mart= mart)
    y_genes <- result[which(result$chromosome_name=="Y"),]

    if (length(y_genes$ensembl_gene_id)!=0){
    #add the found y genes to the y_list var
    y_list <- c(y_list,paste(as.character(unlist(y_genes[1])),col[1],sep = " "))
    }
    check <- gsub( "\\..*", "", col[-1]) %in% y_genes[,1]
  }
  #unlist y_list
  y_gene_mat <- unlist(y_list)
  #create two components (gene & project name) from y_gene_mat vector two use to create a structured data frame
  gene <- sapply(y_gene_mat,function (l) strsplit(l,split = " ")[[1]][1])
  project_name <- sapply(y_gene_mat,function (l) strsplit(l,split = " ")[[1:2]][1])
  #create df
  y_gene_df <- cbind(gene,project_name)
  #aggregate (create the frequency table of ) y_gene_df
  final <- aggregate(project_name~gene,data = y_gene_df,FUN = paste,collapse = ",")
  write.table(final,file = "y_gene_frequencies.txt",sep = "\t",row.names = F)
}

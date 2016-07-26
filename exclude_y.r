#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_data <- read.table("counts_files.merged_data.txt",sep = "\t",header = F,stringsAsFactors = F)
ensembl_genes<- gsub( "\\..*", "", gene_data[,1] )

result <- getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id","chromosome_name"),
  values= ensembl_genes,
  mart= mart)

y_genes <- result[which(result$chromosome_name=="Y"),]
check <- gsub( "\\..*", "", gene_data[,1]) %in% y_genes[,1]
gene_data_without_y_genes <- gene_data[!check,]
write.table(gene_data_without_y_genes,"counts_files.merged_data.without_y.txt",sep = "\t",col.names = F,row.names = F,quote = F)
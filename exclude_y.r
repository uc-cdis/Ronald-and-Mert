#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)

mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_data <- read.table("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Mann-Whitney-unpaired-Wilcoxon.hits.demographic.gender.1.STATS_RESULTS.txt",sep = "\t",header = F,stringsAsFactors = F)
#take the first 10 percent
gene_data <- gene_data[1:(nrow(gene_data)*0.1),]
#apply fdr cut-off
gene_data <- gene_data[which(as.numeric(gene_data[,ncol(gene_data)])<0.000001),]
ensembl_genes<- gsub( "\\..*", "", gene_data[,1] )
#map gene ID to chr
result <- getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id","chromosome_name"),
  values= ensembl_genes,
  mart= mart)
#remove the y chr genes
y_genes <- result[which(result$chromosome_name=="Y"),]
check <- gsub( "\\..*", "", gene_data[,1]) %in% y_genes[,1]
gene_data_without_y_genes <- gene_data[!check,]
write.table(gene_data_without_y_genes,"counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Mann-Whitney-unpaired-Wilcoxon.hits.demographic.gender.1.STATS_RESULTS.10_percent.without_y.txt",sep = "\t",col.names = F,row.names = F,quote = F)

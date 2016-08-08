FilterStats <- function(filename, percent = .01, exclude.y = F, fdr.cutoff = F) {
  gene_data <- read.table(filename,
                          sep = "\t",
                          header = T,
                          stringsAsFactors = F,
                          check.names = F)
  if(percent) {
    #take the first xx percent
    gene_data <- gene_data[1:floor(nrow(gene_data) * percent),]
  }
  if(fdr.cutoff) {
    #apply fdr cut-off
    gene_data <- gene_data[which(as.numeric(gene_data[,ncol(gene_data)]) < fdr.cutoff),]
  }
  if(exclude.y) {
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
    gene_data <- gene_data[!check,]
  }
  write.table(gene_data,
              paste0(file_path_sans_ext(filename), "filtered", file_ext(filename)),
              sep = "\t",
              col.names = T,
              row.names = F,
              quote = F)
}

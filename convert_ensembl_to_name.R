library(biomaRt)

ConvertEnsemblToName <- function(ensembl.genes) {
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  gene.names <- getBM(
    filters= "ensembl_gene_id",
    attributes= c("ensembl_gene_id", "external_gene_name"),
    values= ensembl.genes,
    mart= mart)
  return(gene.names)
}

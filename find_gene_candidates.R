FindGeneCandidates <- function(foldername = ".") {
  all.files <- list.files(foldername)
  genes.list <- list()
  counts.vec <- c()
  for (file in all.files) {
    if (!grepl("significant", file)) {
      next
    }
    file.lines <- readLines(file)
    top.cancer.project <- gsub("\"", "", strsplit(file.lines[3], "\t")[[1]][1])
    top.cancer.count <- strsplit(file.lines[3], "\t")[[1]][2]
    gene.ensemblid <- strsplit(file.lines[1], "\"")[[1]][2]
    genes.list[[gene.ensemblid]] <- c(gene.ensemblid, top.cancer.project, top.cancer.count)
    counts.vec <- c(counts.vec, top.cancer.count)
    if(length(counts.vec) == 578) {
      browser()
    }
  }
  ordered.genes.list <- genes.list[order(counts.vec, decreasing = T)]
  close( file("sorted_genes_list.txt", open="w" ) )
  lapply(ordered.genes.list, write, "sorted_genes_list.txt", append=TRUE, ncolumns=1000)
}

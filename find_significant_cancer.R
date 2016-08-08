library(dunn.test)

FindSignificantCancer <- function(data.filename, metadata.filename, 
                                  gene, metadata.col, 
                                  alpha.level = .001, convert = F) {
  data.table <- read.table(file=data.filename,
                           row.names=1,
                           header=TRUE,
                           sep="\t",
                           colClasses = "character", 
                           check.names=FALSE,
                           comment.char = "",
                           quote="", 
                           fill=TRUE, 
                           blank.lines.skip=FALSE)
  metadata.table <- read.table(file=metadata.filename,
                               row.names=1,
                               header=TRUE,
                               sep="\t",
                               colClasses = "character", 
                               check.names=FALSE,
                               comment.char = "",
                               quote="", 
                               fill=TRUE, 
                               blank.lines.skip=FALSE)
  if (gene == "ALL") {
    gene <- (1:nrow(data.table))
  }
  ensembl.ids <- c()
  for (g in gene) {
    if (suppressWarnings(is.na(as.numeric(g)))) {
      ensembl.ids <- g
    } else {
      ensembl.ids <- rownames(data.table)[g]
    }
  }
  if (convert) {
    source("convert_ensembl_to_name.R")
    gene.name.table <- ConvertEnsemblToName(ensembl.ids)
  }
  for (g in ensembl.ids) {
    dunn.test.results <- dunn.test(as.numeric(data.table[g,]), 
                                   metadata.table[[metadata.col]],
                                   method = "bonferroni",
                                   table = F)
    p.values.ordered <- dunn.test.results$P.adjusted
    names(p.values.ordered) <- dunn.test.results$comparisons
    p.values.ordered <- p.values.ordered[order(p.values.ordered)]
    filtered.p.values <- p.values.ordered[p.values.ordered < alpha.level]
    names.vec <- c()
    for (name in names(lowest.p.values)) {
      names.vec <- c(names.vec, unlist(strsplit(name, " - ")))
    }
    names.count.table <- as.data.frame(table(names.vec))
    names.count.table <- names.count.table[order(names.count.table[[2]], decreasing = T), ]
    if(convert) {
      write(gene.name.table, paste(gene.name, "significant.cancers", sep = "."), sep = "\t")
    } else {
      write(g, paste(gene.name, "significant.cancers", sep = "."))
    }
    write.table(names.count.table, 
                file = paste(gene.name, "significant.cancers", sep = "."),
                sep = "\t",
                row.names = F,
                append = T)
  }
}

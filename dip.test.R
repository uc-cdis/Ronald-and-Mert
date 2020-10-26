source("~/git/Ronald-and-Mert/convert_ensembl_to_name.R")
mytable <- read.table(file="counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Kruskal-Wallis.hits.project.project_id.1.STATS_REMOVED.filtered.byMedian.txt",row.names=1,header=TRUE,sep="\t",colClasses = "character", check.names=FALSE,comment.char = "", quote="", fill=TRUE, blank.lines.skip=FALSE)
p.vector <- c()
for (row in rownames(mytable)) {
  row <- strsplit(row, "\\.")[[1]][1]
  result <- dip.test(as.numeric(mytable[row,]))
  p.vector <- c(p.vector, result$p.value)
  names(p.vector)[length(p.vector)] <- row
}
ordered.p.vector <- sort(p.vector)
unordered.gene.ids <- ConvertEnsemblToName(names(p.vector))
row.names(unordered.gene.ids) <- unordered.gene.ids[[1]]
p.table <- unordered.gene.ids[names(ordered.p.vector),]
p.table <- cbind(p.table, ordered.p.vector)
write.table(p.table, file = "dip.test.sorted", col.names = F, row.names = F, sep = "\t")

CompareCounts <- function(counts.file.name, reference.file.name, tissue.col, project) {
  counts <- read.table(counts.file.name, row.names=1,
                            header=TRUE, sep="\t", colClasses = "character", 
                            check.names=FALSE,comment.char = "", 
                            fill=TRUE, blank.lines.skip=FALSE)
  reference <- read.table(reference.file.name, row.names=1,
                            header=TRUE, sep="\t", colClasses = "character", 
                            check.names=FALSE,comment.char = "", 
                            fill=TRUE, blank.lines.skip=FALSE)
  rownames(counts) <- unlist(strsplit(rownames(counts), "\\."))[(1:nrow(counts))*2-1]
  result.frame <- data.frame(matrix(NA, nrow = nrow(reference), 3))
  rownames(result.frame) <- rownames(reference)
  colnames(result.frame) <- c("reference", "estimate", "p-value")
  for (r in rownames(reference)) {
    if (!(r %in% rownames(counts))) {
      next
    }
    result <- t.test(as.numeric(counts[r,]), mu = as.numeric(reference[r, tissue.col]))
    result.frame[r, 3] <- result$p.value
    result.frame[r, 2] <- result$estimate
    result.frame[r, 1] <- reference[r, tissue.col]
  }
  ordered.result.frame <- result.frame[order(result.frame[["p-value"]], 
                                             -abs(as.numeric(result.frame[["reference"]]) - as.numeric(result.frame[["estimate"]]))), ]
  write.table(ordered.result.frame, file = paste(project, "different_genes.txt", sep = "."), 
              sep = "\t", quote = F, col.names = NA, row.names = T)
}

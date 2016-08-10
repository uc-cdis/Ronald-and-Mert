library(ggplot2)

FPKMVsCountsPlot <- function(fpkm.merged.filename, counts.merged.filename) {
  fpkm.table <- read.table(file = fpkm.merged.filename,
                           row.names = 1,header = TRUE,sep = "\t",
                           colClasses = "character", check.names = FALSE,comment.char = "", 
                           fill = TRUE, blank.lines.skip = FALSE)
  counts.table <- read.table(file = counts.merged.filename,
                             row.names = 1,header = TRUE,sep = "\t",
                             colClasses = "character", check.names = FALSE,comment.char = "", 
                             fill = TRUE, blank.lines.skip = FALSE)
  counts <- c()
  fpkm <- c()
  for (r in row.names(fpkm.table)[1:1000]) {
    counts <- c(counts, as.numeric(counts.table[r, ]))
    fpkm <- c(fpkm, as.numeric(fpkm.table[r, ]))
  }
  plot <- qplot(x=counts, y = fpkm)
  ggsave(plot, file = "fpkm_vs_counts.jpg")
}

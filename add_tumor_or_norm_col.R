AddTumorOrNormCol <- function(metadata.file = "", normal.filenames = "", tumor.filenames = "") {
  metadata.table <- read.table(file = "metadata.file",
      row.names = 1,header = TRUE,sep = "\t",
      colClasses = "character", check.names = FALSE,comment.char = "",
      fill = TRUE, blank.lines.skip = FALSE)
  normal.files <- readLines(normal.filenames)
  tumor.files <- readLines(tumor.filenames)
  normal.ids <- c()
  tumor.ids <- c()
  for (x in normal.files) {
    normal.ids <- c(normal.ids, strsplit(x, "\\.")[[1]][1])
  }
  for (x in tumor.files) {
    tumor.ids <- c(tumor.ids, strsplit(x, "\\.")[[1]][1])
  }
  normal.or.tumor.col <- c()
  for (x in 1:nrow(metadata.table)) {
    if (rownames(metadata.table)[x] %in% normal.ids) {
      normal.or.tumor.col <- c(normal.or.tumor.col, "normal")
    } else if (rownames(metadata.table)[x] %in% tumor.ids) {
      normal.or.tumor.col <- c(normal.or.tumor.col, "tumor")
    } else {
      normal.or.tumor.col <- c(normal.or.tumor.col, "ERROR")
    }
  }
  metadata.table <- cbind(metadata.table, normal.or.tumor.col)
  write.table(metadata.table, file = paste0(metadata.file, ".sampleidcol"),
      sep = "\t", quote = F, col.names = NA, row.names = T)
}

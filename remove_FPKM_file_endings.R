library(tools)

RemoveFPKMColEndings <- function(merged.filename) {
  merged.table <- read.table(file = merged.filename,
                             row.names = 1,header = TRUE,sep = "\t",
                             colClasses = "character", check.names = FALSE,comment.char = "",
                             fill = TRUE, blank.lines.skip = FALSE)
  colnames(merged.table) <- file_path_sans_ext(file_path_sans_ext(colnames(merged.table)))
  write.table(merged.table, file = merged.filename,
              sep = "\t", quote = F, col.names = NA, row.names = T)
}

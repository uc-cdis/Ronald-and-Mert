library(tools)
AppendMeanAndMedian <- function(filename) {
  mytable <- read.table(file=filename,
                        row.names=1,
                        header=TRUE,
                        sep="\t",
                        colClasses = "character",
                        check.names=FALSE,
                        comment.char = "",
                        quote="",
                        fill=TRUE,
                        blank.lines.skip=FALSE,
                        stringsAsFactors = F)
  newtable <- cbind(mytable, 
                    rowMeans(data.matrix(mytable), na.rm = T), 
                    apply(mytable, 1, median, na.rm = T), 
                    stringsAsFactors = F)
  newtable <- cbind(newtable, as.numeric(newtable[[ncol(newtable) - 1]]) - as.numeric(newtable[[ncol(newtable)]]), 
                    stringsAsFactors = F)
  colnames(newtable) <- c(colnames(newtable)[1:(ncol(newtable)-3)], "mean", "median", "mean-median")
  byMedian <- newtable[order(newtable[["median"]], decreasing = T), ]
  byMean <- newtable[order(newtable[["mean"]], decreasing = T), ]
  byDifference <- newtable[order(abs(newtable["mean-median"]), decreasing = T), ]

  write.table(byMedian,
              file = paste(file_path_sans_ext(filename), "byMedian", file_ext(filename), sep = "."),
              sep = "\t",
              quote = F,
              col.names = NA,
              row.names = T)
  write.table(byMean,
              file = paste(file_path_sans_ext(filename), "byMean", file_ext(filename), sep = "."),
              sep = "\t",
              quote = F,
              col.names = NA,
              row.names = T)
  write.table(byDifference,
              file = paste(file_path_sans_ext(filename), "byMeanMedianDiff", file_ext(filename), sep = "."),
              sep = "\t",
              quote = F,
              col.names = NA,
              row.names = T)
}

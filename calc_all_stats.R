source("~/git/Ronald-and-Mert/calc_stats.r")

calc_all_stats <- function(yourdir = ".", metadata_column) {
  for(folder in list.dirs(path=yourdir, recursive = F)) {
    countsFile <- file.path(folder, "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
    metadataFilename <- file.path(folder, list.files(path = folder,
                                                     pattern = ".*METADATA.*"))
    calc_stats(data_table = countsFile, 
               metadata_table = metadataFilename,
               metadata_column = metadata_column,
               stat_test = "Mann-Whitney-unpaired-Wilcoxon")
  }
}

extract_top_percent <- function(yourfile, percent = 10) {
  yourtable <- read.table(file = orderedFile, header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
  subcolNum <- ncol(yourtable) * (percent/100)
  topGenes <- colnames(yourtable)[1:subcolNum]
}

create_top_percent_table <- function(yourdir, percent = 10, metadataCol = hits.demographic.gender.1) {
  projList <- list.dirs(path=yourdir, recursive = F)
  numProj <- length(projList)
  first <- T
  topTable <- matrix()
  for(folder in projList) {
    countsFile <- file.path(folder, paste0("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Mann-Whitney-unpaired-Wilcoxon.", 
                                           metadataCol, 
                                           ".STATS_RESULTS.txt"))
    if(!file.exists(countsFile)) {
      calc_all_stats(yourdir, metadataCol)
    }
    top_percent <- extract_top_percent(countsFile, percent = percent)
    if(first) {
      topTable <- matrix(NA, nrow = numProj, ncol = length(top_percent))
      first <- F
    }
    topTable[,folder] <- top_percent
  }
  namesVec <- c()
  for(proj in projList) {
    append(namesVec, substr(proj, 3, len(proj)))
  }
  colnames(topTable) <- namesVec
}
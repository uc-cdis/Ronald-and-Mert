source("~/git/Ronald-and-Mert/calc_stats.r")

calc_all_stats <- function(yourdir = ".", metadata_column, start = 1) {
  projlist <- list.dirs(path = yourdir, recursive = F)
  for(folder in projlist[start:length(projlist)]) {
    print(paste0("calculating stats for ", folder))
	if(file.exists(file.path(folder, paste0("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Mann-Whitney-unpaired-Wilcoxon.", metadata_column, ".STATS_RESULTS.txt")))) {
      next
    }
    countsFile <- file.path(folder, "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
    metadataFilename <- file.path(folder, list.files(path = folder,
                                                     pattern = ".*METADATA.*"))
    tryCatch({
      calc_stats(data_table = countsFile, 
                 metadata_table = metadataFilename,
                 metadata_column = metadata_column,
                 stat_test = "Mann-Whitney-unpaired-Wilcoxon")
    }, error = function(e){})
  }
}

extract_top_percent <- function(yourfile, percent = 10) {
  yourtable <- read.table(file = yourfile, header = T, sep = "\t", row.names = 1, stringsAsFactors = F, fill = T)
  subcolNum <- floor(nrow(yourtable) * (percent/100))
  topGenes <- rownames(yourtable)[1:subcolNum]
}

create_top_percent_table <- function(yourdir = '.', percent = 10, metadataCol = "hits.demographic.gender.1") {
  projList <- list.dirs(path=yourdir, recursive = F)
  numProj <- length(projList)
  first <- T
  topTable <- matrix()
  for(folder in projList) {
    print(paste("creating table for", folder))
    countsFile <- file.path(folder, paste0("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Mann-Whitney-unpaired-Wilcoxon.", 
                                           metadataCol, 
                                           ".STATS_RESULTS.txt")) # name of file change here
    top_percent <- extract_top_percent(countsFile, percent = percent)
    if(first) {
      topTable <- matrix(NA, ncol = numProj, nrow = length(top_percent) + 1000)
      colnames(topTable) <- projList
      first <- F
    }
    topTable[1:length(top_percent),folder] <- top_percent
  }
  namesVec <- c()
  for(proj in projList) {
    append(namesVec, substr(proj, 3, len(proj)))
  }
  colnames(topTable) <- namesVec
  topTable[rowSums(is.na(topTable)) != ncol(topTable),]
  write.table(topTable, file = paste0("top", percent, "percent.", metadataCol, ".tsv"), sep = "\t")
}

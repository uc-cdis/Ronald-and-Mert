source("~/git/Ronald-and-Mert/calc_stats.r")

calc_all_stats <- function(yourdir = ".", metadata_column, start = 1) {
    # iterates through all directories in your dir and calculates stats
    # using Mann-Whitney-unpaired-Wilcoxon test
    #
    # Args:
    #   yourdir: The directory that is to be iterated through and contains
    #     all the individual directories for each project.
    #   metadata_column: Column of metadata with which to do the analysis.
    #     Can be number or string.
    #   start: Starting project index (used to run many analyses concurrently)
    #
    # Returns:
    #   (None)
    #   produces a file that ends with STATS_RESULTS.txt
    projlist <- list.dirs(path = yourdir, recursive = F)
  for(folder in projlist[start:length(projlist)]) {
        print(paste0("calculating stats for ", folder))
      if(file.exists(file.path(folder, paste0("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Mann-Whitney-unpaired-Wilcoxon.", 
                                                                                          metadata_column, 
                                                                                                                                      ".STATS_RESULTS.txt")))) {
              # if the STATS_RESULTS.txt file exists already, skip this folder
              next
          }
          countsFile <- file.path(folder, "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
          metadataFilename <- file.path(folder, list.files(path = folder,
                                                                                                                pattern = ".*METADATA.*"))
              tryCatch({
                      # tryCatch to stop everything from stopping if running into an error
                      # should probably add error logging at a later time
                      calc_stats(data_table = countsFile, 
                                                  metadata_table = metadataFilename,
                                                                   metadata_column = metadata_column,
                                                                                    stat_test = "Mann-Whitney-unpaired-Wilcoxon")
                            }, error = function(e){})
            }
}

create_top_percent_table <- function(yourdir = '.') {
    # Exports a summary tsv for every project using gene expression files 
    #
    # Args:
    #   yourdir: directory that contains all subdirectories of each individual
    #            project. 
    #
    # Returns:
    #   (None)
    #   Exports table with project ids as column names and the genes that
    #   are associated with said project in its column. 
    projList <- list.dirs(path=yourdir, recursive = F)
  numProj <- length(projList)
    first <- T
    topTable <- matrix()
      for(folder in projList) {
            # CHANGE FILE NAME HERE FOR AFTER REMOVING Y CHROMOSOME GENES
            countsFile <- file.path(folder, paste0("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Mann-Whitney-unpaired-Wilcoxon.hits.demographic.gender.1.STATS_RESULTS.edit.txt"))
        if(!file.exists(countsFile)) {
                calc_all_stats(yourdir, metadataCol)
            }
            gene.names = row.names(read.table(file = countsFile, header = T, sep = '\t', row.names = 1, stringsAsFactors = T))
            
            if(first) {
                    # create large matrix so no need to cbind for efficiency 
                    #   (not actually sure if this helps at all)
                    # hard coded in 1000 because I don't know a better way to correct for variable
                    # number of genes
                    topTable <- matrix(NA, ncol = numProj, nrow = length(gene.names) + 1000)
                  colnames(topTable) <- projList
                        first <- F
                      }
                topTable[1:length(gene.names),folder] <- gene.names
              }
      # add colnames to topTable
      namesVec <- c()
        for(proj in projList) {
              browser()
          namesVec <- c(namesVec, substr(proj, 3, nchar(proj)))
            }
        colnames(topTable) <- namesVec
          topTable <- topTable[rowSums(is.na(topTable)) != ncol(topTable),]
          write.table(topTable, file = "summary_table_for_important_genes.tsv", sep = "\t", row.names = F)
}


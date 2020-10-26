extract_top3_comp <- function(orderedFile, metadataCol)
{
  orderedTable <- read.table(file = orderedFile, header = T, sep = "\t", row.names = NULL, stringsAsFactors = F)
  firstComp <- orderedTable[1, metadataCol]
  secondComp <- orderedTable[2, metadataCol]
  thirdComp <- orderedTable[3, metadataCol]
  firstComp <- sep_compRP(firstComp)
  secondComp <- sep_compRP(secondComp)
  thirdComp <- sep_compRP(thirdComp)
  return(c(firstComp, secondComp, thirdComp))
}

sep_compRP <- function(compRP)
{
  vecCompRP <- strsplit(compRP, split = ':')[[1]]
  vecCompRP[1] <- as.numeric(strsplit(vecCompRP[1], split = '\\.')[[1]][2])
}

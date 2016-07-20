order_corr_matrix <- function(corr_matrix,p_matrix)
{
  corr <<- read.table(file = corr_matrix, header = T, sep = "\t",row.names = 1)
  # replace all NA with 0 so absolute value function works
  p <<- read.table(file = p_matrix, header = T, sep = "\t", row.names = 1)
  for (y in 1:length(corr)){
    for (x in 1:length(corr[,1])){
      corr [x,y] <- paste(corr[x,y],p[x,y],sep = ":")
    }
  }
  corr <<- corr
  
  a <<- apply(corr, 1, function(row) {
    #row <- gsub( " .*$", "", row )
    rowWithName <- paste(names(row), row, sep=":")
    #print(row)
    rowWithName[order(abs(as.numeric(gsub( ":.*$", "", row ))), decreasing = T)]
    })
  
  
  write.table(a, file = paste0(corr_matrix, ".ordered"), sep = "\t", row.names = F)
}

order_all_matrices <- function()
{
	for(dir in list.dirs())
	{
		if(file.exists(file.path(dir, "cor_summary.txt")))
		{
			print(dir)
			order_corr_matrix(file.path(dir,"cor_summary_pearson.txt"),p_matrix = file.path(dir,"cor_summary_pearson_test.txt"))
			order_corr_matrix(file.path(dir,"cor_summary_spearman.txt"),p_matrix = file.path(dir,"cor_summary_spearman_test.txt"))
		}
	}
}
# source("order_corr_matrix.r")
# order_corr_matrix("TCGA-CHOL/cor_summary.txt")

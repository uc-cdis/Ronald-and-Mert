order_corr_matrix <- function(corr_matrix)
{
	corr <- read.table(file = corr_matrix, header = T, sep = "\t", row.names = 1)
	# replace all NA with 0 so absolute value function works

	a <- apply(corr, 1, function(row) {
			   rowWithName <- paste(names(row), row, sep=":")
			   rowWithName[order(abs(as.numeric(row)), decreasing = T)]
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
			order_corr_matrix(file.path(dir,"cor_summary_pearson.txt"))
			order_corr_matrix(file.path(dir,"cor_summary_spearman.txt"))
		}
	}
}
# source("order_corr_matrix.r")
# order_corr_matrix("TCGA-CHOL/cor_summary.txt")

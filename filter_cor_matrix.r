filter_cor_matrix <- function(corr_matrix, p_matrix)
{
	corr <- read.table(file = corr_matrix, header = T, sep = "\t", row.names = 1)
	p <- read.table(file = p_matrix, header = T, sep = "\t", row.names = 1)
	corr.list <- setNames(split(corr, seq(nrow(corr))), rownames(corr))
	p.list <- setNames(split(p, seq(nrow(p))), rownames(p))

	b <- mapply(function(fromc, fromp) {
				fromc[is.na(fromc)] <- 0
				fromp[is.na(fromp)] <- 1
				combinedVec <- paste(names(fromc), fromc, fromp, sep=":")
			   	pos_c <- abs(as.numeric(fromc))
				TFvec <- (pos_c > 0.7) & (fromp < 0.0001)
				combinedVecFiltered <- combinedVec[TFvec]
				return(combinedVecFiltered[order(pos_c[TFvec])])
				},corr.list, p.list)

	if(file.exists(paste0(corr_matrix, ".filtered")))
	{
		file.remove(paste0(corr_matrix, ".filtered"))
		file.create(paste0(corr_matrix, ".filtered"))
	}
	namesList <- names(b)
	index <- 1
	lapply(b, function(elem) {
		   write(append(elem, namesList[index], 0), file = paste0(corr_matrix, ".filtered"), sep = "\t", ncolumns=1000, append=T)
		   index <<- index + 1
			})

}

filter_all_matrices <- function()
{
	for(dir in list.dirs())
	{
		if(file.exists(file.path(dir, "cor_summary.txt")))
		{
			print(dir)
			filter_cor_matrix(file.path(dir,"cor_summary_pearson.txt"), file.path(dir,"cor_summary_pearson_test.txt"))
			filter_cor_matrix(file.path(dir,"cor_summary_spearman.txt"), file.path(dir,"cor_summary_spearman_test.txt"))
		}
	}
}

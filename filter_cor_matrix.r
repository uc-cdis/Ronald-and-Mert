filter_cor_matrix <- function(corr_matrix, p_matrix,r_thresh = 0.7,p_thresh = 0.0001,pco_thresh = 15,file_name =".filtered")
{
	corr <- read.table(file = corr_matrix, header = T, sep = "\t", row.names = 1)
	p <- read.table(file = p_matrix, header = T, sep = "\t", row.names = 1)
	if (pco_thresh < length(corr[1,])){
		corr <- corr[,1:pco_thresh]
		p <- p[,1:pco_thresh]
	}
	corr.list <- setNames(split(corr, seq(nrow(corr))), rownames(corr))
	p.list <- setNames(split(p, seq(nrow(p))), rownames(p))

	b <<- mapply(function(fromc, fromp) {
				 fromc[is.na(fromc)] <- 0
				 fromp[is.na(fromp)] <- 1
				 #print(fromc)
				 combinedVec <- paste(names(fromc), fromc, fromp, sep=":")
				 pos_c <- abs(as.numeric(fromc))
				 TFvec <- (pos_c > r_thresh) & (fromp < p_thresh)
				 combinedVecFiltered <- combinedVec[TFvec]
				 return(combinedVecFiltered[order(pos_c[TFvec])])
},corr.list, p.list)
	ls <- sapply(b,length)
	l <- length(ls[(ls!=0)])
	if(file.exists(paste0(corr_matrix, file_name,".pass_num:",l)))
	{
		file.remove(paste0(corr_matrix, file_name,".pass_num:",l))
		file.create(paste0(corr_matrix, file_name,".pass_num:",l))
	}
	namesList <- names(b)
	index <- 1
	lapply(b, function(elem) {
		   write(append(elem, namesList[index], 0), file = paste0(corr_matrix, file_name,".pass_num:",l), sep = "\t", ncolumns=1000, append=T)
		   index <<- index + 1
})

}

filter_all_matrices <- function(pco_test = c(15),p_test = c(0.0001),r_test = c(0.7))
{
	for(dir in list.dirs())
	{
		if(file.exists(file.path(dir, "cor_summary.txt")))
		{
			print(dir)
			for (pco in pco_test){
				for (p in p_test){
					for (r in r_test){
						filter_cor_matrix(file.path(dir,"cor_summary_pearson.txt"), file.path(dir,"cor_summary_pearson_test.txt"),pco_thresh = pco,r_thresh = r,p_thresh = p,file_name = paste(".filtered.","PCO:",pco,".r:",r,".p:",p,sep = ""))                                                   
						filter_cor_matrix(file.path(dir,"cor_summary_spearman.txt"), file.path(dir,"cor_summary_spearman_test.txt"),pco_thresh = pco,r_thresh = r,p_thresh = p,file_name = paste(".filtered.","PCO:",pco,".r:",r,".p:",p,sep = ""))                                               
					}
				}
			}

			#filter_cor_matrix(file.path(dir,"cor_summary_pearson.txt"), file.path(dir,"cor_summary_pearson_test.txt"))
			#filter_cor_matrix(file.path(dir,"cor_summary_spearman.txt"), file.path(dir,"cor_summary_spearman_test.txt"))
		}
	}
}

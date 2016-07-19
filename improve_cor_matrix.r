improve_cor_matrix <- function(cor_matrix, metadata_matrix)
{
	cor <- read.table(file = cor_matrix, header = F, sep = "\t")
	meta_header <- read.table(file = metadata_matrix, nrows=1, header = F)
	print(unlist(meta_header)[1])
	rownames(cor) <- unlist(meta_header)
	for(i in 1:ncol(cor))
	{
		colnames(cor)[i] = paste("PCO", i)
	}
	cor[,"max"] <- apply(cor, 1, max)
	max.component <- apply(cor, 1, which.max)
	cbind(cor, max_component = max.component)
	cor[,"min"] <- apply(cor, 1, min)
	cor[,"min component"] <- apply(cor, 1, which.min)
	cor[,"total max"] <- ""
   	#cor[,"total max"][1] <- max(cor[,"max"], na.rm=T)
   	#cor[,"total max"][2] <- which.max(cor[,"max"], na.rm=T)
	cor[,"total min"] <- ""
	#cor[,"total min"][1] <- min(cor[,"min"], na.rm=T)
   	#cor[,"total min"][2] <- which.min(cor[,"min"], na.rm=T)

	write.table(cor, file = paste(cor_matrix, ".annotated", sep = ""), sep = "\t", quote = F, col.names = NA, row.names = T)
}

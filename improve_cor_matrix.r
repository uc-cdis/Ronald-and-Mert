improve_corr_matrix <- function(corr_matrix, metadata_matrix)
{
  # Adds row and column names to the corr_matrix
	corr <- read.table(file = corr_matrix, header = F, sep = "\t")
	meta_header <- read.table(file = metadata_matrix, nrows=1, header = F)
	rownames(corr) <- unlist(meta_header)
	for(i in 1:ncol(corr))
	{
		colnames(corr)[i] = paste("PCO", i)
	}

	write.table(corr, file = corr_matrix, sep = "\t", quote = F, col.names = NA, row.names = T)
}

edit_all_matrices <- function()
{
	for (dir in list.dirs())
	{
		if(file.exists(file.path(dir, "cor_summary.txt")))
		{
			print(dir)
			metadata_file <- list.files(path = dir, pattern = "\\.GDC_METADATA\\.txt$")
			metadata_file <- unlist(metadata_file)
			improve_corr_matrix(file.path(dir, "cor_summary_pearson_test.txt"), file.path(dir, metadata_file))
			improve_corr_matrix(file.path(dir, "cor_summary_spearman_test.txt"), file.path(dir, metadata_file))
		}
	}
}

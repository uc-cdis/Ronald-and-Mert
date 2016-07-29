filter_cor_p_matrix <- function(corr_matrix, 
                              p_matrix,
                              r_thresh = 0.7,
                              p_thresh = 0.0001,
                              pco_thresh = 15,
                              file_name =".filtered") {
  # Filters a correlation and p-value matrix pair
  #
  # Args:
  #   corr_matrix: File name for correlation matrix
  #   p_matrix: File name for p-value matrix
  #   r_thresh: Correlation threshold to filter for (lower bound)
  #   p_thresh: P-value threshold to filter for (upper bound)
  #   pco_thresh: Principal component to filter for (upper bound)
  #   file_name: Identifier that is added to end of exported file
  #
  # Returns:
  #   (None)
  #   Exports a table with metadata as rows and information about the data
  #   that passed the filter: 
  #   Metadata  PCO.#:r-value:p-value
	corr <- read.table(file = corr_matrix, header = T, sep = "\t", row.names = 1)
	p <- read.table(file = p_matrix, header = T, sep = "\t", row.names = 1)
  
  # Filter for principal component by chopping off all principal compnents
  # that are too large
	if (pco_thresh < length(corr[1,])){
		corr <- corr[,1:pco_thresh]
		p <- p[,1:pco_thresh]
	}

	# Turn the table into a list of data frames named with metadata
	corr.list <- setNames(split(corr, seq(nrow(corr))), rownames(corr))
	p.list <- setNames(split(p, seq(nrow(p))), rownames(p))

	b <- mapply(function(fromc, fromp) {
         # Takes in the data frame from corr.list and p.list for a piece of metadata
	       # and then spits out a filtered vector that's ordered based on the magnitude 
	       # of the correlation value
	       # mapply returns these vectors as a list
				 fromc[is.na(fromc)] <- 0
				 fromp[is.na(fromp)] <- 1
				 combinedVec <- paste(names(fromc), fromc, fromp, sep=":")
				 pos_c <- abs(as.numeric(fromc))
				 TFvec <- (pos_c > r_thresh) & (fromp < p_thresh)
				 combinedVecFiltered <- combinedVec[TFvec]
				 return(combinedVecFiltered[order(pos_c[TFvec])])
    },corr.list, p.list)
	browser()
	ls <- sapply(b,length)
	l <- length(ls[(ls!=0)])
	# make sure file is empty
	if(file.exists(paste0(corr_matrix, file_name,".pass_num:",l)))
	{
		file.remove(paste0(corr_matrix, file_name,".pass_num:",l))
		file.create(paste0(corr_matrix, file_name,".pass_num:",l))
	}
	namesList <- names(b)
	index <- 1
	lapply(b, function(elem) {
	     # write elements of b line by line
		   write(append(elem, namesList[index], 0), 
		         file = paste0(corr_matrix, file_name,".pass_num:",l), 
		         sep = "\t", ncolumns=1000, append=T)
		   index <<- index + 1
		   })
}

filter_all_matrices <- function(pco_test = c(15),
                                p_test = c(0.0001),
                                r_test = c(0.7)) {
  # Filters the correlation and p-value matrices for every combination of 
  # filters specified
  # 
  # Args: 
  #   pco_test: Vector of all desired principal component filters
  #   p_test: Vector of all desired p-value filters
  #   r_test: Vector of all desired r-value filters
  #
  # Returns:
  #   (None)
  #   Exports two (pearson and spearman) files each with a table that contains 
  #   information about which pieces of metadata passed the filter in a file
  #   ending with ".filtered.PCO:#.r:#.p:#.pass_num:#
	for(dir in list.dirs()) {
		if(file.exists(file.path(dir, "cor_summary.txt"))) {
			print(dir)
			for (pco in pco_test) {
				for (p in p_test) {
					for (r in r_test) {
						filter_cor_p_matrix(corr_matrix = file.path(dir,"cor_summary_pearson.txt"), 
						                    p_matrix = file.path(dir,"cor_summary_pearson_test.txt"),
						                    pco_thresh = pco,
						                    r_thresh = r,
						                    p_thresh = p,
						                    file_name = paste0(".filtered.",
						                                      "PCO:",pco,
						                                      ".r:",r,
						                                      ".p:",p))                                                   
						filter_cor_p_matrix(corr_matrix = file.path(dir,"cor_summary_spearman.txt"), 
						                    p_matrix = file.path(dir,"cor_summary_spearman_test.txt"),
						                    pco_thresh = pco,
						                    r_thresh = r,
						                    p_thresh = p,
						                    file_name = paste0(".filtered.",
						                                      "PCO:",pco,
						                                      ".r:",r,
						                                      ".p:",p))
					}
				}
			}
		}
	}
}

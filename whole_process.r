whole_process <- function(projects, skip_download = F)
{	
	#free | grep Mem | awk '{print $3/$2 * 100.0}'

	time <- Sys.time()
	time <- gsub(" ", "", time, fixed = TRUE)

	for(i in seq_len(sink.number())){
		sink(NULL)
	}
	log <- file(paste(paste(projects,collapse="."),".full.log",sep=""), open="wt")
	sink(log)

	to_log("START")

	source("~/git/Ronald-and-Mert/install_r_prereqs.r")       # from Kevin_R_scripts
	source("~/git/Ronald-and-Mert/preprocessing_tool.r")      # "                  "
	source("~/git/Ronald-and-Mert/heatmap_dendrogram.r")      # "                  "
	source("~/git/Ronald-and-Mert/calculate_pco.r")           # "                  "
	source("~/git/Ronald-and-Mert/render_calculated_pcoa.r")  # "                  "
	source("~/git/Ronald-and-Mert/GDC_raw_count_merge.vRonald.r") # from ronald_r_scripts
	source("~/git/Ronald-and-Mert/get_listof_UUIDs.r")            # "                   "
	source("~/git/Ronald-and-Mert/GDC_metadata_download.RandM.r")
	library(DESeq)

	to_log("Downloading data")
	if(!skip_download)
	{
		if((any(file.exists(projects)) == F) == F)
		{
			list_of_projects <- vector()
				for(i in projects)
				{
					list_of_projects <- scan(i, what = "character")
				}
			for(i in list_of_projects)
			{
				download_project_data(i)
			}
		}
		else
		{
			for(i in projects)
			{
			  print(i)
				readline()
			  download_project_data(i)
			}
		}

	to_log("Download Completed")
	# unzips data into .counts files

	to_log("Unzipping")
	system("gunzip *.gz")
	}
	system("ls | grep .counts$ > counts_files")
	to_log("Unzipping completed")

	# merges abundance data together
	to_log("Merging data")
	GDC_raw_count_merge(id_list="counts_files")
	to_log("Merge completed")
	# gets UUIDs of data merged together
	to_log("Exporting UUID's")
	export_listof_UUIDs(tsv = "counts_files.merged_data.txt")
	to_log("Export Completed")
	# gets metadata file
	to_log("Getting Metadata")
	get_GDC_metadata("counts_files.merged_data_file_UUIDs", my_rot = "yes")
	files <- list.files()[grep("GDC_METADATA.txt", list.files())]
	details <- file.info(list.files()[grep("GDC_METADATA.txt", list.files())])
	details <- details[with(details, order(as.POSIXct(mtime))), ]
	metadata_filename <- rownames(details)[length(files)] 
	print(metadata_filename)
	to_log("Metadata Download Completed")

	# normalizing the data
	to_log("Preprocessing")
	preprocessing_tool(data_in = "counts_files.merged_data.txt", produce_boxplots = TRUE)
	to_log("Preprocessing Completed")
	# visualize the data
	#heatmap_dendrogram(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")

	to_log("Calculating PCoA")
	calculate_pco(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
	print(metadata_filename)
	to_log("Calculated PCoA")
	to_log("Rendering PCoA")
	render_calcualted_pcoa("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", 
			metadata_table = metadata_filename, 
			use_all_metadata_columns = T)
	to_log("Rendering Completed")
	to_log("END")
	sink()

}

to_log <- function(message){
	print(message)
	print(Sys.time())
	print(gc(verbose = T,reset = F))
}

library(RJSONIO)
library(RCurl)

whole_process <- function(projects,
							skip_to_merging = F,
							skip_to_metadata = F, 
							skip_rendering = F, 
							skip_to_corr_matrix = F)
{	
	if(skip_to_metadata == T)
	{
		skip_to_merging = T
	}
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
	source("~/git/Ronald-and-Mert/correlate_pcoa.r")
	library(DESeq)

	if(!skip_to_merging)
	{
		to_log("Downloading data")
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
			  download_project_data(i)
			}
		}

	to_log("Download Completed")
	# unzips data into .counts files

	to_log("Unzipping")
	system("gunzip *.gz")
	}
	if (!skip_to_metadata)
	{
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
	}
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
	if(skip_rendering == F)
	{
		to_log("Rendering PCoA")
		render_calcualted_pcoa("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", 
				metadata_table = metadata_filename, 
				use_all_metadata_columns = T)
		to_log("Rendering Completed")
	}	
	to_log("Creating Correlation Matrix")
	correlate_pcoa("counts_files.merged_data_file_UUIDs.GDC_METADATA.txt","counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA")
	to_log("Created Correlation Matrix")
	to_log("END")
	sink()

}

to_log <- function(message){
	print(message)
	print(Sys.time())
	print(gc(verbose = T,reset = F))
}

individual_analysis <- function()
{
	# API call to get list of projects
	project_list <- fromJSON(getURL("https://gdc-api.nci.nih.gov/projects?fields=project_id&from=1&size=50&sort=project.project_id:asc&pretty=true"))$data$hits
	for(p in project_list)
	{
		#sink(file=NULL)
		print(p)
		# check if folder exists
		if(!dir.exists(file.path("/mnt/single_projects", p)))
		{
			# create folder if it doesn't exist
			dir.create(file.path("/mnt/single_projects", p))
			# go into folder and run whole_process()
			setwd(file.path("/mnt/single_projects", p))
			whole_process(p, skip_rendering = T)
			# return back to /mnt
			setwd("/mnt/single_projects")
		}
		# check if it has already downloaded all the data
		else if(file.exists(file.path("/mnt/single_projects", p, "counts_files.merged_data_file_UUIDs")) & 
				!file.exists(file.path("/mnt/single_projects", p, "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA")))
		{
			setwd(file.path("/mnt/single_projects", p))
			whole_process(p, skip_rendering = T, skip_to_metadata = T)
		}
		temp <- list.files()
		metadataFile <- temp[grepl("GDC_METADATA.txt$", temp)]
		if (!is.null(metadataFile) & 
			 	file.exists(file.path("/mnt/single_projects", p, "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA")) &
			 	!file.exists(file.path("/mnt/single_projects", p, "cor_summary.txt")))
		{
			for(i in metadataFile)
			# something's still wrong with this part, need to test more
			{
				source("~/git/Ronald-and-Mert/correlate_pcoa.r")
				setwd(file.path("/mnt/single_projects", p))
				correlate_pcoa(i,"counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA")
			}
		}
	}
}

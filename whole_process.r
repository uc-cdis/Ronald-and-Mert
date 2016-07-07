whole_process <- function(projects, skip_download = F)
{

	time <- Sys.time()
	time <- gsub(" ", "", time, fixed = TRUE)

	for(i in seq_len(sink.number())){
		sink(NULL)
	}
	error_log <- file(paste(paste(projects,collapse="."),"error.log",sep="."), open="wt")
	sink(error_log,type="message")
	
	system(paste("echo Starting > ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	source("~/git/Ronald-and-Mert/install_r_prereqs.r")       # from Kevin_R_scripts
	source("~/git/Ronald-and-Mert/preprocessing_tool.r")      # "                  "
	source("~/git/Ronald-and-Mert/heatmap_dendrogram.r")      # "                  "
	source("~/git/Ronald-and-Mert/calculate_pco.r")           # "                  "
	source("~/git/Ronald-and-Mert/render_calculated_pcoa.r")  # "                  "
	source("~/git/Ronald-and-Mert/GDC_raw_count_merge.vRonald.r") # from ronald_r_scripts
	source("~/git/Ronald-and-Mert/get_listof_UUIDs.r")            # "                   "
	source("~/git/Ronald-and-Mert/GDC_metadata_download.RandM.r")
	library(DESeq)
	
	free | grep Mem | awk '{print $3/$2 * 100.0}'

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
				download_project_data(i)
			}
		}
	}

	# unzips data into .counts files
	system(paste("echo 'Unzipping files' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	system("gunzip *.gz")

	system("ls | grep .counts$ > counts_files")

	system(paste("echo 'Finished Unzipping' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	# merges abundance data together
	system(paste("echo 'Merging Data Files' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	GDC_raw_count_merge(id_list="counts_files")

	system(paste("echo 'Merge Completed' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	# gets UUIDs of data merged together
	system(paste("echo 'Exporting UUID's >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	export_listof_UUIDs(tsv = "counts_files.merged_data.txt")

	system(paste("echo 'Export Completed >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	# gets metadata file
	system(paste("echo 'Getting Metadata' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	get_GDC_metadata("counts_files.merged_data_file_UUIDs", my_rot = "yes")
	files <- list.files()[grep("GDC_METADATA.txt", list.files())]
	details <- file.info(list.files()[grep("GDC_METADATA.txt", list.files())])
	details <- details[with(details, order(as.POSIXct(mtime))), ]
	metadata_filename <- rownames(details)[length(files)] 
	print(metadata_filename)

	system(paste("echo 'Metadata Downloaded' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	# normalizing the data
	system(paste("echo 'Preprocessing' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	preprocessing_tool(data_in = "counts_files.merged_data.txt", produce_boxplots = TRUE)

	system(paste("echo 'Done Preprocessing' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	# visualize the data
	#heatmap_dendrogram(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
	system(paste("echo 'Calculating PCoA' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	calculate_pco(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
	print(metadata_filename)

	system(paste("echo 'PCoA Calculated' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))



	system(paste("echo 'Rendering PCoA' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	render_calcualted_pcoa("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", 
							metadata_table = metadata_filename, 
							use_all_metadata_columns = T)

	system(paste("echo 'Rendered all PCoA's >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	system(paste("echo 'END' >> ",paste(projects,collapse="."),time,".log",sep=""))
	system(paste("echo ",Sys.time()," >> ",paste(projects,collapse="."),time,".log",sep=""))

	sink()
}

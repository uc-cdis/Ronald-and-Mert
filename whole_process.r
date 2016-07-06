whole_process <- function(project1, project2="", skip_download = F)
{
  
  time <- Sys.time()
  time <- gsub(" ", "", time, fixed = TRUE)
  
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
  error_log <- file(paste(project1,".",project2,time,".error.log",sep=""), open="wt")
  sink(error_log,type="message")
  
  system(paste("echo Starting > ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  source("/Users/mertbozfakioglu/Documents/git/Ronald-and-Mert/install_r_prereqs.r")       # from Kevin_R_scripts
  source("/Users/mertbozfakioglu/Documents/git/Ronald-and-Mert/preprocessing_tool.r")      # "                  "
  source("/Users/mertbozfakioglu/Documents/git/Ronald-and-Mert/heatmap_dendrogram.r")      # "                  "
  source("/Users/mertbozfakioglu/Documents/git/Ronald-and-Mert/calculate_pco.r")           # "                  "
  source("/Users/mertbozfakioglu/Documents/git/Ronald-and-Mert/render_calculated_pcoa.r")  # "                  "
  source("/Users/mertbozfakioglu/Documents/git/Ronald-and-Mert/GDC_raw_count_merge.vRonald.r") # from ronald_r_scripts
  source("/Users/mertbozfakioglu/Documents/git/Ronald-and-Mert/get_listof_UUIDs.r")            # "                   "
  source("/Users/mertbozfakioglu/Documents/git/Ronald-and-Mert/GDC_metadata_download.RandM.r")
  library(DESeq)
  
  if(!skip_download)
  {
  	# downloads data
  	download_project_data(project1)
  	if (project2 != "")
  	{
 		download_project_data(project2)
  	}
  
  # unzips data into .counts files
  system(paste("echo 'Unzipping files' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))

  system("gunzip *.gz")
  }
  system("ls | grep .counts$ > counts_files")
  
  system(paste("echo 'Finished Unzipping' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  # merges abundance data together
  system(paste("echo 'Merging Data Files' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))

  GDC_raw_count_merge(id_list="counts_files")
  
  system(paste("echo 'Merge Completed' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  # gets UUIDs of data merged together
  system(paste("echo 'Exporting UUID's >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  export_listof_UUIDs(tsv = "counts_files.merged_data.txt")
  
  system(paste("echo 'Export Completed >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  # gets metadata file
  system(paste("echo 'Getting Metadata' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  get_GDC_metadata("counts_files.merged_data_file_UUIDs", my_rot = "yes")
  files <- list.files()[grep("GDC_METADATA.txt", list.files())]
  details <- file.info(list.files()[grep("GDC_METADATA.txt", list.files())])
  details <- details[with(details, order(as.POSIXct(mtime))), ]
  metadata_filename <- rownames(details)[length(files)] 
  print(metadata_filename)

  system(paste("echo 'Metadata Downloaded' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))

  # normalizing the data
  system(paste("echo 'Preprocessing' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  preprocessing_tool(data_in = "counts_files.merged_data.txt", produce_boxplots = TRUE)
  
  system(paste("echo 'Done Preprocessing' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  # visualize the data
  #heatmap_dendrogram(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
  system(paste("echo 'Calculating PCoA' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  calculate_pco(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
  print(metadata_filename)

  #dir_path <- paste("/mnt/", project1, "_", project2, sep="")
  #system(paste("mkdir", dir_path))
  render_calcualted_pcoa("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", 
                          metadata_table = metadata_filename,  
                          use_all_metadata_columns = T, 
                          mv_to_mount = dir_path)
  
  system(paste("echo 'PCoA Calculated' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  
  
  system(paste("echo 'Rendering PCoA' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  render_calcualted_pcoa("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", metadata_table = metadata_filename, use_all_metadata_columns = T, mv_to_mount = T)
  
  system(paste("echo 'Rendered all PCoA's >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  system(paste("echo 'END' >> ",project1,".",project2,time,".log",sep=""))
  system(paste("echo ",Sys.time()," >> ",project1,".",project2,time,".log",sep=""))
  
  sink()
}

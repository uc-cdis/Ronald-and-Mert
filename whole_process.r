whole_process <- function(project1, project2="", skip_download = F)
{
  source("~/git/Ronald-and-Mert/install_r_prereqs.r")       # from Kevin_R_scripts
  source("~/git/Ronald-and-Mert/preprocessing_tool.r")      # "                  "
  source("~/git/Ronald-and-Mert/heatmap_dendrogram.r")      # "                  "
  source("~/git/Ronald-and-Mert/calculate_pco.r")           # "                  "
  source("~/git/Ronald-and-Mert/render_calculated_pcoa.r")  # "                  "
  source("~/git/Ronald-and-Mert/GDC_raw_count_merge.vRonald.r") # from ronald_r_scripts
  source("~/git/Ronald-and-Mert/get_listof_UUIDs.r")            # "                   "
  source("~/git/Ronald-and-Mert/GDC_metadata_download.RandM.r")
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
  system("gunzip *.gz")
  }
  system("ls | grep .counts$ > counts_files")
  
  # merges abundance data together
  GDC_raw_count_merge(id_list="counts_files")
  
  # gets UUIDs of data merged together
  export_listof_UUIDs(tsv = "counts_files.merged_data.txt")
  
  # gets metadata file
  get_GDC_metadata("counts_files.merged_data_file_UUIDs", my_rot = "yes")
  files <- list.files()[grep("GDC_METADATA.txt", list.files())]
  details <- file.info(list.files()[grep("GDC_METADATA.txt", list.files())])
  details <- details[with(details, order(as.POSIXct(mtime))), ]
  metadata_filename <- rownames(details)[length(files)] 
  print(metadata_filename)

  # normalizing the data
  preprocessing_tool(data_in = "counts_files.merged_data.txt", produce_boxplots = TRUE)
  
  # visualize the data
  #heatmap_dendrogram(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
  calculate_pco(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
  print(metadata_filename)
  render_calcualted_pcoa("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", metadata_table = metadata_filename, use_all_metadata_columns = T, mv_to_mount = T)
}

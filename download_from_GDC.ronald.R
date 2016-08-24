library(RJSONIO)
library(RCurl)

source("~/git/Ronald-and-Mert/preprocessing_tool.r")
source("~/git/Ronald-and-Mert/GDC_raw_count_merge.vRonald.r")
source("~/git/Ronald-and-Mert/heatmap_dendrogram.r")
source("~/git/Ronald-and-Mert/GDC_metadata_download.RandM.r")
library(DESeq)

CompareNormalCancerPipeline <- function(projects) {
  if (any(file.exists(projects))) {
    vec_of_projects <- vector()
    # loop through every filename and read in projects to vec_of_projects
    for(i in projects) {
      vec_of_projects <- c(vec_of_projects, scan(i, what = "character", sep = "\n"))
    }
    # loop through vec_of_projects and download all HTSeq counts files
    for(i in vec_of_projects) {
      if (size_lim) {
        download_normal_from_GDC(i, size_lim = size_lim)
      } else {
        download_normal_from_GDC(i)
      }
    }
  } else {
    # else if user inputed a vector of project_ids
    for(i in projects) {
      if (size_lim) {
        download_normal_from_GDC(i, size_lim = size_lim)
      } else {
        download_normal_from_GDC(i)
      }
    }
  }
  system("ls | grep .counts$ > normal_filenames")
  vec.normal.filenames <- scan("normal_filenames", what = "character", sep = "\n")
  for (f in vec.normal.filenames) {
    my_call <- paste0("https://gdc-api.nci.nih.gov/files?fields=file_id&size=99999&pretty=true&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.sample_type%22%2C%22value%22%3A%5B%22Blood+Derived+Normal%22%2C%22Solid+Tissue+Normal%22%2C%22Bone+Marrow+Normal%22%2C%22Fibroblasts+from+Bone+Marrow+Normal%22%2C%22Buccal+Cell+Normal%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22analysis.workflow_type%22%2C%22value%22%3A%5B%22HTSeq+-+Counts%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22TXT%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22", p, "%22%5D%7D%7D%5D%7D")
    my_call.json <- fromJSON(getURL(my_call))
    UUID.list <- unlist(my_call.json$data$hits)
    for(j in UUID.list) {
      print(paste0(j, ": ", j))
      system(paste0("curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/data/", 
                    j, "'"))
    }
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

  # gets metadata_filename 
  # (obsolete now but useful for interacting with old metadata filenames)
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
}
download_normal_from_GDC <- function(projects) {
# download normal sample HTSeq Counts data for a given project
# 
# Args:
#   project: vector of files that contain project_ids or 
  #            a vector of projects
  #
  # Returns:
  #   (None)
  #   Downloads all HTSeq Counts files as *.gz, need to unzip
  for (p in projects) {
    print(p)
    my_call <- paste0("https://gdc-api.nci.nih.gov/files?fields=file_id&size=99999&pretty=true&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.sample_type%22%2C%22value%22%3A%5B%22Blood+Derived+Normal%22%2C%22Solid+Tissue+Normal%22%2C%22Bone+Marrow+Normal%22%2C%22Fibroblasts+from+Bone+Marrow+Normal%22%2C%22Buccal+Cell+Normal%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22analysis.workflow_type%22%2C%22value%22%3A%5B%22HTSeq+-+Counts%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22TXT%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22", p, "%22%5D%7D%7D%5D%7D")
    my_call.json <- fromJSON(getURL(my_call))
    UUID.list <- unlist(my_call.json$data$hits)
    for(j in UUID.list) {
      print(paste0(j, ": ", j))
      system(paste0("curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/data/", 
                  j, "'"))
    }
  }
}

download_cancer_from_GDC <- function(projects) {
  # download cancer sample HTSeq Counts data for a given project
  # 
  # Args:
  #   project: vector of files that contain project_ids or 
  #            a vector of projects
  #
  # Returns:
  #   (None)
  #   Downloads all HTSeq Counts files as *.gz, need to unzip
  for (p in projects) {
    my_call <- paste0("https://gdc-api.nci.nih.gov/files?fields=file_id&size=99999&pretty=true&filters=%7B%0D%0A%09%22op%22%3A%22and%22%2C%0D%0A%09%22content%22%3A%5B%7B%0D%0A%09%09%22op%22%3A%22in%22%2C%0D%0A%09%09%22content%22%3A%7B%0D%0A%09%09%09%22field%22%3A%22cases.samples.sample_type%22%2C%0D%0A%09%09%09%22value%22%3A%5B%22Primary+Tumor%22%5D%0D%0A%09%09%09%7D%0D%0A%09%09%7D%2C%7B%0D%0A%09%09%22op%22%3A%22in%22%2C%0D%0A%09%09%22content%22%3A%7B%0D%0A%09%09%09%22field%22%3A%22analysis.workflow_type%22%2C%0D%0A%09%09%09%22value%22%3A%5B%22HTSeq+-+Counts%22%5D%0D%0A%09%09%09%7D%0D%0A%09%09%7D%2C%7B%0D%0A%09%09%22op%22%3A%22in%22%2C%0D%0A%09%09%22content%22%3A%7B%0D%0A%09%09%09%22field%22%3A%22files.data_format%22%2C%0D%0A%09%09%09%22value%22%3A%5B%22TXT%22%5D%0D%0A%09%09%09%7D%0D%0A%09%09%7D%2C%7B%0D%0A%09%09%22op%22%3A%22%3D%22%2C%0D%0A%09++++%22content%22%3A%7B%0D%0A%09++++%09%22field%22%3A%22cases.project.project_id%22%2C%0D%0A%09++++%09%22value%22%3A%5B%22", p, "%22%5D%0D%0A%09++++%7D%0D%0A%09%7D%5D%0D%0A%7D")
    my_call.json <- fromJSON(getURL(my_call))
    UUID.list <- unlist(my_call.json$data$hits)
    for(j in UUID.list) {
      print(paste0(j, ": ", j))
      system(paste("curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/data/", 
                   j,
                   "'",
                   sep=""))
    }
  }
}

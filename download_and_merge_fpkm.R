library(RJSONIO)
library(RCurl)

DownloadAndMergeFPKM <- function(projects, # use files with project_ids or a vector
                          # of project_ids
                          # parameters to skip certain parts of script
                          skip_to_merging = F,
                          # size_lim for every project: only downloads XX
                          # number of files per project if XX files are
                          # available
                          size_lim = F) {
  # from Kevin_R_scripts
  source("~/git/Ronald-and-Mert/install_r_prereqs.r")
  # from Ronald-and-Mert
  source("~/git/Ronald-and-Mert/GDC_raw_count_merge.stripped.R")
  source("~/git/Ronald-and-Mert/remove_FPKM_file_endings.R")

  if (!skip_to_merging) {
    # if user inputed a file name with project_ids
    if (any(file.exists(projects))) {
      vec_of_projects <- vector()
      # loop through every filename and read in projects to vec_of_projects
      for(i in projects) {
        vec_of_projects <- c(vec_of_projects, scan(i, what = "character"))
      }
      # loop through vec_of_projects and download all HTSeq counts files
      for(i in vec_of_projects) {
        if (size_lim) {
          download_project_data(i, size_lim = size_lim)
        } else {
          download_project_data(i)
        }
      }
    } else {
      # else if user inputed a vector of project_ids
      for(i in projects) {
        if (size_lim) {
          download_project_data(i, size_lim = size_lim)
        } else {
          download_project_data(i)
        }
      }
    }


    system("gunzip *.gz")
  }
  # spits all of .counts file names in counts_files
  system("ls | grep FPKM\\.txt$ > counts_files")

  # merges abundance data together
  GDC_raw_count_merge(id_list="counts_files")

  RemoveFPKMColEndings("counts_files.merged_data.txt")
}

download_project_data <- function(project, size_lim = F) {
  # download HTSeq Counts data for a given project
  #
  # Args:
  #   project: vector of files that contain project_ids or
  #            a vector of projects
  #
  # Returns:
  #   (None)
  #   Downloads all HTSeq Counts files as *.gz, need to unzip
  vector_of_files <- unlist(get_UUIDs_from_project(project))
  if(size_lim && size_lim < length(vector_of_files)) {
    vector_of_files <- vector_of_files[sample(1:length(vector_of_files), size_lim)]
  }
  for(j in 1:length(vector_of_files)) {
    print(paste0(j, ": ", vector_of_files[j]))
    system(paste("curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/data/",
                 vector_of_files[j],
                 "'",
                 sep=""))
  }
}

get_UUIDs_from_project <- function(project) {
  # Gets a list of UUIDs from a project id
  #
  # Args:
  #   project: project id
  #
  # Returns:
  #   Vector of all UUIDs
  before_id="https://gdc-api.nci.nih.gov/files?fields=file_id&size=99999&pretty=true&filters=%7B%0D%0A++%22op%22%3A%22and%22%2C%0D%0A++%22content%22%3A%5B%7B%0D%0A++++%22op%22%3A%22%3D%22%2C%0D%0A++++%22content%22%3A%7B%0D%0A++++++%22field%22%3A%22analysis.workflow_type%22%2C%0D%0A++++++%22value%22%3A%5B%0D%0A++++++++%22HTSeq+-+FPKM%22%0D%0A++++++++%5D%0D%0A++++%7D%0D%0A++%7D%2C+%7B%0D%0A++++%22op%22%3A%22%3D%22%2C%0D%0A++++%22content%22%3A%7B%0D%0A++++++%22field%22%3A%22cases.project.project_id%22%2C%0D%0A++++++%22value%22%3A%5B%0D%0A++++++++%22"
  after_id="%22%0D%0A++++++++%5D%0D%0A++++%7D%0D%0A++%7D%5D%0D%0A%7D"
  my_call <- paste0(before_id, project, after_id)
  my_call.json <- fromJSON(getURL(my_call))
  return(my_call.json$data$hits)
}


DownloadAllFPKM <- function(overall.path = "/mnt/FPKM") {
  # perform the whole_process on every project
  #
  # Args:
  #   None
  # Returns:
  #   Produces a folder named the project name in the current directory

  # API call to get list of projects
  project_list <- fromJSON(getURL("https://gdc-api.nci.nih.gov/projects?fields=project_id&from=1&size=50&sort=project.project_id:asc&pretty=true"))$data$hits
  for(p in project_list) {
    print(p)
    # check if folder exists
    if (!dir.exists(file.path(overall.path, p))) {
      # create folder if it doesn't exist
      dir.create(file.path(overall.path, p))
      # go into folder and run function
      setwd(file.path(overall.path, p))
      DownloadAndMergeFPKM(p)
      # return back to /mnt
      setwd(overall.path)
    }
  }
}

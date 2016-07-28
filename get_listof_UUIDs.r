pkgTest <- function(x) {
  if (!require(x,character.only = TRUE)) {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only=TRUE)) {
      stop("Package not found")
    }
  }
}
pkgTest("RJSONIO")
pkgTest("RCurl")
pkgTest("tools")
library(RJSONIO)
library(RCurl)
library(tools)

export_listof_UUIDs <- function(project = "", names = "", tsv = "", size_lim = F) {
  # creates a list of UUIDs
  # 
  # Args:
  #   One of the following:
  #     project: either a vector of files names that contain unique project
  #              a names or a vector of project names
  #     names: vector of file names (not UUIDs)
  #     tsv: a tsv file with file names as column names
  #   size_lim: limit for number of files to download from each project
  # 
  # Returns:
  #   List of UUIDs in a file named *filename*_file_UUIDs
  if (project != "") {
    vector_of_files <- c()
    if (size) {
      vector_of_files <- unlist(get_UUIDs_from_project(project, size))
    } else {
      vector_of_files <- unlist(get_UUIDs_from_project(project))
    }
    if (length(vector_of_files) == 0) {
      # error message if no UUIDs found
      stop(paste("project name:", 
                 paste('"', project, '"', sep=""), 
                 "not found found!"))
    }
    write(sort(vector_of_files), 
          paste(file_path_sans_ext(project),"file_UUIDs",sep="_"), 
          sep = "\n")
  } else if (names != "") {
    names_list <- flatten_list(as.list(scan(file=names, what="character")))
    write(sort(unlist(get_UUIDs_from_names(names_list))), 
          paste(file_path_sans_ext(names),"file_UUIDs",sep="_"), sep = "\n")
  } else if (tsv != "") {
    names_list <- get_names_list_from_tsv(tsv)
    write(unlist(get_UUIDs_from_names(names_list)), 
          paste(file_path_sans_ext(tsv),"file_UUIDs",sep="_"), 
          sep = "\n")
  } else {
    print("Error: no argument given to function")
  }
}

get_UUIDs_from_project <- function(project, size_lim = 99999) { 
  # Gets a list of UUIDs from a project id
  #
  # Args: 
  #   project: project id
  #   size_lim: limit for number of files to download from each project
  #
  # Returns: 
  #   Vector of all UUIDs
  before_id=paste0("https://gdc-api.nci.nih.gov/files?fields=file_id&size=", 
                   size, 
                   "&pretty=true&filters=%7B%0A%09%22op%22%3A%22and%22%2C%0A%09%22content%22%3A%5B%7B%0A%09%09%22op%22%3A%22%3D%22%2C%0A%09%09%22content%22%3A%7B%0A%09%09%09%22field%22%3A%22analysis.workflow_type%22%2C%0A%09%09%09%22value%22%3A%5B%0A%09%09%09%09%22HTSeq%20-%20Counts%22%0A%09%09%09%5D%0A%09%09%7D%0A%09%7D%2C%20%7B%0A%09%09%22op%22%3A%22%3D%22%2C%0A%09%09%22content%22%3A%7B%0A%09%09%09%22field%22%3A%22cases.project.project_id%22%2C%0A%09%09%09%22value%22%3A%5B%0A%09%09%09%22")
  after_id="%22%0A%09%09%09%5D%0A%09%09%7D%0A%09%7D%5D%0A%7D"
  my_call <- gsub(" ", "", paste(before_id, project, after_id))
  my_call.json <- fromJSON(getURL(my_call))
  return(my_call.json$data$hits)
}

# FILTER:
#{
#  "op":"and",
#  "content":[{
#    "op":"=",
#    "content":{
#      "field":"analysis.workflow_type",
#      "value":[
#        "HTSeq - Counts"
#        ]
#    }
#  }, {
#    "op":"=",
#    "content":{
#      "field":"cases.project.project_id",
#      "value":[
#        "YOUR_PROJECT_HERE"
#        ]
#    }
#  }]
#}


get_UUIDs_from_names <- function(names_list) { 
  # gets the UUIDs froma vector of names
  #
  # Args:
  #   names_list: a vector of file_names
  # 
  # Returns:
  #   Vector of UUIDs corresponding to each file_name
  before_id="https://gdc-api.nci.nih.gov/files?pretty=true&fields=file_id&filters=%7B%0A%20%20%20%22op%22%20%3A%20%22%3D%22%20%2C%0A%20%20%20%22content%22%20%3A%20%7B%0A%20%20%20%20%20%20%20%22field%22%20%3A%20%22file_name%22%20%2C%0A%20%20%20%20%20%20%20%22value%22%20%3A%20%5B%20%22"
  after_id=".htseq.counts.gz%22%20%5D%0A%20%20%20%7D%0A%7D"
  UUID_list <- vector()
  for (i in 1:length(names_list)) {
  	print(i)
    my_call <- gsub(" ", "", paste(before_id, names_list[i], after_id))
    my_call.json <- fromJSON(getURL(my_call))
    UUID_list <- c(UUID_list, my_call.json$data$hits)

    if(length(my_call.json$data$hits) == 0) {
      # Throw error if API call gives nothing
      stop(paste("no file_UUID found for", names_list[i], "at name #", i))
    }
  }
  
  return(UUID_list)
}

get_names_list_from_tsv <- function(file) {
  my.table <- read.table(file, row.names=1, header=TRUE, sep="\t", comment.char="", quote="", check.names=FALSE)
  return(colnames(my.table))
}

download_project_data <- function(project, size_lim = 99999) {
  # download HTSeq Counts data for a given project
  # 
  # Args:
  #   project: vector of files that contain project_ids or 
  #            a vector of projects
  #
  # Returns:
  #   (None)
  #   Downloads all HTSeq Counts files as *.gz, need to unzip
  vector_of_files <- unlist(get_UUIDs_from_project(project, size_lim))
  for(j in 1:length(vector_of_files)) {
    print(j)
    system(paste("curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/data/", 
                 vector_of_files[j],
                 "'",
                 sep=""))
  }
}


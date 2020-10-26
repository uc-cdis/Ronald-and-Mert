download_from_GDC <- function(list.of.UUIDs = F,
                              project_id = F,
                              analysis.workflow_type = "HTSeq - Counts"
                              UUIDs.file = F,
                              cases.samples.sample_type = "Primary Tumor",
                              size_lim = F) {
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

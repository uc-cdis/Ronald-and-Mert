source("~/git/Ronald-and-Mert/render_calculated_pcoa.r")
source("~/git/Ronald-and-Mert/extract_top3_comp.R")

# hits.demographic.gender.1

render_all_metadata_pcoa <- function(yourdir, metadataCol) {
  for(folder in list.dirs(path = yourdir, recursive = F)) {
    orderedFilename <- file.path(folder, "cor_summary_spearman.txt.ordered")
    PCoAFilename <- file.path(folder, "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA")
    metadataFilename <- file.path(folder, list.files(path = folder,
                                                  pattern = ".*METADATA.*"))
    top3Comp <- extract_top3_comp(orderedFilename, metadataCol)
    render_calcualted_pcoa(PCoA_in = PCoAFilename, components = top3Comp,
                           metadata_table = metadataFilename,
                           metadata_column_index = metadataCol)
    print(paste(folder, "done!"))

  }
}

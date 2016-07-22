#determine_unique_and_common_metadata_fields
create_frequency_table_for_significant_metadata_fields <- function(sig_meta_fields_file)
{
  fields <- c()
  for(dir in list.dirs())
  {
    if(file.exists(file.path(dir, sig_meta_fields_file)))
    {
      print(dir)
      if(file.size(file.path(dir,sig_meta_fields_file))>0){
        sig_file <- unlist(read.table(file.path(dir,sig_meta_fields_file),stringsAsFactors = F))
        fields <- c(fields,sig_file)
      }
    }
  }
  tab <- (table(fields))
  write.table(tab,"significant_metadata_frequencies")
}

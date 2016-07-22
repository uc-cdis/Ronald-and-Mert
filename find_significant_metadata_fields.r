find_significant_metadata_fields <- function(cor_file){
  #read filtered data
  fil <- scan(cor_file,what = "charachter")
  #places of p/r values
  b <- which(substr(fil,start = 1,stop = 3)=="PCO")
  unprocessed_fields <- fil[b-1]
  c <- which(substr(unprocessed_fields,start = 1,stop = 3)=="PCO")
  fields <- unprocessed_fields[-c]
  fields <- as.vector(fields)
  colnames(fields) <- NULL
  write.table(fields,file = paste(cor_file,".significant_metadata_fields",sep = ""),sep = "\t",row.names = F,col.names = F)
}

find_sig_meta_fields_in_all_projects <- function(filtered_files)
{
  for(dir in list.dirs())
  {
    if(file.exists(file.path(dir, "cor_summary.txt")))
    {
      print(dir)
      for (file in filtered_files){
        #file.path(dir,"cor_summary_pearson.txt")
        find_significant_metadata_fields(file.path(dir,file))
      }
    }
  }
}

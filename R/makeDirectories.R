makeDirectories <- function( path = NULL,
                             twoSex = T){
  ##-- Check if the directory already exists 
  if(is.null(path)){path <- getwd()}
  
  if(dir.exists(path)){
    split <- unlist(strsplit(path, split = "/"))
    message(paste0("A folder named '", split[length(split)], "' already exists in the specified directory:"))
    message(paste(split[-length(split)], collapse = "/"))
    } 

  ##-- Set-up directory structure
  if(twoSex){
  dir.create( file.path(path, "nimbleInFiles/Hann"), showWarnings = F)
  dir.create( file.path(path, "nimbleInFiles/Hunn"), showWarnings = F)
  dir.create( file.path(path, "nimbleOutFiles/Hann"), showWarnings = F)
  dir.create( file.path(path, "nimbleOutFiles/Hunn"), showWarnings = F)
  } else{
    dir.create( file.path(path, "nimbleInFiles"), showWarnings = F)
    dir.create( file.path(path, "nimbleOutFiles"), showWarnings = F)  
  }
  dir.create( file.path(path, "data"), showWarnings = F)
  dir.create( file.path(path, "figures"), showWarnings = F)
  dir.create( file.path(path, "tables"), showWarnings = F)
  dir.create( file.path(path, "rasters"), showWarnings = F)
  dir.create( file.path(path, "reports"), showWarnings = F)
}

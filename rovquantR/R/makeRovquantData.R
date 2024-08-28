makeRovquantData <- function(
    ##-- paths
  species = c("bear","wolf","wolverine"),
  data_dir = "./Data",
  working_dir = NULL,
  
  ##-- data
  years = NULL,
  sex = c("Hann","Hunn"),
  aug.factor = 2,
  sampling.months = list(4,5,6,7,8,9,10,11),
  
  ##-- habitat
  habitat.res = 20000, 
  buffer.size = 50000,
  max.move.dist = 250000,
  
  ##-- detectors
  detector.res = 5000,
  subdetector.res = 1000,
  max.det.dist = 70000,
  resize.factor = 1,
  
  ##-- miscellanious
  plot.check = FALSE,
  print.report = TRUE
) {
  
  ##-- Check species and set corresponding sampling period
  if(sum(grep("bear", species, ignore.case = T))>0|sum(grep("bjorn", species, ignore.case = T))>0){
    out <- makeRovquantData_bear(
      ##-- paths
      data_dir,
      working_dir,
      
      ##-- data
      years,
      sex,
      aug.factor,
      sampling.months,
      
      ##-- habitat
      habitat.res, 
      buffer.size,
      max.move.dist,
      
      ##-- detectors
      detector.res = 5000,
      subdetector.res = 1000,
      max.det.dist = 70000,
      resize.factor = 1,
      
      ##-- miscellanious
      plot.check = FALSE,
      print.report = FALSE)
  }
  # if(sum(grep("wolf", species, ignore.case = T))>0|sum(grep("ulv", species, ignore.case = T))>0){
  #   out <- makeRovquantData_wolf(...)
  # }
  # if(sum(grep("wolverine", species, ignore.case = T))>0|sum(grep("jerv", species, ignore.case = T))>0){
  #   out <- makeRovquantData_wolverine(...)
  # }
  
  return(out)
}
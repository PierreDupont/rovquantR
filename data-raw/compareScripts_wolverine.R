load(file.path(working.dir,"data/CleanData_wolverine_2025-07-17.RData"))
dim(alive)
table(alive$Year)
dim(dead.recovery)


load("C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Test.0.1/data/myFilteredData.sp.RData")
dim(myFilteredData.sp$alive)
table(myFilteredData.sp$alive$Year)
dim(myFilteredData.sp$dead.recovery)


load("C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Test/_myFilteredData.sp.RData")
load("C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Test/myFilteredData_original.RData")
dim(myFilteredData.sp$alive)
table(myFilteredData.sp$alive$Year)
dim(myFilteredData.sp$dead.recovery)


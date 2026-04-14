
##-- Identify individuals and years when in pairs or group 


INDIVIDUAL_ID %>% 
  filter(Status %in% c("Family group","Pair")) %>%
  rename(Year = "ReprodYear (May 1 year y - Apr 30 y+1)") %>%
  select(Id,Year)



INDIVIDUAL_ID$STATUS_Numeric <- 1
INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% "Juvenile"] <- 2
INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% c("Pair" )] <- 3
INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% c("Family group")] <- 4

## FOR THE LAST YEAR GET THROUGH THE PAIR BASED-FILE FROM LINN.
indIDRovBase <- unlist(lapply(strsplit(y.ar$Id.vector," "), function(x) x[1]))

ALLIDS <- c(unique(indIDRovBase))

y.obsALL <- matrix(1, nrow = length(ALLIDS), ncol = dim(y.ar.ALIVE)[3]+1)
yrs <- c(years[1]-1, years)
dimnames(y.obsALL) <- list(ALLIDS, yrs)
for(i in 1:dim(y.obsALL)[1]){
  for(t in 1:(nYears+1)){
    tmp <- unique(INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$ReprodYear..May.1.year.y...Apr.30.y.1. == yrs[t] & 
                                                 INDIVIDUAL_ID$ROVBASE_IndividID == ALLIDS[i]])
    if(length(tmp) > 0){
      if(length(tmp) > 1){
        print(tmp)
        tmp <- tmp[1]
      }
      y.obsALL[i,t] <- tmp
    }
  }#t
}#i

## FOR THE LAST YEAR GET TRHOUGH THE PAIR BASED-FILE FROM LINN.
for(i in 1:dim(y.obsALL)[1]){
  # for(t in 1:(nYears+1)){
  t <- nYears
  ## linn's 2023 file
  tmp <-  Pack_ID2023[Pack_ID2023$Rovbase.ID %in% ALLIDS[i],]
  if(nrow(tmp) > 0){ y.obsALL[i,t-1] <- 3 }
  ## linn's 2024 file
  tmp <-  Pack_ID2024[Pack_ID2024$RovbaseID %in% ALLIDS[i],]
  if(nrow(tmp) > 0){ y.obsALL[i,t] <- 3 }
  ## linn's 2025 file
  tmp <-  Pack_ID2025[Pack_ID2025$IndividID %in% ALLIDS[i],]
  if(nrow(tmp) > 0){ y.obsALL[i,t+1] <- 3 }
}#i

##-- subset y.obs for the individuals present in y.ar
indID <- unlist(lapply(strsplit(y.ar$Id.vector, " "), function(x)x[1]))
y.obs <- y.obsALL[indID, ]
y.obs[y.obs == 4] <- 3
y.obs <- y.obs[ ,as.character(years)]


##-- KERNEL OF INDIVIDUALS IN PAIRS 
kern <- list()
habDens <- matrix(NA, nrow = numHabWindows, ncol = nYears)
IDS <- unlist(lapply(strsplit(as.character(myFullData.sp$alive$Id), " "), function(x)x[1])) 
for(t in 1:nYears){
  id.fam <- which(y.obsALL[ ,as.character(years[t]-1)] %in% c(3,4), arr.ind = T)
  
  m.xy <- matrix(NA, nrow = length(id.fam), ncol = 2)
  colnames(m.xy) <- c("x","y")
  for(i in 1:length(id.fam)){
    tmp <- myFullData.sp$alive[IDS == row.names(y.obsALL)[i], ]
    m.xy[i, ]<- colMeans(st_coordinates(tmp))
  }#i
  if(sum(is.na(m.xy[ ,1])) > 0){
    m.xy <- m.xy[!is.na(m.xy[ ,1]), ]
  }
  locationsFamily <- st_as_sf( as.data.frame(m.xy),
                               coords = c("x","y"),
                               crs = st_crs(habitat$habitat.sp))
  locationsFamily$id <- rep(1, nrow(locationsFamily))
  kern[[t]] <- raster(estUDm2spixdf(kernelUD( as(locationsFamily[ ,"id"], "Spatial"), 
                                              h = 15000,
                                              grid = as(habitat$habitat.r, 'SpatialPixels'))))
  habDens[ ,t] <- scale(kern[[t]][habitat$habitat.r[ ] == 1])
}

## PLOT CHECK 
for(t in 1:nYears){
  plot(kern[[t]], main = years[t])
  plot(habitat$habitat.poly$geometry, add = T, col = NA)
}#t



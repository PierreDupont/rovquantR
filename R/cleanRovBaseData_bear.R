
cleanRovbaseData_bear <- function( years = NULL, 
                              data_dir = "./Data",
                              output_dir = "./Data",
                              Rmd_template = NULL,
                              overwrite = FALSE)
{

engSpecies <- "bear"
norSpecies <- "Bjørn"

months = c("January","February","March","April","May","June",
           "July","August","September","October","November","December")

##-- Load pre-processed habitat shapefiles
data(COUNTRIESWaterHumans, envir = environment()) 

COUNTRIES <- COUNTRIESWaterHumans[COUNTRIESWaterHumans$ISO %in% c("SWE","NOR"), ] %>%
  dplyr::group_by(ISO) %>%
  dplyr::summarise()


##-- Load the most recent .csv file with the focal species name
DNA <- readMostRecent( 
  path = dir.in,
  extension = ".csv",
  pattern = paste0("dna_",engSpecies,".csv")) %>%
  ##-- Filter out samples from other species
  dplyr::filter(., Species == norSpecies)


##-- Load the most recent .csv dead recovery file
DR <- readMostRecent(
  path = dir.in,
  extension = ".csv",
  pattern = paste0("dead_",engSpecies,".csv")) %>%
  ##-- Filter out dead recoveries of other species
  dplyr::filter(., Species == norSpecies)
  

##-- Merge DNA and dead recoveries files using all shared names columns
DATA <- merge( DR, DNA, 
               by = c("Id","RovbaseID","DNAID","Species", "Sex","Date","East_UTM33","North_UTM33", "County"),
               all = TRUE) %>%
  ##-- Add some columns
  dplyr::mutate( 
    ##-- Add "Country" column
    Country_sample = substrRight(County, 3),
    ##-- Change date format
    Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
    ##-- Extract year
    Year = as.numeric(format(Date,"%Y")),
    ##-- Extract month
    Month = as.numeric(format(Date,"%m")))


##-- Process dates :
##-- For sampling periods spanning over two calendar years (wolf & wolverine)
##-- Set all months in given sampling period to the same year
index <- DATA$Month < unlist(samplingMonths)[1]
index[is.na(index)] <- FALSE
DATA$Year[index] <- DATA$Year[index] - 1

##-- Fix unknown "Id"
DATA$Id[DATA$Id == ""] <- NA

##-- Determine Death and Birth Years
DATA$Age <- suppressWarnings(as.numeric(as.character(DATA$Age))) 
DATA$RovbaseID <- as.character(DATA$RovbaseID)
DATA$Death <- NA
DATA$Death[substr(DATA$RovbaseID,1,1) %in% "M"] <- DATA$Year[substr(DATA$RovbaseID,1,1) %in% "M"]
DATA$Birth <- DATA$Death - DATA$Age

##-- Extract useful numbers
noID <- sum(is.na(DATA$Id))              ## number of samples without ID
noDate <- sum(is.na(DATA$Year))          ## number of samples without Date
noCoords <- sum(is.na(DATA$East_UTM33))  ## number of samples without Coords  

##-- Filter out unusable samples
DATA <- DATA %>%
  dplyr::filter(., 
         ##-- Filter out samples with no ID
         !is.na(Id),
         ##-- Filter out samples with no Coordinates
         !is.na(East_UTM33),
         ##-- Filter out samples with no dates  
         !is.na(Year),
         ##-- Filter out samples with 
         Year %in% years) %>%
  droplevels(.)


ID <- unique(as.character(DATA$Id))
DATA$Sex <- as.character(DATA$Sex)
doubleSexID <- IdDoubleSex <- NULL
counter <- 1
for(i in 1:length(ID)){
  ##-- Subset data to individual i
  tmp <- DATA$Sex[DATA$Id == ID[i]]
  
  ##-- Number of times individual i was assigned to each sex
  tab <- table(tmp[tmp %in% c("Hunn","Hann")])
  
  ##-- If conflicting sexes (ID identified as both "Hunn" and "Hann")
  if(length(tab) == 2){
    ##-- If ID assigned the same number of times to the 2 sexes, assign to Ukjent
    if(tab[1] == tab[2]){
      DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"
    } else {
      ##-- Otherwise pick the most common sex
      DATA$Sex[DATA$Id == ID[i]] <- names(tab)[which(tab == max(tab))]
    }
    # print(paste("Warnings!!!", "Individuals", ID[i], "assigned to both sexes. Now assigned to", names(tab)[which(tab == max(tab))])) 
    IdDoubleSex[counter] <- ID[i]
    counter <- counter + 1
  }
  
  ##-- If only one of "Hunn" or "Hann" registered
  if(length(tab) == 1){DATA$Sex[DATA$Id == ID[i]] <- names(tab)}
  
  ##-- If anything else registered : "Ukjent"
  if(length(tab) == 0){DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"}
  
  doubleSexID[i] <- length(tab)
}#i


##-- Split DATA into alive and dead.recovery datasets
alive <- DATA[is.na(DATA$Death), ]
dead.recovery <- DATA[!is.na(DATA$Death), ]

##-- Add earlier detection index
alive$detected.earlier <-
  unlist(lapply(1:nrow(alive),
                function(i){
                  this.id <- alive[i,"Id"]
                  this.date <- alive[i,"Date"]
                  any(alive$Id %in% this.id & alive$Date < this.date)
                }))

dead.recovery$detected.earlier <-
  unlist(lapply(1:nrow(dead.recovery),
                function(i){
                  this.id <- dead.recovery[i,"Id"]
                  this.date <- dead.recovery[i,"Date"]
                  any(alive$Id %in% this.id & alive$Date < this.date)
                }))


##-- Load most recent "flagged" file from HB
flagged <- readMostRecent( 
  path = dir.in,
  extension = ".csv",
  pattern = "dna_bear_to_remove", 
  fileEncoding = "Latin1") 

##-- Remove flagged samples 
remove.alive <- !alive$Barcode_sample %in% flagged$Strekkode
alive <- alive[remove.alive, ]
remove.dead <- !dead.recovery$Barcode_sample %in% flagged$Strekkode
dead.recovery <- dead.recovery[remove.dead, ]

dead.recovery$Missing <- NA
dead.recovery$Individ <- NA



##-- Turn into sf points dataframe
alive <- sf::st_as_sf( x = alive,
                       coords = c("East_UTM33","North_UTM33")) %>%
  sf::st_set_crs(.,sf::st_crs(32633)) 

##-- Intersect and extract country name
alive$Country_sf <- COUNTRIES$ISO[as.numeric(st_intersects(alive, COUNTRIES))]


##-- Turn into sf points dataframe
dead.recovery <- sf::st_as_sf( x = dead.recovery,
                               coords = c("East_UTM33","North_UTM33")) %>%
  sf::st_set_crs(.,sf::st_crs(32633))

##-- Intersect and extract country name
dead.recovery$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(dead.recovery, COUNTRIES))]



##-- Number of NGS samples
samples <- table(alive$Country_sample, alive$Year)
samples2 <- table(alive$Country_sf, alive$Year)
samples <- rbind(samples, "Total" = colSums(samples))
samples <- cbind(samples, "Total" = rowSums(samples))


##-- Number of individuals detected alive
ids <- apply(table(alive$Country_sample, alive$Year, alive$Id),
             c(1,2),
             function(x)sum(x>0))
ids <- rbind(ids,
             "Total" = apply(table(alive$Year,alive$Id),
                             1,
                             function(x)sum(x>0)))
ids <- cbind(ids,
             "Total" = c(apply(table(alive$Country_sample,alive$Id),
                               1,
                               function(x)sum(x>0)),length(unique(alive$Id))))


##-- Number of DR samples
deadSamples <- table(dead.recovery$Country_sample, dead.recovery$Year)
deadSamples <- rbind(deadSamples, "Total" = colSums(deadSamples))
deadSamples <- cbind(deadSamples, "Total" = rowSums(deadSamples))


##-- Number of individuals recovered
deadIds <- apply(table(dead.recovery$Country_sample, dead.recovery$Year, dead.recovery$Id),
                 c(1,2),
                 function(x)sum(x>0))
deadIds <- rbind(deadIds,
                 "Total" = apply(table(dead.recovery$Year,dead.recovery$Id),
                                 1,
                                 function(x)sum(x>0)))
deadIds <- cbind(deadIds,
                 "Total" = c(apply(table(dead.recovery$Country_sample,dead.recovery$Id),
                                   1,
                                   function(x)sum(x>0)),length(unique(dead.recovery$Id))))


# The first steps of the cleaning process consist in :
#   
# -   translating Scandinavian characters  
# -   renaming columns to match the NGS and dead recovery files  
# -   merging both files into one  
# -   extracting year and month from POSIX dates  
# -   filtering out unusable samples :  
  
#-   samples without ID 
noID
#-   samples without spatial coordinates 
noCoords
#-   samples without dates 
noDate



paste0("-   removing ", sum(!remove.alive), " bear DNA samples and ", sum(!remove.dead), " dead recoveries flagged by H.Brøseth")


## Number of samples

#After this initial clean-up, we are left with 
paste0(samples["Total","Total"],"(NOR = ",samples["(N)","Total"]," ; SWE = ",samples["(S)","Total"], ")")
#NGS samples and `r paste0(deadSamples["Total","Total"], "(NOR = ",deadSamples["(N)","Total"]," ; SWE = ",deadSamples["(S)","Total"], ")")` dead recoveries.

##-- Number of samples
dat.alive <- alive %>%
  dplyr::group_by(Date) %>%
  dplyr::summarise(n = n())
dat.alive$type = "NGS"

dat.dead <- dead.recovery %>%
  dplyr::group_by(Date) %>%
  dplyr::summarise(n = -n())
dat.dead$type = "dead.recovery"

dat <- rbind(dat.alive, dat.dead)

ggplot(dat) +
  geom_col(aes(x = Date, y = n, fill = type)) +
  ylab("Number of samples") +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank(),
        legend.position.inside = c(0.1,0.9)) 

##-- Number of NGS samples
kable(samples, align = "lc",
      caption = "Number of NGS samples per year and country") %>%
  kable_styling(full_width = F)

##-- Number of dead recoveries
kable(deadSamples, align = "lc",
      caption = "Number of DNA samples from dead animals per year and country") %>% kable_styling(full_width = F)


## Number of individuals

#In terms of individuals identified, these correspond to  
paste0(ids["Total","Total"], "(NOR = ",ids["(N)","Total"]," ; SWE = ",ids["(S)","Total"], ")")
#individuals detected alive and 
paste0(deadIds["Total","Total"], "(NOR = ",deadIds["(N)","Total"]," ; SWE = ",deadIds["(S)","Total"], ")")` 
#individuals recovered.

##-- Number of IDs
dat.alive <- alive %>% 
  dplyr::group_by(Year) %>% 
  dplyr::summarise(n = length(unique(Id)))
dat.alive$type="NGS"

dat.dead <- dead.recovery %>% 
  dplyr::group_by(Year) %>% 
  dplyr::summarise(n = -length(unique(Id)))
dat.dead$type = "dead.recovery"

dat <- rbind(dat.alive,dat.dead)

ggplot(dat) +
  geom_col(aes(x = Year, y = n, fill = type)) +
  ylab("Number of individuals") +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.1,0.9)) +
  scale_x_continuous(breaks = years, labels = years)


##-- Number of ID detected alive
kable(ids, align = "lc",
      caption = "Number of individuals detected through NGS per year and country") %>%
  kable_styling(full_width = F)

##-- Number of ID detected alive
kable(deadIds, align = "lc",
      caption = "Number of identified dead animals per year and country") %>% 
  kable_styling(full_width = F)

dat.alive <- alive %>%
  dplyr::group_by(Year, detected.earlier) %>%
  dplyr::summarise(n = length(unique(Id)))
dat.alive$type = "NGS"

dat.dead <- dead.recovery %>%
  dplyr::group_by(Year, detected.earlier) %>%
  dplyr::summarise(n = length(unique(Id)))
dat.dead$type = "dead.recovery"

ggplot() +
  geom_col(data = dat.alive,
           aes(x = Year, y = n, fill = detected.earlier)) +
  ylab("Number of individuals detected through NGS") +
  theme(legend.position = c(0.1,0.9)) +
  scale_x_continuous(breaks = years, labels = years)

ggplot() +
  geom_col(data = dat.dead,
           aes(x = Year, y = n, fill = detected.earlier)) +
  ylab("Number of dead individuals") +
  theme(legend.position = c(0.1,0.9)) +
  scale_x_continuous(breaks = years, labels = years)




------------------------------------------------------------------------
  
# Issues
#We can now start digging in the data, looking for potential issues.


## Sex assignment

sexTab <- cbind.data.frame(
  "problems" = c("Unknown sex", "both 'Hunn' and 'Hann'"),
  "number of individuals" = as.numeric(table(doubleSexID)[c(1,3)]))

kable(sexTab, align = "lc") %>%
  kable_styling(full_width = F)


## Multiple deaths

##-- Identify and count individuals dead "more than once"
ID <- names(table(dead.recovery$Id))[table(dead.recovery$Id)>1]
multiDeathDate <- multiDeathYear <- multiDeathLocs <-  NULL
for(i in 1:length(ID)){
  tmp <- dead.recovery[dead.recovery$Id == ID[i], ] 
  ##-- Multiple death dates
  if(length(unique(tmp$Date)) > 1){
    multiDeathDate <- c(multiDeathDate, ID[i])
  }
  ##-- Multiple death years
  if(length(unique(tmp$Year)) > 1){
    multiDeathYear <- c(multiDeathYear, ID[i])
  }
  ##-- Multiple death locations
  if(length(unique(tmp$East)) > 1 | length(unique(tmp$North)) > 1){
    multiDeathLocs <- c(multiDeathLocs, ID[i])
  }
}#i

##-- Remove individuals that died more than once
dead.recovery$Id <- as.character(dead.recovery$Id)
IdDoubleDead <- names(table(dead.recovery$Id))[table(dead.recovery$Id) > 1]
if(length(IdDoubleDead) > 0){
  for(i in IdDoubleDead){
    ##-- Identify repeated deaths
    tmp <- which(dead.recovery$Id %in% i) 
    ##-- Try to keep death with known death cause 
    tmp2 <- which(!is.na(dead.recovery$DeathCause_2[tmp]))[1]
    if(length(tmp2) == 0){tmp <- tmp[-1]} else {tmp <- tmp[!tmp %in% tmp2]}
    ##-- Remove repeated deaths
    dead.recovery <- dead.recovery[-tmp, ]
  }#i
}#if
# dead.recovery <- droplevels(dead.recovery)

#There are 
length(multiDeathDate)
#individuals with multiple death dates, of which 
length(multiDeathYear)
#are recorded dead in different years.

#There are also 
length(multiDeathLocs)
#recorded with different death locations.

## Detections after death

id.list <- unique(c(as.character(dead.recovery$Id), as.character(alive$Id)))
ghosts <- unlist(lapply(id.list, function(id){
  out <- NULL
  try({
    if(id %in% dead.recovery$Id){
      mort.year <- min(dead.recovery$Year[dead.recovery$Id == id])
      this.alive <- alive[alive$Id == id, ]
      ## Was it detected alive in any season after death?
      temp <- this.alive[this.alive$Year > mort.year, ]
      if(length(temp) > 0){
        out <- rownames(temp)
        names(out) <- id
      }## FLAG THOSE FOR HENDRIK
    }
  }, silent = TRUE)
  return(out)
}))
samples.to.remove <- unlist(ghosts)

##-- Remove flagged NGS detections after dead recovery
alive <- alive[!rownames(alive) %in% samples.to.remove, ]


#There are 
length(ghosts)
#individuals identified with NGS samples detected after their supposed death




##------------------------------------------------------------------------
  
# Maps
# We can also add maps of the NGS samples collected year:
  
##-- NGS map
numRows <- ceiling(length(years)/5)
numCols <- 5
NGS_map <- ggplot(data = alive) +
  geom_sf(data = COUNTRIES,
          aes(fill = ISO),
          alpha = 0.3,
          color = NA) +
  geom_sf(color = "black",
          alpha = 0.3, size = 0.8, pch = 3) +
  facet_wrap(~Year, nrow = numRows, ncol = numCols) +
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
NGS_map

##-- Save maps as .png
grDevices::png(filename = file.path(dir.out, 
                                    paste0(species, "_NGS_", years[1]," to ", years[length(years)], ".png")),
               width = 8, height = 6, units = "in", res = 300)
NGS_map
graphics.off()

#and a series of maps for the dead recoveries each year:
  
dead_map <- ggplot(data = dead.recovery) +
  geom_sf(data = COUNTRIES, 
          aes(fill = ISO),
          alpha = 0.3,
          color = NA) + 
  geom_sf(color = "black", alpha = 0.5, size = 0.8, pch = 3) +
  facet_wrap(~Year, nrow = numRows, ncol = numCols) +
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
dead_map

##-- Save maps as .png
grDevices::png(filename = file.path(dir.out, 
                                    paste0(engSpecies, "_DEAD_", years[1]," to ", years[length(years)], ".png")),
               width = 8, height = 6, units = "in", res = 300)
dead_map
graphics.off()




##------------------------------------------------------------------------
  
  # Save clean data
  
fileName <-  paste0("Data_", engSpecies, "_",params$modDate,".RData")

save(alive, 
     dead.recovery,
     IdDoubleSex,
     file = file.path(dir.out, fileName))

# Finally, we save the cleaned **alive** and **dead.recovery** sf objects as a .RData file with name **`r fileName`** located in the `r engSpecies`-specific folder (`r dir.out`).

}



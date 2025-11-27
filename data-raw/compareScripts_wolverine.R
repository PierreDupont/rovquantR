##===== LOAD RAW DATA ====

DNA <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/dna_wolverines.csv"),
                 fileEncoding = "latin1") 
colnames(DNA) <- translateForeignCharacters(dat=colnames(DNA), dir.translation = dir.analysis)
DNA <- DNA[,-which(colnames(DNA)%in% "Kjoenn..Individ.")]
dim(DNA)

## Remove DEAD entries from the DNA data [HB]
DNA <- DNA[substr(DNA$RovbaseID..Proeve.,1,1) != "M", ]
dim(DNA)


DEAD <- read.csv( file.path(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/dead_carnivores.csv"),
                  fileEncoding = "latin1")
colnames(DEAD) <- translateForeignCharacters(dat=colnames(DEAD), dir.translation = dir.analysis)
dim(DEAD)
## Remove un-verified dead recoveries [HB]
## ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
DEAD <- DEAD[!grepl(pattern = "Påskutt", x = as.character(DEAD$Utfall)), ]
dim(DEAD)



### ====    1.1. CLEAN NGS & DEAD RECOVERY DATA ====

# myCleanedData.sp <- CleanDataNew2sf( 
#   dna_samples = DNA,
#   dead_recoveries = DEAD,
#   species_id = myVars$DATA$species,
#   country_polygon = COUNTRIES,
#   threshold_month = unlist(myVars$DATA$samplingMonths)[1],
#   keep_dead = T,
#   age.label.lookup = age.lookup.table)

## Merge DNA and dead recoveries data## Rename Dead Recoveries
names(DEAD)[grep(pattern = "Individ",x = names(DEAD))] <- "Id"
names(DEAD)[grep(pattern = "RovbaseID",x = names(DEAD),fixed = TRUE)] <- "RovBaseId"
names(DEAD)[grep(pattern = "DNAID..Proeve",x = names(DEAD),fixed = TRUE)] <- "DNAID"
names(DEAD)[grep(pattern = "Art",x = names(DEAD))] <- "Species"
names(DEAD)[grep(pattern = "Kjoenn",x = names(DEAD))] <- "Sex"
names(DEAD)[grep(pattern = "Doedsdato",x = names(DEAD))] <- "Date"
names(DEAD)[grep(pattern = "Oest..UTM33",x = names(DEAD))] <- "East"
names(DEAD)[grep(pattern = "Nord..UTM33",x = names(DEAD))] <- "North"
names(DEAD)[grep(pattern = "^Alder..verifisert$",x = names(DEAD))] <- "Age"
names(DEAD)[grep(pattern = "Bakgrunn.arsak.metode",x = names(DEAD))] <- "DeathCause_2"
names(DEAD)[grep(pattern = "Bakgrunn.arsak.formal",x = names(DEAD))] <- "DeathCause_3"
names(DEAD)[grep(pattern = "Bakgrunn.arsak",x = names(DEAD))] <- "DeathCause"

## Rename DNA Samples
names(DNA)[grep(pattern = "Individ",x = names(DNA))] <- "Id"
names(DNA)[grep(pattern = "RovbaseID..Proeve",x = names(DNA))] <- "RovBaseId"
names(DNA)[grep(pattern = "DNAID..Proeve",x = names(DNA))] <- "DNAID"
names(DNA)[grep(pattern = "Art..Analyse",x = names(DNA))] <- "Species"
names(DNA)[grep(pattern = "Kjoenn",x = names(DNA))] <- "Sex"
names(DNA)[grep(pattern = "Funnetdato",x = names(DNA))] <- "Date"
names(DNA)[grep(pattern = "Oest..UTM33",x = names(DNA))] <- "East"
names(DNA)[grep(pattern = "Nord..UTM33",x = names(DNA))] <- "North"
names(DNA)[grep(pattern = "DNA.proeveleverandoer",x = names(DNA))] <- "Origin"

## MERGE DATA
DATA <- merge(DEAD, DNA, by=c("Id","RovBaseId","DNAID","Species","Sex","Date","East","North"), all = TRUE)

## Clean the data
DATA <- DATA[DATA$Id!="",]                            ## Delete unknown individuals
DATA <- DATA[!is.na(DATA$Id), ]                       ## Delete NA individuals
DATA <- DATA[!is.na(DATA$East), ]                     ## Delete NA locations
DATA <- DATA[!is.na(DATA$Date), ]                     ## Delete NA dates
DATA <- DATA[DATA$Species %in% myVars$DATA$species, ] ## Delete unwanted species

## Convert dates to biological years
DATA$Date <- as.POSIXct(strptime(DATA$Date, "%d.%m.%Y"))
DATA$Year <- as.numeric(format(DATA$Date,"%Y"))
DATA$Month <- as.numeric(format(DATA$Date,"%m"))
DATA <- DATA[!is.na(DATA$Year), ]                     ## Delete NA dates

## GET THE SEASON YEAR
DATA$Year[DATA$Month<unlist(myVars$DATA$samplingMonths)[1]] <- DATA$Year[DATA$Month<unlist(myVars$DATA$samplingMonths)[1]] - 1

## Determine Death and Birth Years
DATA$Age <- suppressWarnings(as.numeric(as.character(DATA$Age))) 
DATA$RovBaseId <- as.character(DATA$RovBaseId)
DATA$Death <- NA
DATA$Death[substr(DATA$RovBaseId,1,1)=="M"] <- DATA$Year[substr(DATA$RovBaseId,1,1)=="M"]
DATA$Birth <- DATA$Death-DATA$Age

## Reconstruct minimal & maximal ages
DATA$Age.orig <- DATA$Age
# if(!is.null(age.label.lookup)){
#   temp <- temp1 <- as.character(levels(DATA$Age.orig))  ## list age levels
#   temp <- toupper(temp)                                 ## Upper case all
#   temp <- gsub("\\s", "", temp)                         ## Remove blank spaces
#   DATA$Age.orig2 <- DATA$Age.orig
#   levels(DATA$Age.orig2) <- temp
#   DATA <- merge( DATA, age.label.lookup[ ,-1],
#                  by.x = "Age.orig2",
#                  by.y = "age.label",
#                  all.x = TRUE)                          ## Merge with info from lookup table
#   
#   ## FILL IN THE REST OF THE AGES FROM FOR NUMERIC RECORDS
#   numeric.age.records <- which(!is.na(as.numeric(as.character(DATA$Age.orig2))) & !is.na(DATA$Age.orig2))
#   DATA[numeric.age.records, c("min.age","max.age","age")] <- floor(as.numeric(as.character(DATA$Age.orig2[numeric.age.records])))
# }
dim(DATA)
table(DATA$Age, useNA = "always")

## Convert samples coordinates to the correct spatial projection
DATA <- st_as_sf(DATA, coords = c("East", "North"))
st_crs(DATA) <- st_crs(COUNTRIES)
# Overlay with SpatialPolygons to determine the countries 
if(!is.null(COUNTRIES)){
  DATA$Country <- NA
  DATA$Country[!is.na(as.numeric(st_intersects(DATA, COUNTRIES[which(COUNTRIES$ISO %in% c("FIN")),] )))] <- "F"
  DATA$Country[!is.na(as.numeric(st_intersects(DATA, COUNTRIES[which(COUNTRIES$ISO %in% c("RUS")),] )))] <- "R"
  DATA$Country[!is.na(as.numeric(st_intersects(DATA, COUNTRIES[which(COUNTRIES$ISO %in% c("GOT")),] )))] <- "G"
  DATA$Country[!is.na(as.numeric(st_intersects(DATA, COUNTRIES[which(COUNTRIES$ISO %in% c("NOR")),] )))] <- "N"
  DATA$Country[!is.na(as.numeric(st_intersects(DATA, COUNTRIES[which(COUNTRIES$ISO %in% c("SWE")),] )))] <- "S"
}#if

## the "Factor" trick...needed to suppress unused factor levels in Id
DATA$Id <- factor(as.character(DATA$Id), levels = unique(as.character(DATA$Id)))
dim(DATA)

DATA$Country_sample <- substrRight(DATA$Fylke.y, 3)
DATA$Country_sample2 <- substrRight(DATA$Fylke.x, 3)
DATA$Country_sample[is.na(DATA$Country_sample)] <- DATA$Country_sample2[is.na(DATA$Country_sample)] 
table(DATA$Country_sample,useNA="always")



### ====    1.2. FILTER DATA ====

# myFullData.sp <- FilterDatasf(
#   myData = myCleanedData.sp,
#   poly = myStudyArea,
#   dead.recovery = T ,
#   sex = c("Hann","Hunn"), # do the sex selection at the last moment
#   setSex = T)

# SET THE SEX OF INDIVIDUALS BASED ON ALL INFORMATION AVAILABLE
table(DATA$Sex)


ID <- unique(as.character(DATA$Id))
# initialize the vector of IDs with conflicting sexes
IdDoubleSex <- 0
counter <- 1
for(i in 1:length(ID)){
  # subset data to individual i
  tmp <- DATA$Sex[DATA$Id == ID[i]] 
  # create a table of the number of times individual i was assigned to each sex
  tab <- table(tmp[tmp %in% c("Hunn","Hann")])
  # If conflicting sexes (ID identified as both "Hunn" and "Hann")
  
  if(length(tab) == 2){
    # If ID assigned the same number of times to the 2 sexes, assign to Ukjent
    if(tab[1] == tab[2]){
      DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"
    } else {
      # Otherwise pick the most common sex
      DATA$Sex[DATA$Id == ID[i]] <- names(tab)[which(tab == max(tab))]
    }
    # In any case, print a warning
    print(paste("Warnings!!!", "Individuals", ID[i], "assigned to both sexes. Now assigned to", names(tab)[which(tab == max(tab))])) 
    IdDoubleSex[counter] <- ID[i]
    counter <- counter + 1
  }
  
  # If only one of "Hunn" or "Hann" registered
  if(length(tab) == 1){DATA$Sex[DATA$Id == ID[i]] <- names(tab)}
  
  # If anything else registered : "Ukjent"
  if(length(tab) == 0){DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"}
}#i
table(DATA$Sex)


## Remove all samples outside the polygon of interest
if(!is.null(poly)){DATA <- DATA[!is.na(as.numeric(st_intersects(DATA, myStudyArea))), ]}
#if(!is.null(poly)){DATA <- DATA[!is.na(over(DATA, as(poly,"SpatialPolygons"))), ]}


dead.recovery <- DATA[!is.na(DATA$Death), ]
alive <- DATA[is.na(DATA$Death), ]

dead.recovery$Id <- droplevels(dead.recovery$Id)
alive$Id <- droplevels(alive$Id)

IdDoubleDead <- dead.recovery$Id[duplicated(dead.recovery$Id)]


## REMOVE SUSPECT SAMPLES ACCORDING TO HENRIK
alive$DNAID <- as.character(alive$DNAID)
dead.recovery$DNAID <- as.character(dead.recovery$DNAID)
dim(alive)
dim(dead.recovery)

alive <- alive[!(alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
dead.recovery <- dead.recovery[!(dead.recovery$RovBaseId %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]
dim(alive)
dim(dead.recovery)

##dead recoveries flagged by Henrik that should always be removed (email from the 18/12/2024)
dead.recovery <- dead.recovery[!dead.recovery$RovBaseId %in% c("M495994","M524051","M524052","M524053"), ]
dim(dead.recovery)


## Remove individuals that died twice# [CM] TO BE CHECKED BECAUSE "length(IdDoubleDead) < 0" and so it was desactivated
dead.recovery$Id <- as.character(dead.recovery$Id)
IdDoubleDead <- dead.recovery$Id[duplicated(dead.recovery$Id)]

if(length(IdDoubleDead) > 0){
  duplicatedDeath <- NULL
  for(i in IdDoubleDead){
    tmp  <- which(dead.recovery$Id == i & is.na(dead.recovery$DeathCause_2))
    if(length(tmp)==0){tmp  <- which(dead.recovery$Id == i)[-2]}##[CM] remove the second record.
    duplicatedDeath <- c(duplicatedDeath, tmp)
  }#i
  dead.recovery <- dead.recovery[-duplicatedDeath, ]
}#if

unique(dead.recovery$DeathCause)
unique(dead.recovery$DeathCause_2)



## Remove pups killed before recruitment based on weight (cf. Henrik)
## 1) remove individuals that are "Ja" in column "Doedt.individ..Unge" and recovered dead between March and November
sum(dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
      dead.recovery$Month > 2 &
      dead.recovery$Month < 12)

dead.recovery <- dead.recovery[-which(dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
                                                                    dead.recovery$Month > 2 &
                                                                    dead.recovery$Month < 12),]
dim(dead.recovery)

# 2) remove individuals that have a weight >0 and <4 between March and November
# format the weight correctly
dead.recovery$Helvekt <- as.character(dead.recovery$Helvekt)
dead.recovery$Slaktevekt <- as.character(dead.recovery$Slaktevekt)

#convert to decimals
dead.recovery$Helvekt <- as.numeric(gsub(",", ".", dead.recovery$Helvekt))
dead.recovery$Slaktevekt <- as.numeric(gsub(",", ".", dead.recovery$Slaktevekt))
# get the two weight columns together.
dead.recovery$weight <- ifelse(!is.na(dead.recovery$Helvekt),
                                             dead.recovery$Helvekt,
                                             dead.recovery$Slaktevekt)
#assign negative values to nas to avoid issues
dead.recovery$weight[is.na(dead.recovery$weight)] <- -999
dim(dead.recovery)

# check with Henrik
# this step does not remove dead recoveries on id with weight==0 should it?
# WEIGTH DISTRIBUTION
par(mfrow=c(4,3))
for(t in 1:12){
  hist(dead.recovery$weight[(dead.recovery$weight >-1) &
                                            dead.recovery$Month%in% t],
       breaks=c(0:20), main=t,xlab="Weight")
}
# AGE DISTRIBUTION
par(mfrow=c(4,3))
for(t in 1:12){
  hist(dead.recovery$Age[(dead.recovery$Age >-1) &
                                         dead.recovery$Month%in% t],breaks=seq(-0.01,20.99,by=1),
       main=t,xlab="Age")
}

# check how many dead reco we remove and remove if more than 0
if(sum(dead.recovery$weight > 0 &
       dead.recovery$weight < 4 &
       dead.recovery$Month < 12 &
       dead.recovery$Month > 2)>0){
  dead.recovery <- dead.recovery[-which(dead.recovery$weight > 0 &
                                                                      dead.recovery$weight < 4 &
                                                                      dead.recovery$Month < 12 &
                                                                      dead.recovery$Month > 2),]
}
dim(dead.recovery)

# check how many dead reco with a weight of 0 kg and recovered between march and november
if(sum(dead.recovery$Age %in% 0 &
       dead.recovery$Month < 12 &
       dead.recovery$Month > 2)>0){
  dead.recovery[dead.recovery$Age %in% 0 &
                                dead.recovery$Month < 12 &
                                dead.recovery$Month > 2,  ]
}


dim(alive)
table(alive$Year)

dim(dead.recovery)
table(dead.recovery$Year)

table(dead.recovery$DeathCause)
table(dead.recovery$DeathCause_2)

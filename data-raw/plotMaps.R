library(dplyr)
library(colorspace)

###
#NewCountySwe <- readOGR(file.path(dir.dropbox,"/DATA/GISData/scandinavian_border/rk_lan_07_WGS84.shp"))
NewCountySwe <- st_read(file.path(dir.dropbox,"/DATA/GISData/scandinavian_border/rk_lan_07_WGS84.shp"))
NewCountySwe1 <- st_read(file.path(dir.dropbox,"/DATA/GISData/new_scandinavian_border/alla_lan.shp"))

NorwayManagementRegion <- st_read(file.path(dir.dropbox,"/DATA/GISData/NorwegianManagementRegions/rovviltregioner2024.shp"))


plot(myHabitat.list$habitat.r, legend=F)
plot(NewCountySwe$geometry, add=T, border="red")

par(mar=c(0,0,0,0))
plot(NewCountySwe$geometry, col=as.factor(NewCountySwe$LANSKOD))
plot(NewCountySwe1$geometry,border="red",col=NA,add=T,lwd=2)
## only select the norwegian counties
COMMUNESNOR <- COMMUNES[COMMUNES$NAME_0=="Norway",]
NORWAY <- aggregate(x = COMMUNESNOR, by = "NAME_1")
NORWAY <- COMMUNESNOR %>% group_by(NAME_1) %>% summarize()

plot(NORWAY$geometry,col="red")
NORWAY1 <- st_read(file.path(dir.dropbox,"/DATA/GISData/new_scandinavian_border/fylker-2024.shp"))
plot(NORWAY1$geometry,add=T,border="blue")
plot(NorwayManagementRegion$geometry,add=T,border="green")
habitatRasters$Regions


NAME_1 <- as.character(NORWAY$NAME_1)
df.CountiesRegions <- matrix(c(
  "Finnmark",         8,
  "Troms",            8,
  "Nordland",         7,
  NAME_1[14],         6,
  NAME_1[9],         6,
  NAME_1[8],         6,
  "Hedmark",          5,
  "Oppland",          3,
  NAME_1[1],          4,
  "Oslo",             4,
  "Akershus",         4,
  "Sogn og Fjordane", 1,
  "Hordaland",        1,
  "Rogaland" ,        1,
  "Vest-Agder",       1,
  "Aust-Agder",       2,
  "Telemark" ,        2,
  "Buskerud" ,        2,
  "Vestfold",         2), byrow=T,ncol=2)


NORWAY$NAME_1 <- as.character(NORWAY$NAME_1) 
for(i in 1:nrow(df.CountiesRegions)){
  NORWAY$NAME_1[NORWAY$NAME_1 %in% df.CountiesRegions[i,1]] <- df.CountiesRegions[i,2]
}
NORWAY1 <- aggregate(NORWAY,by="NAME_1")
plot(NORWAY1)



#RENAME THE FIELDS SO THEY MATCH BETWEEN NORWEGIAN AND SWEDISH LAYERS
NewCountySwe <- NewCountySwe[,"LANSNAMN"]
colnames(NewCountySwe@data) <- "NAME_1"
NewCountySwe$Country <- "SWE"
NORWAY1$Country <- "NOR"

# MERGE THE 2 LAYERS.
# THERE IS SOME SPACE BETWEEN THE TWO LAYERS, BUT IT DOESNT MATTER. 
# IF A CELL IS NOT ASSIGNED TO ANY COUNTY THEN WE ASSIGN IT THE CLOSEST COUNTY BELOW. 
COUNTIESsimp <- rbind(NORWAY1,NewCountySwe)#, NewCountySwe, makeUniqueIDs = TRUE) 



country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
col <- c("firebrick2", "deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway", "Sweden")
border.col <- NA

pdf(file = file.path(myVars$WDTables, "RegionMaps.pdf"),
    width = 9, height = 13, pointsize = 12)
par(mar=c(0,0,0,0))
plot(gSimplify(COUNTIESsimp, tol=5000), col=NA, border=NA)

#NORWAY
# CARNIVORE.REGIONS <- COUNTIESsimpCarnivoreRegions[COUNTIESsimpCarnivoreRegions$Country%in% "NOR",]
# CARNIVORE.REGIONS1 <- gSimplify(CARNIVORE.REGIONS, tol=200, topologyPreserve = T)
# CARNIVORE.REGIONS1$Region <- CARNIVORE.REGIONS$NAME_1
# region <- gUnaryUnion(CARNIVORE.REGIONS1, id = CARNIVORE.REGIONS1$Region)
region <- gSimplify(NORWAY1, tol=500)


#---SPECIFY COLOR PALETTE

this.col <- sequential_hcl(1+length(unique(NORWAY1$NAME_1)), "Reds 3")
this.col <- this.col[-length(this.col)]
set.seed(100)
this.col <- sample(this.col)
plot(region,add=T, col=this.col, border=border.col, lwd=1)


# SWEDEN
#swedenCounties1 <- COUNTIESsimpCarnivoreRegions[COUNTIESsimpCarnivoreRegions$Country%in% "SWE",]
NewCountySwe$NAME_1 <- as.character(NewCountySwe$NAME_1)
#swedenCounties1 <- swedenCounties1[-which(swedenCounties1$NAME_1%in%"Gotlands län"),]

NewCountySwe1 <- gSimplify(NewCountySwe,tol=200,topologyPreserve = T)
swedenCounties2 <- NewCountySwe1
swedenCounties2$NAME_1 <- as.character(NewCountySwe$NAME_1)

swedenCounties2 <- swedenCounties2[!(swedenCounties2$NAME_1%in% "Gotlands lÃ¤n"), ]

#---SPECIFY COLOR PALETTE
this.col <- sequential_hcl(25+length(unique(swedenCounties2$NAME_1)), "Blues 3")
this.col <- this.col[1:length(unique(swedenCounties2$NAME_1))]
set.seed(100)
this.col <- sample(this.col)

plot(RemoveHolesSp(swedenCounties2),add=T, col=this.col, border=border.col, lwd=1)

#---LABELS
CARNIVORE.REGIONS2 <- aggregate(CARNIVORE.REGIONS1, by="Region")
#raster::text(CARNIVORE.REGIONS2,labels = paste("Region ",CARNIVORE.REGIONS2$Region,sep=""),col=grey(0.05),cex=1.3)
raster::text(NORWAY1,labels = NORWAY1$NAME_1, cex=1.3,
             col=ifelse(NORWAY1$NAME_1 %in% c(5), grey(0.7), grey(0)))

swedenCounties2$NAME_1 <- as.character(unlist(lapply( strsplit(swedenCounties2$NAME_1, " "),
                                                      function(x) x[1])))
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="GÃ¤vleborgs"] <- "Gävleborg"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤rmlands"] <- "Värmland"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Ã–stergÃ¶tlands"] <- "Östergötland"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="JÃ¤mtlands"] <- "Jämtland"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="JÃ¶nkÃ¶pings"] <- "Jönköping" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="SkÃ¥ne"] <- "Skåne" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤stmanlands"] <- "Västmanland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤stra"] <- "Västra Götaland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤sternorrlands"] <- "Västernorrland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤sterbottens"] <- "Västerbotten" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="SÃ¶dermanlands"] <- "Södermanland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Ã–rebro"] <- "Örebro" 

swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Norrbottens"] <- "Norrbotten" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Hallands"] <- "Halland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Kronobergs"] <- "Kronoberg" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Stockholms"] <- "Stockholm" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Dalarnas"] <- "Dalarna" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Dalarnas"] <- "Dalarna" 

####
swedenCounties2$Region <- c("Mellersta",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Norra",
                            "Norra",
                            "Norra",
                            "Norra",
                            "Mellersta")
swedenCounties2Regions <- aggregate(swedenCounties2,by="Region")

plot(RemoveHolesSp(swedenCounties2Regions),add=T, border=grey(0.0), lwd=3)

polygonsLabel(swedenCounties2, swedenCounties2$NAME_1,
              method = "buffer", cex=1.1,#gridpoints = 50000,
              col=ifelse(swedenCounties2$NAME_1%in%c("Jämtland"),grey(0.7),grey(0.7)))

dev.off()


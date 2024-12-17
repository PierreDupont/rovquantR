##-- Plots

## Weigth distribution
par(mfrow = c(4,3))
for(t in 1:12){
  hist(myFullData.sp$dead.recovery$weight[(myFullData.sp$dead.recovery$weight > -1) &
                                            myFullData.sp$dead.recovery$Month %in% t],
       breaks = c(0:30), main = t, xlab = "Weight")
}


## AGE DISTRIBUTION
par(mfrow = c(4,3))
for(t in 1:12){
  hist(myFullData.sp$dead.recovery$Age[(myFullData.sp$dead.recovery$Age >-1) &
                                         myFullData.sp$dead.recovery$Month%in% t],
       breaks = seq(-0.01,20.99, by = 1), main = t, xlab = "Age")
}


## Plot check filtered data
for(t in 1:nYears){
  plot( st_geometry(studyArea))
  plot( st_geometry(COUNTIESNorrbotten), add = T, col = "blue")
  plot( st_geometry(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t], ]),
        col = "red", add = T, pch = 16)
}


par(mfrow = c(1,3))
for(t in 1:nYears){
  ## DEAD RECOVERIES TOTAL
  tempTotal <-  myFilteredData.sp$dead.recovery[ myFilteredData.sp$dead.recovery$Year == years[t] & myFilteredData.sp$dead.recovery$Sex %in% DATA$sex, ]
  NGS_TabTotal <- table(tempTotal$Country)
  ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country), 2, function(x) sum(x>0))
  ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD
  ## PLOT NGS SAMPLES
  plot(st_geometry(GLOBALMAP), col="gray80")
  plot(st_geometry(studyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
  plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
  plot(st_geometry(tempTotal), pch = 21, bg = "darkred",add=T)
  ## ADD NUMBER OF NGS samples and IDs per COUNTRY
  graphics::text(x = 100000, y = 7250000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
  graphics::text(x = 820000, y = 6820000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
  ## ADD OVERALL NUMBERS
  mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
}#t



## ------      1.4.3. PLOT CHECKS 
pdf(file = paste(WD,"/",modelName,"/",modelName,"DetectionsStructuredOppBarplot",".pdf", sep="" ))
par(mfrow=c(2,1),mar=c(4,4,3,2))
barplot(rbind(table(myFilteredData.spStructured$Year),
              table(myFilteredData.spOthers$Year)),beside=T,ylim=c(0,2000),col=c(grey(0.2),grey(0.8)),ylab="Number of samples")
abline(h=seq(0,2000,by=500),lty=2,col=grey(0.8))
title(main="500m threshold")
legend("topleft",fill=c(grey(0.2),grey(0.8)),legend=c("Structured","Other"))

###
distanceThreshold1 <- 2000
whichStructured2000 <- myFilteredData.sp$alive$Proeveleverandoer %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
  !is.na(myFilteredData.sp$alive$TrackRovbsID) &
  myFilteredData.sp$alive$TrackDist <= distanceThreshold1
myFilteredData.spStructured2000 <- myFilteredData.sp$alive[whichStructured2000,]
myFilteredData.spOthers2000 <- myFilteredData.sp$alive[!whichStructured2000,]

barplot(rbind(table(myFilteredData.spStructured2000$Year),
              table(myFilteredData.spOthers2000$Year)),beside=T,ylim=c(0,2000),col=c(grey(0.2),grey(0.8)),ylab="Number of samples")
abline(h=seq(0,2000,by=500),lty=2,col=grey(0.8))
title(main="2000m threshold")
legend("topleft",fill=c(grey(0.2),grey(0.8)),legend=c("Structured","Other"))
dev.off()

## CONSTRAIN TO SAMPLES COLLECTED "Fylkesmannen","SNO" 
tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Proeveleverandoer %in% 
                                 c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),]
tab <- table(tmp$Year,tmp$TrackRovbsID,useNA ="always" )


## plot check 
pdf(file = paste(WD,"/",modelName,"/",modelName,"DetectionsStructuredOpp",".pdf", sep="" ))
for(t in 1:nYears){
  par(mar=c(0,0,3,0),mfrow=c(1,3))
  tmp1 <- tmp[tmp$Year%in% years[t],]
  tmpNoTracks <-  tmp1[is.na(tmp1$TrackRovbsID), ]
  tmpTracks <-  tmp1[!is.na(tmp1$TrackRovbsID), ]
  
  plot(studyArea, main="Structured with track")
  plot(st_geometry(tmpTracks), pch=21, col="black", cex=1,bg="red",add=T)
  
  plot(studyArea, main="Structured without track")
  plot(st_geometry(tmpNoTracks), pch=21, col="black", cex=1,bg="blue",add=T)
  
  tmpOpp <- myFilteredData.sp$alive[!myFilteredData.sp$alive$Proeveleverandoer %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),]
  tmpOpp <- tmpOpp[tmpOpp$Year%in% years[t],]
  
  plot(studyArea, main="Other samples")
  plot(st_geometry(tmpOpp), pch=21, col="black", cex=1,bg="green",add=T)
  
  mtext(years[t],adj = -0.8,padj = 1)
  
}
barplot(tab[,which(is.na(colnames(tab)))]/rowSums(tab),main="% of samples from Statsforvalteren and \nSNO that cannot be assigned to a track") 

dev.off()

### plot check 
pdf(file = paste(WD,"/",modelName,"/",modelName,"OverallDetectionsDeadRecoveries",".pdf", sep="" ))
plot(st_geometry(GLOBALMAP))
plot(st_geometry(studyArea),add=T)
plot(st_geometry(myFullData.sp$alive),pch=16, col="red", cex=0.3,add=T)
plot(st_geometry(myFullData.sp$dead.recovery),pch=16, col="blue", cex=0.3,add=T)
mtext(paste("Live detections", length(myFullData.sp$alive),
            "; ID:", length(unique(myFullData.sp$alive$Id))
),line = +1)
mtext(paste("Dead recovery:",length(myFullData.sp$dead.recovery)))
dev.off()




## CHECK THAT ALL RESULTS FROM HAIR TRAPS WERE ASSIGNED TO OPPORTUNISTIC
whichHair <- which(myFilteredData.sp$alive$DNAID%in% HairTrapSamples$DNAID)
plot(COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten",]$geometry)
plot(myFilteredData.sp$alive[whichHair, ]$geometry, add = T, col = "red", pch = 16)



## ------    1.5.SEPARATE MORTALITY CAUSES ------ 
par(mfrow = c(1,3))
for(t in 1:nYears){
  ## DEAD RECOVERIES TOTAL
  tempTotal <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]
  NGS_TabTotal <- table(tempTotal$Country)
  ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country), 2, function(x) sum(x>0))
  ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD
  tempIn <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]
  NGS_TabIn <- table(tempIn$Country)
  ID_TabIn <- apply(table(tempIn$Id, tempIn$Country), 2, function(x) sum(x>0))
  ## PLOT NGS SAMPLES
  plot(st_geometry(GLOBALMAP), col="gray80")
  plot(st_geometry(studyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
  plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
  # points(tempTotal, pch = 21, bg = "darkred")
  plot(st_geometry(tempIn), pch = 21, bg = "blue",add=T)
  ## ADD NUMBER OF NGS samples and IDs per COUNTRY
  graphics::text(x = 100000, y = 7250000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
  graphics::text(x = 820000, y = 6820000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
  ## ADD OVERALL NUMBERS
  mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
  mtext(text = paste(sum(NGS_TabIn), "Dead Recoveries /", sum(ID_TabIn), "IDs IN"), side = 3, line = 0)
  # mtext(text = paste(sum(NGS_TabTotal), "Recoveries /", sum(ID_TabTotal)-sum(ID_TabIn), "IDs OUT"), side = 3, line = -1)
}#t



## ------    PLOT TREND DETECTIONS AND DEAD RECOVERIES OVER TIME AND SPACE ------ 
#DETECTIONS
pdf(file=file.path(WD, modelName, paste(modelName,"_TRENDDetections.pdf",sep="")))
temp <- unique(myFilteredData.sp$alive[,c("Year","Country","DNAID")])
tab_Country.Year <- table(temp$Year, temp$Country)
country.colors <- c("goldenrod1","goldenrod3")

par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)), ylab="N Detections", xlab="Years")
lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("bottomright",c("N","S"), fill=country.colors)

#ID DETECTED
temp <- table(myFilteredData.sp$alive$Year,myFilteredData.sp$alive$Country,myFilteredData.sp$alive$Id)
tab_Country.Year1 <- apply(temp,c(1,2),function(x) sum(x>0))
country.colors <- c("goldenrod1","goldenrod3")
par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year1))), ylim=c(0,max(tab_Country.Year1)), ylab="N Id detected", xlab="Years")
lines(tab_Country.Year1[,"N"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year1[,"S"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("bottomright",c("N","S"), fill=country.colors)

## Average number of detection per detected ID  #[CM]
tab_Country.Year2 <- tab_Country.Year/tab_Country.Year1
par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year2))), ylim=c(0,max(tab_Country.Year2)),
     ylab="Average Number of detections", xlab="Years")
lines(tab_Country.Year2[,"N"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year2[,"S"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("bottomright",c("N","S"), fill=country.colors)

## deadrecovery #[CM]
temp <- unique(myFilteredData.sp$dead.recovery[,c("Year","Country","Id")])
tab_Country.Year <- table(temp$Year, temp$Country)
country.colors <- c("goldenrod1","goldenrod3")

par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)),
     ylab="N Id Dead recovered", xlab="Years")
lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("topright",c("N","S"), fill=country.colors)
dev.off()



## ------    PLOT CHECK ----- 
par(mar=c(0,0,0,0))
plot(st_geometry(studyArea))
plot(st_geometry(myFilteredData.spAllSex$alive), pch=21, bg="red", cex=0.5,add=T)
plot(st_geometry(myFilteredData.spAllSex$dead.recovery), pch=21, bg="blue", cex=0.5,add=T)
plot(st_geometry(BuffDead),add=T, border="blue")
plot(st_geometry(myBufferedArea), border="red", add=T)
plot(st_geometry(studyArea), border="grey", add=T)

plot(habitat$habitat.r)
plot(st_geometry(studyArea), add = T, col = rgb(150/250,150/250,150/250, alpha = 0.75))
plot(st_geometry(GLOBALMAP), add = T)
plot(st_geometry(habitat$buffered.habitat.poly), add = T)
plot(st_geometry(myFilteredData.sp$alive), pch=21, bg="red", cex=0.5,add=T)
plot(st_geometry(myFilteredData.sp$dead.recovery),pch=21, bg="blue", cex=0.5,add=T)



## ------ PLOT CHECK TIME BETWEEN DEAD RECO AND MONITORING SEASON ------
pdf(file=file.path(WD, modelName, paste(modelName,"Prop id detected_Time available.pdf",sep="")))
plot(ndet ~ timeDiff, ylab="Total number of detections", xlab="Number of days between dec 1 and dead recovery")
hh <- hist(timeDiff[ndet>0], breaks = seq(0,400,by=25))
hh1 <- hist(timeDiff[ndet==0], breaks = seq(0,400,by=25))
barplot(rbind(hh$counts/(hh$counts+hh1$counts),
              hh1$counts/(hh$counts+hh1$counts)),names.arg=hh$breaks[1:(length(hh$breaks)-1)],
        xlab="number of days between dead reco and start monitoring",
        ylab="%"
)
legend("topright",fill=c(grey(0.2),grey(0.8)),legend = c("detected","notDetected"))
dev.off()



## ------ 7.1 calculate realized phi########
# get location of individuals 
sxy.initscaled <- scaleCoordsToHabitatGrid(coordsData = sxy.init,
                                           coordsHabitatGridCenter = habitat$habitat.xy,
                                           scaleToGrid =F )$coordsDataScaled

####
#initialize objects 
recruit <- 0
z_caculate <- nimData$z#[,1:(nYears-1)]
z_caculate[is.na(z_caculate)] <- 0
lev <- levels(habitatRasterResolution$'2.5km'$Countries)
#EXTRACT LOCATION BASED ON INITIAL AC
countryId <- list()
for(t in 1:dim(z_caculate)[2]){
  tmp <- st_as_sf(data.frame(sxy.initscaled[,,t]), coords = c("x", "y"))
  countryId[[t]] <- raster::extract( habitatRasterResolution$'2.5km'$Countries ,tmp,sparse = F)
}
phi <- phiind1 <- culled <- recruit <- recruitnb<- matrix(0,nrow=nYears-1,ncol=length(lev[[1]]$ID))
colnames(phi) <- colnames(phiind1) <- colnames(culled) <- 
  colnames(recruitnb) <-  colnames(recruit)  <- factorValues(habitatRasterResolution$'2.5km'$Countries,lev[[1]]$ID)[,1]
for(c in 1:ncol(phi)){
  for(t in 2:dim(z_caculate)[2]){
    #phi
    alivet <- which(z_caculate[,t-1] %in% c(2) & countryId[[t-1]] %in% c )
    phi[t-1, c] <- sum(z_caculate[alivet,t] %in% c(2))/length(alivet)
    #culled
    #culled[t-1] <- sum(z_caculate[alivet,t] %in% c(3))/length(alivet)
    #recru
    notentered <- which(z_caculate[,t-1] == 0)
    recruitnb[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
                              countryId[[t]][notentered] %in% c)#/sum(z_caculate[,t-1] %in% c(2))
    recruit[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
                            countryId[[t]][notentered] %in% c)/sum(z_caculate[,t-1] %in% c(2) & 
                                                                     countryId[[t-1]] %in% c)
  }
}
phi <- phi[,c(2,4)]        # overall phi
recruit <- recruit[,c(2,4)]        # overall phi
recruitnb <- recruitnb[,c(2,4)]        # overall phi


###
pdf(file = file.path( WD,
                      modelName,paste(modelName,"realizedPhiCountry",".pdf",sep="")
))
par(mfrow = c(1,1))
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z")
axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)]+1)
yr <- c(1:(nYears-1))
for(c in 1:ncol(phi)){
  points(phi[,c]~yr,pch=16,type="b", col=c)
}
legend("bottomright",colnames(phi),col=c(1:4),pch=16)


par(mfrow = c(1,1))
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized recruit from z")
axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)]+1)
yr <- c(1:(nYears-1))
for(c in 1:ncol(recruit)){
  points(recruit[,c]~yr,pch=16,type="b", col=c)
}
legend("topright",colnames(recruit),col=c(1:4),pch=16)
dev.off()



###
lev <- levels(habitatRasterResolution$'10km'$Counties)
countryId <- list()
for(t in 1:dim(z)[2]){
  tmp <- st_as_sf(data.frame(sxy.initscaled[,,t]), coords = c("x", "y"))
  countryId[[t]] <- raster::extract( habitatRasterResolution$'10km'$Counties ,tmp,sparse = F)
}

phi <- phiind1 <- culled <- recruit <- recruitnb<- matrix(0,nrow=nYears-1,ncol=length(lev[[1]]$ID))
colnames(phi) <- colnames(phiind1) <- colnames(culled) <- 
  colnames(recruitnb) <-  colnames(recruit)  <- factorValues(habitatRasterResolution$'10km'$Counties,lev[[1]]$ID)[,1]

for(c in 1:ncol(phi)){
  for(t in 2:dim(z)[2]){
    #phi
    alivet <- which(z_caculate[,t-1] %in% c(2) & countryId[[t-1]] %in% c )
    phi[t-1, c] <- sum(z_caculate[alivet,t] %in% c(2))/length(alivet)
    #culled
    #culled[t-1] <- sum(z_caculate[alivet,t] %in% c(3))/length(alivet)
    #recru
    notentered <- which(z_caculate[,t-1] == 0)
    recruitnb[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
                              countryId[[t]][notentered] %in% c)#/sum(z_caculate[,t-1] %in% c(2))
    recruit[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
                            countryId[[t]][notentered] %in% c)/sum(z_caculate[,t-1] %in% c(2) & 
                                                                     countryId[[t-1]] %in% c)
  }
}
phi <- phi[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]        # overall phi

phiNOR <- phi[,c(8,9,10,11,12,13)]
phiSWE <- phi[,-c(8,9,10,11,12,13,14,15,16)]

recruitnb <- recruitnb[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]        # overall phi
recruitnbNOR <- recruitnb[,c(8,9,10,11,12,13)]
recruitnbSWE <- recruitnb[,-c(8,9,10,11,12,13,14,15,16)]

recruit <- recruit[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]        # overall phi
recruitNOR <- recruit[,c(8,9,10,11,12,13)]
recruitSWE <- recruit[,-c(8,9,10,11,12,13,14,15,16)]


###
pdf(file = file.path( WD,
                      modelName,paste(modelName,"realizedPhiCounties",".pdf",sep="")
),width = 11,height=6)
# PHI
#NORWAY
par(mfrow = c(1,2))
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z", main="Norway")
axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
yr <- c(1:(nYears-1))
for(c in 1:ncol(phiNOR)){
  points(phiNOR[,c]~yr,pch=16,type="b", col=c)
}
legend("bottomleft",colnames(phiNOR),col=c(1:ncol(phiNOR)),pch=16)

###SWEDEN
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z", main="Sweden")
axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
yr <- c(1:(nYears-1))
for(c in 1:ncol(phiSWE)){
  points(phiSWE[,c]~yr,pch=16,type="b", col=c)
}
legend("bottomleft",colnames(phiSWE),col=c(1:ncol(phiSWE)),pch=16)


# RECRUITS
#NORWAY
par(mfrow = c(1,2))
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", 
     ylab = "Realized recruitment from z", main="Norway")
axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
yr <- c(1:(nYears-1))
for(c in 1:ncol(recruitNOR)){
  points(recruitNOR[,c]~yr,pch=16,type="b", col=c)
}
legend("topleft",colnames(phiNOR),col=c(1:ncol(phiNOR)),pch=16)

###SWEDEN
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized recruitment from z", main="Sweden")
axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
yr <- c(1:(nYears-1))
for(c in 1:ncol(recruitSWE)){
  points(recruitSWE[,c]~yr,pch=16,type="b", col=c)
}
legend("topleft",colnames(recruitSWE),col=c(1:ncol(phiSWE)),pch=16)
dev.off()

##prop detected vs Alive in z
propDet <- 0
for(t in 1:nYears){
  whichdets <- unique(c(which(nimData$nbDetections[,t]>0),
                        which(nimData$nbDetectionsOth[,t]>0)))
  whichAlive <- which(nimData$z[,t]%in%2)
  propDet[t] <- length(whichdets)/length(whichAlive)
}



for(c in 1:4){
  ## ------    3.4. LIST NIMBLE INITS ------
  nimInits <- list( "sxy" = sxy.init,
                    "dmean" = runif(1,0,10),
                    "z" = z.init,
                    "omeg1" = c(0.5,0.5),
                    "gamma" = runif(dim(y.alive)[3]-1,0,1),
                    "p01" = array(runif(18,0,0.2), c(nimConstants$n.counties,dim(y.alive)[3])),
                    "p01Oth" = array(runif(18,0,0.2), c(nimConstants$n.countries+1,dim(y.alive)[3])),
                    "sigma" = runif(nYears,1,4),
                    "betaDens" = runif(1,-0.1,0.1),#[CM]#0,
                    "betaCovs" = array( runif(dim(detCovs)[3],-0.1,0.1),c(dim(detCovsOth)[3],nYears)),#[CM]rep(0,dim(detCovs)[3]),
                    "betaCovsOth" = array( runif(dim(detCovsOth)[3],-0.1,0.1),c(dim(detCovsOth)[3],nYears)),#[CM]rep(0,dim(detCovs)[3]),
                    "betaResponseOth" = runif(dim(y.alive)[3], -0.1, 0.1),#[CM]#0,
                    "betaResponse" = runif(dim(y.alive)[3], -0.1, 0.1),#[CM]#0,
                    "detResponse" = InitsDetResponse,
                    "pResponse"  = runif(1, 0.4, 0.5),#[CM]#0,
                    # "h" = runif(dim(y.alive)[3]-1,0.1,0.3), #  runif(dim(y.alive)[3]-1,0.2,0.4),
                    "phi" = runif(dim(y.alive)[3]-1,0.1,0.3)) #,runif(dim(y.alive)[3]-1,0.2,0.4))
  
  # SXY
  # nimInits$sxy <- UTMToGrid(data.sxy = sxy.init,
  #                           grid.sp = habitat$habitat.sp)$data.scaled.xy
  # dimnames(sxy.init)[[2]] <- c("x","y")
  # nimInits$sxy <- scaleCoordsToHabitatGrid(coordsData = sxy.init,
  #                                             coordsHabitatGridCenter = habitat$habitat.xy,
  #                                             scaleToGrid =F )$coordsDataScaled
  # 
  # nimData$sxy <- UTMToGrid(data.sxy = nimData$sxy,
  #                          grid.sp = habitat$habitat.sp)$data.scaled.xy
  
  ### TEST IF THE LESS RESTRICTION ON DETECTORS WILL WORK 
  ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
  i=586
  t=10
  idDEtected <- which(!rownames(z) %in%"Augmented")
  
  for(i in 1:length(idDEtected)){
    for(t in 1:nimConstants$n.years){
      if(!is.na(nimInits$sxy[i,1,t])){
        SXY <- nimInits$sxy[i,,t]  
      }else{SXY <- nimData$sxy[i,,t]}
      sxyID <- nimData$habitatID[trunc(SXY[2]/nimConstants$ResizeFactor)+1, trunc(SXY[1]/nimConstants$ResizeFactor)+1]
      DETECTIndexdetectorIndex <- nimData$detectorIndex[1:nimConstants$n.cellsSparse, 
                                                        1:nimConstants$maxNBDets] 
      DETECTLESS <- nimData$nDetectorsLESS[1:nimConstants$n.cellsSparse]
      index <- DETECTIndexdetectorIndex [sxyID,1:DETECTLESS[sxyID]]
      
      #table(detectorIndex)
      ## GET NECESSARY INFO 
      n.detectors <- length(index)
      #maxDist_squared <- maxDist*maxDist
      
      YDET <- nimData$yDets[i,1:nimConstants$nMaxDetectors, t]
      YDETOth <- nimData$yDetsOth[i,1:nimConstants$nMaxDetectorsOth, t]
      
      ## RECREATE Y
      if(nimData$nbDetections[i, t] > 0){
        for(j in 1:nimData$nbDetections[i, t]){
          ## check if a detection is out of the "detection window"
          if(sum(YDET[j]==index)==0){
            print(paste("id",i,"t",t,"j",j))
          }
        }
      }
    }}
  
  plot(nimData$detector.xy[,2]~nimData$detector.xy[,1])
  points(nimData$detector.xy[index,2]~nimData$detector.xy[index,1], col="red")
  points(SXY[2]~SXY[1], col="blue", pch=16)
  # points( s[i,2,t]~ s[i,1,t], col="green", pch=16)
  # 
  # 
  # 
  
  points(nimData$detector.xy[YDET[1:nimData$nbDetections[i, t]],2]~
           nimData$detector.xy[YDET[1:nimData$nbDetections[i, t]],1], col="green", pch=16)
  points(nimData$detector.xy[YDETOth[1:nimData$nbDetectionsOth[i, t]],2]~
           nimData$detector.xy[YDETOth[1:nimData$nbDetectionsOth[i, t]],1], col="purple", pch=16)
  
  # 
  # i=1307
  # i=1486
  # 
  # myData.dead[myData.dead$Id %in% row.names(y.ar.ALIVE)[i],]
  plot(st_geometry(COUNTRIES))
  tmp <- myData.alive$myData.sp[myData.alive$myData.sp$Id %in% row.names(y.ar.ALIVE)[i] & myData.alive$myData.sp$Year %in% years[t],]
  # 
  # tmp <- myData.aliveOthers$myData.sp[myData.aliveOthers$myData.sp$Id %in% row.names(y.ar.ALIVE)[i] &
  #                                          myData.aliveOthers$myData.sp$Year %in% years[t],]
  # tmp <- myData.aliveStruc$myData.sp[myData.aliveStruc$myData.sp$Id %in% row.names(y.ar.ALIVE)[i] &
  #                                      myData.aliveStruc$myData.sp$Year %in% years[t],]
  # 
  plot(st_geometry(tmp),col="red",add=T)
  # 
  ## PLOT CHECK 
  
  nimInits$sxy <- round(nimInits$sxy, 5)#---an extreme number of decimals may cause a number to appear as an integer to Nimble, and then coincide with habitat window boundaries
  
  
  ##CHECK WHERE IS NORRBOTTEN. IT IS ON THE 5TH INDEX
  plot(st_geometry(GLOBALMAP), col = "gray80", main = "Aggregated Counties")
  plot(st_geometry(studyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  plot(st_geometry(detectors$main.detector.sp[detCounties%in% c(1),]), col = myCol[5], pch = 16, cex = 0.8,add=T)
  
  nimConstants$countyToggle <- nimInits$p01
  nimConstants$countyToggle[] <- 1
  
  yearsNotSampled <- which(!years%in% yearsSampledNorrb)
  for(t in yearsNotSampled){
    nimConstants$countyToggle[1,t] <- 0
  }
  
  ## add another category to detcountry if in norrbotten, to turnoff detection to 0 there. 
  detCountriesNorb <- matrix(NA, nrow=length(detCountries),ncol=nYears)
  detCountries1 <- detCountries
  detCountries1[detCounties %in% 1] <- 3
  for(t in 1:nYears){
    if(t %in% yearsNotSampled){
      detCountriesNorb[,t] <- detCountries1
    }else{
      detCountriesNorb[,t] <- detCountries
    }
  }  
  
  nimData$detCountries <-  detCountriesNorb
  
  nimConstants$countyToggleOth <- nimInits$p01Oth
  nimConstants$countyToggleOth[] <- 1
  yearsNotSampled <- which(!years%in% yearsSampledNorrb)
  for(t in yearsNotSampled){
    nimConstants$countyToggleOth[3,t] <- 0
  }
  
  ## ------ 7. SAVE NIMBLE INPUT ------
  save(nimData,
       nimConstants,
       y.dead,
       nimParams,
       nimParams2,
       modelCode,
       nimInits,
       file = file.path(WD, modelName,
                        paste(modelName,"Chain", c, ".RData", sep = "")))
  #####
}#c



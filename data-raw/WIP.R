
## ------       2.2.4. EXTRACT DISTANCE TO ROADS ------ 

DistAllRoads <- raster::stack(file.path(data.dir, "Roads/MinDistAllRoads1km.tif"))
##-- Fasterize to remove values that fall in the sea
r <- fasterize::fasterize(sf::st_as_sf(REGIONS), DistAllRoads$MinDistAllRoads1km)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r
DistAllRoads <- raster::crop(DistAllRoads, detectors$grid)
rm(list = c("r"))
##-- Aggregate to match the detectors resolution
DistAllRoads <- raster::aggregate( 
  x = DistAllRoads,
  fact = detectors$resolution/raster::res(DistAllRoads),
  fun = mean)
##-- Extract distance to roads for each detector
detRoads <- raster::extract(DistAllRoads, detectors$main.detector.sp)



##-- ALTERNATIVE WITH TERRA
detRoads2 <- terra::rast(file.path(data.dir, "Roads/MinDistAllRoads1km.tif")) %>%
  terra::crop(., terra::ext(detectors$grid)) %>%
  terra::aggregate(.,
    fact = detectors$resolution/terra::res(.),
    fun = mean) %>%
  terra::extract(., st_coordinates(detectors$main.detector.sp))



##------------------------------------------------------------------------------
## ------       2.2.5. EXTRACT DAYS OF SNOW ------ 

##-- Load raster stack of snow cover
SNOW <- stack(file.path(data.dir, "Snow/AverageSnowCoverModisSeason2014_2025_Wolf.tif"))
##-- Select snow data corresponding to the monitoring period
SNOW <- SNOW[[paste("X", years, "_", years + 1, sep = "")]]
SNOW <- raster::crop(SNOW, c(0,40,55,75))
##-- Extract snow 
detSnow <- matrix(0, nrow = dim(detectors$main.detector.sp)[1], ncol = n.years)
det.sptransf <- st_transform(detectors$main.detector.sp, st_crs(SNOW))
detSnow[ ,1:n.years] <- raster::extract(SNOW, det.sptransf)
##-- if NA returns the average value of the cells within 15km 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract( SNOW, det.sptransf[isna, ],
                        buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:n.years] <- tmp
##-- still some NA... Increase buffer again 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract( SNOW, det.sptransf[isna, ],
                        buffer = 35000, fun = mean, na.rm = T)
detSnow[isna,1:n.years] <- tmp
##-- Put into "nimble2SCR" format
colnames(detSnow) <- paste0("snow.", years)
detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detSnow)



##-- ALTERNATIVE WITH TERRA
detSnow <- terra::rast(file.path(data.dir, "Snow/AverageSnowCoverModisSeason2014_2025_Wolf.tif")) %>%
  terra::crop(., c(0,40,55,75)) %>%
  terra::extract(., st_transform(detectors$main.detector.sp, st_crs(.))) %>%
  rename(., )
  



library(dplyr)
library(colorspace)
library(rovquantR)
library(ggsflabel)

##-- Get polygons of regions in Norway
REGIONS_NOR <- REGIONS %>%
  filter(country == "NOR") %>%
  group_by(region) %>%
  summarise() %>%
  mutate(num_region = 1:nrow(.))


##-- Get polygons of counties in Norway
COUNTY_NOR <- REGIONS %>%
  filter(country == "NOR") %>%
  group_by(county) %>%
  summarise()


##-- Specify color palette
this.col <- sequential_hcl(1+length(unique(REGIONS_NOR$region)), "Reds 3")
this.col <- this.col[c(2,7,6,3,5,8,1,4)]
# this.col <- c("#990E1F", "#FFBEC1", "#FF999E", "#CC1C2F", "#FF7078", "#FFDEE0", "#69000C", "#F53C4B")


##-- Map of large carnivore management regions and counties in Norway
pdf(file = file.path(working.dir, "RegionMaps.pdf"),
    width = 9, height = 13, pointsize = 12)

par(mar=c(0,0,0,0))

ggplot( data = REGIONS_NOR) +
  ## Plot management regions in red
  geom_sf( fill = this.col, col = NA) + 
  ## Add county borders in grey
  geom_sf( data = COUNTY_NOR, fill = NA, col = "gray20",lwd = 0.3) + 
  ## Add county labels
  geom_sf_label_repel( data = COUNTY_NOR,
                       aes(label = county),
                       force = 10,
                       nudge_y = c(-150000,-10000,50000,120000,60000,
                                   80000,-50000,-50000,-50000,-200000,
                                   90000,0,-100000,0,-50000),
                       nudge_x = c( -20000,200000,350000,-200000,250000,
                                    -250000,200000,200000,-180000,65000,
                                    -200000,200000,100000,-200000,100000)) +
  ## Add region numbers
  geom_sf_text( data = REGIONS_NOR,
                aes(label = num_region),
                cex = 8, fontface = "bold",
                nudge_y = c(80000,rep(0,7)),
                nudge_x = rep(0,8),
                col = c("white", "black", "black", "black",
                        "black", "black", "white", "black" )) +
  theme_void() +
  ## Remove background and axes
  theme( legend.position = "none", 
         panel.grid=element_blank(), 
         axis.text = element_blank(), 
         axis.ticks = element_blank())

dev.off()

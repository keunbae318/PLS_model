##Map

install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel",
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))

library("ggplot2")
theme_set(theme_bw())
library("sf")
library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf()

ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$NAME)), " countries)"))

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-10, 30), ylim = c(45, 60), expand = FALSE)

ggplot(data = world) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-10, 30), ylim = c(45, 60), expand = FALSE)

library("sf")
world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))

str(Variables.dat)
Undrained.dat= subset(Variables.dat,Treatment=="U"&Depth=="T")
ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  geom_text(data= Undrained.dat,aes(x= Longitude, y=Latitude, label=Site),
            color = "darkblue", fontface = "bold", check_overlap = FALSE,size=2) +
  # annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", 
  #          fontface = "italic", color = "grey22", size = 6) +
  coord_sf(xlim = c(-7, 26), ylim = c(46, 58), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))



library(tidyverse)
library(OpenStreetMap)

str(Variables.dat)
Zwarte_Beek.dat= subset(Variables.dat,Site=="Zwarte Beek"&Depth=="T")
Zwarte_Beek.dat$Treatment = factor(Zwarte_Beek.dat$Treatment, levels = c("U","R","D"),
                                   labels= c("Undrained","Rewetted","Drained"))
Zwarte_Beek.dat$Latitude <-Zwarte_Beek.dat$Latitude + 0.008 
lat1 <- 51.07; lat2 <- 51.095; lon1 <- 5.27; lon2 <- 5.33
#install.packages("OpenStreetMap")
sa_map <- openmap(c(lat2, lon1), c(lat1, lon2), 
                  type = "osm", mergeTiles = TRUE)
sa_map2 <- openproj(sa_map)

sa_map2_plt <- OpenStreetMap::autoplot.OpenStreetMap(sa_map2) + 
  # annotate("text", label = "Atlantic\nOcean", 
  #          x = 18.2, y = -33.8, size = 5.0, angle = -60, colour = "navy") +
  geom_point(data = Zwarte_Beek.dat,
             aes(x = Longitude  + 0.002, y = Latitude - 0.007, color=Treatment ,shape =Treatment), # slightly shift the points
              size =  6) +
  scale_color_manual(values=c("#00AFBB","#E7B800","#FC4E07"))+
  geom_text(data = Zwarte_Beek.dat, # Choose dataframe
            aes(Longitude+0.006, Latitude-0.009, label = Treatment), # Set aesthetics
            hjust = 1.15, vjust = 0.5, # Adjust vertical and horizontal
            size = 3, colour ="black")+
  xlab("Longitude (°E)") + ylab("Latitude (°S)")
sa_map2_plt
annotation_scale()

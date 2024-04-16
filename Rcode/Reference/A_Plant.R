library("sf")

library("ggplot2")
theme_set(theme_bw())
library("sf")
library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")

Plant.dat=read.csv("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/Data/Plant_data.csv")
Plant.dat= subset(Plant.dat,DrainStatus=="natural")
ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  geom_text(data= Plant.dat,aes(x= East, y=North, label=ID),
            color = "darkblue", fontface = "bold", check_overlap = FALSE,size=2) +
  # annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", 
  #          fontface = "italic", color = "grey22", size = 6) +
  coord_sf(xlim = c(-7, 26), ylim = c(46, 58), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))



ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  geom_text(data= Plant.dat,aes(x= East, y=North, label=ID),
            color = "darkblue",  check_overlap = TRUE ,size=3) +
  # annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", 
  #          fontface = "italic", color = "grey22", size = 6) +
  coord_sf(xlim = c(5, 24), ylim = c(53.5, 54.5), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))


str(Subset_map)
Subset_map <- subset(Plant.dat, ID >= 300 & ID <= 318)
Subset_map[,1:5]
View(Subset_map)

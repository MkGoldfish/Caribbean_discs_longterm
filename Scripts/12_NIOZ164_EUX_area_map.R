##%#####################################################################################%##
#NIOZ164 Statia Discs - Area maps                                             #####                                            
#Author: Maaike Goudriaan, NIOZ, MMB                                                     #
#                                                                                         #
# Purpose:                                                                                # 
# NIOZ164 is an illumina sequencing lane consisting of different samplesets.              # 
# Here we focus on analysis of 16S amplicon data, from amplified DNA                      #
# extracted from discs covered in polymer films that were incubated at 2 different        #
# locations in 2 different waterdepths (seafloor and watercolumn) in                      #
# coastal Caribbean waters close to the island of St. EUstatius.                          #
# The tested polymers were PE, PP, PS, PEt and Nylon, PE-13C and PP-13C.                  #
# Of each polymere there was an UV pretreated and a non-treated version.                  #
##%#####################################################################################%##

# Date: 2023 - 08 - 23
# R-version: 4.3.1 

## Set working directory ------------------------------------------------------------------
setwd("C:/Users/mgoudriaan/Documents/Github/Carribean_foils_1week/Scripts")

# Load libraries -------------------------------------------------------------------------
library("ggspatial") 
library("sf")
library("ggpp")
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")
library(cowplot)
library(broom)
library(svglite)
library(tidyverse)

# Create area map of the Caribbean ------------------------------------------------------
## Import geographical data from Natural Earth database
carrib <- ne_countries(scale = 10, returnclass = "sf")
world_points<- st_centroid(carrib)
world_points <- cbind(carrib, st_coordinates(st_centroid(carrib$geometry)))

## Provide coordinates of EUX to indicate position in map
st.eux <- data.frame(longitude = -62.9666628, latitude = 17.499998)

## Colors to use
pal.loc <- c("#FF6DB6FF" , "#004949FF",  "#66A61E")

## Plot the Caribbean overview w EUX indicated. 
# Also put labels for the different Carribean islands
carrib.p <- ggplot(data=carrib) + 
  geom_sf(fill= "grey20")+
  coord_sf(xlim = c(-58, -86), ylim = c(9,22)) +
  annotation_scale(bar_cols = c("black", "white"), location = "bl", width_hint = 0.4, pad_x = unit(1.4, "cm"), pad_y = unit(0.45, "cm")) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0, "cm"), pad_y = unit(0, "cm"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(data = st.eux, aes(x = longitude, y = latitude),  
             size = 6, shape = "diamond filled", fill = "#332288", stroke = 2)+
  annotate(geom = "text", x = -71, y = 15.5, label = "Carribean Sea", 
           fontface = "italic", color = "grey20", size = 10) +
  annotate(geom = "text", x = -69.4, y = 20.5, label = "Greater Antilles", 
           fontface = "bold.italic", color = "black", size = 7, angle = 350) +
  annotate(geom = "text", x = -61.3, y = 18.0, label = "Leeward Islands", 
           fontface = "bold.italic", color = "black", size = 5.5, angle = 310) +
  annotate(geom = "text", x = -59.2, y = 13.2, label = "Windward Islands", 
           fontface = "bold.italic", color = "black", size = 5.5, angle = 70) +
  annotate(geom = "text", x = -67.6, y = 12.6, label = "Leeward Antilles", 
           fontface = "bold.italic", color = "black", size = 5.5, angle = 350) +
  theme(panel.grid.major = element_line(color = gray(.4), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "white"),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12, colour = "black"),
        plot.title = element_text(size = 16)) +
   xlab("Longitude") + ylab("Latitude")

carrib.p

# Create map of Sint Eustatius with Samplesites ------------------------------------------------------
# Import datafiles, doawnloaded from Dutch Caribbean Biodiversity Database
# https://www.dcbd.nl/document/topography-steustatius
steux <- sf::st_read("topographyStatia/Contours_20n.shp")
roads <- readOGR("../topographyStatia/main_road_20n.shp") %>% broom::tidy(region = "NAME")
sites <- readOGR("../topographyStatia/EUX_DiveSites_Moorings_from_coords.shp")  %>% data.frame() %>% filter(Depth_m > 6)
map <-  system.file("topographyStatia/referenceStatia_AerialPhotos.shp", package="sf")
map <- sf::st_read("./topographyStatia/referenceStatia_AerialPhotos.shp")

unique(sites$Name)
summary(map@data)

# Filter/create a df with divesites we used for incubations
sites.f <- sites %>% filter(Name %in% c("Crook's Castle", "The Charles L Brown"))
sites.d  <- data.frame(Name = c("Crook's Castle", "The Charles L Brown", "Zeelandia"), longitude = c(-62.98757, -62.99413, -62.9812 ), latitude = c(17.47192, 17.46400, 17.5068 ))

# Plot contour map of the island with dive/collection sites indicated. 
EUX <- ggplot() + 
  geom_polygon(data=map,aes(x = long, y = lat, group = group),fill= "grey30", color = "black") +
  geom_point(data = sites.d, aes(x = longitude, y = latitude),  
             size = 6, shape = "diamond filled", fill = pal.loc, stroke = 2) +
  annotate(geom = "text", x = -62.9975, y = 17.474, label = "Crook's\nCastle",
           fontface = "bold", color = "black", size = 6) +
  annotate(geom = "text", x = -62.97, y = 17.464, label = "Charles Brown",
           fontface = "bold", color = "black", size = 6) +
  annotate(geom = "text", x = -62.965, y = 17.508 , label = "Zeelandia",
           fontface = "bold", color = "black", size = 6) +
  theme_classic()+ 
  theme(panel.grid.major = element_line(color = gray(.4), linetype = "dashed", linewidth = 0.5),
        panel.background = element_rect(fill = "white", color = "grey10", linewidth = 2),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(face = "bold", size = 18, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 18, colour = "black"),
        plot.title = element_text(size = 16),
        plot.background = element_rect(fill = NA, color = NA)) +
  xlab("") + ylab("")

EUX

eux.tb <- tibble(x = -64, y =16, plot = list(EUX))

# Make embedded map of the island within the region
ggdraw(carrib.p) +
  draw_plot(EUX, 
            x = 0.1,
            y = 0.3,
            width = 0.3,
            height = 0.45)

# Other way is to plot grid with both maps
plot_grid(carrib.p, EUX,   
          ncol=2,  
          rel_widths =c(1,0.3),
          align = 'h',
          axis = "tb")

ggsave("Genus_bubble.tiff", 
       width = 25, height  = 19, unit = "cm", 
       dpi = 500, bg = 'white')

# ggplot(data=carrib) + 
#   geom_sf(fill= "grey30")+
#   coord_sf(xlim = c(-61, -71), ylim = c(9,19)) +
#   annotation_scale(bar_cols = c("black", "white"), location = "bl", width_hint = 0.4, pad_x = unit(1.4, "cm"), pad_y = unit(0.45, "cm")) +
#   annotation_north_arrow(location = "bl", which_north = "true", 
#                          pad_x = unit(0, "cm"), pad_y = unit(0, "cm"),
#                          style = north_arrow_fancy_orienteering) +
#   geom_point(data = st.eux, aes(x = longitude, y = latitude),  
#              size = 5, shape = "diamond filled", fill = "firebrick4", stroke = 2)+
#   annotate(geom = "text", x = -69, y = 17, label = "Carribean Sea", 
#            fontface = "italic", color = "grey10", size = 8) +
#   geom_plot(data = eux.tb, aes(x,y, label = plot)) +
#   theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
#         panel.background = element_rect(fill = "grey95"),
#         axis.text.y = element_text(colour = "black", size = 25),
#         axis.text.x = element_text(colour = "black", size = 25),
#         axis.title.x = element_text(face = "bold", size = 30, colour = "black"),
#         axis.title.y = element_text(face = "bold", size = 30, colour = "black"),
#         plot.title = element_text(size = 16)) +
#   xlab("Longitude") + ylab("Latitude") 
#          
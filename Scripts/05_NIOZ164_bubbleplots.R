##%#####################################################################################%##
#NIOZ164 Statia Discs - Alpha Diversity                                              #####                                            
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
setwd("C:/Users/mgoudriaan/Documents/GitHub/Caribbean_discs_longterm/Scripts")
set.seed(42)

## Load libraries -------------------------------------------------------------------------
library(devtools)
library(tidyverse)
library(ggpubr)
library(ggh4x)
library(cowplot)
library (emojifont)
library(ggthemes)
library(stringr)
library(cowplot)

## Import data ----------------------------------------------------------------------------
tt <- read.csv("../Processed-Data/NIOZ164_EUX_discs_tidy_RA_data_filtered.csv", row.names = NULL)
head(tt)
tt$X <- NULL
tt.1<- tt %>%  mutate(Polymer_Isotope = if_else(Polymer %in% c("PE", "PP"), paste(Polymer, Isotope, sep = "-"), Polymer))
head(tt.1)

tt.inc <- tt.1 %>% filter(Phase == "Disc" ) %>% filter(!Treatment == "NA ")
unique(tt.inc$Location)
unique(tt.inc$Treatment)
unique(tt.inc$Polymer)

tt.wild <- tt.1 %>% filter(Method == "Collection") 
unique(tt.wild$Location)
unique(tt.wild$Method)
unique(tt.wild$Polymer)

## Colors for plotting --------------------------------------------------------------------
pal_isme <- c("#006d77", "#ffddd2", "#00C49A", "#e29578", "#83c5be")
Pal.plast <- c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )
pal.time <- c("#44AA99", "#882255")
pal.uv <- c("#999933", "#CC6677")

colors_M1 <- c("#004e64", "#ecc8af", "#F2AF29", "#436436", "#00a5cf", 
               "#c18c5d", "#5f0f40", "#DC602E", "#495867", "#A29F15", 
               "#570000", "#FFF5B2", "#20221B", "#9fffcb", "#c08497", 
               "#8D6346", "#FF4B3E", "#149911", "#472d30")

colors_by_Maaike <- c("#004e64", "#ecc8af", "#F2AF29", "#436436", "#00a5cf", 
                      "#c18c5d", "#5f0f40", "#DC602E", "#495867", "#A29F15", 
                      "#570000", "#FFF5B2", "#20221B", "#9fffcb", "#c08497",
                      "#8D6346", "#FF4B3E", "#149911", "#472d30", "#ce796b",
                      "#25a18e", "#BC412B", "#95D9DA", "#B10F2E", "#0E273C",
                      "#E3FDD8", "#353535", "#e7ad99", "#0F8B8D", "#7ae582",
                      "#F2AF29", "#606c38", "#3d405b", "#94d2bd", "#772e25",
                      "#344e41", "#0047E0", "#6c584c", "#5f0f40", "#D7F171", 
                      "#c89f9c", "#339989", "#faf3dd", "#04724d", "#98B9AB",
                      "#b09e99", "#AD343E", "#F2AF29", "#362C28", "#5171A5",
                      "#F7FE72", "#F4978E", "#7A9B76", "#8A7E72", "#143642", 
                      "#662C91")

colors_M2 <- c("#20221B", "#9fffcb", "#c08497", 
               "#8D6346", "#FF4B3E", "#149911", "#472d30", "#ce796b", 
               "#25a18e", "#BC412B", "#95D9DA", "#B10F2E", "#0E273C",
               "#E3FDD8", "#353535", "#e7ad99", "#0F8B8D", "#7ae582",
               "#F2AF29", "#606c38", "#3d405b", "#94d2bd", "#772e25",
               "#344e41", "#0047E0", "#6c584c", "#5f0f40", "#D7F171", "#c89f9c" )


colors_M3 <- c(   "#339989", "#faf3dd", "#04724d", "#98B9AB",
                  "#b09e99", "#AD343E", "#F2AF29", "#362C28", "#5171A5",
                  "#F7FE72", "#F4978E", "#7A9B76", "#8A7E72", "#143642", 
                  "#662C91")

pal.cb.tol.1 <- c('#77AADD', '#EE8866', '#FFAABB', '#AAAA00', '#99DDFF', 
                  '#44BB99', '#BBCC33', '#222255', '#225555', '#225522', 
                  '#666633', '#663333', '#555555', '#EE7733', '#0077BB',  
                  '#EE3377', '#33BBEE', '#CC3311', '#009988' )

pal.cb.tol.2 <- c('#125A56', '#238F9D', '#60BCE9', '#9DCCEF', 
                  '#FFB954', '#FD9A44', '#E94C1F', '#A01813',
                  '#762A83', '#9970AB', '#5AAE61', '#1B7837',
                  '#EC7014', '#993404', '#662506' )

# Generating bubbleplots --------------------------------------------------------------------

## Genus Relative Abundance --------------------------------------------------------------------
Genus <- tt.inc  %>%  select(Description, Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, Treatment, 
                             Phylum, Order, Genus, Genus_rel_abund_Sample)%>% 
  distinct() 

Genus.pel <- Genus %>% filter(Habitat == "Pelagic")
Genus.bent <- Genus %>% filter(Habitat == "Benthic")

#Cout how many unique genera we have. 
Genus.pel %>% select(Genus) %>%  unique() %>% count()
Genus.bent %>% select(Genus) %>%  unique() %>% count()

### Incubation -----------------------------------------------------------------
#### Pelagic -------------------------------------------------------------------
# select top_n genera per sample Genera for plotting
Genus_top_pel <- Genus.pel %>% dplyr::select(Description, Genus, Genus_rel_abund_Sample) %>% 
  filter ( !Genus %in% c("NA", "unassigned")) %>%  #First remove unassigned genera
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Genus_rel_abund_Sample, n = 3) %>% ungroup()

# #Check how much/which genera we find
# unique(Genus_top$Genus)
# # and with an RA above 1%
# Genus_top %>% filter(Genus_rel_abund_Sample > 0.01) %>% select(Genus) %>%  unique()

#Filter genera with RA>1%, to avoid bubbles with 0 value
top_genus_pel= Genus.pel %>% filter(Genus%in%unique((c(Genus_top_pel$Genus)))) %>%  
  filter(Genus_rel_abund_Sample > 0.01) 
top_genus %>% select(Genus) %>% unique() 

top_genus_pel$Genus <- factor(top_genus_pel$Genus, levels=rev(sort(unique(top_genus$Genus))))
top_genus_pel <- top_genus_pel %>% arrange(Order)

Genus_bubble <- ggplot(top_genus_pel,aes(x=interaction(Polymer_Isotope,Backbone),y= Genus)) +
  geom_point(aes(size=Genus_rel_abund_Sample, fill = factor(Polymer)), shape = "circle filled", stroke = 1, colour = "black") +
  scale_fill_manual(values = Pal.plast) +
  scale_size(range = c(2,9.5))+
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 10))) +
  ylab("") +
  xlab("") +  
  facet_nested(Phylum  ~ Location + Treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins',
               nest_line = element_line(),
               strip = strip_nested(
                 background_y =  elem_list_rect(color = "grey25",
                                                fill = "white", linewidth = 1),
                 text_y = elem_list_text(size = 12, angle = 0, 
                                         color = "grey25", by_layer_y = F),
                 background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=13),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Genus_bubble

#### Benthic -------------------------------------------------------------------
# select top_n genera per sample Genera for plotting
Genus_top_bent <- Genus.bent %>% dplyr::select(Description, Genus, Genus_rel_abund_Sample) %>% 
  filter ( !Genus %in% c("NA", "unassigned")) %>%  #First remove unassigned genera
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Genus_rel_abund_Sample, n = 3) %>% ungroup()

# #Check how much/which genera we find
# unique(Genus_top$Genus)
# # and with an RA above 1%
# Genus_top %>% filter(Genus_rel_abund_Sample > 0.01) %>% select(Genus) %>%  unique()

#Filter genera with RA>1%, to avoid bubbles with 0 value
top_genus_bent= Genus.bent %>% filter(Genus%in%unique((c(Genus_top_bent$Genus)))) %>%  
  filter(Genus_rel_abund_Sample > 0.01) 
top_genus %>% select(Genus) %>% unique() 

top_genus_bent$Genus <- factor(top_genus_bent$Genus, levels=rev(sort(unique(top_genus$Genus))))
top_genus_bent <- top_genus_bent %>% arrange(Order)

Genus_bubble <- ggplot(top_genus_bent,aes(x=interaction(Polymer_Isotope,Backbone),y= Genus)) +
  geom_point(aes(size=Genus_rel_abund_Sample, fill = factor(Polymer)), shape = "circle filled", stroke = 1, colour = "black") +
  scale_fill_manual(values = Pal.plast) +
  scale_size(range = c(2,9.5))+
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 10))) +
  ylab("") +
  xlab("") +  
  facet_nested(Phylum  ~ Location + Treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins',
               nest_line = element_line(),
               strip = strip_nested(
                 background_y =  elem_list_rect(color = "grey25",
                                                fill = "white", linewidth = 1),
                 text_y = elem_list_text(size = 12, angle = 0, 
                                         color = "grey25", by_layer_y = F),
                 background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=13),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Genus_bubble

#### Wild plastic-----------------------------------------------------------------
Genus.wild <- tt.wild  %>%  select(Location, Habitat, Description, Phylum, Order, Genus, Genus_rel_abund_Sample)%>% 
  distinct() 

# select top_n genera per sample Genera for plotting
Genus_top_wild <- Genus.wild %>% dplyr::select(Description, Genus, Genus_rel_abund_Sample) %>% 
  filter ( !Genus %in% c("NA", "unassigned", " ")) %>%  #First remove unassigned genera
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Genus_rel_abund_Sample, n = 3) %>% ungroup()

#Check how much/which genera we find
unique(Genus_top_wild$Genus)
# and with an RA above 1%
Genus_top_wild %>% filter(Genus_rel_abund_Sample > 0.01) %>% select(Genus) %>%  unique()

#Filter genera with RA>1%, to avoid bubbles with 0 value
top_genus_wild= Genus.wild %>% filter(Genus%in%unique((c(Genus_top_wild$Genus)))) %>%  
  filter(Genus_rel_abund_Sample > 0.01) 
top_genus_wild %>% select(Genus) %>% unique() 

top_genus_wild$Genus <- factor(top_genus_wild$Genus, levels=rev(sort(unique(top_genus_wild$Genus))))
top_genus_wild <- top_genus_wild %>% arrange(Order)

Genus_bubble <- ggplot(top_genus_wild,aes(x=Description, y= Genus)) +
  geom_point(aes(size=Genus_rel_abund_Sample, fill = factor(Phylum)), shape = "circle filled", stroke = 1, colour = "black") +
  scale_fill_manual(values = rev(colors_M3)) +
  scale_size(range = c(2,9.5))+
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 10))) +
  ylab("") +
  xlab("") +  
  facet_nested(Phylum + Order ~ Location, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins',
               nest_line = element_line(),
               strip = strip_nested(
                 background_y =  elem_list_rect(color = "grey25",
                                                fill = "white", linewidth = 1),
                 text_y = elem_list_text(size = 12, angle = 0, 
                                         color = "grey25", by_layer_y = F),
                 background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=13),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Genus_bubble

## Family Relative Abundance --------------------------------------------------------------------
Family <- tt.inc  %>%  select(Description, Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, Treatment, 
                             Phylum, Class, Family, Family_rel_abund_Sample)%>% 
  distinct() 

Family.pel <- Family %>% filter(Habitat == "Pelagic")
Family.bent <- Family %>% filter(Habitat == "Benthic")

#Cout how many unique genera we have. 
Family.pel %>% select(Family) %>%  unique() %>% count()
Family.bent %>% select(Family) %>%  unique() %>% count()

### Incubation -----------------------------------------------------------------
#### Pelagic -------------------------------------------------------------------
# select top_n genera per sample Genera for plotting
Family_top_pel <- Family.pel %>% dplyr::select(Description, Family, Family_rel_abund_Sample) %>% 
  filter ( !Family %in% c("NA", "unassigned")) %>%  #First remove unassigned genera
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Family_rel_abund_Sample, n = 5) %>% ungroup()

# #Check how much/which genera we find
# unique(Family_top$Family)
# # and with an RA above 1%
# Family_top %>% filter(Family_rel_abund_Sample > 0.01) %>% select(Family) %>%  unique()

#Filter genera with RA>1%, to avoid bubbles with 0 value
top_Family_pel= Family.pel %>% filter(Family%in%unique((c(Family_top_pel$Family)))) %>%  
  filter(Family_rel_abund_Sample > 0.02) 
top_Family %>% select(Family) %>% unique() 

top_Family_pel$Family <- factor(top_Family_pel$Family, levels=rev(sort(unique(top_Family$Family))))
top_Family_pel <- top_Family_pel %>% arrange(Class)

Family_bubble <- ggplot(top_Family_pel,aes(x=interaction(Polymer_Isotope,Backbone),y= Family)) +
  geom_point(aes(size=Family_rel_abund_Sample, fill = factor(Polymer)), shape = "circle filled", stroke = 1, colour = "black") +
  scale_fill_manual(values = Pal.plast) +
  scale_size(range = c(2,9.5))+
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 10))) +
  ylab("") +
  xlab("") +  
  facet_nested(Phylum + Class ~ Location + Treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins',
               nest_line = element_line(),
               strip = strip_nested(
                 background_y =  elem_list_rect(color = "grey25",
                                                fill = "white", linewidth = 1),
                 text_y = elem_list_text(size = 12, angle = 0, 
                                         color = "grey25", by_layer_y = F),
                 background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=13),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Family_bubble

#### Benthic -------------------------------------------------------------------
# select top_n genera per sample Genera for plotting
Family_top_bent <- Family.bent %>% dplyr::select(Description, Family, Family_rel_abund_Sample) %>% 
  filter ( !Family %in% c("NA", "unassigned")) %>%  #First remove unassigned genera
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Family_rel_abund_Sample, n = 5) %>% ungroup()

# #Check how much/which genera we find
# unique(Family_top$Family)
# # and with an RA above 1%
# Family_top %>% filter(Family_rel_abund_Sample > 0.01) %>% select(Family) %>%  unique()

#Filter genera with RA>1%, to avoid bubbles with 0 value
top_Family_bent= Family.bent %>% filter(Family%in%unique((c(Family_top_bent$Family)))) %>%  
  filter(Family_rel_abund_Sample > 0.02) 
top_Family_bent %>% select(Family) %>% unique() 

top_Family_bent$Family <- factor(top_Family_bent$Family, levels=rev(sort(unique(top_Family_bent$Family))))
top_Family_bent <- top_Family_bent %>% arrange(Class)

Family_bubble <- ggplot(top_Family_bent,aes(x=interaction(Polymer_Isotope,Backbone),y= Family)) +
  geom_point(aes(size=Family_rel_abund_Sample, fill = factor(Polymer)), shape = "circle filled", stroke = 1, colour = "black") +
  scale_fill_manual(values = Pal.plast) +
  scale_size(range = c(2,9.5))+
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 10))) +
  ylab("") +
  xlab("") +  
  facet_nested(Phylum + Class ~ Location + Treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins',
               nest_line = element_line(),
               strip = strip_nested(
                 background_y =  elem_list_rect(color = "grey25",
                                                fill = "white", linewidth = 1),
                 text_y = elem_list_text(size = 12, angle = 0, 
                                         color = "grey25", by_layer_y = F),
                 background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=13),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Family_bubble

#### Wild plastic-----------------------------------------------------------------
Family.wild <- tt.wild  %>%  select(Location, Habitat, Description, Phylum, Class, Family, Family_rel_abund_Sample)%>% 
  distinct() 

# select top_n genera per sample Genera for plotting
Family_top_wild <- Family.wild %>% dplyr::select(Description, Family, Family_rel_abund_Sample) %>% 
  filter ( !Family %in% c("NA", "unassigned")) %>%  #First remove unassigned genera
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Family_rel_abund_Sample, n = 8) %>% ungroup()

# #Check how much/which genera we find
# unique(Family_top$Family)
# # and with an RA above 1%
# Family_top %>% filter(Family_rel_abund_Sample > 0.01) %>% select(Family) %>%  unique()

#Filter genera with RA>1%, to avoid bubbles with 0 value
top_Family_wild= Family.wild %>% filter(Family%in%unique((c(Family_top_wild$Family)))) %>%  
  filter(Family_rel_abund_Sample > 0.01) 
top_Family_wild %>% select(Family) %>% unique() 

top_Family_wild$Family <- factor(top_Family_wild$Family, levels=rev(sort(unique(top_Family_wild$Family))))
top_Family_wild <- top_Family_wild %>% arrange(Class)

Family_bubble <- ggplot(top_Family_wild,aes(x=Description,y= Family)) +
  geom_point(aes(size=Family_rel_abund_Sample, fill = factor(Phylum)), shape = "circle filled", stroke = 1, colour = "black") +
  scale_fill_manual(values = colors_M3) +
  scale_size(range = c(2,9.5))+
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 10))) +
  ylab("") +
  xlab("") +  
  facet_nested(Phylum + Class ~ Location, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins',
               nest_line = element_line(),
               strip = strip_nested(
                 background_y =  elem_list_rect(color = "grey25",
                                                fill = "white", linewidth = 1),
                 text_y = elem_list_text(size = 12, angle = 0, 
                                         color = "grey25", by_layer_y = F),
                 background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=13),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Family_bubble




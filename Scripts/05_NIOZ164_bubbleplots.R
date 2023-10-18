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
library(ggpubr)
library(ggh4x)
library(cowplot)
library(emojifont)
library(ggthemes)
library(stringr)
library(tidyverse)


## Import data ----------------------------------------------------------------------------
tt <- read.csv('../Processed-Data/NIOZ164_EUX_discs_RA_tidy_data_decontamed_tax.correct_pruned.csv', na.strings = c(""))
head(tt)
tt$X <- NULL
# # Remove the _ from the location names and treatments for plotting
# tt.1 <- tt %>%  mutate(Location = ifelse(Location == "Crooks_Castle", "Crooks Castle",
#                                    ifelse(Location == "Charles_Brown", "Charles Brown", Location)))
# tt.1 <- tt.1 %>%  mutate(Treatment = ifelse(Treatment == "no_UV", "no UV", Treatment))
# head(tt.1)

tt.inc <- tt %>% filter(Phase == "Disc" ) 
unique(tt.inc$Location)
unique(tt.inc$Treatment)
unique(tt.inc$Polymer)

tt.wild <- tt %>% filter(Method == "Collection") 
unique(tt.wild$Location)
unique(tt.wild$Method)
unique(tt.wild$Polymer)

# tt.filt <- tt %>% filter(Phase == "Filter" ) %>% filter(!Treatment == "NA ")
# unique(tt.inc$Location)
# unique(tt.inc$Treatment)
# unique(tt.inc$Polymer)

## Colors for plotting --------------------------------------------------------------------
pal_isme <- c("#006d77", "#ffddd2", "#00C49A", "#e29578", "#83c5be")

pal.loc <- c("#FF6DB6FF" , "#004949FF",  "#66A61E")
# CB, CC, Zeelandia
pal.habs<- c("#CC5800FF", "#51C3CCFF")
# Benthic, Pelagic
pal.loc.hab <- c("#FF6DB6FF","#FFB6DBFF", "#004949FF", "#009292FF", "#66A61E")
# CB_P, CB_B, CC_P, CB_B, Zeelandia
pal.pols.isotop <- c("#E31A1C", "#7570B3", "#1F78B4","#A6CEE3", "#E6AB02","#A6761D", "#E5C494","#1B9E77")
# Blanco;Nylon;PE;PE-13C;PET;PP;PP-13C;PS;
pal.uv <- c("#DDCC77","#332288") 
# UV, noUV
colors_M1 <- c("#004e64", "#ecc8af", "#F2AF29", "#436436", "#00a5cf", 
               "#c18c5d", "#5f0f40", "#DC602E", "#495867", "#A29F15", 
               "#570000", "#FFF5B2", "#20221B", "#9fffcb", "#c08497", 
               "#8D6346", "#FF4B3E", "#149911", "#472d30")

# Generating bubbleplots --------------------------------------------------------------------

## Genus Relative Abundance --------------------------------------------------------------------
Genus <- tt.inc  %>%  select(Description, Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, Treatment, 
                             Phylum, Genus, Genus_rel_abund_Sample, Sample_rel_abund, Sample_st_dev) %>% 
  distinct() 

Genus.pel <- Genus %>% filter(Habitat == "Pelagic")  

Genus.bent <- Genus %>% filter(Habitat == "Benthic")

#Count how many unique genera we have. 
Genus.pel %>% select(Genus) %>%  unique() %>% count()
Genus.bent %>% select(Genus) %>%  unique() %>% count()

Genus %>%  group_by(Description) %>%  filter(Genus_rel_abund_Sample > 0.01) %>% 
  ungroup() %>% select(Genus) %>% unique() 

### Incubation -----------------------------------------------------------------
#### Pelagic -------------------------------------------------------------------
# select top_n genera per sample Genera for plotting
Genus_top_pel <- Genus.pel %>% dplyr::select(Description, Genus, Genus_rel_abund_Sample) %>% 
  filter ( !Genus %in% c("NA", "unassigned")) %>%  #First remove unassigned genera
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Genus_rel_abund_Sample, n = 4) %>% ungroup()

# #Check how much/which genera we find
# unique(Genus_top$Genus)
# # and with an RA above 1%
# Genus_top %>% filter(Genus_rel_abund_Sample > 0.01) %>% select(Genus) %>%  unique()

#Filter genera with RA>1%, to avoid bubbles with 0 value
top_genus_pel= Genus.pel %>% filter(Genus%in%unique((c(Genus_top_pel$Genus)))) %>%  
  filter(Genus_rel_abund_Sample > 0.005) %>% distinct()
top.gen.pel <- top_genus_pel %>% select(Genus) %>% unique() 

top_genus_pel$Genus <- factor(top_genus_pel$Genus, levels=rev(sort(unique(top_genus_pel$Genus))))
# top_genus_pel <- top_genus_pel %>% arrange(Order)

Genus_bubble_pel <- ggplot(top_genus_pel,aes(x=interaction(Polymer_Isotope,Backbone),y= Genus)) +
  geom_point(aes(size=Genus_rel_abund_Sample, fill = factor(Polymer_Isotope)), shape = "circle filled", stroke = NA, alpha = 0.5) +
  scale_fill_manual(values = pal.pols.isotop) +
  scale_size(range = c(3,9))+
  
  ylab("") +
  xlab("") +  
  facet_nested(Habitat + Phylum  ~ Location + Treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins', as.table = F, 
               nest_line = element_line(),
               strip = strip_nested( size = "variable",
                 background_y =  elem_list_rect(color = c("#51C3CCFF", rep_len("grey25",7)),
                                                fill = "white", linewidth = 1, by_layer_y = F),
                 text_y = elem_list_text(size = c(18,rep_len(13,7)), angle = c(270,rep_len(0,7)), 
                                         color = c("#51C3CCFF",rep_len("grey25",7)), by_layer_y = F),
                 background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
    axis.ticks = element_line(color = "grey50"),
    legend.text=element_text(size = 12),
    legend.title = element_text(size=14),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "right") +
  guides( x = "axis_nested",
  fill = FALSE) +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Genus_bubble_pel

#### Benthic -------------------------------------------------------------------
# select top_n genera per sample Genera for plotting
Genus_top_bent <- Genus.bent %>% dplyr::select(Description, Genus, Genus_rel_abund_Sample) %>% 
  filter ( !Genus %in% c("NA", "unassigned")) %>%  #First remove unassigned genera
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Genus_rel_abund_Sample, n = 4) %>% ungroup()

# #Check how much/which genera we find
# unique(Genus_top$Genus)
# # and with an RA above 1%
# Genus_top %>% filter(Genus_rel_abund_Sample > 0.01) %>% select(Genus) %>%  unique()

#Filter genera with RA>1%, to avoid bubbles with 0 value
top_genus_bent= Genus.bent %>% filter(Genus%in%unique((c(Genus_top_bent$Genus)))) %>%  
  filter(Genus_rel_abund_Sample > 0.005) 
top.gen.bent <- top_genus_bent %>% select(Genus) %>% unique() 

top_genus_bent$Genus <- factor(top_genus_bent$Genus, levels=rev(sort(unique(top_genus$Genus))))
# top_genus_bent <- top_genus_bent %>% arrange(Order)

Genus_bubble_bent <- ggplot(top_genus_bent,aes(x=interaction(Polymer_Isotope,Backbone),y= Genus)) +
  geom_point(aes(size=Genus_rel_abund_Sample, fill = factor(Polymer_Isotope)), shape = "circle filled", stroke = NA, alpha = 0.5) +
  scale_fill_manual(values = pal.pols.isotop) +
  scale_size(range = c(2,7))+
  ylab("") +
  xlab("") +  
  facet_nested(Habitat + Phylum  ~ Location + Treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins', as.table = F, 
               nest_line = element_line(),
               strip = strip_nested( size = "variable",
                                     background_y =  elem_list_rect(color = c("#CC5800FF", rep_len("grey25",10)),
                                                                    fill = "white", linewidth = 1, by_layer_y = F),
                                     text_y = elem_list_text(size = c(18,rep_len(13,10)), angle = c(270,rep_len(0,10)), 
                                                             color = c("#CC5800FF",rep_len("grey25",10)), by_layer_y = F),
                                     background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
    axis.ticks = element_line(color = "grey50"),
    legend.text=element_text(size = 12),
    legend.title = element_text(size=14),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "right") +
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 10, color = NA)), 
          size  = guide_legend(order = 1),
         color = guide_legend(order = 2)) +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Genus_bubble_bent

#### Combine incubations with cowplot -----------------------------------------------------------------
# legend.a <- get_legend(Genus_bubble_bent+
#                          theme(legend.direction = "horizontal",
#                                legend.title.align = 0.5))

plot_grid(Genus_bubble_pel + theme( axis.text.x = element_blank(),
                                   axis.title.x = element_blank(),ggh4x.axis.nestline.x = element_blank(),
                                   ggh4x.axis.nesttext.x =element_blank(), plot.margin = unit(c(0,0,-1,0), "cm")),
          Genus_bubble_bent + theme(plot.margin = unit(c(0,0,0,0), "cm")),
          ncol = 1,
          nrow = 2,
          align = 'v',
          axis = "ltbr",
          rel_heights = c(1,1.5))

#### Calculate sum_percentage of the top genera per incubation habitat ----------------------------
top.gen.pel <- Genus %>% filter(Genus %in% Genus_top_pel$Genus)
top.gen.pel.avg <- top.gen.pel %>% group_by(Description) %>% 
            summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

top.gen.pel.sd <- top.gen.pel %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

top.gen.pel.avg
top.gen.pel.sd


top.gen.bent <- Genus %>% filter(Genus %in% Genus_top_bent$Genus)
top.gen.bent.avg <- top.gen.bent %>% group_by(Description) %>%
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

top.gen.bent.sd <- top.gen.ben %>% group_by(Description) %>%
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

top.gen.bent.avg
top.gen.bent.sd

### Wild plastic-----------------------------------------------------------------
Genus.wild <- tt.wild  %>%  select(Location, Habitat, Description, Phylum, Order, Genus, Genus_rel_abund_Sample, Sample_rel_abund)%>% 
  distinct() 

# select top_n genera per sample Genera for plotting
Genus_top_wild <- Genus.wild %>% dplyr::select(Description, Genus, Genus_rel_abund_Sample) %>% 
  filter ( !Genus %in% c("NA", "unassigned", " ")) %>%  #First remove unassigned genera
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Genus_rel_abund_Sample, n = 5) %>% ungroup()

#Check how much/which genera we find
unique(Genus_top_wild$Genus)
# and with an RA above 1%
Genus_top_wild %>% filter(Genus_rel_abund_Sample > 0.01) %>% select(Genus) %>%  unique()

#Filter genera with RA>1%, to avoid bubbles with 0 value
top_genus_wild= Genus.wild %>% filter(Genus%in%unique((c(Genus_top_wild$Genus)))) %>%  
  filter(Genus_rel_abund_Sample > 0.01) 
top.gen.wild <- top_genus_wild %>% select(Genus) %>% unique() 

top_genus_wild$Genus <- factor(top_genus_wild$Genus, levels=rev(sort(unique(top_genus_wild$Genus))))
top_genus_wild <- top_genus_wild %>% arrange(Order)

Genus_bubble <- ggplot(top_genus_wild,aes(x=Description, y= Genus)) +
  geom_point(aes(size=Genus_rel_abund_Sample, fill = factor(Phylum)), shape = "circle filled", stroke = NA, alpha = 0.5) +
  scale_fill_manual(values = colors_M1) +
  scale_size(range = c(3,10))+
  guides( x = "axis_nested",  fill = FALSE) +
  ylab("") +
  xlab("") +  
  facet_nested(Phylum  ~ Location, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins', as.table = F,
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
    legend.text=element_text(size = 11),
    legend.title = element_text(size=12),
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

#### Calculate sum_percentage of the top genera on wild plastics ----------------------------
top.gen.wild <- Genus.wild %>% filter(Genus %in% Genus_top_wild$Genus)
top.gen.wild.avg <- top.gen.wild %>% group_by(Description) %>% 
  summarise(sum = sum(Genus_rel_abund_Sample)) %>% summarise(mean = mean(sum)) *100

top.gen.wild.sd <- top.gen.wild %>% group_by(Description) %>% 
  summarise(sd = sd(Genus_rel_abund_Sample)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) * 100 

            
top.gen.wild.avg
top.gen.wild.sd

# Family Relative Abundance --------------------------------------------------------------------
Family <- tt.inc  %>%  select(Description, Location, Habitat, Polymer_Isotope, Backbone, Treatment, 
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
  filter(Family_rel_abund_Sample > 0.02) %>% distinct()
top_Family %>% select(Family) %>% unique() 

top_Family_pel$Family <- factor(top_Family_pel$Family, levels=rev(sort(unique(top_Family_pel$Family))))
# top_Family_pel <- top_Family_pel %>% arrange(Class)

Family_bubble_pel <- ggplot(top_Family_pel,aes(x=interaction(Polymer_Isotope,Backbone),y= Family)) +
  geom_point(aes(size=Family_rel_abund_Sample, fill = factor(Polymer_Isotope)), shape = "circle filled", stroke = NA, alpha = 0.75) +
  scale_fill_manual(values = pal.pols.isotop) +
  scale_size(range = c(3,10))+
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 10, color = NA))) +
  ylab("") +
  xlab("") +  
  facet_nested(Habitat + Phylum  ~ Location + Treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins', as.table = F, 
               nest_line = element_line(),
               strip = strip_nested( size = "variable",
                                     background_y =  elem_list_rect(color = c("#51C3CCFF", rep_len("grey25",10)),
                                                                    fill = "white", linewidth = 1, by_layer_y = F),
                                     text_y = elem_list_text(size = c(18,rep_len(13,10)), angle = c(270,rep_len(0,10)), 
                                                             color = c("#51C3CCFF",rep_len("grey25",10)), by_layer_y = F),
                                     background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=14),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Family_bubble_pel

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
#top_Family_bent <- top_Family_bent %>% arrange(Class)

Family_bubble_bent<- ggplot(top_Family_bent,aes(x=interaction(Polymer_Isotope,Backbone),y= Family)) +
  geom_point(aes(size=Family_rel_abund_Sample, fill = factor(Polymer_Isotope)), shape = "circle filled", stroke = NA, alpha  = 0.75) +
  scale_fill_manual(values = pal.pols.isotop) +
  scale_size(range = c(3,10))+
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 10, color = NA))) +
  ylab("") +
  xlab("") +  
  facet_nested(Habitat + Phylum  ~ Location + Treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins', as.table = F, 
               nest_line = element_line(),
               strip = strip_nested( size = "variable",
                                     background_y =  elem_list_rect(color = c("#CC5800FF", rep_len("grey25",10)),
                                                                    fill = "white", linewidth = 1, by_layer_y = F),
                                     text_y = elem_list_text(size = c(18,rep_len(13,10)), angle = c(270,rep_len(0,10)), 
                                                             color = c("#CC5800FF",rep_len("grey25",10)), by_layer_y = F),
                                     background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=14),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Family_bubble_bent

#### Combine incubations with cowplot -----------------------------------------------------------------
legend.a <- get_legend(Family_bubble_bent+
                         theme(legend.direction = "horizontal",
                               legend.title.align = 0.5))

plot_grid(Family_bubble_pel + theme(legend.position ="none", axis.text.x = element_blank(), 
                                   axis.title.x = element_blank(),ggh4x.axis.nestline.x = element_blank(),
                                   ggh4x.axis.nesttext.x =element_blank(), plot.margin = unit(c(0,0,-1,0), "cm")),
          Family_bubble_bent + theme(legend.position ="none", plot.margin = unit(c(0,0,0,0), "cm")),
          legend.a,
          ncol = 1,
          nrow = 3,
          align = 'v',
          axis = "ltbr",
          rel_heights = c(1,1.15,0.1))

### Wild plastic-----------------------------------------------------------------
Family.wild <- tt.wild  %>%  select(Location, Habitat, Description, Phylum, Family, Family_rel_abund_Sample)%>% 
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

Family_bubble_wild <- ggplot(top_Family_wild,aes(x=Description,y= Family)) +
  geom_point(aes(size=Family_rel_abund_Sample, fill = factor(Phylum)), shape = "circle filled", stroke = NA, alpha = 0.5) +
  scale_fill_manual(values = colors_M1) +
  scale_size(range = c(3,10))+
  guides( x = "axis_nested",  fill = FALSE) +
  ylab("") +
  xlab("") +  
  facet_nested(Phylum  ~ Location, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins', as.table = F,
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
    legend.text=element_text(size = 11),
    legend.title = element_text(size=12),
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

Family_bubble_wild

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


# Set working directory ------------------------------------------------------------------
setwd("C:/Users/mgoudriaan/Documents/GitHub/Caribbean_discs_longterm/Scripts")

# Load libraries -------------------------------------------------------------------------
library(devtools)
library(phyloseq)
library(vegan)
library(microbiome)
library("ggh4x")
library(ggpubr)
library("ggsci")
library("cowplot")
library(rvg)
library(officer)
library(tidyverse)

# Import data ----------------------------------------------------------------------------
physeq_object <- readRDS("../Analysis/NIOZ164_physeq_object_decontamed.rds")

# Colors for plotting --------------------------------------------------------------------
pal.iso <- c( "#00C49A", "#e29578",)

pal.loc <- c("#FF6DB6FF" , "#004949FF",  "#66A61E")
# CB, CC, Zeelandia
pal.habs<- c("#FF8E32FF", "#51C3CCFF")
# Benthic, Pelagic
pal.loc.hab <- c("#FF6DB6FF","#FFB6DBFF", "#004949FF", "#009292FF", "#66A61E")
# CB_P, CB_B, CC_P, CB_B, Zeelandia
pal.pols.isotop <- c("#A6CEE3", "#1F78B4","#E5C494","#A6761D", "#7570B3" , "#E31A1C", "#E6AB02", "#1B9E77")
# PE;PE-13C;PP;PP-13C;PS;PET;Nylon;Blanco
pal.uv <- c("#DDCC77","#332288") 
# UV, noUV
pal.iso <- c( "#00C49A", "#e29578")

# Calculate alpha diversity indices -----------------------------------------------------
Alpha.Div <- phyloseq::estimate_richness(physeq_object, measures=c("Observed", "Simpson", "Shannon", "Chao1") )
head(Alpha.Div)

# Add metadata to the dataframe for plotting
Alpha.Discs   <- Alpha.Div %>% cbind(data.frame(physeq_object@sam_data)) 

# add new factors combined column
Alpha.Discs$Method_Phase <- str_c(Alpha.Discs$Method, "_", Alpha.Discs$Phase)
Alpha.Discs$Backbone_Isotope <- str_c(Alpha.Discs$Backbone, "_", Alpha.Discs$Isotope)
Alpha.Discs <- Alpha.Discs %>%  mutate(Polymer_Isotope = if_else(Polymer %in% c("PE", "PP"), paste(Polymer, Isotope, sep = "-"), Polymer))

#remove unecessary columns and rows
Alpha.Discs$InputFileName <- NULL
#We want to focus only on the discs themselfs and the collected polymers for the alpha diversity
Alpha.Discs <- Alpha.Discs %>%  filter( !Polymer %in% c("NA") & Habitat != "Unknown")
unique(Alpha.Discs$Polymer)
unique(Alpha.Discs$Habitat)
unique(Alpha.Discs$Method)

head(Alpha.Discs)

# Make exploratory raincloud plots to store in pptx---------------------------------------------------
## Chao1 -----------------------------------------------------------------------------------
Chao1.Discs <- read_pptx()
Chao1.Discs  <- read_pptx("../Reports/Chao1.NIOZ164.pptx")

Chao1 <- ggplot(Alpha.Discs,          #Pick data to plot
                  aes(x = Treatment, y = Chao1,  fill = Treatment)) + #Pick factors to use
  geom_boxplot(aes(alpha = 0.5), # Make boxplot
               width = 0.1, outlier.colour = NULL) +
  geom_half_violin(aes(alpha = 0.5),side = "r", nudge = 0.15, scale = "count") + # Make half-violin
  geom_half_point(aes(color = Treatment),                                       # Chose factor for points
                  side = 'l', position = position_nudge(x=-0.1), range_scale = 0.3) +
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  guides(alpha = "none") +
 facet_grid( ~ Location) +
  scale_colour_manual(values = colors_M1) +
  scale_fill_manual(values = colors_M1) 

Chao1

editable_graph <- dml(ggobj = Chao1)
Chao1.Discs  <- add_slide(Chao1.Discs ) 
Chao1.Discs <- ph_with(x = Chao1.Discs, editable_graph,location = ph_location_type(type = "body") )
print(Chao1.Discs , target = "../Reports/Chao1.NIOZ164.pptx")

## Observed Features -----------------------------------------------------------------------------------
Observed.Discs <- read_pptx()
Observed.Discs  <- read_pptx("../Reports/Observed.Features.NIOZ164.pptx")

Obs.feat <- ggplot(Alpha.Discs,          #Pick data to plot
                aes(x = Treatment, y = Observed,  fill = Treatment)) + #Pick factors to use
  geom_boxplot(aes(alpha = 0.5), # Make boxplot
               width = 0.1, outlier.colour = NULL) +
  geom_half_violin(aes(alpha = 0.5),side = "r", nudge = 0.15, scale = "count") + # Make half-violin
  geom_half_point(aes(color = Treatment),                                       # Chose factor for points
                  side = 'l', position = position_nudge(x=-0.1), range_scale = 0.3) +
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  guides(alpha = "none") +
  facet_grid( ~ Location) +
  scale_colour_manual(values = colors_M1) +
  scale_fill_manual(values = colors_M1) 

Obs.feat

editable_graph <- dml(ggobj = Obs.feat)
Observed.Discs  <- add_slide(Observed.Discs ) 
Observed.Discs <- ph_with(x = Observed.Discs, editable_graph,location = ph_location_type(type = "body") )
print(Chao1.Discs , target = "../Reports/Observed.Features.NIOZ164.pptx")

## Simpson -----------------------------------------------------------------------------------
Simpson.Discs <- read_pptx()
Simpson.Discs  <- read_pptx("../Reports/Simpson.NIOZ164.pptx")

Simpson <- ggplot(Alpha.Discs,          #Pick data to plot
                   aes(x = Treatment, y = Simpson,  fill = Treatment)) + #Pick factors to use
  geom_boxplot(aes(alpha = 0.5), # Make boxplot
               width = 0.1, outlier.colour = NULL) +
  geom_half_violin(aes(alpha = 0.5),side = "r", nudge = 0.15, scale = "count") + # Make half-violin
  geom_half_point(aes(color = Treatment),                                       # Chose factor for points
                  side = 'l', position = position_nudge(x=-0.1), range_scale = 0.3) +
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  guides(alpha = "none") +
  facet_grid( ~ Location) +
  scale_colour_manual(values = colors_M1) +
  scale_fill_manual(values = colors_M1) 

Simpson

editable_graph <- dml(ggobj = Simpson)
Simpson.Discs  <- add_slide(Simpson.Discs ) 
Simpson.Discs<- ph_with(x = Simpson.Discs, editable_graph,location = ph_location_type(type = "body") )
print(Chao1.Discs , target = "../Reports/Simpson.NIOZ164.pptx")

## Shannon -----------------------------------------------------------------------------------
Shannon.Discs <- read_pptx()
Shannon.Discs  <- read_pptx("../Reports/Shannon.NIOZ164.pptx")

Shannon <- ggplot(Alpha.Discs,          #Pick data to plot
                  aes(x = Treatment, y = Shannon,  fill = Treatment)) + #Pick factors to use
  geom_boxplot(aes(alpha = 0.5), # Make boxplot
               width = 0.1, outlier.colour = NULL) +
  geom_half_violin(aes(alpha = 0.5),side = "r", nudge = 0.15, scale = "count") + # Make half-violin
  geom_half_point(aes(color = Treatment),                                       # Chose factor for points
                  side = 'l', position = position_nudge(x=-0.1), range_scale = 0.3) +
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  guides(alpha = "none") +
  facet_grid( ~ Location) +
  scale_colour_manual(values = colors_M1) +
  scale_fill_manual(values = colors_M1) 

Shannon

editable_graph <- dml(ggobj = Shannon)
Shannon.Discs  <- add_slide(Shannon.Discs ) 
Shannon.Discs<- ph_with(x = Shannon.Discs, editable_graph,location = ph_location_type(type = "body") )
print(Chao1.Discs , target = "../Reports/Shannon.NIOZ164.pptx")

# ## Summarize alpha diversity indices in average and SD----------------------------------
# # First pito dataframe into long format
# Alpha.Discs.long <- Alpha.Discs  %>% select(!se.chao1) %>%  pivot_longer(cols = c("Observed", "Chao1", "Shannon", "Simpson"), names_to = "Diversity_Index", 
#                                                                          values_to = "Value")
# head(Alpha.Discs.long)
# #load needed formula
# source("summarySE.R")
# 
# # Summarize and calculate avg and SD
# # Define groups per which you want to have the mean and SD and you want to keep for plotting
# Alpha.Discs.long.c <-  summarySE(Alpha.Discs.long, measurevar="Value", groupvars=c("Location", "Habitat","Polymer", "Isotope", "Backbone", "Treatment", 
#                                                                                    "Location_Habitat", "Diversity_Index"))
# # Select richness and diversity seperately for plotting
# Richness.long.c <- Alpha.Discs.long.c %>% filter(Diversity_Index %in% c("Observed", "Chao1"))
# Diversity.long.c <- Alpha.Discs.long.c %>% filter(Diversity_Index %in% c("Shannon", "Simpson"))

# Make final boxplot plots -------------------------------------------------------------
### 1. Alpha diversity for habitats per location                                    ----
#%#----------------------------------------------------------------------------------#%#

## Observed ----------------------------------------------------------------------------
Observed <- ggplot(Alpha.Discs,         #Pick data to plot
                aes(x= Location, y = Observed, color = Location)) + #Pick factors to use
 geom_boxplot(stat = "boxplot", outlier.colour =  NULL, linewidth = 1) +
  geom_point( position = "jitter", size = 3, alpha = 0.8) +
  theme_pubclean()+
  scale_y_continuous(position = 'right', sec.axis = dup_axis()) +
  theme_pubclean()+
  theme(axis.text.x = element_text(size = 13, angle = 60, hjust = 1), 
        axis.text.y.right=element_text(size= 13), 
        axis.text.y.left=element_blank(),
        legend.text=element_text(size = 13),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y.left = element_text(size=16),
        axis.title.y.right =  element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  scale_colour_manual(values = pal.loc) +
  scale_fill_manual(values = pal.loc) +
  xlab("") +
  labs( title = "", color = "Location")

Observed

## Chao1 -----------------------------------------------------------------------------
Chao1 <- ggplot(Alpha.Discs,         #Pick data to plot
                aes(x= Location, y = Chao1, color = Location)) + #Pick factors to use
  geom_boxplot(stat = "boxplot", outlier.colour =  NULL, linewidth = 1) +
  geom_point( position = "jitter", size = 3, alpha = 0.8) +
  theme_pubclean()+
  scale_y_continuous(position = 'right', sec.axis = dup_axis()) +
  theme_pubclean()+
  theme(axis.text.x = element_text(size = 13, angle = 60, hjust = 1), 
        axis.text.y.right=element_text(size= 13), 
        axis.text.y.left=element_blank(),
        legend.text=element_text(size = 13),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y.left = element_text(size=16),
        axis.title.y.right =  element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  scale_colour_manual(values = pal.loc) +
  scale_fill_manual(values = pal.loc) +
  xlab("") +
  labs( title = "", color = "Location")

Chao1

## Simpson --------------------------------------------------------------------------
Simpson <- ggplot(Alpha.Discs,         #Pick data to plot
                  aes(x= Location, y = Simpson, color = Location)) + #Pick factors to use
  geom_boxplot(stat = "boxplot", outlier.colour =  NULL, linewidth = 1) +
  geom_point( position = "jitter", size = 3, alpha = 0.8) +
  scale_y_continuous(position = 'right', sec.axis = dup_axis()) +
  theme_pubclean()+
  theme(axis.text.x = element_text(size = 13, angle = 60, hjust = 1), 
        axis.text.y.right=element_text(size= 13), 
        axis.text.y.left=element_blank(),
        legend.text=element_text(size = 13),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y.left = element_text(size = 16),
        axis.title.y.right =  element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  scale_colour_manual(values = pal.loc) +
  scale_fill_manual(values = pal.loc) +
  xlab("") +
  labs( title = "", color = "Location")

Simpson

## Shannon -------------------------------------------------------------------------

Shannon <- ggplot(Alpha.Discs,         #Pick data to plot
                  aes(x= Location, y = Shannon, color = Location)) + #Pick factors to use
  geom_boxplot(stat = "boxplot", outlier.colour =  NULL, linewidth = 1) +
  geom_point( position = "jitter", size = 3, alpha = 0.8) +
  scale_y_continuous(position = 'right', sec.axis = dup_axis()) +
  theme_pubclean()+
  theme(axis.text.x = element_text(size = 13, angle = 60, hjust = 1), 
        axis.text.y.right=element_text(size= 13), 
        axis.text.y.left=element_blank(),
        legend.text=element_text(size = 13),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y.left = element_text(size=16),
        axis.title.y.right =  element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.key = element_rect(fill = NA)) +
  scale_colour_manual(values = pal.loc) +
  scale_fill_manual(values = pal.loc) +
  xlab("") +
  labs( title = "", color = "Location")

Shannon

## Combine panels with cowplot -------------------------------------------------- 
legend.a <- get_legend(Shannon+
                         theme(legend.direction = "horizontal"))

plot_grid(# Observed + theme(legend.position ="none"),
          Chao1 + theme(legend.position ="none"),
          Simpson + theme(legend.position ="none"),
          Shannon + theme(legend.position ="none"),
          NULL,
          legend.a,
          ncol = 3,
          nrow = 2,
          align = 'v',
          axis = "l",
          rel_heights = c(1,0.1),
          rel_widths = c(1,1))


# Make final boxplot plots -------------------------------------------------------------
### 2. Incubation only - Alpha diversity for treatments per habitat/location        ----
### 3. Incubation only - Alpha diversity for 12C/13C PE/PP per habitat and location ----                          
#%#----------------------------------------------------------------------------------#%#

#Here we only want to focus on the discs
Alpha.Discs.inc <- Alpha.Discs %>%  filter( Habitat %in% c("Benthic", "Pelagic")) %>%  filter( Polymer != "B")
unique(Alpha.Discs.inc$Polymer)
unique(Alpha.Discs.inc$Habitat)
unique(Alpha.Discs.inc$Method)
colnames(Alpha.Discs.inc)

Alpha.Discs.isotopes <- Alpha.Discs %>%  filter( Polymer %in% c("PE", "PP"))
colnames(Alpha.Discs.isotopes)

## Observed ----------------------------------------------------------------------------
Observed <- ggplot(Alpha.Discs.isotopes,          #Pick data to plot
                aes(x = Isotope, y = Observed,  color = Isotope)) + #Pick factors to use
  geom_boxplot(stat = "boxplot", outlier.colour =  NULL, linewidth = 1) +
  geom_point( position = "jitter", size = 3, alpha = 0.8) +
  facet_grid(fct_relevel(Habitat, "Pelagic", "Benthic") ~ Location, 
             switch = 'y') +
  scale_y_continuous(position = 'right') +
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  guides(alpha = "none") +
  scale_colour_manual(values = pal.iso) +
  scale_fill_manual(values = pal.iso) 

Observed

## Chao1 ----------------------------------------------------------------------------
Chao1 <- ggplot(Alpha.Discs.isotopes,          #Pick data to plot
                   aes(x = Isotope, y = Chao1,  color =  Isotope)) + #Pick factors to use
  geom_boxplot(stat = "boxplot", outlier.colour =  NULL, linewidth = 1) +
  geom_point( position = "jitter", size = 3, alpha = 0.8) +
  facet_grid(fct_relevel(Habitat, "Pelagic", "Benthic") ~ Location, 
             switch = 'y') +
  scale_y_continuous(position = 'right') +
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  guides(alpha = "none") +
  scale_colour_manual(values = pal.iso) +
  scale_fill_manual(values = pal.iso) 

Chao1

## Simpson ----------------------------------------------------------------------------
Simpson <- ggplot(Alpha.Discs.isotopes,          #Pick data to plot
                aes(x =  Isotope, y = Simpson,  color =  Isotope)) + #Pick factors to use
  geom_boxplot(stat = "boxplot", outlier.colour =  NULL, linewidth = 1) +
  geom_point( position = "jitter", size = 3, alpha = 0.8) +
  facet_grid(fct_relevel(Habitat, "Pelagic", "Benthic") ~ Location, 
             switch = 'y') +
  scale_y_continuous(position = 'right') +
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  guides(alpha = "none") +
  scale_colour_manual(values = pal.iso) +
  scale_fill_manual(values = pal.iso) 

Simpson

## Shannon ----------------------------------------------------------------------------
Shannon <- ggplot(Alpha.Discs.isotopes,          #Pick data to plot
                  aes(x =  Isotope, y = Shannon,  color =  Isotope)) + #Pick factors to use
  geom_boxplot(stat = "boxplot", outlier.colour =  NULL, linewidth = 1) +
  geom_point( position = "jitter", size = 3, alpha = 0.8) +
  facet_grid(fct_relevel(Habitat, "Pelagic", "Benthic") ~ Location, 
             switch = 'y') +
  scale_y_continuous(position = 'right') +
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  guides(alpha = "none") +
  scale_colour_manual(values = pal.iso) +
  scale_fill_manual(values = pal.iso) 

Shannon

## Combine panels with cowplot -------------------------------------------------- 
legend.a <- get_legend(Shannon+
                         theme(legend.direction = "horizontal",
                               legend.title.align = 0.5))

plot_grid(#Observed + theme(legend.position ="none"),
          Chao1 + theme(legend.position ="none"),
          Simpson + theme(legend.position ="none"),
          Shannon + theme(legend.position ="none"),
          NULL,
          legend.a,
          ncol = 3,
          nrow = 2,
          align = 'v',
          axis = "l",
          rel_heights = c(1,0.1),
          rel_widths = c(1,1))


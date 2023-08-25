##%#####################################################################################%##
# NIOZ164 Statia Discs - FA SIP Fingerprint Plots                                     #####                                            
# Author: Maaike Goudriaan, NIOZ, MMB                                                     #
#                                                                                         #
# Purpose:                                                                                # 
# NIOZ164 is an illumina sequencing lane consisting of different samplesets.              # 
#                                                                                         #
# Discs covered in polymer films that were incubated at 2 different                       #
# locations in 2 different waterdepths (seafloor and watercolumn) in                      #
# coastal Caribbean waters close to the island of St. EUstatius.                          #
# The tested polymers were PE, PP, PS, PEt and Nylon, PE-13C and PP-13C.                  #
# Of each polymere there was an UV pretreated and a non-treated version.                  #
# In addition to amplicon analysis, we performed SIP by extracting membrane FA's          #
# and analyzing the extracts with FID, MS and IRMS. This script is used to make plots of  #
# the resulting fingerprints                                                              #
##%#####################################################################################%##

# Date: 2023 - 08 - 24
# R-version: 4.3.1 

# Set working directory ------------------------------------------------------------------
setwd("C:/Users/mgoudriaan/Documents/GitHub/Caribbean_discs_longterm/Scripts")

# Load libraries -------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library("ggh4x")
library(ggpubr)
library("ggsci")
library("cowplot")

# Import data ----------------------------------------------------------------------------
Areas <- read.delim("../Data/Lipid_areas_corrected_signals.txt", sep = '\t', dec = ".")
d13C <- read.delim("../Data/Lipid_d13C_corrected_signals.txt", sep = '\t', na.strings = c(" "), dec = ".")
metadata  <- read.delim("../Data/Metadata_lipids_R_plots.txt", sep = '\t')

# Create one tibble for easy plotting ----------------------------------------------------
# Remove the row of the standard and pivot
Areas.long <- Areas %>% filter(Fatty.Acid != "C19:0 (standard!)") %>% pivot_longer( -Fatty.Acid, names_to = "Sample", values_to = "Area")
Areas.long.ra <- Areas.long %>%  group_by(Sample) %>%  mutate(Rel.Abund = Area/sum(Area))

d13C.long <- d13C %>% filter(Fatty.Acid != "C19:0 (standard!)") %>% pivot_longer( -Fatty.Acid, names_to = "Sample", values_to = "d13C")

# Combine tibbles and add metadata
df <- inner_join(Areas.long.ra, d13C.long)
df.c <- df %>% select(-Area) %>%  pivot_longer(Rel.Abund:d13C,
                                               names_to = "Variable",
                                               values_to = "Value")

df.m <-inner_join(df.c, metadata, by = "Sample") 
colnames(df.m)
unique(df.m$Location)

# Plot Charles Brown ---------------------------------------------------------------------
df.cb <- df.m %>% filter(Location == "Charles Brown")

CB <- ggplot(df.cb) +         #Pick data to plot
  geom_bar(data = subset(df.cb, Variable =="Rel.Abund"), aes(x = Fatty.Acid, y = Value),
           stat="identity", position="stack")+
  geom_point(data = subset(df.cb, Variable =="d13C"), aes(x = Fatty.Acid, y = Value))+
  scale_y_continuous(position = "right") +
  facet_nested(fct_relevel(Habitat, 'Pelagic', 'Benthic')  + Variable ~ Polymer + fct_relevel(Treatment, "UV", "noUV"),
                drop = T, scale = "free_y",
                axes = 'margins', switch = "y", 
               as.table = T, 
                nest_line = element_line()) +
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
        panel.grid.major.x = element_blank()) 
CB

# #ratio for logarithmic d13C y-axis
# ylim.prim = c(0,0.5)
# ylim.sec = c(-50, 5000)
# b <- diff(ylim.prim)/diff(ylim.sec)
# a <- ylim.prim[1] - b*ylim.sec[1]

# CB <- ggplot(df.cb) +         #Pick data to plot
#   geom_bar(aes(x = Fatty.Acid, y = Rel.Abund),
#            stat="identity", position="stack")+
#   geom_point(aes(x = Fatty.Acid, y = d13C)) +
#   scale_y_continuous(
#     name = "Relative Abundance",
#     sec.axis = sec_axis(~ ./0.1, name = "d13C")
#   ) +
#   facet_nested(fct_relevel(Habitat, 'Pelagic', 'Benthic') + fct_relevel(Treatment, "UV", "noUV")  ~ Polymer,
#                drop = T, scale = "free_y",
#                axes = 'margins',switch = "y", 
#                nest_line = element_line()) +
#   
#   theme_pubclean()+
#   theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
#         axis.text.y=element_text(size= 12), 
#         legend.text=element_text(size = 12),
#         legend.title = element_text(size=15, face = "bold"),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         strip.text.x = element_text(size = 13),
#         plot.title = element_blank(),
#         panel.border = element_rect(color = "grey90", fill = NA),
#         panel.grid.major.y = element_line(color = "grey90", linetype = 3),
#         panel.grid.major.x = element_blank()) 
# CB

# Plot Crooks Castle ---------------------------------------------------------------------
df.cc <- df.m %>% filter(Location == "Crooks Castle ")

CC <- ggplot(df.cc) +         #Pick data to plot
  geom_bar(data = subset(df.cc, Variable =="Rel.Abund"), aes(x = Fatty.Acid, y = Value),
           stat="identity", position="stack")+
  geom_point(data = subset(df.cc, Variable =="d13C"), aes(x = Fatty.Acid, y = Value))+
  scale_y_continuous(position = "right") +
  facet_nested(fct_relevel(Habitat, 'Pelagic', 'Benthic')  + Variable ~ Polymer + fct_relevel(Treatment, "UV", "noUV"),
               drop = T, scale = "free_y",
               axes = 'margins', switch = "y", 
               as.table = T, 
               nest_line = element_line()) +
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
        panel.grid.major.x = element_blank()) 
CC

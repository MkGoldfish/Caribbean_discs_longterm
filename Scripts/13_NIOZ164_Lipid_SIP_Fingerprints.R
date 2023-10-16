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
library("ggh4x")
library(ggpubr)
library("ggsci")
library("cowplot")
library(stringr)
library(tidyverse)

# Import data ----------------------------------------------------------------------------
Areas <- read.delim("../Data/Lipid_areas_corrected_signals.txt", sep = '\t', dec = ".")
d13C <- read.delim("../Data/Lipid_d13C_corrected_signals.txt", sep = '\t', na.strings = c(" "), dec = ".")
metadata  <- read.delim("../Data/Metadata_lipids_R_plots.txt", sep = '\t')
Biomass_average <- read.delim("../Data/Abundance_weigthed_average_13c_lipid.txt", sep = '\t')

# Create one tibble for easy plotting ----------------------------------------------------
# Remove the row of the standard and pivot
Areas.long <- Areas %>% filter(Fatty.Acid != "C19:0 (standard!)") %>% pivot_longer( -Fatty.Acid, names_to = "Sample", values_to = "Area")
Areas.long.ra <- Areas.long %>%  group_by(Sample) %>%  mutate(Rel.Abund = Area/sum(Area))

d13C.long <- d13C %>% filter(Fatty.Acid != "C19:0 (standard!)") %>% pivot_longer( -Fatty.Acid, names_to = "Sample", values_to = "d13C")


# Combine tibbles and add metadata
df <- inner_join(Areas.long.ra, d13C.long)
df.c <- df  %>%  pivot_longer(Rel.Abund:d13C,
                                               names_to = "Variable",
                                               values_to = "Value")

df.m <-inner_join(df.c, metadata, by = "Sample") 
colnames(df.m)
df.m$Fatty.Acid <- gsub( "w", "\u{03C9}", df.m$Fatty.Acid)
# df.m$Variable <- gsub( "Rel.Abund", "Relative Abundance", df.m$Variable)
unique(df.m$Polymer)


# Plot Charles Brown ---------------------------------------------------------------------
df.cb <- df.m %>% filter(Location == "Charles Brown")

var.labs <- as_labeller(c(Rel.Abund = "Relative~Abundance", d13C = "\u{03b4}^13~C", 
                          Pelagic = "Pelagic", Benthic = "Benthic",
                          UV = "UV", noUV = "noUV",
                          PE = "phantom()^12~C -PE", "PE-13C" = "phantom()^13~C -PE", Si = "Si~disc",
                          PP = "phantom()^12~C -PP", "PP-13C" = "phantom()^13~C -PP"), default = label_parsed)


CB <- ggplot(df.cb) +         #Pick data to plot
  geom_bar(data = subset(df.cb, Variable == "Rel.Abund"), aes(x = Fatty.Acid, y = Value),
           stat="identity", position="stack")+
  geom_point(data = subset(df.cb, Variable =="d13C"), aes(x = Fatty.Acid, y = Value))+
  scale_y_continuous(position = "right") +
  facet_nested(fct_relevel(Habitat, 'Pelagic', 'Benthic')  + Variable ~ Polymer + fct_relevel(Treatment, "UV", "noUV"),
                drop = T, scale = "free_y",
                axes = 'margins', switch = "y", 
               as.table = T, 
                nest_line = element_line(),
               labeller = var.labs)+
                 theme_pubclean()+
   theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 13), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  labs (y = " ", x = "Fatty Acid")

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
  geom_bar(data = subset(df.cc, Variable == "Rel.Abund"), aes(x = Fatty.Acid, y = Value),
           stat="identity", position="stack")+
  geom_point(data = subset(df.cc, Variable =="d13C"), aes(x = Fatty.Acid, y = Value))+
  scale_y_continuous(position = "right") +
  facet_nested(fct_relevel(Habitat, 'Pelagic', 'Benthic')  + Variable ~ Polymer + fct_relevel(Treatment, "UV", "noUV"),
               drop = T, scale = "free_y",
               axes = 'margins', switch = "y", 
               as.table = T, 
               nest_line = element_line(),
               labeller = var.labs)+
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 13), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  labs (y = " ", x = "Fatty Acid")

CC


# Plot average biomass value ---------------------------------------------------------------------
df.bm.avg <-  inner_join(Biomass_average, metadata, by = "Sample", suffix = c("",".y")) 
colnames(df.bm.avg)
df.bm.avg$Description.y <-NULL
df.bm.avg$Location.y <-NULL

df.bm.avg$Polymer_Treatment <- str_c(df.bm.avg$Polymer, "_", df.bm.avg$Treatment)
# 
pol.treat <- unique(df.bm.avg$Polymer_Treatment) %>% as.character()
df.bm.avg$Polymer_Treatment <- factor(df.bm.avg$Polymer_Treatment,levels=pol.treat)
levels(df.bm.avg$Polymer_Treatment) <- c("PE_noUV", "PE_UV", "PE^13*CnoUV", "PE^13*C_UV",
                                         "PP_noUV", "PP_UV", "PP^13*C_noUV", "PP^13*C_UV", "Si_noUV" )


ggplot(df.bm.avg) +
  geom_col(aes(x = Polymer_Treatment, y = d13C.FA)) +
  geom_hline(yintercept = 0,linewidth = 1, colour = "black") +
    coord_flip() +
    facet_nested(Location + fct_relevel(Habitat, 'Pelagic', 'Benthic') ~., drop = T,
             axes = 'margins',
             scale = "free",
             as.table = T) +
  scale_x_discrete(labels = parse(text = levels(df.bm.avg$Polymer_Treatment))) +
  theme_pubclean() +
  labs(y = expression("Abundance weighted average \u{03b4}"^13* "C fatty acids (\u2030)"), x = " ") +
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 13), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) 

  

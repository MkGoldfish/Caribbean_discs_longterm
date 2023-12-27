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
dim(Biomass_average)

showtext::showtext_opts(dpi=500)

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
unique(df.m$Fatty.Acid)

# Plot Charles Brown ---------------------------------------------------------------------
df.cb <- df.m %>% filter(Location == "Charles Brown")

var.labs <- as_labeller(c(Rel.Abund = "Fractional~Abundance", d13C = "\u{03b4}^13~C", 
                          Pelagic = "Pelagic", Benthic = "Benthic",
                          UV = "UV", noUV = "noUV",
                          PE = "phantom()^12~C -PE", "PE-13C" = "phantom()^13~C -PE", Si = "Si~disc",
                          PP = "phantom()^12~C -PP", "PP-13C" = "phantom()^13~C -PP"), default = label_parsed)


CB <- ggplot(df.cb) +         #Pick data to plot
  geom_bar(data = subset(df.cb, Variable == "Rel.Abund"), aes(x = Fatty.Acid, y = Value),
           stat="identity", position="stack")+
  geom_point(data = subset(df.cb, Variable =="d13C"), aes(x = Fatty.Acid, y = Value))+
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = c("C14:1", "C14:0", "C15-ω5c",  "C15:1","C15:0", "C16:1", "C16:1-ω7c",
                              "C16:0", "C17:1", "C17:0", "C18:1-ω9c", "C18:1-ω9t", "C18:0 ")) +
  facet_nested(fct_relevel(Habitat, 'Pelagic', 'Benthic')  + Variable ~ Polymer + fct_relevel(Treatment, "UV", "noUV"),
                drop = T, scale = "free_y",
                axes = 'margins', switch = "y", 
               as.table = T, 
                nest_line = element_line(),
               labeller = var.labs)+
                 theme_pubclean()+
   theme(axis.text.x=element_text(size = 8, angle = 60, hjust = 1), 
        axis.text.y=element_text(size = 10), 
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        strip.text.x = element_text(size = 10, margin = margin(0.075,0.1,0.075,0.1, "cm")),
        strip.text.y = element_text(size = 10, margin = margin(0.075,0.1,0.075,0.1, "cm")),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(1, 'pt')) +
  labs (y = " ", x = "Fatty Acid")

CB

ggsave("CharlesBrown_lipid_fingerprint.pdf", 
       width = 35, height  = 18, unit = "cm", 
       dpi = 500, bg = 'white')


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
  scale_x_discrete(limits = c("C14:1", "C14:0", "C15-ω5c",  "C15:1","C15:0", "C16:1", "C16:1-ω7c",
                             "C16:0", "C17:1", "C17:0", "C18:1-ω9c", "C18:1-ω9t", "C18:0 ")) +
  facet_nested(fct_relevel(Habitat, 'Pelagic', 'Benthic')  + Variable ~ Polymer + fct_relevel(Treatment, "UV", "noUV"),
               drop = T, scale = "free_y",
               axes = 'margins', switch = "y", 
               as.table = T, 
               nest_line = element_line(),
               labeller = var.labs)+
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 8, angle = 60, hjust = 1), 
        axis.text.y=element_text(size = 10), 
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        strip.text.x = element_text(size = 10, margin = margin(0.075,0.1,0.075,0.1, "cm")),
        strip.text.y = element_text(size = 10, margin = margin(0.075,0.1,0.075,0.1, "cm")),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(1, 'pt')) +
  labs (y = " ", x = "Fatty Acid")

CC

ggsave("CrooksCastle_lipid_fingerprint.pdf", 
       width = 35, height  = 18, unit = "cm", 
       dpi = 500, bg = 'white')


# Plot average biomass value ---------------------------------------------------------------------
df.bm.avg <-  inner_join(Biomass_average, metadata, by = "Sample", suffix = c("",".y")) 

summary(df.bm.avg$Description == df.bm.avg$Description.y)
summary(df.bm.avg$Location == df.bm.avg$Location.y)
colnames(df.bm.avg)
df.bm.avg$Description.y <-NULL
df.bm.avg$Location.y <-NULL

df.bm.avg$Polymer_Treatment <- str_c(df.bm.avg$Polymer, "_", df.bm.avg$Treatment)
# df.bm.avg$Polymer_Treatment <- gsub( "-13C", "^13*C", df.bm.avg$Polymer_Treatment)
# 
# 
# pol.treat <- unique(df.bm.avg$Polymer_Treatment) %>% as.character()
# df.bm.avg$Polymer_Treatment <- factor(df.bm.avg$Polymer_Treatment,levels=pol.treat)
# levels(df.bm.avg$Polymer_Treatment) <- c("PE_noUV", "PE_UV", "PE^13*C_noUV", "PE^13*C_UV",
#                                          "PP_noUV", "PP_UV", "PP^13*C_noUV", "PP^13*C_UV", "Si_noUV" )

df.bm.avg$Polymer_Treatment <- factor(df.bm.avg$Polymer_Treatment, levels=rev(sort(unique(df.bm.avg$Polymer_Treatment))))

ggplot(df.bm.avg) +
  geom_col(aes(x = Polymer_Treatment, y = d13C.FA)) +
  geom_hline(yintercept = 0,linewidth = 0.5, colour = "black") +
    coord_flip() +
    facet_nested(Location + fct_relevel(Habitat, 'Pelagic', 'Benthic') ~., drop = T,
             axes = 'margins',
             scale = "free",
             as.table = T) +
  
  theme_pubclean() +
  labs(y = expression("Abundance weighted average \u{03b4}"^13* "C fatty acids (\u2030)"), x = " ") +
  theme(axis.text.x=element_text(size = 10, angle = 60, hjust = 1), 
        axis.text.y=element_text(size = 9), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.x = element_text(size = 10, margin = margin(0.075,0.1,0.075,0.1, "cm")),
        strip.text.y = element_text(size = 10, margin = margin(0.075,0.1,0.075,0.1, "cm")),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(1, 'pt')) 

ggsave("Avg_d13C_biomass.eps", 
       width = 10, height  = 16, unit = "cm", 
       dpi = 500, bg = 'white')
  

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
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(ggh4x)



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
                      "#662C91","#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )

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

## Import data ----------------------------------------------------------------------------
tt <- read.csv('../Processed-Data/NIOZ164_EUX_discs_tidy_data_decontamed_pruned_RA.csv', row.names = NULL)
head(tt)
tt$X <- NULL
head(tt)

tt.inc <- tt %>% filter(Phase == "Disc" ) %>% filter(!Treatment == "NA ")
unique(tt.inc$Location)
unique(tt.inc$Treatment)
unique(tt.inc$Polymer)

tt.wild <- tt %>% filter(Method == "Collection") 
unique(tt.wild$Location)
unique(tt.wild$Method)
unique(tt.wild$Polymer)

tt.filt <- tt %>% filter(Phase == "Filter" ) %>% filter(!Treatment == "NA ")
unique(tt.inc$Location)
unique(tt.inc$Treatment)
unique(tt.inc$Polymer)

# Generating barplots --------------------------------------------------------------------

## Kingdom Relative Abundance --------------------------------------------------------------------
### Incubation -----------------------------------------------------------------
Kingdom <- tt.inc  %>%  select(Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, Treatment, Kingdom, Kingdom_rel_abund_Sample)%>% 
  distinct() 

gg_King <- ggplot(Kingdom, aes(x=interaction(Polymer_Isotope, Backbone), y= Kingdom_rel_abund_Sample, 
                         fill=forcats::fct_reorder(Kingdom,Kingdom_rel_abund_Sample )))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(pal_isme)) + 
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
  facet_nested(fct_relevel(Habitat, "Pelagic", "Benthic") ~ Location + Treatment, drop = T, 
                axes = 'margins',
                nest_line = element_line(),
                strip = strip_nested(
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                )) + 
  
  theme_classic2() + 
  theme(
    axis.text.x=element_text( size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text( size= 12), 
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
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Material") + labs(fill = "Kingdom")

gg_King

### Wild plastic-----------------------------------------------------------------
Kingdom <- tt.wild  %>%  select(Location, Habitat,Description, Kingdom, Kingdom_rel_abund_Sample)%>% 
  distinct() 

gg_King <- ggplot(Kingdom, aes(x=Description, y= Kingdom_rel_abund_Sample, 
                               fill=forcats::fct_reorder(Kingdom,Kingdom_rel_abund_Sample )))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(pal_isme)) + 
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
  facet_nested(  ~ Location + Habitat, drop = T,
                axes = 'margins',
                nest_line = element_line(),
                strip = strip_nested(
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                )) +
  
  theme_classic2() + 
  theme(
    axis.text.x=element_text( size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text( size= 12), 
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
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Material") + labs(fill = "Kingdom")

gg_King

### Wild plastic-----------------------------------------------------------------
Kingdom <- tt.filt  %>%  select(Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, Treatment, Kingdom, Kingdom_rel_abund_Sample)%>%
  distinct() 

gg_King <- ggplot(Kingdom, aes(x=interaction(Polymer_Isotope, Backbone), y= Kingdom_rel_abund_Sample, 
                               fill=forcats::fct_reorder(Kingdom,Kingdom_rel_abund_Sample ))) +
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(pal_isme)) + 
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
  facet_nested( . ~ Location + Treatment, drop = T, 
                axes = 'margins',
                nest_line = element_line(),
                strip = strip_nested(
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                )) + 
  
  theme_classic2() + 
  theme(
    axis.text.x=element_text( size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text( size= 12), 
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
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Material") + labs(fill = "Kingdom")

gg_King

## Phylum Relative Abundance --------------------------------------------------------------------
### Incubation -----------------------------------------------------------------
Phyla <- tt.inc  %>% filter(Kingdom != "Eukaryota") %>% select(Description, Location, Habitat, Polymer_Isotope, 
                                                               Polymer, Backbone, Treatment, Phylum, Phylum_rel_abund_Sample) %>% 
  distinct() 

Phyla_filtAb <- Phyla %>% mutate(Phylum = ifelse(Phylum_rel_abund_Sample<0.02, "others<2.0%", Phylum))

gg_Phylum<- ggplot(Phyla_filtAb, aes(x=interaction(Polymer_Isotope, Backbone), y= Phylum_rel_abund_Sample, 
                               fill=forcats::fct_reorder(Phylum, Phylum_rel_abund_Sample )))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = colors_M1) + 
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
  facet_nested( fct_relevel(Habitat, "Pelagic", "Benthic") ~ Location + Treatment, drop = T, 
                axes = 'margins',
                nest_line = element_line(),
                strip = strip_nested(
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                )) + 
  theme_classic2() + 
  theme(
    axis.text.x=element_text( size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text( size= 12), 
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
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Polymer") + labs(fill = "Phylum")

gg_Phylum

### Wild plastic-----------------------------------------------------------------
Phyla <- tt.wild %>% filter(Kingdom != "EUkaryota") %>% select(Description, Location, Habitat, Polymer_Isotope, 
                                                               Polymer, Backbone, Treatment, Phylum, Phylum_rel_abund_Sample) %>% 
 
  distinct() 

Phyla_filtAb <- Phyla %>% mutate(Phylum = ifelse(Phylum_rel_abund_Sample<0.02, "others<2.0%", Phylum))

gg_Phylum<- ggplot(Phyla_filtAb, aes(x=Description, y= Phylum_rel_abund_Sample, 
                                    fill=forcats::fct_reorder(Phylum, Phylum_rel_abund_Sample )))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = colors_M1) + 
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
  facet_nested(  ~ Location + Treatment, drop = T, 
                axes = 'margins',
                nest_line = element_line(),
                strip = strip_nested(
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                )) + 
  theme_classic2() + 
  theme(
    axis.text.x=element_text( size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text( size= 12), 
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
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Material") + labs(fill = "Phylum")

gg_Phylum

### Filters -----------------------------------------------------------------
Phyla <- tt.filt %>% filter(Kingdom != "EUkaryota") %>% select(Location, Habitat, Polymer, Isotope, 
                                                               Polymer_Isotope, Backbone, Treatment, Phylum, Phylum_rel_abund_Sample)%>% 
  distinct() 

Phyla_filtAb <- Phyla %>% mutate(Phylum = ifelse(Phylum_rel_abund_Sample<0.02, "others<2.0%", Phylum))

gg_Phylum <- gg_Phylum<- ggplot(Phyla_filtAb, aes(x=interaction(Polymer_Isotope, Backbone), y= Phylum_rel_abund_Sample, 
                                                  fill=forcats::fct_reorder(Phylum, Phylum_rel_abund_Sample )))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = colors_M2) + 
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
  facet_nested( . ~ Location + Treatment, drop = T, 
                axes = 'margins',
                nest_line = element_line(),
                strip = strip_nested(
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                )) + 
  
  theme_classic2() + 
  theme(
    axis.text.x=element_text( size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text( size= 12), 
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
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Material") + labs(fill = "Phylum")

gg_Phylum

## Order Relative Abundance --------------------------------------------------------------------
### Incubation -----------------------------------------------------------------
Orders <- tt.inc  %>% filter(Kingdom != "EUkaryota") %>%  select(Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, 
                                                                 Treatment, Order, Order_rel_abund_Sample)%>% 
  distinct() 

Order_filtAb <- Orders %>% mutate(Order = ifelse(Order_rel_abund_Sample<0.03, "others<3%", Order))

gg_Order<- ggplot(Order_filtAb, aes(x=interaction(Polymer_Isotope, Backbone), y= Order_rel_abund_Sample, 
                                     fill=forcats::fct_reorder(Order, Order_rel_abund_Sample )))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(colors_by_Maaike)) + 
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
  facet_nested( fct_relevel(Habitat, "Pelagic", "Benthic") ~ Location + Treatment, drop = T, 
                axes = 'margins',
                nest_line = element_line(),
                strip = strip_nested(
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                )) + 
  theme_classic2() + 
  theme(
    axis.text.x=element_text( size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text( size= 12), 
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
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Polymer") + labs(fill = "Order")

gg_Order

### Wild plastic-----------------------------------------------------------------
Orders <- tt.wild %>% filter(Kingdom != "Eukaryota") %>%   select(Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, 
                                                                Treatment, Order, Order_rel_abund_Sample) %>% 
  dictinct()

Order_filtAb <- Orders %>% mutate(Order = ifelse(Order_rel_abund_Sample<0.03, "others<3%", Order))

gg_Order<- ggplot(Order_filtAb, aes(x=Description, y= Order_rel_abund_Sample, 
                                      fill=forcats::fct_reorder(Order, Order_rel_abund_Sample)))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(colors_by_Maaike)) + 
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
  facet_nested(  ~ Location + Habitat, drop = T,
                 axes = 'margins',
                 nest_line = element_line(),
                 strip = strip_nested(
                   background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                 )) +
  
  theme_classic2() + 
  theme(
    axis.text.x=element_text( size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text( size= 12), 
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
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Material") + labs(fill = "Order")

gg_Order

### Filters -----------------------------------------------------------------
Orders <- tt.filt %>% filter(Kingdom != "EUkaryota")  %>%   select(Description, Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, 
                                                                   Treatment, Order, Order_rel_abund_Sample) %>% 
  distinct() 

Order_filtAb <- Orders %>% mutate(Order = ifelse(Order_rel_abund_Sample<0.03, "others<3%", Order))

gg_Order<- ggplot(Order_filtAb, aes(x=interaction(Polymer_Isotope, Backbone), y= Order_rel_abund_Sample, 
                                    fill=forcats::fct_reorder(Order, Order_rel_abund_Sample)))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(colors_by_Maaike)) + 
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
  facet_nested( . ~ Location + Treatment, drop = T, 
                axes = 'margins',
                nest_line = element_line(),
                strip = strip_nested(
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                )) + 
  
  theme_classic2() + 
  theme(
    axis.text.x=element_text( size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text( size= 12), 
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
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Material") + labs(fill = "Order")

gg_Order


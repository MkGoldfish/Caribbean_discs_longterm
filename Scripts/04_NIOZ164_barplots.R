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
pal.phylum.i <- c('#44AA99', '#999933', '#332288', '#88CCEE', '#CC6677',
                  '#DDCC77', '#117733', '#882255', '#EE8866', '#663333', 
                  '#AA4499', '#0077BB', '#EE3377', '#E6AB02FF', '#CC3311', 
                  '#CAB2D6FF', "#E5C494",'#B2DF8AFF','#225555','#666633', '#555555')

pal.phylum.w <- c('#44AA99', '#332288', '#CC6677', '#CC3311', "#E5C494",
                  '#88CCEE', '#DDCC77', '#AA4499', '#117733', '#0077BB',
                  '#663333', '#E6AB02FF', '#999933', '#555555')

pal.order.i <- c("#cb7c77", "#68d359", "#6b42c8", "#c9d73d", "#c555cb", 
                 "#aed688", "#502e71", "#c49a3f", "#6a7dc9", "#d7652d", 
                 "#7cd5c8", "#c5383c", "#507d41", "#cf4c8b", "#5d8d9c", 
                 "#722e41", "#c8b693", "#33333c", "#674c2a", "#F1BB7B",
                 "#5B1A18", "#AF4E24", "#E6A0C4", "#C6CDF7", "#D8A499",
                 "#149911", "#046C9A", "#D69C4E", "#175149", "#CCC591", 
                 "#D3DDDC", "#FF0000", "#0B775E", "#FAD510", "#F98400", 
                 "#0047E0", "#354823", "#F3DF6C", "#5F5647", "#9A8822",
                 "#4E2A1E", "#9B110E", "#0C1707", "#FF4B3E","#79402E", 
                 "#74A089", "#8D8680")



BottleRocket2 = c("#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E")
Rushmore1 = c("#E1BD6D", "#EABE94", "#0B775E", "#35274A" ,"#F2300F")
Royal2 = c("#9A8822", "#F5CDB4", "#F8AFA8", "#FDDDA0", "#74A089")
Darjeeling1 = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")
Chevalier1 = c("#446455", "#FDD262", "#D3DDDC", "#C7B19C")
Darjeeling2 = c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000")
BottleRocket1 = c("#A42820", "#5F5647", "#9B110E", "#3F5151", "#4E2A1E", "#550307", "#0C1707")
Moonrise1 = c("#F3DF6C", "#CEAB07", "#D5D5D3", "#24281A")
Moonrise2 = c("#798E87", "#C27D38", "#CCC591", "#29211F")
Moonrise3 = c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B")
Cavalcanti1 = c("#D8B70A", "#02401B", "#A2A475", "#81A88D", "#972D15")
GrandBudapest1 = c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
GrandBudapest2 = c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
IsleofDogs1 = c("#9986A5", "#79402E", "#CCBA72", "#0F0D0E", "#D9D0D3", "#8D8680")
IsleofDogs2 = c("#EAD3BF", "#AA9486", "#B6854D", "#39312F", "#1C1718")
FrenchDispatch = c("#90D4CC", "#BD3027", "#B0AFA2", "#7FC0C6", "#9D9C85")
AsteroidCity1 = c("#0A9F9D", "#CEB175", "#E54E21", "#6C8645", "#C18748")
AsteroidCity2 = c("#C52E19", "#AC9765", "#54D8B1", "#b67c3b", "#175149", "#AF4E24")
AsteroidCity3 = c("#FBA72A", "#D3D4D8", "#CB7A5C", "#5785C1")

pal.carto <- cartography::carto.pal(pal1 = "multi.pal", n1 = 20)
pal.carto <- c("#cb7c77", "#68d359", "#6b42c8", "#c9d73d", "#c555cb", "#aed688", "#502e71", 
  "#c49a3f", "#6a7dc9", "#d7652d", "#7cd5c8", "#c5383c", "#507d41", "#cf4c8b", 
  "#5d8d9c", "#722e41", "#c8b693", "#33333c", "#c6a5cc", "#674c2a", '#DDDDDD')

pal.igv <- pal_igv("default")(55)

[1] "#5050FFFF" "#CE3D32FF" "#749B58FF" "#F0E685FF" "#466983FF" "#BA6338FF" "#5DB1DDFF" "#802268FF" "#6BD76BFF" 
"#D595A7FF" "#924822FF" "#837B8DFF" "#C75127FF" "#D58F5CFF" "#7A65A5FF" "#E4AF69FF" "#3B1B53FF"
[18] "#CDDEB7FF" "#612A79FF" "#AE1F63FF" "#E7C76FFF" "#5A655EFF" "#CC9900FF" "#99CC00FF" "#A9A9A9FF" 
"#CC9900FF" "#99CC00FF" "#33CC00FF" "#00CC33FF" "#00CC99FF" "#0099CCFF" "#0A47FFFF" "#4775FFFF" "#FFC20AFF"
[35] "#FFD147FF" "#990033FF" "#991A00FF" "#996600FF" "#809900FF" "#339900FF"

## Import data ----------------------------------------------------------------------------
tt <- read.csv('../Processed-Data/NIOZ164_EUX_discs_RA_tidy_data_decontamed_tax.correct_pruned.csv', row.names = NULL)
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
                                                               Polymer, Backbone, Treatment, Phylum, Phylum_rel_abund_Sample) %>% distinct()

Phyla_filtAb <- Phyla %>% mutate(Phylum = ifelse(Phylum_rel_abund_Sample<0.02, "others<2.0%", Phylum)) 

length(unique(Phyla_filtAb$Phylum))

gg_Phylum<- ggplot(Phyla_filtAb, aes(x=interaction(Polymer_Isotope, Backbone), y= Phylum_rel_abund_Sample, 
                               fill=forcats::fct_reorder(Phylum, Phylum_rel_abund_Sample )))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(pal.phylum.i)) + 
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
    ggh4x.axis.nesttext.x = element_blank(),
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
  scale_fill_manual(values = rev(pal.phylum.w)) + 
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
    ggh4x.axis.nesttext.x = element_blank(),
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

Order_filtAb <- Orders %>% mutate(Order = ifelse(Order_rel_abund_Sample<0.04, "others<4%", Order))
length(unique(Order_filtAb$Order))

gg_Order<- ggplot(Order_filtAb, aes(x=interaction(Polymer_Isotope, Backbone), y= Order_rel_abund_Sample, 
                                     fill=forcats::fct_reorder(Order, Order_rel_abund_Sample )))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(pal.order.i)) + 
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
    ggh4x.axis.nesttext.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Polymer") + labs(fill = "Order")

gg_Order

### Wild plastic-----------------------------------------------------------------
Orders <- tt.wild %>% filter(Kingdom != "Eukaryota") %>%  select(Description, Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, 
                                                                Treatment, Order, Order_rel_abund_Sample) %>% 
  distinct()

Order_filtAb <- Orders %>% mutate(Order = ifelse(Order_rel_abund_Sample<0.04, "others<4%", Order))

gg_Order<- ggplot(Order_filtAb, aes(x=Description, y= Order_rel_abund_Sample, 
                                      fill=forcats::fct_reorder(Order, Order_rel_abund_Sample)))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(pal.order.i)) + 
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
    ggh4x.axis.nesttext.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Material") + labs(fill = "Order")

gg_Order



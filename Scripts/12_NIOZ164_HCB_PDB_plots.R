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


## Import sequencing data ----------------------------------------------------------------------------
tt <- read.csv('../Processed-Data/NIOZ164_EUX_discs_tidy_data_decontamed_pruned_RA.csv', na.strings = c(""))
head(tt)
tt$X <- NULL

# Remove the _ from the location names and treatments for plotting
tt.1 <- tt %>%  mutate(Location = ifelse(Location == "Crooks_Castle", "Crooks Castle",
                                   ifelse(Location == "Charles_Brown", "Charles Brown", Location)))
tt.1 <- tt.1 %>%  mutate(Treatment = ifelse(Treatment == "no_UV", "no UV", Treatment))
head(tt.1)

tt.inc <- tt.1 %>% filter(Phase == "Disc" ) 
unique(tt.inc$Location)
unique(tt.inc$Treatment)
unique(tt.inc$Polymer)

tt.wild <- tt.1 %>% filter(Method == "Collection") 
unique(tt.wild$Location)
unique(tt.wild$Method)
unique(tt.wild$Polymer)

tt.filt <- tt.1 %>% filter(Phase == "Filter" ) %>% filter(!Treatment == "NA ")
unique(tt.inc$Location)
unique(tt.inc$Treatment)
unique(tt.inc$Polymer)

## PDB data ---------------------------------------------------------------------
PDB_all <- read.delim("../Data/PlasticDB_all_data_20230904.txt", sep = '\t', dec = ".")

#Select only plastic and taxonomy column, and separate the taxonomy
PDB.sep <- PDB_all %>% select(Plastic, lineage) %>% 
  separate(lineage, into = c( "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ",") %>% 
  separate(col = Kingdom, into = c("A", "Kingdom"), sep = ":") %>% 
  separate(col = Phylum, into = c("B", "Phylum"), sep = ":") %>% 
  separate(col = Class, into = c("C", "Class"), sep = ":") %>% 
  separate(col =  Order, into = c("D", "Order"), sep = ":") %>% 
  separate(col = Family, into = c("E", "Family"), sep = ":") %>% 
  separate(col = Genus , into = c("F", "Genus"), sep = ":") %>% 
  separate(col =Species, into = c("G", "Species"), sep = ":") 

# Remove eukaryotes and unculterd bacteria  
PDB_clean <- PDB.sep %>%  select(Plastic, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  filter(Kingdom == "Bacteria") %>% filter(Phylum != "uncultured bacterium")

# Check! All good
unique(PDB_clean$Kingdom)
unique(PDB_clean$Phylum)
unique(PDB_clean$Class)

# What plastics can be found in PlasticDB
Plastics <- PDB$Plastic %>% unique() %>% sort()

# What genera can be found in PlasticDB
pdb.gens <- PDB$Genus %>% unique()
length(pdb.gens)

## HCB data ---------------------------------------------------------------------
HCB <- as.character(read_lines("../data/Hydrocarbon_degraders_sorted_22_08.txt"))
HCB_only <- setdiff(HCB, pdb.gens)
PDB_only <- setdiff(pdb.gens, HCB)
plast.hcb.intersect <- intersect(pdb.gens, HCB)
plast.hcb <- c(HCB_only, pdb.gens)

Genus_plastic <- Genus %>%  filter(Genus%in%PDB$Genus)
Genus_oil <-  Genus %>%  filter(Genus%in%HCB)
Genus_plastoil <- Genus %>%  filter(Genus%in%plast.HCB) 
Genus_HCBonly <-  Genus %>%  filter(Genus%in%HCB_only)


## Colors for plotting --------------------------------------------------------------------
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

# Genera sequencing data incubations --------------------------------------------------------------------
Genus <- tt.inc  %>%  select(Description, Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, Treatment, 
                             Phylum, Genus, Genus_rel_abund_Sample)%>% 
  distinct() 

Genus.pel <- Genus %>% filter(Habitat == "Pelagic")
Genus.bent <- Genus %>% filter(Habitat == "Benthic")

#Count how many unique genera we have. 
Genus.pel %>% select(Genus) %>%  unique() %>% count()
Genus.bent %>% select(Genus) %>%  unique() %>% count()

# Create a unique vector of all genera we found
found.genera <- Genus$Genus %>% unique() %>%  sort()
str(found.genera)

# Whih of our found genera are in PlasticDB?
gen.in.pdb <- PDB_clean %>% filter(Genus %in% found.genera)
length(unique(gen.in.pdb$Genus))

#What non-PE Plastics that have the PE string do we find?
PDB_clean %>% filter(str_detect(Plastic,"PE")) %>% select(Plastic) %>% unique()

# Group plastics per category
pdb.gen.filt <-gen.in.pdb %>% mutate(Plastic_group = if_else(str_detect(Plastic, "PET"), "PET",
                                                       (if_else(str_detect(Plastic, "PS"), "PS", 
                                                                (if_else(str_detect(Plastic, "PP"), "PP", 
                                                                         (if_else(Plastic %in% c("PEA", "PEC", "PEF", "PEG", "PES" ), "Other", 
                                                                         (if_else(str_detect(Plastic, "PE"), "PE", 
                                                                                  (if_else(str_detect(Plastic, "Nylon"), "Nylon", 
                                                                                           "Other"))))))))))))

# Select the genera from these plastic groups
PE <- pdb.gen.filt %>% filter(Plastic_group == "PE") %>% pull(Genus) %>% unique()
PP <- pdb.gen.filt %>% filter(Plastic_group == "PP") %>% pull(Genus) %>% unique()
PS <- pdb.gen.filt %>% filter(Plastic_group== "PS") %>% pull(Genus) %>% unique()
PET <- pdb.gen.filt %>% filter(Plastic_group == "PET") %>% pull(Genus) %>% unique()
Nylon <- pdb.gen.filt %>% filter(Plastic_group == "Nylon") %>% pull(Genus) %>% unique()

## PA/map of genera in PlasticDB ------------------------------------------------



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
  filter(Genus_rel_abund_Sample > 0.01) %>% distinct()
top_genus_pel %>% select(Genus) %>% unique() 

top_genus_pel$Genus <- factor(top_genus_pel$Genus, levels=rev(sort(unique(top_genus_pel$Genus))))
top_genus_pel <- top_genus_pel %>% arrange(Order)

Genus_bubble_pel <- ggplot(top_genus_pel,aes(x=interaction(Polymer_Isotope,Backbone),y= Genus)) +
  geom_point(aes(size=Genus_rel_abund_Sample, fill = factor(Polymer_Isotope)), shape = "circle filled", stroke = NA, alpha = 0.7) +
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

Genus_bubble_pel

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

Genus_bubble_bent <- ggplot(top_genus_bent,aes(x=interaction(Polymer_Isotope,Backbone),y= Genus)) +
  geom_point(aes(size=Genus_rel_abund_Sample, fill = factor(Polymer_Isotope)), shape = "circle filled", stroke = NA, alpha = 0.7) +
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

Genus_bubble_bent

#### Combine incubations with cowplot -----------------------------------------------------------------
legend.a <- get_legend(Genus_bubble_bent+
                         theme(legend.direction = "horizontal",
                               legend.title.align = 0.5))

plot_grid(Genus_bubble_pel + theme(legend.position ="none", axis.text.x = element_blank(), 
                                   axis.title.x = element_blank(),ggh4x.axis.nestline.x = element_blank(),
                                   ggh4x.axis.nesttext.x =element_blank(), plot.margin = unit(c(0,0,-1,0), "cm")),
          Genus_bubble_bent + theme(legend.position ="none", plot.margin = unit(c(0,0,0,0), "cm")),
          legend.a,
          ncol = 1,
          nrow = 3,
          align = 'v',
          axis = "ltbr",
          rel_heights = c(1,1.25,0.1))


### Wild plastic-----------------------------------------------------------------
Genus.wild <- tt.wild  %>%  select(Location, Habitat, Description, Phylum, Order, Genus, Genus_rel_abund_Sample)%>% 
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
top_genus_wild %>% select(Genus) %>% unique() 

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


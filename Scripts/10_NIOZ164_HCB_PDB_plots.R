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
library (emojifont)
library(ggthemes)
library(stringr)
library(tidyverse)

## Import sequencing data ----------------------------------------------------------------------------
tt.1 <- read.csv('../Processed-Data/NIOZ164_EUX_discs_RA_tidy_data_decontamed_tax.correct_pruned.csv', na.strings = c(""))
head(tt.1)
tt.1$X <- NULL

unique(tt.1$Location)
unique(tt.1$Method)
unique(tt.1$Polymer)
unique(tt.1$Phase)

tt.inc <- tt.1 %>% filter(Phase == "Disc" ) 
unique(tt.inc$Location)
unique(tt.inc$Habitat)
unique(tt.inc$Treatment)
unique(tt.inc$Polymer)

tt.wild <- tt.1 %>% filter(Method == "Collection") 
unique(tt.wild$Location)
unique(tt.wild$Habitat)
unique(tt.wild$Method)
unique(tt.wild$Polymer)

tt.i.w. <- tt.1 # %>% filter(Phase != "Filter" ) 
# unique(tt.i.w.$Location)
# unique(tt.i.w.$Treatment)
# unique(tt.i.w.$Polymer)

# tt.filt <- tt.1 %>% filter(Phase == "Filter" ) %>% filter(!Treatment == "NA ")
# unique(tt.inc$Location)
# unique(tt.inc$Treatment)
# unique(tt.inc$Polymer)

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
Plastics <- PDB_all$Plastic %>% unique() %>% sort()

# What genera can be found in PlasticDB
pdb.gens <- PDB_clean$Genus %>% unique()
length(pdb.gens)

## HCB data ---------------------------------------------------------------------
HCB <- as.character(read_lines("../data/Hydrocarbon_degraders_sorted_22_08.txt"))
HCB_only <- setdiff(HCB, pdb.gens)
PDB_only <- setdiff(pdb.gens, HCB)
plast.hcb.intersect <- intersect(pdb.gens, HCB)
plast.hcb <- c(HCB_only, pdb.gens)

## Colors for plotting --------------------------------------------------------------------
pal.loc <- c("#6A3D9A" , "#33A02C", "#51C3CCFF")
# CB, CC, Zeelandia
pal.habs<- c("#FF8E32FF", "#009292FF")
# Benthic, Pelagic
pal.loc.hab <- c("#6A3D9A", "#CAB2D6", "#33A02C","#B2DF8A", "#51C3CCFF")
# CB_P, CB_B, CC_P, CB_B, Zeelandia
pal.pols.isotop <- c("#E31A1C","#E7298A", "#1F78B4","#A6CEE3","#E6AB02", "#A6761D", "#E5C494","#1B9E77")
# PE;PE-13C;PP;PP-13C;PS;PET;Nylon;Blanco
pal.uv <- c("#DDCC77","#332288") 
# UV, noUV

colors_M1 <- c("#004e64", "#ecc8af", "#F2AF29", "#436436", "#00a5cf", 
               "#c18c5d", "#5f0f40", "#DC602E", "#495867", "#A29F15", 
               "#570000", "#FFF5B2", "#20221B", "#9fffcb", "#c08497", 
               "#8D6346", "#FF4B3E", "#149911", "#472d30")

showtext::showtext_opts(dpi=500)

# PA/map of genera in PlasticDB ------------------------------------------------
## Combined data of wild and incubations ---------------------------------------
Genus <- tt.i.w.  %>%  select(Description, Location, Habitat, Polymer, Polymer_Isotope, Backbone, Treatment, 
                             Phylum, Genus, Genus_rel_abund_Sample, Sample_rel_abund) %>% 
  distinct() 

# Create a unique vector of all genera we found
found.genera <- Genus$Genus %>% unique() %>%  sort()
str(found.genera)

# Which of our found genera are in PlasticDB?
gen.in.pdb <- PDB_clean %>% filter(Genus %in% found.genera)
length(unique(gen.in.pdb$Genus))

hcb.gen <- HCB_only[HCB_only %in% found.genera] 
length(unique(hcb.gen))

# Which of our found genera are in PlasticDB and HCB w RA> 0.5%
Genus.filt <- Genus  %>%  filter(Genus_rel_abund_Sample > 0.005)%>% 
  distinct() 

found.gen.filt <- Genus.filt$Genus %>% unique() %>%  sort()

gen.in.pdb <- PDB_clean %>% filter(Genus %in% found.gen.filt) 
length(unique(gen.in.pdb$Genus))

hcb.gen <- HCB_only[HCB_only %in% found.gen.filt] 
length(unique(hcb.gen))


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

# Filter for HCB and PDB in incubations, and for minimum RA
# PDB/HCB in our dataset
Genus_plastic <- Genus %>%  filter(Genus%in%pdb.gens) %>% filter(Genus_rel_abund_Sample > 0)
unique(Genus_plastic$Polymer)

## add symbols matching PDB polymers
PDB.Gens.pols <- Genus_plastic  %>% mutate(Genus.pol = 
                                                 if_else(Genus %in% PE, paste(Genus, "\u{25CF}", sep = " "), Genus)) 
PDB.Gens.pols <- PDB.Gens.pols %>% mutate(Genus.pol =                                                 
                                                 if_else(Genus %in% PP, paste(Genus.pol, "\u{25B2}", sep = " "), Genus.pol))
PDB.Gens.pols <- PDB.Gens.pols %>% mutate(Genus.pol =                                                 
                                                 if_else(Genus %in% PS, paste(Genus.pol, "\u{2217}", sep = " "), Genus.pol))
PDB.Gens.pols <- PDB.Gens.pols %>% mutate(Genus.pol =                                                 
                                                 if_else(Genus %in% PET, paste(Genus.pol, "\u{25A0}", sep = " "), Genus.pol))
PDB.Gens.pols <- PDB.Gens.pols %>% mutate(Genus.pol =                                                 
                                                 if_else(Genus %in% Nylon, paste(Genus.pol, "\u{25C6}", sep = " "), Genus.pol))

# Put genera in alphabetcial order
PDB.Gens.pols$Genus.pol <- factor(PDB.Gens.pols$Genus.pol, levels=rev(sort(unique(PDB.Gens.pols$Genus.pol))))

# set colors and shape
shapes = c(19,17,8,15,18,4,6)

#### Plot genera from our dataset present in plasticDB, and the polymers upon which they are found in plasticDB ####
plot.pdb <- ggplot(PDB.Gens.pols, aes(x = Polymer, y = Genus.pol)) +
  geom_point(aes(shape = Polymer), size = 3) +
  scale_shape_manual(values = shapes, limits = c("PE", "PP", "PS", "PET", "Nyl", "Unknown", "B"), name = "Polymer") +
  scale_x_discrete(limits = c("PE", "PP", "PS", "PET", "Nyl", "Unknown", "B")) +
  theme_pubclean() +
  theme(
    axis.text.x=element_text(size = 11, angle = 60, hjust = 0.9), 
    axis.text.y=element_text(size= 11), 
    legend.text=element_text(size = 11),
    legend.title = element_text(size=12),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    plot.title = element_text(size = 14),
    plot.subtitle = element_text(),
    panel.border = element_rect(color = "grey", fill = NA),
    legend.key = element_rect(fill = NA), 
    legend.position = "right") +
  xlab("Plastics in PlasticDB")+ 
  ylab("Genera in dataset found in PlasticDB") + 
  labs(title = "")

plot.pdb


# Bubbleplots -----------------------------------------------------------------
## Genera sequencing data incubations --------------------------------------------------------------------
Genus <- tt.inc  %>%  select(Description, Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, Treatment, 
                             Phylum, Genus, Genus_rel_abund_Sample, Sample_rel_abund) %>% 
                             filter(Genus_rel_abund_Sample > 0.005) %>% 
  distinct() 

# Create a unique vector of all genera we found
found.genera.inc <- Genus$Genus %>% unique() %>%  sort()
str(found.genera)

# Whih of our found genera are in PlasticDB?
gen.in.pdb.inc <- PDB_clean %>% filter(Genus %in% found.genera.inc)
length(unique(gen.in.pdb.inc$Genus))

# Which of our found genera are HCB?
gen.hcb.only.inc<-HCB_only[HCB_only %in% found.genera.inc]
unique(gen.hcb.only.inc) #12

gen.in.hcb.inc<-HCB[HCB%in% found.genera.inc]
unique(gen.in.hcb.inc) #23


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

# Filter for HCB and PDB in incubations, and for minimum RA
# PDB/HCB in our dataset
Genus_plastic <- Genus %>%  filter(Genus%in%PDB_clean$Genus)
Genus_oil <-  Genus %>%  filter(Genus%in%HCB)
Genus_plastoil <- Genus %>%  filter(Genus%in%plast.hcb) 
Genus_HCBonly <-  Genus %>%  filter(Genus%in%HCB_only)
Genus_PDBonly <-  Genus %>%  filter(Genus%in%PDB_only)


HCB.PDB.Genera <- Genus_plastoil %>% dplyr::select(Description, Genus, Genus_rel_abund_Sample, Sample_rel_abund,
                                                   Location, Habitat, Polymer_Isotope, Backbone, Treatment) %>%
  filter(Genus_rel_abund_Sample > 0.005) %>% distinct()

unique(HCB.PDB.Genera$Genus) %>% length()

# Add degradation category including intersect for facetting
HCB.PDB.Genera.subs <- HCB.PDB.Genera  %>% mutate(Degrading = case_when(
  Genus %in% HCB_only ~ "HCB",
  Genus %in% PDB_only ~ "PlasticDB genus",
  Genus %in% plast.hcb.intersect ~ "HCB in PlasticDB",
) )

# Add symbols to Genus names based on the Plastic they are detected on in PDB
## add symbols matching PDB polymers
HCB.PDB.Genera.pdbpols <- HCB.PDB.Genera.subs  %>% mutate(Genus.pol = 
                                                            if_else(Genus %in% PE, paste(Genus, "\u{25CF}", sep = " "), Genus)) 
HCB.PDB.Genera.pdbpols <- HCB.PDB.Genera.pdbpols %>% mutate(Genus.pol =                                                 
                                                              if_else(Genus %in% PP, paste(Genus.pol, "\u{25B2}", sep = " "), Genus.pol))
HCB.PDB.Genera.pdbpols <- HCB.PDB.Genera.pdbpols %>% mutate(Genus.pol =                                                 
                                                              if_else(Genus %in% PS, paste(Genus.pol, "*", sep = " "), Genus.pol))
HCB.PDB.Genera.pdbpols <- HCB.PDB.Genera.pdbpols %>% mutate(Genus.pol =                                                 
                                                              if_else(Genus %in% PET, paste(Genus.pol, "\u{25A0}", sep = " "), Genus.pol))
HCB.PDB.Genera.pdbpols <- HCB.PDB.Genera.pdbpols %>% mutate(Genus.pol =                                                 
                                                              if_else(Genus %in% Nylon, paste(Genus.pol, "\u{25C6}", sep = " "), Genus.pol))

### Calculate abundanof HCB and PDB ------------------------
#### Pelagic HCB reads
Samp.HCB.pel.avg <- Genus_HCBonly %>% filter(Habitat == "Pelagic") %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean = mean(sum)) *100

Samp.HCB.pel.sd <- Genus_HCBonly %>% filter(Habitat == "Pelagic")  %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Samp.HCB.pel.avg
Samp.HCB.pel.sd

#### Average reads per Genus to determine dominant Genera
Genus.HCB.pel.avg <- Genus_HCBonly %>% filter(Habitat == "Pelagic") %>% select(Description, Genus, Genus_rel_abund_Sample) %>% group_by(Genus) %>%
  mutate(avg = 100*mean(Genus_rel_abund_Sample)) %>% 
  mutate(sd = sd(Genus_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

Genus.HCB.pel.avg 


#### Benthic HCB reads
Samp.HCB.ben.avg <- Genus_HCBonly %>% filter(Habitat == "Benthic") %>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean = mean(sum)) *100

Samp.HCB.ben.sd <- Genus_HCBonly %>% filter(Habitat == "Benthic")  %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Samp.HCB.ben.avg
Samp.HCB.ben.sd

#### Average reads per Genus to determine dominant Genera
Genus.HCB.ben.avg <- Genus_HCBonly %>% filter(Habitat == "Benthic") %>% select(Description, Genus, Genus_rel_abund_Sample) %>% group_by(Genus) %>%
  mutate(avg = 100*mean(Genus_rel_abund_Sample)) %>% 
  mutate(sd = sd(Genus_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

Genus.HCB.ben.avg 


#### Pelagic PDB reads
Gen.plast.pel.avg <- Genus_plastic %>% filter(Habitat == "Pelagic")%>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

Gen.plast.pel.sd <- Genus_plastic %>% filter(Habitat == "Pelagic") %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Gen.plast.pel.avg
Gen.plast.pel.sd

#### Average reads per Genus to determine dominant Genera
Genus.PDB.pel.avg <- Genus_plastic %>% filter(Habitat == "Pelagic") %>% select(Description, Genus, Genus_rel_abund_Sample) %>% group_by(Genus) %>%
  mutate(avg = 100*mean(Genus_rel_abund_Sample)) %>% 
  mutate(sd = sd(Genus_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

Genus.PDB.pel.avg 


#### Benthic PDB reads
Gen.plast.ben.avg <- Genus_plastic %>% filter(Habitat == "Benthic")%>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

Gen.plast.ben.sd <- Genus_plastic %>% filter(Habitat == "Benthic") %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Gen.plast.ben.avg
Gen.plast.ben.sd

#### Average reads per Genus to determine dominant Genera
Genus.PDB.ben.avg <- Genus_plastic %>% filter(Habitat == "Benthic") %>% select(Description, Genus, Genus_rel_abund_Sample) %>% group_by(Genus) %>%
  mutate(avg = 100*mean(Genus_rel_abund_Sample)) %>% 
  mutate(sd = sd(Genus_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

Genus.PDB.ben.avg 


## Bubbleplot incubations  -------------------------------------------------------------------                                                                
HCB.PDB.Genera.pdbpols$Genus.pol <- factor(HCB.PDB.Genera.pdbpols$Genus.pol, levels=rev(sort(unique(HCB.PDB.Genera.pdbpols$Genus.pol))))

HCB.PDB.Bubble <- ggplot(HCB.PDB.Genera.pdbpols,aes(x=interaction(Polymer_Isotope,Backbone),y= Genus.pol)) +
  geom_point(aes(size=Genus_rel_abund_Sample, fill = factor(Polymer_Isotope)), shape = "circle filled", stroke = NA) +
  scale_fill_manual(values = pal.pols.isotop) +
  scale_size(range = c(1.5,6))+
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 6, color = NA))) +
  ylab("") +
  xlab("") +  
  facet_nested(Habitat + Degrading  ~ Location + Treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins', as.table = F, 
               nest_line = element_line(),
               strip = strip_nested( size = "variable",
                 background_y =  elem_list_rect(color = c("#009292FF", "#FF8E32FF",
                                                          "grey10", "grey50", "grey30",
                                                          "grey10", "grey50", "grey30" ),
                                                fill = "white", linewidth = 1, by_layer_y = F),
                 text_y = elem_list_text(size = c(10,10,rep_len(9,7)), angle = c(270,270,rep_len(0,7)), 
                                         color = c("#009292FF", "#FF8E32FF",
                                                          "grey10", "grey50", "grey30",
                                                          "grey10", "grey50", "grey30" ), by_layer_y = F),
                 background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 9, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 9, color = "black"), 
    legend.text=element_text(size = 9),
    legend.title = element_text(size=10),
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10),
    strip.text.x = element_text(size = 9),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1.5, 'pt'),
    legend.position = "bottom") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

HCB.PDB.Bubble


ggsave("NIOZ164_PDB_HCB_incubations_bubble.tiff", 
       width = 19.5, height  = 26, unit = "cm", 
       dpi = 500, bg='white')

## Genera sequencing Wild plastic --------------------------------------------------------------------
Genus <- tt.wild %>%  select(Description, Location, Habitat, Polymer, Isotope, Polymer_Isotope, Backbone, Treatment, 
                             Phylum, Genus, Genus_rel_abund_Sample, Sample_rel_abund) %>%  
  filter(Genus_rel_abund_Sample > 0.005) %>% 
  distinct() 

# Create a unique vector of all genera we found
found.genera.w <- Genus$Genus %>% unique() %>%  sort()
str(found.genera.w)

# Which of our found genera are in PlasticDB?
gen.in.pdb.w <- PDB_clean %>% filter(Genus %in% found.genera.w)
length(unique(gen.in.pdb.w$Genus))

intersect(gen.in.pdb.w$Genus, gen.in.pdb.inc$Genus)

# Which of our found genera are HCB?
gen.hcb.only.w<-HCB_only[HCB_only %in% found.genera.w]
unique(gen.hcb.only.w) #12

gen.in.hcb.w<-HCB[HCB%in% found.genera.w]
unique(gen.in.hcb.w) #21

intersect(gen.in.hcb.w, gen.in.hcb.inc) #15
setdiff(gen.in.hcb.w, gen.in.hcb.inc) #6
setdiff( gen.in.hcb.inc, gen.in.hcb.w) #8


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

# Filter for HCB and PDB in incubations, and for minimum RA
# PDB/HCB in our dataset
Genus_plastic <- Genus %>%  filter(Genus%in%PDB_clean$Genus)
Genus_oil <-  Genus %>%  filter(Genus%in%HCB)
Genus_plastoil <- Genus %>%  filter(Genus%in%plast.hcb) 
Genus_HCBonly <-  Genus %>%  filter(Genus%in%HCB_only)
Genus_PDBonly <-  Genus %>%  filter(Genus%in%PDB_only)


HCB.PDB.Genera <- Genus_plastoil  %>% dplyr::select(Description, Genus, Genus_rel_abund_Sample, 
                                                   Location, Habitat, Polymer_Isotope, Backbone, Treatment) %>%
  filter(Genus_rel_abund_Sample > 0.005) %>% distinct()

unique(HCB.PDB.Genera$Genus) %>% length()

### Calculate abundance of HCB and PDB ------------------------
#### Beach HCB reads
Samp.HCB.beach.avg <- Genus_HCBonly  %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean = mean(sum)) *100

Samp.HCB.beach.sd <- Genus_HCBonly  %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Samp.HCB.beach.avg
Samp.HCB.beach.sd

#### Average reads per Genus to determine dominant Genera
Genus.HCB.beach.avg <- Genus_HCBonly %>%  select(Description, Genus, Genus_rel_abund_Sample) %>% group_by(Genus) %>%
  mutate(avg = 100*mean(Genus_rel_abund_Sample)) %>% 
  mutate(sd = sd(Genus_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

Genus.HCB.beach.avg 

#### Beach PDB reads
Gen.plast.beach.avg <- Genus_plastic %>%  select(Description,Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

Gen.plast.beach.sd <- Genus_plastic %>%  select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Gen.plast.beach.avg
Gen.plast.beach.sd

#### Average reads per Genus to determine dominant Genera
Genus.PDB.beach.avg <- Genus_plastic  %>% select(Description, Genus, Genus_rel_abund_Sample) %>% group_by(Genus) %>%
  mutate(avg = 100*mean(Genus_rel_abund_Sample)) %>% 
  mutate(sd = sd(Genus_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

Genus.PDB.beach.avg 


## Bubbleplot wild plastics  ------------------------------------------------------------------- 
# Add degradation category including intersect for facetting
HCB.PDB.Genera.subs <- HCB.PDB.Genera  %>% mutate(Degrading = case_when(
  Genus %in% HCB_only ~ "HCB",
  Genus %in% PDB_only ~ "PlasticDB genus",
  Genus %in% plast.hcb.intersect ~ "HCB in PlasticDB",
))

# Add symbols to Genus names based on the Plastic they are detected on in PDB
## add symbols matching PDB polymers
HCB.PDB.Genera.pdbpols <- HCB.PDB.Genera.subs  %>% mutate(Genus.pol = 
                                                            if_else(Genus %in% PE, paste(Genus, "\u{25CF}", sep = " "), Genus)) 
HCB.PDB.Genera.pdbpols <- HCB.PDB.Genera.pdbpols %>% mutate(Genus.pol =                                                 
                                                              if_else(Genus %in% PP, paste(Genus.pol, "\u{25B2}", sep = " "), Genus.pol))
HCB.PDB.Genera.pdbpols <- HCB.PDB.Genera.pdbpols %>% mutate(Genus.pol =                                                 
                                                              if_else(Genus %in% PS, paste(Genus.pol, "*", sep = " "), Genus.pol))
HCB.PDB.Genera.pdbpols <- HCB.PDB.Genera.pdbpols %>% mutate(Genus.pol =                                                 
                                                              if_else(Genus %in% PET, paste(Genus.pol, "\u{25A0}", sep = " "), Genus.pol))
HCB.PDB.Genera.pdbpols <- HCB.PDB.Genera.pdbpols %>% mutate(Genus.pol =                                                 
                                                              if_else(Genus %in% Nylon, paste(Genus.pol, "\u{25C6}", sep = " "), Genus.pol))


HCB.PDB.Genera.pdbpols$Genus.pol <- factor(HCB.PDB.Genera.pdbpols$Genus.pol, levels = rev(sort(unique(HCB.PDB.Genera.pdbpols$Genus.pol))))

HCB.PDB.Bubble <- ggplot(HCB.PDB.Genera.pdbpols,aes(x=Description,y= Genus.pol)) +
  geom_point(aes(size=Genus_rel_abund_Sample, fill = factor(Degrading)), shape = "circle filled", stroke = NA, alpha = 0.7) +
  scale_fill_manual(values = c("#149911","#A29F15","#436436")) +
  scale_size(range = c(1.5,6))+
  guides( x = "axis_nested",  fill = FALSE) +
  ylab("") +
  xlab("") +  
  facet_nested(Degrading  ~ Location + Treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins', as.table = F, 
               nest_line = element_line(),
               strip = strip_nested( size = "variable",
                                     background_y =  elem_list_rect(color = c("grey10", "grey50", "grey30"),
                                                                    fill = "white", linewidth = 1, by_layer_y = F),
                                     text_y = elem_list_text(size = rep_len(8,7), angle = rep_len(270,7), 
                                                             color = c("grey10", "grey50", "grey30"), by_layer_y = F),
                                     background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 9, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 9, color = "black"), 
    legend.text=element_text(size = 9),
    legend.title = element_text(size=10),
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9),
    strip.text.x = element_text(size = 8),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.border = element_rect(color = "grey80", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1.5, 'pt'),
    legend.position = "right") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

HCB.PDB.Bubble

showtext::showtext_opts(dpi=500)

ggsave("NIOZ164_PDB_HCB_wild_bubble.eps", 
       width = 14, height  = 15, unit = "cm", 
       dpi = 500, bg='white')
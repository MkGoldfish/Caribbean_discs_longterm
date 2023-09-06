##%#####################################################################################%##
#NIOZ164 Statia Discs - Ordinations                                                  #####                                            
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

# Load libraries -------------------------------------------------------------------------
library(devtools)
library(phyloseq)
library(tidyverse)
library(vegan)
library(ggpubr)
library(ggh4x)
library(compositions)
library(ade4)
library(FactoMineR)
library(labdsv)
library(zCompositions)
library(microbiome)
library(ggpubr)

# Colors for plotting --------------------------------------------------------------------
pal_isme <- c("#006d77", "#ffddd2", "#00C49A", "#e29578", "#83c5be")

pal.loc <- c("#FF6DB6FF" , "#004949FF",  "#66A61E")
# CB, CC, Zeelandia
pal.habs<- c("#66A61E","#FF8E32FF", "#51C3CCFF")
# Benthic, Pelagic
pal.loc.hab <- c("#FF6DB6FF","#FFB6DBFF", "#004949FF", "#009292FF", "#66A61E")
# CB_P, CB_B, CC_P, CB_B, Zeelandia

pal.pols.isotop <- c("#E31A1C", "#7570B3", "#1F78B4","#A6CEE3", "#E6AB02","#A6761D", "#E5C494","#1B9E77","#882255")
pal.pols <- c("#1F78B4","#A6761D", "#7570B3" , "#E31A1C", "#E6AB02", "#1B9E77","#882255")
# PE;PE-13C;PP;PP-13C;PS;PET;Nylon;Blanco
pal.uv <- c("#332288","#DDCC77", "grey") 
# UV, noUV

 
# Import data ----------------------------------------------------------------------------
source("basic_info_physeq_object.R")
# Import earlier created physeq object
physeq_object <- readRDS("../Analysis/NIOZ164_physeq_object_decontamed_filtered.rds")
summarize_phyloseq(physeq_object)
basic_info_physeq_object(physeq_object)

# "Lowest readnumber is 7301"
# "Highest readnumber is 539187"
# "Lowest taxa sum is 2"
# "Highest taxa sum is 225095"

## Trimming and transforming data -----------------------------------------------------
# Filter to ASVs with at least 10 reads per sample, present in at least 5% of the samples
physeq_pruned <- filter_taxa(physeq_object, function(x)  sum(x>=10) > (0.05*length(x)), prune = T)

summarize_phyloseq(physeq_pruned)
basic_info_physeq_object(physeq_pruned)
# "Lowest readnumber is 5876"
# "Highest readnumber is 379436"
# "Lowest taxa sum is 78"
# "Highest taxa sum is 225095"

#transform counts to rclr for distances/ordinations
physeq_rclr <- microbiome::transform(physeq_pruned, "rclr")
summarize_phyloseq(physeq_rclr)
basic_info_physeq_object(physeq_rclr)


# Select dataframes for the statistics --------------------------------------------------
# create tiblle to compute distances and for plotting 
source("tidy_tibble_maker.R")
tt_rclr<- tidy_tibble_maker(physeq_rclr)
head(tt_rclr)

# Remove underscores from data to get nicer plots, and remove filters from ordination
tt.1 <- tt_rclr %>%  mutate(Location = ifelse(Location == "Crooks_Castle", "Crooks Castle",
                                         ifelse(Location == "Charles_Brown", "Charles Brown", Location)))
tt.1 <- tt.1 %>%  mutate(Treatment = ifelse(Treatment == "no_UV", "no UV", Treatment))
head(tt.1)
tt.filt <- tt.1 %>% filter(Phase != "Filter" ) 
unique(tt.filt$Phase)

tt.inc <- tt.filt %>% filter(Method == "Incubation" ) 
unique(tt.inc$Location_Habitat)

## Create df for all 3 locations ---------------------------------------------------------
# Select metadata
ASV.metadata <- tt.filt %>% select(Sample,Description, Location, Habitat, Polymer, Isotope, Backbone,
                                   Treatment, Method, Location_Habitat, Polymer_Treatment, Polymer_Isotope) %>% 
              distinct()  %>% arrange(Sample) %>%  column_to_rownames(var = "Sample")

head(ASV.metadata)
colnames(ASV.metadata)
rownames(ASV.metadata)

# Extract RA data from tidy tibble 
# Filter out all rel_abundance = 0 values, to safe computing time w ordination
ASV.rclr <-  tt.filt %>%  select(Sample, OTU, Abundance)%>% mutate(across(c(OTU),factor)) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance)  %>%  
  arrange(Sample)%>% column_to_rownames(var = "Sample") %>% as.data.frame()

rownames(ASV.rclr)

# Are the rownames of both DFs matching?
rownames(ASV.rclr) == rownames(ASV.metadata)
# yes :) 

## Calculate distances for 3 locations --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
ASV.rclr.aitd <- vegdist(ASV.rclr, method = "euclidean")

ASV.rclr.aitd.pcoa <- wcmdscale(ASV.rclr.aitd, eig = T)
# Lets plot the relative eigenvalues of all axes
barplot(as.vector(ASV.rclr.aitd.pcoa$eig)/sum(ASV.rclr.aitd.pcoa$eig))

pct_xplnd <- (as.vector(ASV.rclr.aitd.pcoa$eig)/sum(ASV.rclr.aitd.pcoa$eig)) *100

#How much do the first 2 axes explain in total?
sum((as.vector(ASV.rclr.aitd.pcoa$eig)/sum(ASV.rclr.aitd.pcoa$eig))[1:2]) * 100
# ~21%

#and the first 3
sum((as.vector(ASV.rclr.aitd.pcoa$eig)/sum(ASV.rclr.aitd.pcoa$eig))[1:3]) * 100
# ~27%

# Set-up DF for plotting with GGPLOT ------------------------------------------------
data.scores <- as.data.frame(ASV.rclr.aitd.pcoa$points)[1:2] 
colnames(data.scores) <- c("PCo1", "PCo2")
dim(data.scores)

# Check
rownames(data.scores) == rownames(ASV.metadata)
# True

data.scores$Description = ASV.metadata$Description
data.scores$Location= ASV.metadata$Location
data.scores$Location_Habitat = ASV.metadata$Location_Habitat
data.scores$Polymer = ASV.metadata$Polymer
data.scores$Polymer_Isotope = ASV.metadata$Polymer_Isotope
data.scores$Habitat= ASV.metadata$Habitat
data.scores$Backbone = ASV.metadata$Backbone
data.scores$Treatment = ASV.metadata$Treatment

head(data.scores)

shapes <- c(15,17,18,19,4,8,13)

## Create NMDS ordination plot
PCoA_plot <- 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = PCo1, y = PCo2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 5, aes(color = Treatment, shape = Location_Habitat ,stroke = 1))  +
  theme_pubr() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10, colour ="black"),
        legend.position="right", legend.box = "vertical", 
        axis.title.y = element_text(face = "bold", size = 11),
        axis.title.x = element_text(face = "bold", size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1)) +
  # Set axislabels and title 
  labs(title = "PCoA plastic communities")  + 
 scale_shape_manual(values = shapes) + #<- Specify shapes you want as vector. 
  scale_fill_manual(values = pal.uv) +      #<- Pick the correct colorpalette, based on your fill 
  scale_color_manual(values = pal.uv) +     #<- Pick the correct colorpalette, based on your fill 
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL)),  #<- Costumize the legend. All color swatches are now squares, to avoid confusion with the other shapes
         shape = guide_legend(override.aes = list(fill = "black", colopr = "black")))  #<- Make ure shapes are same as the ones you picked, and are all filled and black
 
PCoA_plot

# Create df for incubations only ---------------------------------------------------------
# Select metadata
ASV.metadata <- tt.inc %>% select(Sample,Description, Location, Habitat, Polymer, Isotope, Backbone,
                                   Treatment, Method, Location_Habitat, Polymer_Treatment, Polymer_Isotope) %>% 
  distinct()  %>% arrange(Sample) %>%  column_to_rownames(var = "Sample")

head(ASV.metadata)
colnames(ASV.metadata)
rownames(ASV.metadata)

# Extract RA data from tidy tibble 
# Filter out all rel_abundance = 0 values, to safe computing time w ordination
ASV.rclr <-  tt.inc %>%  select(Sample, OTU, Abundance)%>% mutate(across(c(OTU),factor)) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance)  %>%  
  arrange(Sample)%>% column_to_rownames(var = "Sample") %>% as.data.frame()

rownames(ASV.rclr)

# Are the rownames of both DFs matching?
rownames(ASV.rclr) == rownames(ASV.metadata)
# yes :) 

## Calculate distances for 3 locations --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
ASV.rclr.aitd <- vegdist(ASV.rclr, method = "euclidean")

ASV.rclr.aitd.pcoa <- wcmdscale(ASV.rclr.aitd, eig = T)
# Lets plot the relative eigenvalues of all axes
barplot(as.vector(ASV.rclr.aitd.pcoa$eig)/sum(ASV.rclr.aitd.pcoa$eig))

pct_xplnd <- (as.vector(ASV.rclr.aitd.pcoa$eig)/sum(ASV.rclr.aitd.pcoa$eig)) *100

#How much do the first 2 axes explain in total?
sum((as.vector(ASV.rclr.aitd.pcoa$eig)/sum(ASV.rclr.aitd.pcoa$eig))[1:2]) * 100
# ~22.7%

#and the first 3
sum((as.vector(ASV.rclr.aitd.pcoa$eig)/sum(ASV.rclr.aitd.pcoa$eig))[1:3]) * 100
# ~29.5%

# Set-up DF for plotting with GGPLOT ------------------------------------------------
data.scores <- as.data.frame(ASV.rclr.aitd.pcoa$points)[1:2] 
colnames(data.scores) <- c("PCo1", "PCo2")
dim(data.scores)

# Check
rownames(data.scores) == rownames(ASV.metadata)
# True

data.scores$Description = ASV.metadata$Description
data.scores$Location= ASV.metadata$Location
data.scores$Location_Habitat = ASV.metadata$Location_Habitat
data.scores$Polymer = ASV.metadata$Polymer
data.scores$Polymer_Isotope = ASV.metadata$Polymer_Isotope
data.scores$Habitat= ASV.metadata$Habitat
data.scores$Backbone = ASV.metadata$Backbone
data.scores$Treatment = ASV.metadata$Treatment

head(data.scores)


shapes <- c(15,17,18,19,4,8,13)

## Create NMDS ordination plot
PCoA_plot <- 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = PCo1, y = PCo2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 5, aes(color = Polymer_Isotope, shape = Treatment ,stroke = 1))  +
  theme_pubr() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10, colour ="black"),
        legend.position="right", legend.box = "vertical", 
        axis.title.y = element_text(face = "bold", size = 11),
        axis.title.x = element_text(face = "bold", size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1)) +
  # Set axislabels and title 
  labs(title = "PCoA plastic communities incubations")  + 
  scale_shape_manual(values = shapes) + #<- Specify shapes you want as vector. 
  scale_fill_manual(values = pal.pols.isotop) +      #<- Pick the correct colorpalette, based on your fill 
  scale_color_manual(values = pal.pols.isotop) +     #<- Pick the correct colorpalette, based on your fill 
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL)),  #<- Costumize the legend. All color swatches are now squares, to avoid confusion with the other shapes
         shape = guide_legend(override.aes = list(fill = "black", colopr = "black")))  #<- Make ure shapes are same as the ones you picked, and are all filled and black

PCoA_plot

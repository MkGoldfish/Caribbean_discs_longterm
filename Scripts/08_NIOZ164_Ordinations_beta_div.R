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
library(zCompositions)
library(microbiome)
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
library(glue)
library(cowplot)

# Colors for plotting --------------------------------------------------------------------
pal_isme <- c("#006d77", "#ffddd2", "#00C49A", "#e29578", "#83c5be")

pal.loc <- c("#FF6DB6FF" , "#004949FF",  "#66A61E")
# CB, CC, Zeelandia"
pal.habs<- c("#66A61E","#FF8E32FF", "#51C3CCFF")
pal.habs.i<- c("#FF8E32FF", "#51C3CCFF")
# Beach, Benthic, Pelagic
pal.loc.hab <- c("#FF6DB6FF","#FFB6DBFF", "#004949FF", "#009292FF", "#66A61E")
# CB_P, CB_B, CC_P, CB_B, Zeelandia
pal.CB <- c("#FF6DB6FF", "#FFB6DBFF")
pal.CC <- c("#004949FF", "#009292FF")
pal.pel <- c("#FF6DB6FF", "#004949FF")
pal.ben <- c("#FFB6DBFF", "#009292FF")

pal.pols.isotop <- c("#E31A1C", "#7570B3", "#1F78B4","#A6CEE3", "#E6AB02","#A6761D", "#E5C494","#1B9E77","#882255")
pal.pols <- c("#E31A1C", "#7570B3" ,"#1F78B4","#E6AB02", "#A6761D", "#1B9E77","#882255")
# PE;PE-13C;PP;PP-13C;PS;PET;Nylon;Blanco
pal.uv <- c("#332288","#DDCC77", "grey") 
# UV, noUV, NA
pal.uv.2 <- c("#332288","#DDCC77") 

 
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

# Remove underscores from data to get nicer plots
meta <- as.data.frame(as.matrix(sample_data(physeq_object))) 
colnames(meta)
meta.1 <- meta %>%  mutate(Location = ifelse(Location == "Crooks_Castle", "Crooks Castle",
                                              ifelse(Location == "Charles_Brown", "Charles Brown", Location)))
meta.2 <- meta.1 %>%  mutate(Treatment = ifelse(Treatment == "no_UV", "no UV", Treatment)) %>% sample_data()

physeq_object.1 <- merge_phyloseq(otu_table(physeq_object), tax_table(physeq_object), meta.2)

# Trimming and transforming data -----------------------------------------------------
# Filter to ASVs with at least 10 reads per sample, present in at least 5% of the samples
physeq_pruned <- filter_taxa(physeq_object.1, function(x)  sum(x>=10) > (0.05*length(x)), prune = T)

summarize_phyloseq(physeq_pruned)
basic_info_physeq_object(physeq_pruned)
# "Lowest readnumber is 5876"
# "Highest readnumber is 379436"
# "Lowest taxa sum is 78"
# "Highest taxa sum is 225095"

## Subset the physeq for different analyses and extract OTU tables and sample data ----------------
### 1. inc and wild ----
physeq_inc_wild<- physeq_pruned %>% subset_samples(Phase != "Filter") 
physeq_inc_wild <- prune_taxa(taxa_sums(physeq_inc_wild)>0, physeq_inc_wild)
basic_info_physeq_object(physeq_inc_wild)
otu_i.w <- as.data.frame(otu_table(physeq_inc_wild)) %>% t() # We need samples per rows and asv per columns
samp_i.w <- as.data.frame(as.matrix(sample_data(physeq_inc_wild)))
colnames(otu_i.w)
# "Total taxa is 2634"
# "Total samples is 68"
# "Lowest readnumber is 5876"
# "Highest readnumber is 379436"
# "Lowest taxa sum is 1"
# "Highest taxa sum is 225095"
# Number of singletons = 87"

### 2. inc only ----
physeq_inc <- subset_samples(physeq_pruned, Phase == "Disc")
physeq_inc <- prune_taxa(taxa_sums(physeq_inc)>0, physeq_inc)
summarize_phyloseq(physeq_inc)
basic_info_physeq_object(physeq_inc)
otu_i <- as.data.frame(otu_table(physeq_inc)) %>% t()
samp_i <- as.data.frame(as.matrix(sample_data(physeq_inc)))
# "Total taxa is 2634"
# "Total samples is 59"
# "Lowest readnumber is 5876"
# "Highest readnumber is 379436"
# "Lowest taxa sum is 0"
# "Highest taxa sum is 225095"
# Number of singletons = 119"

### 3. CC Only----
physeq_inc_CC <- subset_samples(physeq_inc, Location == "Crooks Castle")
physeq_inc_CC <- prune_taxa(taxa_sums(physeq_inc_CC)>0, physeq_inc_CC)
summarize_phyloseq(physeq_inc_CC)
basic_info_physeq_object(physeq_inc_CC)
otu_i.cc<- as.data.frame(otu_table(physeq_inc_CC)) %>% t()
samp_i.cc <- as.data.frame(as.matrix(sample_data(physeq_inc_CC)))
# "Total taxa is 2634"
# "Total samples is 29"
# "Lowest readnumber is 31016"
# "Highest readnumber is 341333"
# "Lowest taxa sum is 0"
# "Highest taxa sum is 225095"
# Number of singletons = 397"

### 4. CB Only   ----
physeq_inc_CB <- subset_samples(physeq_inc, Location == "Charles Brown")
physeq_inc_CB <- prune_taxa(taxa_sums(physeq_inc_CB)>0, physeq_inc_CB)
summarize_phyloseq(physeq_inc_CB)
basic_info_physeq_object(physeq_inc_CB)
otu_i.cb<- as.data.frame(otu_table(physeq_inc_CB)) %>% t()
samp_i.cb <- as.data.frame(as.matrix(sample_data(physeq_inc_CB)))
# "Total taxa is 2634"
# "Total samples is 30"
# "Lowest readnumber is 5876"
# "Highest readnumber is 379436"
# "Lowest taxa sum is 0"
# "Highest taxa sum is 83017"
# Number of singletons = 380"

### 5.Pelagic Only----
physeq_inc_P <- subset_samples(physeq_inc, Habitat == "Pelagic")
physeq_inc_P <- prune_taxa(taxa_sums(physeq_inc_P)>0, physeq_inc_P)
summarize_phyloseq(physeq_inc_P)
basic_info_physeq_object(physeq_inc_P)
otu_p<- as.data.frame(otu_table(physeq_inc_P)) %>% t()
samp_p <- as.data.frame(as.matrix(sample_data(physeq_inc_P)))


### 6. Benthic Only   ----
physeq_inc_B <- subset_samples(physeq_inc, Habitat == "Benthic")
physeq_inc_B <- prune_taxa(taxa_sums(physeq_inc_B)>0, physeq_inc_B)
summarize_phyloseq(physeq_inc_B)
basic_info_physeq_object(physeq_inc_B)
otu_b<- as.data.frame(otu_table(physeq_inc_B)) %>% t()
samp_b <- as.data.frame(as.matrix(sample_data(physeq_inc_B)))


## CLR transformation of ASV tables -----
# Check that samples are rows and ASVs are columns
rownames(otu_i)
zPatterns(otu_i, label = 0)

# Remove zero-count values
otu_i.w.no0 <- cmultRepl(otu_i.w, method = "CZM", output = "p-counts")
otu_i.no0 <- cmultRepl(otu_i, method = "CZM", output = "p-counts")
otu_i.cc.no0 <- cmultRepl(otu_i.cc, method = "CZM", output = "p-counts")
otu_i.cb.no0 <- cmultRepl(otu_i.cb, method = "CZM", output = "p-counts")
otu_p.no0 <- cmultRepl(otu_p, method = "CZM", output = "p-counts")
otu_b.no0 <- cmultRepl(otu_b, method = "CZM", output = "p-counts")

# perform clr 
clr_i.w <- as.data.frame(clr(otu_i.w.no0))
clr_i <- as.data.frame(clr(otu_i.no0))
clr_i_CC <- as.data.frame(clr(otu_i.cc.no0))
clr_i_CB <- as.data.frame(clr(otu_i.cb.no0))
clr_i_P <- as.data.frame(clr(otu_p.no0))
clr_i_B <- as.data.frame(clr(otu_b.no0))

# Check that clr_data and metadata rownames match
summary(rownames(samp_i.w) == rownames(clr_i.w))
summary(rownames(samp_i) == rownames(clr_i))
summary(rownames(samp_i.cc) == rownames(clr_i_CC))
summary(rownames(samp_i.cb) == rownames(clr_i_CB))
summary(rownames(samp_p) == rownames(clr_i_P))
summary(rownames(samp_b) == rownames(clr_i_B))


# 1. Calculate distances for 3 locations --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
aitd_i.w <- vegdist(clr_i.w, method = "euclidean")

pcoa_i.w <- wcmdscale(aitd_i.w, eig = T)
# Lets plot the relative eigenvalues of all axes
barplot(as.vector(pcoa_i.w$eig)/sum(pcoa_i.w$eig))

pct_xplnd <- pcoa_i.w$eig/sum(pcoa_i.w$eig) *100
pretty_px <- format(round(pct_xplnd[1:2], digits = 1), nsmall = 1, trim = TRUE)
labs <- c(glue("PCo1 ({pretty_px[1]}%)"),
          glue("PCo2 ({pretty_px[2]}%)"))

#How much do the first 2 axes explain in total?
sum((as.vector(pcoa_i.w$eig)/sum(pcoa_i.w$eig))[1:2]) * 100
# ~19.2%

#and the first 3
sum((as.vector(pcoa_i.w$eig)/sum(pcoa_i.w$eig))[1:3]) * 100
# ~23.7%

## Plotting with GGPLOT ------------------------------------------------
data.scores <- as.data.frame(pcoa_i.w$points)[1:2] 
colnames(data.scores) <- c("PCo1", "PCo2")
dim(data.scores)

data.scores$Description = samp_i.w$Description
data.scores$Location= samp_i.w$Location
data.scores$Location_Habitat = samp_i.w$Location_Habitat
data.scores$Polymer = samp_i.w$Polymer
data.scores$Polymer_Isotope = samp_i.w$Polymer_Isotope
data.scores$Habitat= samp_i.w$Habitat
data.scores$Backbone = samp_i.w$Backbone
data.scores$Treatment = samp_i.w$Treatment

head(data.scores)

shapes <- c(15,17,18,19,4,8,13)

## Create NMDS ordination plot
PCoA_plot_all <- 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = PCo1, y = PCo2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 5, aes(color = Location_Habitat, shape = Polymer, stroke = 1))  +
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
  labs(title = "PCoA plastic communities",
       x = labs[1], y = labs[2])  + 
 scale_shape_manual(values = shapes) + #<- Specify shapes you want as vector. 
  scale_fill_manual(values = pal.loc.hab) +      #<- Pick the correct colorpalette, based on your fill 
  scale_color_manual(values =pal.loc.hab) +     #<- Pick the correct colorpalette, based on your fill 
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL)),  #<- Costumize the legend. All color swatches are now squares, to avoid confusion with the other shapes
         shape = guide_legend(override.aes = list(fill = "black", colopr = "black")))  #<- Make ure shapes are same as the ones you picked, and are all filled and black
 
PCoA_plot_all

# 2. Calculate distances for 2 locations, incubations only --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
aitd_i <- vegdist(clr_i, method = "euclidean")

pcoa_i <- wcmdscale(aitd_i, eig = T)
# Lets plot the relative eigenvalues of all axes
barplot(as.vector(pcoa_i$eig)/sum(pcoa_i$eig))

pct_xplnd <- pcoa_i.w$eig/sum(pcoa_i.w$eig) *100
pretty_px <- format(round(pct_xplnd[1:2], digits = 1), nsmall = 1, trim = TRUE)
labs <- c(glue("PCo1 ({pretty_px[1]}%)"),
          glue("PCo2 ({pretty_px[2]}%)"))

#How much do the first 2 axes explain in total?
sum((as.vector(pcoa_i$eig)/sum(pcoa_i$eig))[1:2]) * 100
# ~20.6%

#and the first 3
sum((as.vector(pcoa_i$eig)/sum(pcoa_i$eig))[1:3]) * 100
# ~25.3%

## Plotting with GGPLOT ------------------------------------------------
data.scores <- as.data.frame(pcoa_i$points)[1:2] 
colnames(data.scores) <- c("PCo1", "PCo2")
dim(data.scores)

data.scores$Description = samp_i$Description
data.scores$Location= samp_i$Location
data.scores$Location_Habitat = samp_i$Location_Habitat
data.scores$Polymer = samp_i$Polymer
data.scores$Polymer_Isotope = samp_i$Polymer_Isotope
data.scores$Habitat= samp_i$Habitat
data.scores$Backbone = samp_i$Backbone
data.scores$Treatment = samp_i$Treatment

head(data.scores)

shapes <- c(15,17,18,19,4,8,13)

## Create NMDS ordination plot
PCoA_plot_inc <- 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = PCo1, y = PCo2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 5, aes(color = Location_Habitat, shape = Polymer, stroke = 1))  +
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
  labs(x = labs[1], y = labs[2])  + 
  scale_shape_manual(values = shapes) + #<- Specify shapes you want as vector. 
  scale_fill_manual(values = pal.loc.hab) +      #<- Pick the correct colorpalette, based on your fill 
  scale_color_manual(values = pal.loc.hab) +     #<- Pick the correct colorpalette, based on your fill 
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL)),  #<- Costumize the legend. All color swatches are now squares, to avoid confusion with the other shapes
         shape = guide_legend(override.aes = list(fill = "black", colopr = "black")))  #<- Make ure shapes are same as the ones you picked, and are all filled and black

PCoA_plot_inc

# 3. Calculate distances for CC only --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
aitd_i_cc <- vegdist(clr_i_CC, method = "euclidean")

pcoa_cc <- wcmdscale(aitd_i_cc, eig = T)
# Lets plot the relative eigenvalues of all axes
barplot(as.vector(pcoa_cc$eig)/sum(pcoa_cc$eig))

pct_xplnd <- pcoa_cc$eig/sum(pcoa_cc$eig) *100
pretty_px <- format(round(pct_xplnd[1:2], digits = 1), nsmall = 1, trim = TRUE)
labs <- c(glue("PCo1 ({pretty_px[1]}%)"),
          glue("PCo2 ({pretty_px[2]}%)"))

#How much do the first 2 axes explain in total?
sum((as.vector(pcoa_cc$eig)/sum(pcoa_cc$eig))[1:2]) * 100
# ~24%

#and the first 3
sum((as.vector(pcoa_cc$eig)/sum(pcoa_cc$eig))[1:3]) * 100
# ~30.5%

## Plotting with GGPLOT ------------------------------------------------
data.scores <- as.data.frame(pcoa_cc$points)[1:2] 
colnames(data.scores) <- c("PCo1", "PCo2")
dim(data.scores)

data.scores$Description = samp_i.cc$Description
data.scores$Location= samp_i.cc$Location
data.scores$Location_Habitat = samp_i.cc$Location_Habitat
data.scores$Polymer = samp_i.cc$Polymer
data.scores$Polymer_Isotope = samp_i.cc$Polymer_Isotope
data.scores$Habitat= samp_i.cc$Habitat
data.scores$Backbone = samp_i.cc$Backbone
data.scores$Treatment = samp_i.cc$Treatment

head(data.scores)

shapes <- c(15,17,18,19,4,8,13)

## Create NMDS ordination plot
PCoA_plot_CC <- 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = PCo1, y = PCo2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 5, aes(color = Habitat, shape = Polymer, stroke = 1))  +
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
  labs(x = labs[1], y = labs[2])  + 
  scale_shape_manual(values = shapes) + #<- Specify shapes you want as vector. 
  scale_fill_manual(values = pal.habs.i) +      #<- Pick the correct colorpalette, based on your fill 
  scale_color_manual(values = pal.habs.i) +     #<- Pick the correct colorpalette, based on your fill 
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL)),  #<- Costumize the legend. All color swatches are now squares, to avoid confusion with the other shapes
         shape = guide_legend(override.aes = list(fill = "black", colopr = "black")))  #<- Make ure shapes are same as the ones you picked, and are all filled and black

PCoA_plot_CC

# 4. Calculate distances for CB only --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
aitd_i_cb <- vegdist(clr_i_CB, method = "euclidean")

pcoa_cb <- wcmdscale(aitd_i_cb, eig = T)
# Lets plot the relative eigenvalues of all axes
barplot(as.vector(pcoa_cb$eig)/sum(pcoa_cb$eig))

pct_xplnd <- pcoa_cb$eig/sum(pcoa_cb$eig) *100
pretty_px <- format(round(pct_xplnd[1:2], digits = 1), nsmall = 1, trim = TRUE)
labs <- c(glue("PCo1 ({pretty_px[1]}%)"),
          glue("PCo2 ({pretty_px[2]}%)"))

#How much do the first 2 axes explain in total?
sum((as.vector(pcoa_cb$eig)/sum(pcoa_cb$eig))[1:2]) * 100
# ~27.8%

#and the first 3
sum((as.vector(pcoa_cb$eig)/sum(pcoa_cb$eig))[1:3]) * 100
# ~35.11%

## Plotting with GGPLOT ------------------------------------------------
data.scores <- as.data.frame(pcoa_cb$points)[1:2] 
colnames(data.scores) <- c("PCo1", "PCo2")
dim(data.scores)

data.scores$Description = samp_i.cb$Description
data.scores$Location= samp_i.cb$Location
data.scores$Location_Habitat = samp_i.cb$Location_Habitat
data.scores$Polymer = samp_i.cb$Polymer
data.scores$Polymer_Isotope = samp_i.cb$Polymer_Isotope
data.scores$Habitat= samp_i.cb$Habitat
data.scores$Backbone = samp_i.cb$Backbone
data.scores$Treatment = samp_i.cb$Treatment

head(data.scores)

shapes <- c(15,17,18,19,4,8,13)

## Create NMDS ordination plot
PCoA_plot_CB <- 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = PCo1, y = PCo2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 5, aes(color = Habitat, shape = Polymer, stroke = 1))  +
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
  labs(x = labs[1], y = labs[2])  + 
  scale_shape_manual(values = shapes) + #<- Specify shapes you want as vector. 
  scale_fill_manual(values = pal.habs.i) +      #<- Pick the correct colorpalette, based on your fill 
  scale_color_manual(values = pal.habs.i) +     #<- Pick the correct colorpalette, based on your fill 
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL)),  #<- Costumize the legend. All color swatches are now squares, to avoid confusion with the other shapes
         shape = guide_legend(override.aes = list(fill = "black", colopr = "black")))  #<- Make ure shapes are same as the ones you picked, and are all filled and black

PCoA_plot_CB

# Grid of both incubation locations w cowplot ---------------
legend.a <- get_legend(PCoA_plot_CC)

plot_grid(PCoA_plot_CB + theme(legend.position ="none"),
          legend.a,
          PCoA_plot_CC + theme(legend.position ="none"),
          ncol = 3,
          nrow = 1,
          labels = c('A', '', 'B'),
          align = 'v',
          axis = "l",
          rel_widths = c(1,0.2,1))

# 5. Calculate distances for Pelagic only --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
aitd_P <- vegdist(clr_i_P, method = "euclidean")

pcoa_P <- wcmdscale(aitd_P, eig = T)
# Lets plot the relative eigenvalues of all axes
barplot(as.vector(pcoa_P$eig)/sum(pcoa_P$eig))

pct_xplnd <- pcoa_P$eig/sum(pcoa_P$eig) *100
pretty_px <- format(round(pct_xplnd[1:2], digits = 1), nsmall = 1, trim = TRUE)
labs <- c(glue("PCo1 ({pretty_px[1]}%)"),
          glue("PCo2 ({pretty_px[2]}%)"))

#How much do the first 2 axes explain in total?
sum((as.vector(pcoa_P$eig)/sum(pcoa_P$eig))[1:2]) * 100
# ~22%

#and the first 3
sum((as.vector(pcoa_P$eig)/sum(pcoa_P$eig))[1:3]) * 100
# ~30.1%

## Plotting with GGPLOT ------------------------------------------------
data.scores <- as.data.frame(pcoa_P$points)[1:2] 
colnames(data.scores) <- c("PCo1", "PCo2")
dim(data.scores)

data.scores$Description = samp_p$Description
data.scores$Location= samp_p$Location
data.scores$Location_Habitat = samp_p$Location_Habitat
data.scores$Polymer = samp_p$Polymer
data.scores$Polymer_Isotope = samp_p$Polymer_Isotope
data.scores$Habitat= samp_p$Habitat
data.scores$Backbone = samp_p$Backbone
data.scores$Treatment = samp_p$Treatment

head(data.scores)

shapes <- c(15,17,18,19,4,8,13)

## Create NMDS ordination plot
PCoA_plot_P <- 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = PCo1, y = PCo2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 5, aes(color = Location, shape = Polymer, stroke = 1))  +
  theme_pubr() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10, colour ="black"),
        legend.position="bottom", legend.box = "vertical", 
        axis.title.y = element_text(face = "bold", size = 11),
        axis.title.x = element_text(face = "bold", size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1)) +
  # Set axislabels and title 
  labs(x = labs[1], y = labs[2])  + 
  scale_shape_manual(values = shapes) + #<- Specify shapes you want as vector. 
  scale_fill_manual(values = pal.pel) +      #<- Pick the correct colorpalette, based on your fill 
  scale_color_manual(values = pal.pel) +     #<- Pick the correct colorpalette, based on your fill 
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL)),  #<- Costumize the legend. All color swatches are now squares, to avoid confusion with the other shapes
         shape = guide_legend(override.aes = list(fill = "black", colopr = "black")))  #<- Make ure shapes are same as the ones you picked, and are all filled and black

PCoA_plot_P

# 6. Calculate distances for Benthic only --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
aitd_B <- vegdist(clr_i_B, method = "euclidean")

pcoa_B <- wcmdscale(aitd_B, eig = T)
# Lets plot the relative eigenvalues of all axes
barplot(as.vector(pcoa_B$eig)/sum(pcoa_B$eig))

pct_xplnd <- pcoa_B$eig/sum(pcoa_B$eig) *100
pretty_px <- format(round(pct_xplnd[1:2], digits = 1), nsmall = 1, trim = TRUE)
labs <- c(glue("PCo1 ({pretty_px[1]}%)"),
          glue("PCo2 ({pretty_px[2]}%)"))

#How much do the first 2 axes explain in total?
sum((as.vector(pcoa_B$eig)/sum(pcoa_B$eig))[1:2]) * 100
# ~19.2%

#and the first 3
sum((as.vector(pcoa_B$eig)/sum(pcoa_B$eig))[1:3]) * 100
# ~25.8%

## Plotting with GGPLOT ------------------------------------------------
data.scores <- as.data.frame(pcoa_B$points)[1:2] 
colnames(data.scores) <- c("PCo1", "PCo2")
dim(data.scores)

data.scores$Description = samp_b$Description
data.scores$Location= samp_b$Location
data.scores$Location_Habitat = samp_b$Location_Habitat
data.scores$Polymer = samp_b$Polymer
data.scores$Polymer_Isotope = samp_b$Polymer_Isotope
data.scores$Habitat= samp_b$Habitat
data.scores$Backbone = samp_b$Backbone
data.scores$Treatment = samp_b$Treatment

head(data.scores)

shapes <- c(15,17,18,19,4,8,13)

## Create NMDS ordination plot
PCoA_plot_B <- 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = PCo1, y = PCo2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 5, aes(color = Location, shape = Polymer, stroke = 1))  +
  theme_pubr() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10, colour ="black"),
        legend.position="bottom", legend.box = "vertical", 
        axis.title.y = element_text(face = "bold", size = 11),
        axis.title.x = element_text(face = "bold", size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1)) +
  # Set axislabels and title 
  labs(x = labs[1], y = labs[2])  + 
  scale_shape_manual(values = shapes) + #<- Specify shapes you want as vector. 
  scale_fill_manual(values = pal.ben) +      #<- Pick the correct colorpalette, based on your fill 
  scale_color_manual(values = pal.ben) +     #<- Pick the correct colorpalette, based on your fill 
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL)),  #<- Costumize the legend. All color swatches are now squares, to avoid confusion with the other shapes
         shape = guide_legend(override.aes = list(fill = "black", colopr = "black")))  #<- Make ure shapes are same as the ones you picked, and are all filled and black

PCoA_plot_B

# Grid of both incubation locations w cowpolt ---------------
legend.a <- get_legend(PCoA_plot_CC)

plot_grid(PCoA_plot_P,
          PCoA_plot_B,
          ncol = 2,
          nrow = 1,
          labels = c('A', 'B'),
          align = 'v',
          axis = "tb",
          rel_widths = c(1,1))

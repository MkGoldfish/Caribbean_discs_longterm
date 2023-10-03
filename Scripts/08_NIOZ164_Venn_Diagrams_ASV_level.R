##%#####################################################################################%##
#NIOZ164 Statia Discs - Venn Diagrams shared ASVs                                  #####                                            
#Author: Maaike Goudriaan, NIOZ, MMB                                                      #
#                                                                                         #
# Purpose:                                                                                # 
# NIOZ164 is an illumina sequencing lane consisting of different samplesets.              # 
# Here we focus on analysis of 16S amplicon data, from amplified DNA                      #
# extracted from discs covered in polymer films that were incubated at 2 different        #
# locations in 2 different waterdepths (seafloor and watercolumn) in                      #
# coastal Caribbean waters close to the island of St. EUstatius.                          #
# The tested polymers were PE, PP, PS, PEt and Nylon, PE-13C and PP-13C.                  #
# Of each polymere there was an UV pretreated and a non-treated version.                  #
#                                                                                         #
# Determine shared ASVs among samples                                                #
##%#####################################################################################%##

# Date: 2023 - 08 - 23
# R-version: 4.3.1 

## Set working directory ------------------------------------------------------------------
setwd("C:/Users/mgoudriaan/Documents/GitHub/Caribbean_discs_longterm/Scripts")
set.seed(42)

# Load libraries -------------------------------------------------------------------------
library(ggplot2)
library(ggVennDiagram)
library(stringr)
library(tidyverse)
library(rvg)
library(officer)
library(MicEco)
library(eulerr)
library(cowplot)

# Import data ----------------------------------------------------------------------------
source("basic_info_physeq_object.R")
# Import earlier created physeq object
physeq_object <- readRDS("../Analysis/NIOZ164_physeq_object_subset_decontamed_tax.corrected_pruned.rds")
summarize_phyloseq(physeq_object)
basic_info_physeq_object(physeq_object)

# "Lowest readnumber is 7301"
# "Highest readnumber is 549477"
# "Lowest taxa sum is 2"
# "Highest taxa sum is 225095"

# what metadata do we have?
meta <- as.data.frame(as.matrix(sample_data(physeq_object))) 
head(meta)

# Remove underscores from data to get nicer plots
meta.1 <- meta %>%  mutate(Location = ifelse(Location == "Crooks_Castle", "Crooks Castle",
                                             ifelse(Location == "Charles_Brown", "Charles Brown", Location)))
meta.2 <- meta.1 %>%  mutate(Treatment = ifelse(Treatment == "no_UV", "no UV", Treatment)) %>% sample_data()

physeq_object.1 <- merge_phyloseq(otu_table(physeq_object), tax_table(physeq_object), meta.2)

# Import RA tibble to calculate percentages of top taxa
tt <- read.csv('../Processed-Data/NIOZ164_EUX_discs_RA_tidy_data_decontamed_tax.correct_pruned.csv', row.names = NULL)

# Trimming and transforming data -----------------------------------------------------
# Filter to ASVs with at least 10 reads per sample, present in at least 5% of the samples
# physeq_pruned <- filter_taxa(physeq_object.1, function(x)  sum(x>=10) > (0.05*length(x)), prune = T)
# 
# summarize_phyloseq(physeq_pruned)
# basic_info_physeq_object(physeq_pruned)
# # "Lowest readnumber is 5544"
# # "Highest readnumber is  391980"
# # "Lowest taxa sum is 55"
# # "Highest taxa sum is 225095"
# 
# colnames(sample_data(physeq_pruned))
# ntaxa(physeq_pruned)

# Subset the physeq for different locations and habitats, OTU level ----------------
ps_wild <- physeq_object.1 %>% subset_samples(Location == "Zeelandia") 
ps_inc <- physeq_object.1 %>% subset_samples(Phase == "Disc")
ps_CC <- physeq_object.1%>% subset_samples(Location == "Crooks Castle") 
ps_CB <- physeq_object.1 %>% subset_samples(Location == "Charles Brown") 
ps_pel <- physeq_object.1 %>% subset_samples(Habitat == "Pelagic") 
ps_bent <- physeq_object.1 %>% subset_samples(Habitat == "Benthic") 

# Use MicEco to create Euler diagram from physeq objects ----------------------------
pal.loc <- c("#FF6DB6FF" , "#009292FF",  "#66A61E")
# CB, CC, Zeelandia
pal.habs<- c("#CC5800FF", "#51C3CCFF")
# Benthic, Pelagic
pal.loc.hab <- c("#FF6DB6FF","#FFB6DBFF", "#004949FF", "#009292FF", "#66A61E")
# CB_P, CB_B, CC_P, CB_B, Zeelandia
pal.pols.isotop <- c("#E31A1C", "#7570B3", "#1F78B4","#A6CEE3", "#E6AB02","#A6761D", "#E5C494","#1B9E77")
# Blanco;Nylon;PE;PE-13C;PET;PP;PP-13C;PS;
pal.uv <- c("#DDCC77","#332288") 


venn.3locations <- ps_euler(physeq_object.1, "Location", fraction = 0.01, 
                      weight = F, plot = T, relative = F, 
                      fills = list(fill = c("#FF6DB6FF" , "#009292FF",  "#66A61E")),
                      edges = list(col = "white", lwd = 3), 
                      labels = list(fontsize = 18, col = "black"),
                      shape = "ellipse", 
                      quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "black", fontface = "bold"),
                      legend = F)

venn.3locations

venn.inc.loc <- ps_euler(ps_inc, "Location", fraction = 0.01, 
                     weight = F, plot = T, relative = F, 
                     fills = list(fill = c( "#FF6DB6FF" , "#009292FF")),
                     edges = list(col = "white", lwd = 3), 
                     labels = list(fontsize = 18, col = "black"),
                     shape = "ellipse", 
                     quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "black", fontface = "bold"),
                     legend = F)

venn.inc.loc

venn.inc.hab <- ps_euler(ps_inc, "Habitat", fraction = 0.01, 
                         weight = F, plot = T, relative = F, 
                         fills = list(fill = c( "#CC5800FF", "#51C3CCFF")),
                         edges = list(col = "white", lwd = 3), 
                         labels = list(fontsize = 18, col = "black"),
                         shape = "ellipse", 
                         quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "black", fontface = "bold"),
                         legend = F)

venn.inc.hab

venn.pelagic <- ps_euler(ps_pel, "Location", fraction = 0.01, 
                         weight = F, plot = T, relative = F, 
                         fills = list(fill = c("#FFB6DBFF", "#009292FF")),
                         edges = list(col = "white", lwd = 3), 
                         labels = list(fontsize = 18, col = "black"),
                         shape = "ellipse", 
                         quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "black", fontface = "bold"),
                         legend = F)

venn.pelagic

venn.benthic <- ps_euler(ps_bent, "Location", fraction = 0.01, 
                         weight = F, plot = T, relative = F, 
                         fills = list(fill = c("#FF6DB6FF" , "#004949FF")),
                         edges = list(col = "white", lwd = 3), 
                         labels = list(fontsize = 18, col = "black"),
                         shape = "ellipse", 
                         quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "black", fontface = "bold"),
                         legend = F)

venn.benthic

venn.CB <- ps_euler(ps_CB, "Habitat", fraction = 0.01, 
                    weight = F, plot = T, relative = F, 
                    fills = list(fill = c("#FF6DB6FF","#FFB6DBFF")),
                    edges = list(col = "white", lwd = 3), 
                    labels = list(fontsize = 18, col = "black"),
                    shape = "ellipse", 
                    quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "black", fontface = "bold"),
                    legend = F)

venn.CB

venn.CC <- ps_euler(ps_CC, "Habitat", fraction = 0.01, 
                   weight = F, plot = T, relative = F, 
                   fills = list(fill = c("#004949FF", "#009292FF")),
                   edges = list(col = "white", lwd = 3), 
                   labels = list(fontsize = 18, col = "black"),
                   shape = "ellipse", 
                   quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "black", fontface = "bold"),
                   legend = F)

venn.CC

## Plot grid of Venns ----------
plot_grid(venn.3locations,
          NULL,
          venn.inc.loc,
          venn.inc.hab,
          venn.pelagic,
          venn.benthic, 
          venn.CB,
          venn.CC,
          ncol = 2,
          align = 'v',
          axis = "tbrl",
          labels = c("A", " ", "B", "C", "D", "E", "F", "G"),
          label_size = 15,
          rel_heights = c(1.5,1,1,1),
          rel_widths = c(1,1,1))


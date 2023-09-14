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

# Permutations for Permanova and Permdisp
perm <- 9999

# Load libraries -------------------------------------------------------------------------
library(devtools)
library(phyloseq)
library(vegan)
library(compositions)
library(ade4)
library(FactoMineR)
library(labdsv)
library(zCompositions)
library(microbiome)
library(glue)
library(tidyverse)


# # Colors for plotting --------------------------------------------------------------------
# pal_isme <- c("#006d77", "#ffddd2", "#00C49A", "#e29578", "#83c5be")
# 
# pal.loc <- c("#FF6DB6FF" , "#004949FF",  "#66A61E")
# # CB, CC, Zeelandia"
# pal.habs<- c("#66A61E","#FF8E32FF", "#51C3CCFF")
# pal.habs.i<- c("#FF8E32FF", "#51C3CCFF")
# # Beach, Benthic, Pelagic
# pal.loc.hab <- c("#FF6DB6FF","#FFB6DBFF", "#004949FF", "#009292FF", "#66A61E")
# # CB_P, CB_B, CC_P, CB_B, Zeelandia
# 
# pal.pols.isotop <- c("#E31A1C", "#7570B3", "#1F78B4","#A6CEE3", "#E6AB02","#A6761D", "#E5C494","#1B9E77","#882255")
# pal.pols <- c("#E31A1C", "#7570B3" ,"#1F78B4","#E6AB02", "#A6761D", "#1B9E77","#882255")
# # PE;PE-13C;PP;PP-13C;PS;PET;Nylon;Blanco
# pal.uv <- c("#332288","#DDCC77", "grey") 
# # UV, noUV, NA
# pal.uv.2 <- c("#332288","#DDCC77") 

 
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

## CLR transformation of ASV tables -----
# Check that samples are rows and ASVs are columns
rownames(otu_i)
zPatterns(otu_i, label = 0)

# Remove zero-count values
otu_i.w.no0 <- cmultRepl(otu_i.w, method = "CZM", output = "p-counts")
otu_i.no0 <- cmultRepl(otu_i, method = "CZM", output = "p-counts")
otu_i.cc.no0 <- cmultRepl(otu_i.cc, method = "CZM", output = "p-counts")
otu_i.cb.no0 <- cmultRepl(otu_i.cb, method = "CZM", output = "p-counts")

# perform clr 
clr_i.w <- as.data.frame(clr(otu_i.w.no0))
clr_i <- as.data.frame(clr(otu_i.no0))
clr_i_CC <- as.data.frame(clr(otu_i.cc.no0))
clr_i_CB <- as.data.frame(clr(otu_i.cb.no0))

# Check that clr_data and metadata rownames match
summary(rownames(samp_i.w) == rownames(clr_i.w))
summary(rownames(samp_i) == rownames(clr_i))
summary(rownames(samp_i.cc) == rownames(clr_i_CC))
summary(rownames(samp_i.cb) == rownames(clr_i_CB))


# 1. Calculate distances for 3 locations --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
aitd_i.w <- vegdist(clr_i.w, method = "euclidean")

## PERMANOVA i.w. -------------------------------------------------------
Permanova_i.w.1 <- adonis2(aitd_i.w ~ Habitat * Location * Backbone * Treatment, data = samp_i.w, add = T, na.rm = T, permutations = perm )
Permanova_i.w.1 

Permanova_i.w.2 <- adonis2(aitd_i.w ~ Habitat + Location + Backbone + Treatment, data = samp_i.w, add = T, na.rm = T, permutations = perm )
Permanova_i.w.2

# Adjust p-values to q-values, add to df, write table to store results
q_i.w_p1 <- p.adjust(Permanova_i.w.1$`Pr(>F)`, method = "BH") %>% enframe(value = "qvalue")
Permanova_i.w.1.q <- bind_cols(Permanova_i.w.1,q_i.w_p1$qvalue, .name_repair = "minimal")

write.table(Permanova_i.w.1.q,"../Analysis/Permanova/PERM_inc.wild.full.txt" )


q_i.w_p2 <- p.adjust(Permanova_i.w.2$`Pr(>F)`, method = "BH") %>% enframe(value = "qvalue")
Permanova_i.w.2.q <- bind_cols(Permanova_i.w.2,q_i.w_p2$qvalue, .name_repair = "minimal")

write.table(Permanova_i.w.1.q,"../Analysis/Permanova/PERM_inc.wild_additive.txt" )

####
# Permanova_i.w.3 <- adonis2(aitd_i.w ~ Habitat * Location * (Backbone + Treatment), data = samp_i.w, add = T, na.rm = T, permutations = perm )
# Permanova_i.w.3 

# 2. Calculate distances for 2 locations, incubations only --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
aitd_i <- vegdist(clr_i, method = "euclidean")

## PERMANOVA i. -------------------------------------------------------
Permanova_i.1 <- adonis2(aitd_i ~ Habitat * Location * Backbone * Treatment, data = samp_i, add = T, na.rm = T, permutations = perm )
Permanova_i.1 

Permanova_i.2 <- adonis2(aitd_i ~ Habitat + Location + Backbone + Treatment, data = samp_i, add = T, na.rm = T, permutations = perm )
Permanova_i.2 

# Adjust p-values to q-values, add to df, write table to store results
q_i_p1 <- p.adjust(Permanova_i.1$`Pr(>F)`, method = "BH") %>% enframe(value = "qvalue")
Permanova_i.1.q <- bind_cols(Permanova_i.1,q_i_p1$qvalue, .name_repair = "minimal") 

write.table(Permanova_i.1.q,"../Analysis/Permanova/PERM_inc.full.txt" )



q_i_p2 <- p.adjust(Permanova_i.2$`Pr(>F)`, method = "BH") %>% enframe(value = "qvalue")
Permanova_i.2.q <- bind_cols(Permanova_i.2,q_i_p2$qvalue, .name_repair = "minimal")

write.table(Permanova_i.2.q,"../Analysis/Permanova/PERM_inc_additive.txt" )


# 3. Calculate distances for CC only --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
aitd_i_cc <- vegdist(clr_i_CC, method = "euclidean")

## PERMANOVA CC -------------------------------------------------------
Permanova_cc.1<- adonis2(aitd_i_cc ~ Habitat * Backbone * Treatment, data = samp_i.cc, add = T, na.rm = T, permutations = perm )
Permanova_cc.1

Permanova_cc.2 <- adonis2(aitd_i_cc ~ Habitat +  Backbone + Treatment, data = samp_i.cc, add = T, na.rm = T, permutations = perm )
Permanova_cc.2 

# Adjust p-values to q-values, add to df, write table to store results
q_cc_p1 <- p.adjust(Permanova_cc.1$`Pr(>F)`, method = "BH") %>% enframe(value = "qvalue")
Permanova_cc.1.q <- bind_cols(Permanova_cc.1,q_cc_p1$qvalue, .name_repair = "minimal") 

write.table(Permanova_cc.1.q,"../Analysis/Permanova/PERM_CC.full.txt" )



q_cc_p2 <- p.adjust(Permanova_cc.2$`Pr(>F)`, method = "BH") %>% enframe(value = "qvalue")
Permanova_cc.2.q <- bind_cols(Permanova_cc.2,q_cc_p2$qvalue, .name_repair = "minimal")

write.table(Permanova_i.w.1.q,"../Analysis/Permanova/PERM_CC_additive.txt" )


# 4. Calculate distances for CB only --------------------------------------------------------------
# Aitchison distance it the Euclidean distance of (r)clr transformed data
aitd_i_cb <- vegdist(clr_i_CB, method = "euclidean")

## PERMANOVA CB -------------------------------------------------------
Permanova_cb.1 <- adonis2(aitd_i_cb ~ Habitat * Backbone * Treatment, data = samp_i.cb, add = T, na.rm = T, permutations = perm )
Permanova_cb.1 

Permanova_cb.2 <- adonis2(aitd_i_cb ~ Habitat +  Backbone + Treatment, data = samp_i.cb, add = T, na.rm = T, permutations = perm )
Permanova_cb.2 

Permanova_cb.3 <- adonis2(aitd_i_cb ~ Habitat * (Backbone + Treatment), data = samp_i.cb, add = T, na.rm = T, permutations = perm )
Permanova_cb.3 

# Adjust p-values to q-values, add to df, write table to store results
q_cb_p1 <- p.adjust(Permanova_cb1$`Pr(>F)`, method = "BH") %>% enframe(value = "qvalue")
Permanova_cb.1.q <- bind_cols(Permanova_cb.1,q_cb_p1$qvalue, .name_repair = "minimal") 

write.table(Permanova_cc.1.q,"../Analysis/Permanova/PERM_CC.full.txt" )



q_cb_p2 <- p.adjust(Permanova_cb.2$`Pr(>F)`, method = "BH") %>% enframe(value = "qvalue")
Permanova_cb.2.q <- bind_cols(Permanova_cb.2,q_cb_p2$qvalue, .name_repair = "minimal")

write.table(Permanova_i.w.1.q,"../Analysis/Permanova/PERM_CC_additive.txt" )


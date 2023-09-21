##%#####################################################################################%##
#NIOZ164 Statia Discs - PERMANOVA & PERMDISP                                          #####                                            
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
# Test which factor causes significant community differences                              #
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
library(microbiome)
library(tidyverse)

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

# what metadata do we have?
meta <- as.data.frame(as.matrix(sample_data(physeq_object))) 
head(meta)


# Trimming and transforming data -----------------------------------------------------
# Filter to ASVs with at least 10 reads per sample, present in at least 5% of the samples
physeq_pruned.0 <- filter_taxa(physeq_object, function(x)  sum(x>=10) > (0.05*length(x)), prune = T)
physeq_pruned <- subset_taxa(physeq_pruned.0, !Order == "unassigned")

summarize_phyloseq(physeq_pruned)
basic_info_physeq_object(physeq_pruned)
# "Lowest readnumber is 5876"
# "Highest readnumber is 379436"
# "Lowest taxa sum is 78"
# "Highest taxa sum is 225095"

ps.pruned.gen <- aggregate_taxa(physeq_pruned, "Genus")
ps.pruned.gen.rel <-  microbiome::transform(ps.pruned.gen, "compositional")

## Subset the physeq for different locations and habitats ----------------
### 1. wild ----
ps.wild.gen.rel<- ps.pruned.gen.rel %>% subset_samples(Location == "Zeelandia") 
summarize_phyloseq(ps.wild.gen.rel)
basic_info_physeq_object(ps.wild.gen.rel)

### 2. inc only ----
ps.inc.gen.rel <- subset_samples(ps.pruned.gen.rel, Phase == "Disc")
summarize_phyloseq(ps.inc.gen.rel)
basic_info_physeq_object(ps.inc.gen.rel)

### 3. CC Only----
ps.cc.gen.rel <- subset_samples(ps.inc.gen.rel, Location == "Crooks_Castle")
summarize_phyloseq(ps.cc.gen.rel)
basic_info_physeq_object(ps.cc.gen.rel)

### 4. CB Only   ----
ps.cb.gen.rel  <- subset_samples(ps.inc.gen.rel, Location == "Charles_Brown")
summarize_phyloseq(ps.cb.gen.rel)
basic_info_physeq_object(ps.cb.gen.rel)

### 5.Pelagic Only----
ps.p.gen.rel <- subset_samples(ps.inc.gen.rel, Habitat == "Pelagic")
summarize_phyloseq(physeq_P)
basic_info_physeq_object(physeq_P)

### 6. Benthic Only   ----
ps.b.gen.rel <- subset_samples(ps.inc.gen.rel, Habitat == "Benthic")
summarize_phyloseq(physeq_B)
basic_info_physeq_object(physeq_B)

### 7. CC Pelagic----
ps.cc.p.gen.rel <- subset_samples(ps.cc.gen.rel, Habitat == "Pelagic")
summarize_phyloseq(physeq_CC_P)
basic_info_physeq_object(physeq_CC_P)

### 8. CC Benthic----
ps.cc.b.gen.rel<- subset_samples(ps.cc.gen.rel, Habitat == "Benthic")
summarize_phyloseq(physeq_CC_B)
basic_info_physeq_object(physeq_CC_B)

### 9. CB Pelagic----
ps.cb.p.gen.rel <- subset_samples(ps.cc.gen.rel, Habitat == "Pelagic")
summarize_phyloseq(physeq_CB_P)
basic_info_physeq_object(physeq_CB_P)

### 10. CC Benthic----
ps.cb.b.gen.rel <- subset_samples(ps.cc.gen.rel, Habitat == "Benthic")
summarize_phyloseq(physeq_CB_B)
basic_info_physeq_object(physeq_CB_B)


# A. Calculate top genera per data subset --------------------------------------------------------------
top_wild <- top_taxa(ps.wild.gen.rel, n = 15)
top_inc <- top_taxa(ps.inc.gen.rel, n = 15)
top_CB <- top_taxa(ps.cb.gen.rel, n = 15)
top_CC <- top_taxa(ps.cc.gen.rel, n = 15)
top_B <- top_taxa(ps.b.gen.rel, n = 15)
top_P <- top_taxa(ps.p.gen.rel, n = 15)
top_CB_B <- top_taxa(ps.cb.b.gen.rel, n = 15)
top_CB_P <- top_taxa(ps.cb.p.gen.rel, n = 15)
top_CC_B <- top_taxa(ps.cc.b.gen.rel, n = 15)
top_CC_P <- top_taxa(ps.cc.p.gen.rel, n = 15)

# B.Determine core community per dataset --------------------------------------------------------------
# determine core taxa
# detection at least 1.0%, prevalence at least 10 of all the samples samples
# gives the names of the core genera under given conditions
core.wild <- core_members(ps.wild.gen.rel, detection = 0.01, prevalence = 1/4) %>% enframe(value = "genus")
length(core.wild) #12 core taxa

write.table(core.families.taxa , '../core_community/21_core_families_detection0.01_prevalence0.25.txt',
            na = "NA", dec = '.', quote = F, sep = "\t")

core.inc <- core_members(ps.inc.gen.rel, detection = 0.01, prevalence = 1/4) %>% enframe(value = "genus")
length(core.inc) #23 core taxa

core.CB <- core_members(ps.cb.gen.rel, detection = 0.01, prevalence = 1/4) %>% enframe(value = "genus")
length(core.CB) #25 core taxa

core.CC <- core_members(ps.cc.gen.rel, detection = 0.01, prevalence = 1/4) %>% enframe(value = "genus")
length(core.CC) #25 core taxa

core.B <- core_members(ps.b.gen.rel, detection = 0.01, prevalence = 1/4) %>% enframe(value = "genus")
length(core.B) #25 core taxa

core.P <- core_members(ps.p.gen.rel, detection = 0.01, prevalence = 1/4) %>% enframe(value = "genus")
length(core.P) #30 core taxa

core.CB.B <- core_members(ps.cb.b.gen.rel, detection = 0.01, prevalence = 1/4) %>% enframe(value = "genus")
length(core.CB.B) #25 core taxa

core.CB.P <- core_members(ps.cb.p.gen.rel, detection = 0.01, prevalence = 1/4) %>% enframe(value = "genus")
length(core.CB.P) #30 core taxa

core.CC.B <- core_members(ps.cb.b.gen.rel, detection = 0.01, prevalence = 1/4) %>% enframe(value = "genus")
length(core.CB.B) #25 core taxa

core.CC.P <- core_members(ps.cb.p.gen.rel, detection = 0.01, prevalence = 1/4) %>% enframe(value = "genus")
length(core.CB.P) #30 core taxa

core <- bind_rows(Incubation = core.inc, 
                  Wild = core.wild, 
                  .id = "Community" )

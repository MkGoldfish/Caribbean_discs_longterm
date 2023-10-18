##%#####################################################################################%##
#NIOZ164 Statia Discs - Core community & Top genera                                   #####                                            
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
# Determine top and core genera of samples                                                #
##%#####################################################################################%##

# Date: 2023 - 08 - 23
# R-version: 4.3.1 

## Set working directory ------------------------------------------------------------------
setwd("C:/Users/mgoudriaan/Documents/GitHub/Caribbean_discs_longterm/Scripts")
set.seed(42)

# Load libraries -------------------------------------------------------------------------
library(devtools)
library(phyloseq)
library(vegan)
library(microbiome)
library(tidyverse)

# Import data ----------------------------------------------------------------------------
source("basic_info_physeq_object.R")
# Import earlier created physeq object
physeq_object <- readRDS("../Analysis/NIOZ164_physeq_object_subset_decontamed_tax.corrected_pruned.rds")
summarize_phyloseq(physeq_object)
basic_info_physeq_object(physeq_object)

# "Lowest readnumber is 7301"
# "Highest readnumber is 539187"
# "Lowest taxa sum is 2"
# "Highest taxa sum is 225095"

# what metadata do we have?
meta <- as.data.frame(as.matrix(sample_data(physeq_object))) 
head(meta)

# Import RA tibble to calculate percentages of top taxa
tt <- read.csv('../Processed-Data/NIOZ164_EUX_discs_RA_tidy_data_decontamed_tax.correct_pruned.csv', row.names = NULL)
Phyla <- tt %>%  select( Description, Location, Phylum, Phylum_rel_abund_Sample, Sample_rel_abund, Sample_st_dev)  %>% distinct()
Orders<- tt %>%  select( Description, Location, Order, Order_rel_abund_Sample, Sample_rel_abund, Sample_st_dev)  %>% distinct()
Arch <- tt %>% filter(Kingdom == "Archaea")
  

# Trimming and transforming data -----------------------------------------------------
# Filter to ASVs with at least 10 reads per sample, present in at least 5% of the samples
physeq_pruned <- filter_taxa(physeq_object, function(x)  sum(x>=10) > (0.05*length(x)), prune = T)

summarize_phyloseq(physeq_pruned)
basic_info_physeq_object(physeq_pruned)
# "Lowest readnumber is 5876"
# "Highest readnumber is 379436"
# "Lowest taxa sum is 78"
# "Highest taxa sum is 225095"

# 1. Subset the physeq for different wild vs incubations on different tax-levels ----------------
ps.pruned.rel <-  microbiome::transform(physeq_pruned, "compositional")
ps.prune.rel.phyl <- aggregate_taxa(ps.pruned.rel , "Phylum") %>% subset_taxa(Phylum != "unassigned")
ps.prune.rel.ord <- aggregate_taxa(ps.pruned.rel , "Order") %>% subset_taxa(Order != "unassigned")

## 1.1 Detemine amount of Phyla and Orders in all samples together ----------------------
basic_info_physeq_object(ps.prune.rel.phyl)
 # "Total taxa is 31"
 # "Total samples is 68"

basic_info_physeq_object(ps.prune.rel.ord)
# "Total taxa is 138"
# "Total samples is 86"

### 1.1.1 wild only ----
ps.wild.phyl.rel<- ps.prune.rel.phyl %>% subset_samples(Location == "Zeelandia") 
summarize_phyloseq(ps.wild.phyl.rel)
basic_info_physeq_object(ps.wild.phyl.rel)
# "Total taxa is 31"
# "Total samples is 9"

ps.wild.ord.rel <- ps.prune.rel.ord %>% subset_samples(Location == "Zeelandia") 
summarize_phyloseq(ps.wild.ord.rel)
basic_info_physeq_object(ps.wild.ord.rel)
# "Total taxa is 138"
# "Total samples is 9"

### 1.1.2. inc only ----
ps.inc.phyl.rel<- ps.prune.rel.phyl %>% subset_samples(Phase == "Disc") 
summarize_phyloseq(ps.inc.phyl.rel)
basic_info_physeq_object(ps.inc.phyl.rel)
# "Total taxa is 31"
# "Total samples is 59"

ps.inc.ord.rel <- ps.prune.rel.ord %>% subset_samples(Phase == "Disc") 
summarize_phyloseq(ps.inc.ord.rel)
basic_info_physeq_object(ps.inc.ord.rel)
# "Total taxa is 138"
# "Total samples is 59"

## 1.A. Calculate top taxa per data subset --------------------------------------------------------------
### Phyla ----
top_wild <- top_taxa(ps.wild.phyl.rel, n = 7)
#[1] "Proteobacteria"   "Bacteroidota"     "Firmicutes"       "Actinobacteriota" "Cyanobacteria"    "Planctomycetota" "Dadabacteria"  
top_inc <- top_taxa(ps.inc.phyl.rel, n = 7)
# "Proteobacteria"   "Bacteroidota"     "Actinobacteriota" "Desulfobacterota" "Planctomycetota"  "Firmicutes" "Cyanobacteria"   
top_all <-top_taxa(ps.prune.rel.phyl, n = 7)
# "Proteobacteria"   "Bacteroidota"     "Actinobacteriota" "Firmicutes"       "Desulfobacterota" "Planctomycetota" "Cyanobacteria" 

### Calculate the percentage of reads belonging to these taxa in the sample set
top.phyl.wild <- Phyla %>% filter(Location == "Zeelandia") %>% filter(Phylum %in% top_wild) 
top.phyl.wild.avg <- top.phyl.wild %>%
   group_by(Description) %>% summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100
top.phyl.wild.sd <- top.phyl.wild %>% 
  group_by(Description) %>%  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

top.phyl.wild.avg
top.phyl.wild.sd


top.phyl.inc <- Phyla %>% filter(Location %in% c("Charles Brown", "Crooks Castle")) %>% filter(Phylum %in% top_inc) 
top.phyl.inc.avg <- top.phyl.inc %>% 
  group_by(Description) %>% summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100
top.phyl.inc.sd <- top.phyl.inc %>% 
  group_by(Description) %>%  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

top.phyl.inc.avg
top.phyl.inc.sd


top.phyl.all <- Phyla %>% filter(Phylum %in% top_all) 
top.phyl.all.avg <- top.phyl.all %>% 
  group_by(Description) %>% summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100
top.phyl.all.sd <- top.phyl.all %>% 
  group_by(Description) %>%  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

top.phyl.all.avg
top.phyl.all.sd


### Order ----
top_ord_wild <- top_taxa(ps.wild.ord.rel, n = 10)
# [1] "Rhodobacterales"   "Chitinophagales"   "Enterobacterales"  "Flavobacteriales"  "Caulobacterales"   "Staphylococcales"  
# "Rhizobiales"       "Lactobacillales"   "Burkholderiales"   "Pseudomonadales"

top_ord_inc <- top_taxa(ps.inc.ord.rel, n = 10)
# [1] "Rhodobacterales"     "Rhizobiales"         "Enterobacterales"    "Flavobacteriales"    "Chitinophagales"     
# "Corynebacteriales"   "Staphylococcales"    "Desulfuromonadales" "Pseudomonadales"     "Thiotrichales" 

top_ord_all <-top_taxa(ps.prune.rel.ord, n = 14)
# [1] "Rhodobacterales"     "Rhizobiales"         "Enterobacterales"    "Chitinophagales"     "Flavobacteriales"    "Corynebacteriales"   
# "Staphylococcales"    "Caulobacterales"    
# [9] "Pseudomonadales"     "Desulfuromonadales"  "Thiotrichales"       "Propionibacteriales" "Burkholderiales"     "Steroidobacterales"  

top.ord.all.avg <- Orders %>% filter(Order %in% top_ord_all) %>% 
  group_by(Description) %>% summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100 
top.ord.all.sd <- Orders %>% filter(Order %in% top_ord_all) %>% 
  group_by(Description) %>%  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

top.ord.all.avg
top.ord.all.sd

# 2. Subset the physeq for different locations and habitats, Genus level ----------------
physeq_pruned.1 <- subset_taxa(physeq_pruned, !Order == "unassigned")
ps.pruned.gen <- aggregate_taxa(physeq_pruned.1, "Genus")
ps.pruned.gen.rel <-  microbiome::transform(ps.pruned.gen, "compositional")

### 2.1. wild ----
ps.wild.gen.rel<- ps.pruned.gen.rel %>% subset_samples(Location == "Zeelandia") 
summarize_phyloseq(ps.wild.gen.rel)
basic_info_physeq_object(ps.wild.gen.rel)

### 2.2. inc only ----
ps.inc.gen.rel <- subset_samples(ps.pruned.gen.rel, Phase == "Disc")
summarize_phyloseq(ps.inc.gen.rel)
basic_info_physeq_object(ps.inc.gen.rel)

### 2.3. CC Only----
ps.cc.gen.rel <- subset_samples(ps.inc.gen.rel, Location == "Crooks_Castle")
summarize_phyloseq(ps.cc.gen.rel)
basic_info_physeq_object(ps.cc.gen.rel)

### 2.4. CB Only   ----
ps.cb.gen.rel  <- subset_samples(ps.inc.gen.rel, Location == "Charles_Brown")
summarize_phyloseq(ps.cb.gen.rel)
basic_info_physeq_object(ps.cb.gen.rel)

### 2.5.Pelagic Only----
ps.p.gen.rel <- subset_samples(ps.inc.gen.rel, Habitat == "Pelagic")
summarize_phyloseq(ps.p.gen.rel)
basic_info_physeq_object(ps.p.gen.rel)

### 2.6. Benthic Only   ----
ps.b.gen.rel <- subset_samples(ps.inc.gen.rel, Habitat == "Benthic")
summarize_phyloseq(ps.b.gen.rel)
basic_info_physeq_object(ps.b.gen.rel)

### 2.7. CC Pelagic----
ps.cc.p.gen.rel <- subset_samples(ps.cc.gen.rel, Habitat == "Pelagic")
summarize_phyloseq(ps.cc.p.gen.rel)
basic_info_physeq_object(ps.cc.p.gen.rel)

### 2.8. CC Benthic----
ps.cc.b.gen.rel<- subset_samples(ps.cc.gen.rel, Habitat == "Benthic")
summarize_phyloseq(ps.cc.b.gen.rel)
basic_info_physeq_object(ps.cc.b.gen.rel)

### 2.9. CB Pelagic----
ps.cb.p.gen.rel <- subset_samples(ps.cb.gen.rel, Habitat == "Pelagic")
summarize_phyloseq(ps.cb.p.gen.rel)
basic_info_physeq_object(ps.cb.p.gen.rel)

### 2.10. CC Benthic----
ps.cb.b.gen.rel <- subset_samples(ps.cb.gen.rel, Habitat == "Benthic")
summarize_phyloseq(ps.cb.b.gen.rel)
basic_info_physeq_object(ps.cb.b.gen.rel)


# 2.A. Calculate top genera per data subset --------------------------------------------------------------
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

# 2.B.Determine core community per dataset --------------------------------------------------------------
# determine core taxa
# detection at least 1.0%, prevalence at least 10 of all the samples samples
# gives the names of the core genera under given conditions
core.wild <- core_members(ps.wild.gen.rel, detection = 0.0025, prevalence = 1/2) %>% enframe(value = "genus")
dim(core.wild) #15 core taxa
core.abundance.wild <- sample_sums(core(ps.wild.gen.rel, detection = 0.0025, prevalence = 1/2))
mean(core.abundance.wild)
sd(core.abundance.wild)

# write.table(core.families.taxa , '../core_community/21_core_families_detection0.01_prevalence0.25.txt',
#             na = "NA", dec = '.', quote = F, sep = "\t")

core.inc <- core_members(ps.inc.gen.rel, detection = 0.0025, prevalence = 1/2) %>% enframe(value = "genus")
dim(core.inc) #17 core taxa
core.abundance.inc <- sample_sums(core(ps.inc.gen.rel, detection = 0.0025, prevalence = 1/2)) 
mean(core.abundance.inc)
sd(core.abundance.inc)

core.CB <- core_members(ps.cb.gen.rel, detection = 0.0025, prevalence = 1/2) %>% enframe(value = "genus")
dim(core.CB) #17 core taxa
core.abundance.CB <- sample_sums(core(ps.cb.gen.rel, detection = 0.0025, prevalence = 1/2)) 
mean(core.abundance.CB)
sd(core.abundance.CB)

core.CC <- core_members(ps.cc.gen.rel, detection = 0.0025, prevalence = 1/2) %>% enframe(value = "genus")
dim(core.CC) #23 core taxa
core.abundance.CC <- sample_sums(core(ps.cc.gen.rel, detection = 0.0025, prevalence = 1/2)) 
mean(core.abundance.CC)
sd(core.abundance.CC)

core.B <- core_members(ps.b.gen.rel, detection = 0.0025, prevalence = 1/2) %>% enframe(value = "genus")
dim(core.B) #38 core taxa
core.abundance.B <- sample_sums(core(ps.b.gen.rel, detection = 0.0025, prevalence = 1/2)) 
mean(core.abundance.B)
sd(core.abundance.B)

core.P <- core_members(ps.p.gen.rel, detection = 0.0025, prevalence = 1/2) %>% enframe(value = "genus")
dim(core.P) #17 core taxa
core.abundance.P <- sample_sums(core(ps.p.gen.rel, detection = 0.0025, prevalence = 1/2)) 
mean(core.abundance.P)
sd(core.abundance.P)

core.CB.P <- core_members(ps.cb.p.gen.rel, detection = 0.0025, prevalence = 1/2) %>% enframe(value = "genus")
dim(core.CB.P) #24 core taxa
core.abundance.CB.P <- sample_sums(core(ps.cb.p.gen.rel, detection = 0.0025, prevalence = 1/2)) 
mean(core.abundance.CB.P)
sd(core.abundance.CB.P)

core.CC.P <- core_members(ps.cc.p.gen.rel, detection = 0.001, prevalence = 1/2) %>% enframe(value = "genus")
dim(core.CC.P) #46 core taxa
core.abundance.B <- sample_sums(core(ps.b.gen.rel, detection = 0.0025, prevalence = 1/2)) 
mean(core.abundance.B)
sd(core.abundance.B)

core.CB.B <- core_members(ps.cb.b.gen.rel, detection = 0.0025, prevalence = 1/2) %>% enframe(value = "genus")
dim(core.CB.B) #40 core taxa
core.abundance.CB.B <- sample_sums(core(ps.cb.b.gen.rel, detection = 0.0025, prevalence = 1/2)) 
mean(core.abundance.CB.B)
sd(core.abundance.CB.B)

core.CC.B <- core_members(ps.cc.b.gen.rel, detection = 0.0025, prevalence = 1/2) %>% enframe(value = "genus")
dim(core.CC.B) #35 core taxa
core.abundance.CC.B <- sample_sums(core(ps.cc.b.gen.rel, detection = 0.0025, prevalence = 1/2)) 
mean(core.abundance.CC.B)
sd(core.abundance.CC.B)



# Bind into one tibble, and pivot
core <- bind_rows(Incubation = core.inc, 
                  Wild = core.wild,
                  Benthic = core.B,
                  Pelagic = core.P,
                  CharlesBrown = core.CB,
                  CrooksCastle = core.CC,
                  Pelagic.CB = core.CB.P,
                  Pelagic.CC = core.CC.P,
                  Benthic.CB = core.CB.B,
                  Benthic.CC = core.CC.B,
                  .id = "Community" )
head(core)

# Pivot wider to get test per column, to compare core communities
core_wide <- core %>% pivot_wider(names_from = Community, values_from = genus)
head(core_wide)

write.csv(core_wide, "../Analysis/Core_Communites.csv", row.names = F)

## Compare core communities -------------------------------------
inc.vs.wild <- setdiff(core.inc, core.wild) #16
wild.vs.inc <- setdiff(core.wild, core.inc) #14
inc.and.wild <- intersect(core.inc, core.wild) #1 

P.vs.B <- setdiff(core.P, core.B) #17
B.vs.P<- setdiff(core.B, core.P) #38
P.and.B <- intersect(core.P, core.B) #0

CB.vs.CC <- setdiff(core.CB, core.CC) #16
CC.vs.CB<- setdiff(core.CC, core.CB) #22
CB.and.CC <- intersect(core.CB, core.CC) #1

P.CB.vs.CC <- setdiff(core.CB.P, core.CC.P) #22
P.CC.vs.CB<- setdiff(core.CC.P, core.CB.P) #44
P.CB.and.CC <- intersect(core.CB.P, core.CC.P) #2

B.CB.vs.CC <- setdiff(core.CB.B, core.CC.B) #38
B.CC.vs.CB <- setdiff(core.CC.B, core.CB.B) #33
B.CB.and.CC <- intersect(core.CB.B, core.CC.B) #2

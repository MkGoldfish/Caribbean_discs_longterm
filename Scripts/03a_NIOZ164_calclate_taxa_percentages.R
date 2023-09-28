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
# Calculate percentages for taxa                                                          #
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
# Import RA tibble to calculate percentages of top taxa
tt <- read.csv('../Processed-Data/NIOZ164_EUX_discs_RA_tidy_data_decontamed_tax.correct_pruned.csv', row.names = NULL)

# Explore present Archaea and calculate  RAs ------------------------------------------------------------------
Arch <- tt %>% filter(Kingdom == "Archaea")

phyls <- unique(Arch$Phylum)
ords <- unique(Arch$Order)

Arch.pel <- Arch %>% filter(Habitat == "Pelagic") 
Arch.pel.CB <- Arch %>% filter(Location_Habitat == "Charles_Brown_Pelagic") 
Arch.pel.CC <- Arch %>% filter(Location_Habitat == "Crooks_Castle_Pelagic") 

Arch.bent <- Arch %>% filter(Habitat == "Benthic")
Arch.ben.CB <- Arch %>% filter(Location_Habitat == "Charles_Brown_Benthic") 
Arch.ben.CC <- Arch %>% filter(Location_Habitat == "Crooks_Castle_Benthic") 

Arch.wild <- Arch %>% filter(Habitat == "Beach")

## Overall dataset --------------------------------
# Average Archaeal reads per sample 
Arch.sample.avg <- Arch%>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

Arch.sample.sd <- Arch %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Arch.sample.avg
Arch.sample.sd

# Average reads per phylum to determine order of phyla
Arch.phyl.avg <- Arch %>% select(Description, Phylum, Phylum_rel_abund_Sample) %>% group_by(Phylum) %>% 
  summarise(avg = 100*sum(Phylum_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))

# Avreage reads per Order to determine dominant Orders
Arch.ord.avg <- Arch %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>% 
  summarise(avg = 100*(sum(Order_rel_abund_Sample)/length(Arch.pel$Description))) %>% arrange(desc(avg))


## Pelagic Archaea --------------------------------
# Average Archaeal reads per sample in pelagic
Arch.pel.sample.avg <- Arch.pel %>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

Arch.pel.sample.sd <- Arch.pel %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Arch.pel.sample.avg
Arch.pel.sample.sd

# 
# # Average reads per phylum to determine order of phyla 
# Arch.pel.phyl.avg <- Arch.pel %>% select(Description, Phylum, Phylum_rel_abund_Sample) %>% group_by(Phylum) %>% 
#   summarise(avg = 100*sum(Phylum_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))
# 
# # Average reads per Order to determine dominant Orders
# Arch.pel.phyl.avg <- Arch.pel %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>% 
#   summarise(avg = 100*sum(Order_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))

# Not really since the order is the same 

## Pelagic CB

## Pelagic CC

## Benthic Archaea --------------------------------
# Average Archaeal reads per sample in benthic

## Benthic CB

## Benthic CC
## Beach Archaea --------------------------------
# Average Archaeal reads per sample in beach
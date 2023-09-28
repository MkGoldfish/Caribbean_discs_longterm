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
tt$X <- NULL

# Explore present Archaea and calculate  RAs ------------------------------------------------------------------
Arch <- tt %>% filter(Kingdom == "Archaea")

phyls <- unique(Arch$Phylum)
ords <- unique(Arch$Order)

Arch.pel <- Arch %>% filter(Habitat == "Pelagic") 
Arch.pel.CB <- Arch %>% filter(Location_Habitat == "Charles_Brown_Pelagic") 
Arch.pel.CC <- Arch %>% filter(Location_Habitat == "Crooks_Castle_Pelagic") 

Arch.bent <- Arch %>% filter(Habitat == "Benthic")
Arch.bent.CB <- Arch %>% filter(Location_Habitat == "Charles_Brown_Benthic") 
Arch.bent.CC <- Arch %>% filter(Location_Habitat == "Crooks_Castle_Benthic") 

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
  summarize(avg = 100*sum(mean(Phylum_rel_abund_Sample)))  %>% arrange(desc(avg))

Arch.phyl.avg.sd <- Arch %>% select(Description, Phylum, Phylum_rel_abund_Sample) %>% group_by(Phylum) %>% 
 mutate(avg = 100*sum(mean(Phylum_rel_abund_Sample))) %>% 
  mutate(sd = sd(Phylum_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

Arch.phyl.avg.sd

# Avreage reads per Order to determine dominant Orders
Arch.ord.avg <- Arch %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>% 
  summarize(avg = 100*sum(mean(Order_rel_abund_Sample))) %>% arrange(desc(avg))

Arch.ord.avg <- Arch %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>%
  mutate(avg = 100*sum(mean(Order_rel_abund_Sample))) %>% 
  mutate(sd = sd(Order_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

Arch.ord.avg

## Pelagic Archaea --------------------------------
# Average Archaeal reads per sample in pelagic
Arch.pel.sample.avg <- Arch.pel %>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)*100) 

Arch.pel.sample.sd <- Arch.pel %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Arch.pel.sample.avg
Arch.pel.sample.sd


# Average reads per phylum to determine order of phyla
## Same order as overall dataset
Arch.phyl.pel.avg.sd <- Arch.pel %>% select(Description, Phylum, Phylum_rel_abund_Sample) %>% group_by(Phylum) %>% 
  mutate(avg = 100*sum(mean(Phylum_rel_abund_Sample))) %>% 
  mutate(sd = sd(Phylum_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

# Average reads per Order to determine dominant Orders
Arch.pel.ord.avg <- Arch.pel %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>%
  summarise(avg = 100*sum(Order_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))

Arch.pel.ord.avg.sd <- Arch.pel %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>%
  mutate(avg = 100*sum(mean(Order_rel_abund_Sample))) %>% 
  mutate(sd = sd(Order_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

### Pelagic CB------------
Arch.pel.CB.sample.avg <- Arch.pel.CB %>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

Arch.pel.CB.sample.sd <- Arch.pel.CB %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Arch.pel.CB.sample.avg
Arch.pel.CB.sample.sd


# Average reads per phylum to determine order of phyla
Arch.phyl.pel.CB.avg.sd <- Arch.pel.CB %>% select(Description, Phylum, Phylum_rel_abund_Sample) %>% group_by(Phylum) %>% 
  mutate(avg = 100*sum(mean(Phylum_rel_abund_Sample))) %>% 
  mutate(sd = sd(Phylum_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

# Average reads per Order to determine dominant Orders
Arch.pel.CB.ord.avg <- Arch.pel.CB %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>%
  mutate(avg = 100*sum(mean(Order_rel_abund_Sample))) %>% 
  mutate(sd = sd(Order_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))


Nano <- Arch.pel.CB  %>% filter(Phylum ==  "Halobacterota") 
unique(Nano$Order)


### Pelagic CC -------------------------------
Arch.pel.CC.sample.avg <- Arch.pel.CC %>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

Arch.pel.CC.sample.sd <- Arch.pel.CC %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Arch.pel.CC.sample.avg
Arch.pel.CC.sample.sd


# Average reads per phylum to determine order of phyla
Arch.pel.CC.phyl.avg <- Arch.pel.CC %>% select(Description, Phylum, Phylum_rel_abund_Sample) %>% group_by(Phylum) %>%
  summarise(avg = 100*sum(Phylum_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))

# Average reads per Order to determine dominant Orders
Arch.pel.CC.ord.avg <- Arch.pel.CC %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>%
  summarise(avg = 100*sum(Order_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))

## Benthic Archaea --------------------------------
# Average Archaeal reads per sample in benthic
Arch.bent.sample.avg <- Arch.bent %>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

Arch.bent.sample.sd <- Arch.bent %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Arch.bent.sample.avg
Arch.bent.sample.sd


# Average reads per phylum to determine order of phyla
Arch.phyl.bent.avg.sd <- Arch.bent %>% select(Description, Phylum, Phylum_rel_abund_Sample) %>% group_by(Phylum) %>% 
  mutate(avg = 100*sum(mean(Phylum_rel_abund_Sample))) %>% 
  mutate(sd = sd(Phylum_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

# Average reads per Order to determine dominant Orders
Arch.bent.ord.avg <- Arch.bent %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>%
  summarise(avg = 100*sum(Order_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))

### Benthic CB -------------------------------
Arch.bent.CB.sample.avg <- Arch.bent.CB %>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)*100) 

Arch.bent.CB.sample.sd <- Arch.bent.CB %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Arch.bent.CB.sample.avg
Arch.bent.CB.sample.sd


# Average reads per phylum to determine order of phyla
Arch.bent.CB.phyl.avg <- Arch.bent.CB %>% select(Description, Phylum, Phylum_rel_abund_Sample) %>% group_by(Phylum) %>%
  summarise(avg = 100*sum(Phylum_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))

# Average reads per Order to determine dominant Orders
Arch.bent.CB.ord.avg <- Arch.bent.CB %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>%
  summarise(avg = 100*sum(Order_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))


### Benthic CC ------------------
Arch.bent.CC.sample.avg <- Arch.bent.CC %>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

Arch.bent.CC.sample.sd <- Arch.bent.CC %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Arch.bent.CC.sample.avg
Arch.bent.CC.sample.sd


# Average reads per phylum to determine order of phyla
Arch.bent.CC.phyl.avg <- Arch.bent.CC %>% select(Description, Phylum, Phylum_rel_abund_Sample) %>% group_by(Phylum) %>%
  summarise(avg = 100*sum(Phylum_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))

# Average reads per Order to determine dominant Orders
Arch.bent.CC.ord.avg <- Arch.bent.CC %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>%
  summarise(avg = 100*sum(Order_rel_abund_Sample)/length(Arch.pel$Description)) %>% arrange(desc(avg))


## Beach Archaea --------------------------------
# Average Archaeal reads per sample in beach
Arch.wild.sample.avg <- Arch.wild %>% select(Description, Sample_rel_abund)  %>% group_by(Description) %>% 
  summarise(sum = sum(Sample_rel_abund)) %>% summarise(mean= mean(sum)) *100

Arch.wild.sample.sd <- Arch.wild %>% select(Description, Sample_rel_abund) %>% group_by(Description) %>% 
  summarise(sd = sd(Sample_rel_abund)) %>% summarise(sum_sd = sqrt(sum((sd)^2))) *100

Arch.wild.sample.avg
Arch.wild.sample.sd

Arch.phyl.wild.avg.sd <- Arch.wild %>% select(Description, Phylum, Phylum_rel_abund_Sample) %>% group_by(Phylum) %>% 
  mutate(avg = 100*sum(mean(Phylum_rel_abund_Sample))) %>% 
  mutate(sd = sd(Phylum_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

Arch.ord.avg <- Arch.wild %>% select(Description, Order, Order_rel_abund_Sample) %>% group_by(Order) %>%
  mutate(avg = 100*sum(mean(Order_rel_abund_Sample))) %>% 
  mutate(sd = sd(Order_rel_abund_Sample)*100) %>% summarize_at(vars(avg,sd), mean) %>% arrange(desc(avg))

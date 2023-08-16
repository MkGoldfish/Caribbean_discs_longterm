#%############################################################################%#
#%                                                                          ##%#
# Title: Filter, correct and tidy the data                                  ####
# Project: NIOZ164 Discs EUX                                                ####
# author: Maaike Goudriaan, MMB, NIOZ                                       ##%#    
# date: 2023-08-09                                                          ##%#
# Purpose:  Apply filters to the data, correct weird taxonomy               ##%#
#           Transform data to RA and CLR, store as tidy tibble              ##%#
#       ##%#
#%############################################################################%#
  

# Prepare workspace -------------------------------------------------------
## Set working Directory
getwd()
setwd()
print(dir()) 

## Load relevant libraries
library(devtools)
library(phyloseq)
library(grid)
library(tidyverse)
library(vegan)
library(rmarkdown)
library(knitr)
library(microbiome)


# Load data ---------------------------------------------------------------
tax <- as.matrix(read.delim('../Data/representative_seq_set_tax_assignments.txt', 
                            row.names = 1, na.strings = c(" ")))
tax <- tax_table(tax)
otu <- read.delim('../Data/asvTable.txt', row.names = 1)
otu <- otu %>%  select(-taxonomy) %>% as.matrix() 
otu <- otu_table(otu, taxa_are_rows = T)
map <- sample_data(read.delim('../Data/Metadata_EUX_Discs.txt', 
                              row.names = 1, na.strings = c("")))

# Add extra combi columns to the map metadata matirx
Location_Habitat <- str_c(map$Location, "_", map$Habitat)
map <- map %>%
  add_column(Location_Habitat)

# Create and check phyloseq -----------------------------------------------
physeq_object = merge_phyloseq(otu, tax, map)
summarize_phyloseq(physeq_object) 

source("basic_info_physeq_object.R")
basic_info_physeq_object(physeq_object)

# Check and correct taxonomy
# First, pruning to get rid of ASVs with less then 20 reads in total, 
# After, scheck if this leaves samples with only 1 ASV and remove these
# Wonky_tonky_taxonomy removes ASVs unassigned on Kingdom and Phylum level, and removes Chloroplast and mitochondria.
# After, we apply a loop on the taxonomy table to redefine "weird" taxonomy into a single communstring unassigned

physeq_pruned<- prune_taxa(taxa_sums(physeq_object) > 20, physeq_object) 
physeq_pruned <- prune_samples(sample_sums(physeq_pruned) >1, physeq_object)
summarize_phyloseq(physeq_pruned)

source("wonky_tonky_taxonomy.R")
wonky_tonky_taxonomy(physeq_pruned)

# extract tax table as dataframe
taxo <- as.data.frame(physeq_pruned@tax_table)

# Loooop
for (i in 1:nrow(taxo)) {
  for (y in 1:ncol(taxo)) {
    if 
    (any(str_detect(taxo[i,y], c("uncultured","Uncultured","metagenome","Metagenome","unknown","Unknown","NA")))) {taxo[i,y] <- "unassigned" }
  }
} 

# re-define as tax table object
taxo <- tax_table(as.matrix(taxo))
# merge updated taxonomy
physeq_filtered <- merge_phyloseq(physeq_object@otu_table, taxo, map) 


## check physeq object

# check the physeq object
source("basic_info_physeq_object.R")
basic_info_physeq_object(physeq_filtered)
summarize_phyloseq(physeq_filtered)


## Transform data and store as a tidy tibble
# Data is transformed to RA for plotting of commuity, and hellinger and rclr transformed for statistical analysis
source("tidy_psmelt.R")
physeq_filtered_ra <- transform_sample_counts(physeq_filtered,function(x) x / sum(x) )
tidy_RA_filtered <- tidy_tibble_maker(physeq_filtered_ra)
write.csv(tidy_RA_filtered, '../Processed-Data/NIOZ165_EUX_discs_RA_data_filtered.csv')

physeq_filtered_hell <- microbiome::transform(physeq_filtered, 'hellinger' )
tidy_hell_filtered <- tidy_tibble_maker(physeq_filtered_hell)
write.csv(tidy_hell_filtered, '../Processed-Data/NIOZ165_EUX_discs_sqrt_RA_data_filtered.csv')

physeq_filtered_rclr <- microbiome::transform(physeq_filtered, 'rclr' )
tidy_rclr_filtered <- tidy_tibble_maker(physeq_filtered_rclr)
write.csv(tidy_rclr_filtered, '../Processed-Data/NIOZ165_EUX_discs_rclr_data_filtered.csv')




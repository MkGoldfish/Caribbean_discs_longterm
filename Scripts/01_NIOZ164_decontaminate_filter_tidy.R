##%#####################################################################################%##
#NIOZ164 Statia Discs - Decontamination and filtering                                 #####                                            
#Author: Maaike Goudriaan, NIOZ, MMB                                                      #
#                                                                                         #
# Purpose:                                                                                # 
# NIOZ164 is an illumina sequencing lane consisting of different samplesets.              # 
# Here we focus on analysis of 16S amplicon data, from amplified DNA                      #
# extracted from discs covered in polymer films that were incubated at 2 different        #
# locations in 2 different waterdepths (seafloor and watercolumn) in                      #
# coastal Caribbean waters close to the island of St. EUstatius.                          #
# The tested polymers were PE, PP, PS, PET and Nylon, PE-13C and PP-13C.                  #
# Of each polymere there was an UV pretreated and a non-treated version.                  #
#                                                                                         #
# In this script, we use the package microDecon to decontaminate the samples based on the #
# PCR blanks, filter/trim/prune the data with phyloseq, and turn into tidy format         # 
##%#####################################################################################%##

# Date: 2023 - 08 - 29
# R-version: 4.3.1 

# Set working directory ------------------------------------------------------------------
setwd("C:/Users/mgoudriaan/Documents/GitHub/Caribbean_discs_longterm/Scripts")

# Load libraries -------------------------------------------------------------------------
library(tidyverse)
library("microDecon")
library(microbiome)
library(phyloseq)
library(stringr)

# Import data ----------------------------------------------------------------------------
asv <- read.delim('../Data/asvTable.txt', row.names = 1)
# asv.tab <- asv %>%  select(-taxonomy) %>% as.matrix()
# asv.tab <- otu_table(asv.tab, taxa_are_rows = T)
# tax <- as.matrix(read.delim('../Data/representative_seq_set_tax_assignments.txt', row.names = 1, na.strings = c(" ")))
# tax <- tax_table(tax)

# Add extra combi columns to the map metadata matrix
map <- read.delim('../Data/Metadata_EUX_Discs.txt', row.names = 1, na.strings = c(""))
map$InputFileName <-NULL
map$Location_Habitat <- str_c(map$Location, "_", map$Habitat)
map$Polymer_Treatment <- str_c(map$Polymer, "_", map$Treatment)
map1 <- map %>%  mutate(Polymer_Isotope = if_else(Polymer %in% c("PE", "PP"), paste(Polymer, Isotope, sep = "-"), Polymer))
head(map1)
map1 <- sample_data(map1) 

# Decontamination of OTU table w microDecon -------------------------------------------------
## Prep table -----
## Find samplenames for negative controls
meta <- map %>% rownames_to_column("Sample.ID")
nc <- map %>% filter( Location =="NC")

# Determine the number of individuals per subgroup
# microDecon appears to work best when grouping individuals
# This is the order in wich samples are numbered and found in both meta and asv table
CC_pelagic <- meta %>% filter(Location_Habitat == "Crooks_Castle_Pelagic") %>% pull(Sample.ID)
CC_benthic <- meta %>% filter(Location_Habitat == "Crooks_Castle_Benthic") %>% pull(Sample.ID)
CB_pelagic <- meta %>% filter(Location_Habitat == "Charles_Brown_Pelagic") %>% pull(Sample.ID)
CB_benthic <- meta %>% filter(Location_Habitat == "Charles_Brown_Benthic") %>% pull(Sample.ID)
Filters <- meta %>% filter(Phase == "Filter") %>% pull(Sample.ID)
Wild <- meta %>% filter(Location_Habitat == "Zeelandia_Beach") %>% pull(Sample.ID)

Numb.Ind <- c(15, 14, 15, 15, 18, 9)

## We learned from barplot NTC1 has a much higher read abundace than the other 2
## I don't trust this one and will not use it 
nc.v <- c("NIOZ164.126", "NIOZ164.127")

# Negatitve controls need to be in the first two columns
asv.reor <- asv %>% select(!contains(c("NIOZ156","NIOZ164.8","NIOZ164.9"))) %>%          #Remove the 2nd library from the CASCABEL run
  select(-NIOZ164.128) %>%                                    #Remove the NC we will not use for decontamination
  relocate(all_of(nc.v)) %>% rownames_to_column("OTU.ID")     #Move NCs to the first 2 columns

colnames(asv.reor)

## Decontamination -----
# ## This runs the whole package at once
# asv.decon.complete <- decon(asv.reor, numb.blanks = 2, numb.ind = c(rep_len(1,86)), taxa= F, thresh = 0.08, prop.thresh= 1e-5)

# Remove decontaminations from the samples based on two negative controls with microDecon
# Samples are treated as individuals
asv.decon.1 <- remove.cont(asv.reor, numb.blanks = 2,  taxa= T)
ncol(asv.decon.1)

# Remove residual contamination from the output from the previous steps
## Since all samples are so different and we have no replicates, 
# it does not feel right to threat samples as group per location/habitat, so we use them as individuals. 
asv.decon.2 <- remove.thresh(asv.decon.1, taxa = T, numb.ind = Numb.Ind, thresh = 0.08, prop.thresh= 1e-5)

#Create a summary to results of decontamination
decon.summary <- decon.diff(asv.reor, asv.decon.2, numb.blanks = 2, numb.ind = Numb.Ind, taxa = T)

# Acces the 4 different elements of the summary, presenting the results of deconminating
# And store the results
Reads.removed <- decon.summary$reads.removed
Diff.Sum <- decon.summary$sum.per.group
Mean <- decon.summary$mean.per.group
OTU.removed <- decon.summary$OTUs.removed

write.table(Reads.removed, "../Analysis/microDecon/NIOZ164_reads_removed.txt")
write.table(OTU.removed, "../Analysis/microDecon/NIOZ164_ASVs_removed.txt")
write.table(Diff.Sum, "../Analysis/microDecon/NIOZ164_removed_reads_per_group_sum.txt")
write.table(Mean, "../Analysis/microDecon/NIOZ164_removed_reads_per_group_avg.txt")

# After filtering, curation with phyloseq package ----------------------------------------------
## Correct taxonomy table ----------------------------------------------------------------------
# First, extract tax table as dataframe from the decontaminated asv table
# It is there present as 1 string column, were the different taxonomic levels are seperated by  ";"
# This column needs to be split for phyloseq
tax.decon <- asv.decon.2 %>% select(OTU.ID, taxonomy) %>% column_to_rownames("OTU.ID") %>%  
separate(taxonomy, into = c( "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
tax.decon <- tax.decon  %>% tax_table()
head(tax.decon)

## Create new physeq object ----------------
# define the decontaminated asv table as otu.table
# First we need to remove rownames, and the columns of blanks mean and taxonomy, and set OTU.IDs as rownames
rownames(asv.decon.2) <- NULL
asv.decon <- asv.decon.2 %>% select(-Mean.blank, -taxonomy) %>% column_to_rownames("OTU.ID")
asv.decon <- otu_table(asv.decon, taxa_are_rows = T)

# Merge all data into phyloseq
physeq_decon <- merge_phyloseq(asv.decon, tax.decon, map1) 
physeq_object <- physeq_decon

## Store Physeq Decontamed object --------------------
saveRDS(physeq_decon, "../Analysis/NIOZ164_physeq_object_decontamed.rds")

## Correct Physeq Taxonomy -----------------------------
# Wonky_tonky_taxonomy removes ASVs unassigned on Kingdom and Phylum level, and removes Chloroplast and mitochondria.
# After, we apply a loop on the taxonomy table to redefine "weird" taxonomy into a single communstring unassigned
# source("wonky_tonky_taxonomy.R")
# wonky_tonky_taxonomy(physeq_object)

get_taxa_unique(physeq_object, "Kingdom") # unassigned in Kingdom
physeq_object <- subset_taxa(physeq_object, !is.na(Kingdom) & !Kingdom%in% c(" ", "Unassigned", "NA")) #let's eliminate those otus
get_taxa_unique(physeq_object, "Kingdom") # all good now

# get_taxa_unique(physeq_object, "Phylum") # let's check the Phyla, there's "NA"
length(get_taxa_unique(physeq_object,"Phylum"))  
physeq_object <- subset_taxa(physeq_object, !is.na(Phylum) & !Phylum%in% c("NA", " ")) 
get_taxa_unique(physeq_object, "Phylum")
length(get_taxa_unique(physeq_object,"Phylum")) 

length(get_taxa_unique(physeq_object,"Order"))
physeq_object <- subset_taxa(physeq_object, !Order%in% c(" Chloroplast", "Chloroplast", "chloroplast", " chloroplast"))
length(get_taxa_unique(physeq_object,"Order"))

length(get_taxa_unique(physeq_object,"Family"))
physeq_object <- subset_taxa(physeq_object, !Family%in% c("Mitochondria", " Mitochondria"))
length(get_taxa_unique(physeq_object,"Family"))

physeq_decon <- physeq_object

# extract tax table as dataframe
taxo <- as.data.frame(physeq_decon@tax_table)

# Loooop
for (i in 1:nrow(taxo)) {
  for (y in 1:ncol(taxo)) {
    if 
    (any(str_detect(taxo[i,y], c("uncultured","Uncultured","metagenome","Metagenome","unknown","Unknown","NA")))) {taxo[i,y] <- "unassigned" }
  }
} 

# re-define as tax table object
taxo <- tax_table(as.matrix(taxo))

# Merge corrected taxonomy with other elements to create phyloseq object -------------------------------------
physeq_decon.1 <- merge_phyloseq(asv.decon, taxo, map1) 
get_taxa_unique(physeq_decon.1, "Kingdom")

summarize_phyloseq(physeq_decon.1)
#"1] Min. number of reads = 7301"
# "2] Max. number of reads = 539188"
#  "3] Total number of reads = 10279736"
#  "4] Average number of reads = 119531.813953488"
# "5] Median number of reads = 89666.5"
# "7] Sparsity = 0.980317379653689"
#  "6] Any OTU sum to 1 or less? YES"
#  "8] Number of singletons = 32911"

source("basic_info_physeq_object.R")
basic_info_physeq_object(physeq_decon.1)
# "Total taxa is 71573"
# "Total samples is 86"
#  "Lowest readnumber is 7301"
#  "Highest readnumber is 539188"
# "Lowest taxa sum is 0"
# "Highest taxa sum is 225095"

# Prune Physeq-object --------------------
# Remove all samples with less than 10 reads
physeq_pruned <- prune_samples(sample_sums(physeq_decon.1) > 10, physeq_decon)
# Remove all single singletons
physeq_pruned <- prune_taxa(taxa_sums(physeq_pruned) > 1, physeq_decon) 
summarize_phyloseq(physeq_pruned)
basic_info_physeq_object(physeq_pruned)

#"1] Min. number of reads = 7301"
# "2] Max. number of reads = 539188"
#  "3] Total number of reads = 10277947"
# [1] "4] Average number of reads = 119511.011627907""
# [1] "5] Median number of reads = 89630"
# [1] "7] Sparsity = 0.964100619200633"
# [1] "6] Any OTU sum to 1 or less? NO"
# [1] "8] Number of singletons = 0"

# [1] "Total taxa is Total taxa is 38662"
# [1] "Total samples is 86"
# "Lowest readnumber is 7301"
# [1] "Highest readnumber is  539188"
# [1] "Lowest taxa sum is 2"
# [1] "Highest taxa sum is 225095"


## Visualizing the read abundance per sample plot -----
sample_sums <- sample_sums(physeq_pruned)
sample_sums_df <- data.frame(sum = sample_sums) %>% rownames_to_column("Sample_ID") %>% arrange(Sample_ID)
Sample_description <- data.frame(physeq_pruned@sam_data) %>% rownames_to_column("Sample_ID") %>% arrange(Sample_ID) %>%  filter(Location != "NC")
sample_sums_df$Description <- Sample_description$Description
ggplot(sample_sums_df, aes(x = Description, y = sum)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Sum of Counts", title = "Sample Sums after decontamination and pruning")


## Store Physeq object --------------------
saveRDS(physeq_pruned, "../Analysis/NIOZ164_physeq_object_decontamed_filtered.rds")

# Create and store tidy tibble 3 domains ------------
source("tidy_tibble_maker.R")
tidy_decont_pruned <- tidy_tibble_maker(physeq_pruned)

## Calculate the relative abundance per Taxonomic level for each sample
tidy_decont_RA <- tidy_decont_pruned  %>% group_by(Description) %>% mutate(Sample_rel_abund = Abundance / sum(Abundance)) %>% #relative abundance of each otu per sample
  mutate(Sample_st_dev = sd(Sample_rel_abund)) %>% 
  ungroup() %>% 
  
  #_Kingdom_section__#
  group_by(Description, Kingdom) %>% 
  mutate(Kingdom_rel_abund_Sample = sum(Sample_rel_abund)) %>%  #Kingdom relative abundance per sample 
  ungroup() %>% 
  
  #_Phylum_section__#
  group_by(Description, Phylum) %>% 
  mutate(Phylum_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  
  #_Class_section__#
  group_by(Description, Class) %>% 
  mutate(Class_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  
  #_Order_section__#
  group_by(Description, Order) %>% 
  mutate(Order_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  
  #_Family_section__#
  group_by(Description, Family) %>% 
  mutate(Family_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  
  #_Genus_section__#
  group_by(Description, Genus) %>% 
  mutate(Genus_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() 

head(tidy_decont_RA)
write.csv(tidy_decont_RA, '../Processed-Data/NIOZ164_EUX_discs_Proks_tidy_data_decontamed_pruned_RA.csv')

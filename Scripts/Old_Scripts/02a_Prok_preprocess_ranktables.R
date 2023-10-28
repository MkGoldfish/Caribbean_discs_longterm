##%######################################################%##
#                                                          #
####              Exploratory Analysis of               ####
####            16S Amplicon Sequencing Data            ####
#                                                          #
##%######################################################%##

##%######################################################%##
#                                                          #
####        Project: OpenBio Prokaryotic data           ####
####        Written by: Maaike Goudriaan July 2022      ####
####                                                    ####
#                                                          #
##%######################################################%##
####_______________________________________________________________________________________####
####                 Working Directory                 
####_______________________________________________________________________________________####
setwd('C:/Users/mgoudriaan/Documents/R-files/Projects/OpenBio/Scripts')


####_______________________________________________________________________________________####
####                   Load libraries                                                      
####_______________________________________________________________________________________####
library(devtools)
library(phyloseq)
library(grid)
library(tidyverse)
library(vegan)
library(rmarkdown)
library(knitr)
library(patchwork)
library(microbiome)
library("ggsci")

####_______________________________________________________________________________________####
####                  Import Data for phyloseq object                
####_______________________________________________________________________________________####
# Transform data into 3 needed obkects for physeq
tax <- as.matrix(read.delim('../Data/representative_seq_set_tax_assignments.txt', row.names = 1, na.strings = c(" ")))
tax <- tax_table(tax)
otu <- read.delim('../Data/20220225_PE_asvTable_noSingletons_Decon_Filt.txt', row.names = 1)
otu <- otu %>%  select(-taxonomy) %>% as.matrix() 
otu <- otu_table(otu, taxa_are_rows = T)
# map <- read.delim('../Data/OpenBio_Metadata_20230202.txt', row.names = 1, na.strings = c(" "))
map <- read.delim('../Data/OpenBio_Sample_Metadata3.txt', row.names = 1, na.strings = c(" "))
map$X_1 <-NULL
map <- sample_data(map)
colnames(map)
head(map)
head(tax)


####_______________________________________________________________________________________####
####           Create & check physeq object            
####_______________________________________________________________________________________####
physeq_object = merge_phyloseq(otu, tax, map) 
summarize_phyloseq(physeq_object)

# check the physeq object
source("basic_info_physeq_object.R")
basic_info_physeq_object(physeq_object)

Phyla <- phyloseq::tax_glom(physeq_object, taxrank="Genus", NArm = FALSE)
length(get_taxa_unique(Phyla, "Genus"))
Arch <- subset_taxa(Phyla, Kingdom %in% c("Archaea")) 
ntaxa(Phyla, )
ntaxa(Arch)
get_taxa_unique(Arch, "Genus") %>% length()
summarize_phyloseq(Arch)
mean(taxa_sums(Arch))
sd(taxa_sums(Arch))
basic_info_physeq_object(Arch)

####_______________________________________________________________________________________####
####      Check & correct taxonomy assignments         
####_______________________________________________________________________________________####
# Removes ASVs unassigned on Kingdom and Phylum level, and removes Chloroplast and mitochondria in case they are still there
source("wonky_tonky_taxonomy.R")
wonky_tonky_taxonomy(physeq_object)

# Pruning to get rid of ASVs that are only present ones, and samples with only 1 ASV
physeq_object <- prune_taxa(taxa_sums(physeq_object) > 1, physeq_object) #no singletons based on taxasums
physeq_object <- prune_samples(sample_sums(physeq_object) >1, physeq_object) # Remove Singletons from samples

summarize_phyloseq(physeq_object)


# transform into dataframe
taxo <- as.data.frame(physeq_object@tax_table)
taxo$Family %>% unique() %>% sort()
taxo$Order %>% unique() %>% sort()

####_Loop to redefine weird taxonomy to a single common character string "unassigned"_____####
for (i in 1:nrow(taxo)) {
  for (y in 1:ncol(taxo)) {
    if 
    (any(str_detect(taxo[i,y], c("uncultured","Uncultured","metagenome","Metagenome","unknown","Unknown","NA")))) {taxo[i,y] <- "unassigned" }
  }
} 

# For some reason, CASCABEL gave us Alcanivoracaceae and Alcanivoracaceae1, remove the latter
taxo <- taxo %>%  mutate(Family = ifelse(Family == "Alcanivoracaceae1", "Alcanivoracaceae", Family))


# re-define as tax table object
taxo <- tax_table(as.matrix(taxo))

# merge updated taxonomy
physeq_object <- merge_phyloseq(physeq_object@otu_table, taxo, map) 


# check the physeq object
source("basic_info_physeq_object.R")
basic_info_physeq_object(physeq_object)
summarize_phyloseq(physeq_object)

####_______________________________________________________________________________________####
####      Merge samples per surface all replicates together         
####_______________________________________________________________________________________####
# First create a tidy tibble
source("tidy_tibble_maker.R")
tidy_physeq <- tidy_tibble_maker(physeq_object)

# Get rid of columns we no longer need for plotting
tidy_physeq$LinkerPrimerSequence <- NULL
tidy_physeq$ReversePrimerSequence <- NULL
tidy_physeq$InputFileName <- NULL
tidy_physeq$BarcodeSequence <- NULL
tidy_physeq$BarcodeSequence_1 <- NULL


colnames(tidy_physeq)

# Optionally, add extra column with genus+species combined, to deal with 
# tidy_physeq$g_species <- str_c(tidy_physeq$Genus, " ", tidy_physeq$Species)

# source("merge_sample_replicates.R")
# merge_sample_replicates(tidy_physeq, "TixPolxHa")
# ddply(tidy_physeq, merge_sample_replicates())# 
# lapply(tidy_physeq, merge_sample_replicates,  rep_group)
# c("TixPolxHa") %>%  merge_sample_replicates(tidy_physeq, rep_group = "TixPolxHa") # 
# merge_sample_replicates(tidy_physeq, "TixPolxHa") # 
# apply(tidy_physeq, "TixPolxHa", merge_sample_replicates)


####_Calculate relative abundance per grouped replicates_________________________________#### 

## First choose column you want to use for grouping, provide as string 
# SelectColName <- c("TixPolxHa")

## Function to calculate rel.ab. per replicates 
# Custom function to deal with the fact that we do not always have 3 replicates. Takes sums and amount of present replicates

tidy_grouped_replis <- tidy_physeq  %>% group_by(Description) %>% mutate(Sample_rel_abund = Abundance / sum(Abundance)) %>% #relative abundance of each otu per sample
  mutate(Sample_st_dev = sd(Sample_rel_abund)) %>% 
  ungroup() %>% 
  group_by(TixPolxHa) %>% mutate( rep_rel_abund = Sample_rel_abund/ sum(Sample_rel_abund)) %>% #relative abundance of each otu per number of samples in replicates
  mutate(Rep_st_dev = sd(Sample_rel_abund)) %>% 
  ungroup() %>%
  
  #_Kingdom_section__#
  group_by(Description, Kingdom) %>% 
  mutate(Kingdom_rel_abund_Sample = sum(Sample_rel_abund)) %>%  #Kingdom relative abundance per sample 
  ungroup() %>% 
  group_by(TixPolxHa, Kingdom) %>% 
  mutate(st_dev_Kingdom_abund = sd(Kingdom_rel_abund_Sample)) %>% # standard dev of Kingdom relative abundances between replicates of TixPolxHa (ployner_timepoint_treatment)
  mutate(Kingdom_rep_rel_abund = sum(rep_rel_abund)) %>% #Kingdom relative abundance per samples of desc 
  ungroup() %>%
  
  #_Phylum_section__#
  group_by(Description, Phylum) %>% 
  mutate(Phylum_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(TixPolxHa, Phylum) %>% 
  mutate(st_dev_Phylum_abund = sd(Phylum_rel_abund_Sample)) %>%
  mutate(Phylum_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  
  #_Class_section__#
  group_by(Description, Class) %>% 
  mutate(Class_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(TixPolxHa, Class) %>% 
  mutate(st_dev_Class_abund = sd(Class_rel_abund_Sample)) %>%
  mutate(Class_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  
  #_Order_section__#
  group_by(Description, Order) %>% 
  mutate(Order_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(TixPolxHa, Order) %>% 
  mutate(st_dev_Order_abund = sd(Order_rel_abund_Sample)) %>%
  mutate(Order_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  
  #_Family_section__#
  group_by(Description, Family) %>% 
  mutate(Family_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(TixPolxHa, Family) %>% 
  mutate(st_dev_Family_abund = sd(Family_rel_abund_Sample)) %>%
  mutate(Family_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  
  #_Genus_section__#
  group_by(Description, Genus) %>% 
  mutate(Genus_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(TixPolxHa, Genus) %>% 
  mutate(st_dev_Genus_abund = sd(Genus_rel_abund_Sample)) %>%
  mutate(Genus_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() 
  
  # #_Species_Section__#
  # #### Correcting species by combinig genus+species# 
  # #### Easy fix idea merge gen/spec in new col & replace in following paragraph
  # group_by(TixPolxHa, g_species) %>% 
  # mutate(Species_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  # ungroup() %>% 
  # group_by(TixPolxHa, g_species) %>% 
  # mutate(st_dev_Species_abund = sd(Species_rel_abund_Sample)) %>%
  # mutate(Species_rep_rel_abund = sum(rep_rel_abund)) %>%
  # ungroup() 

# Add a column with months to the tibble, with monthnumber as string. 
# this is easier with plotting and analysis
tidy_grouped_replis <- tidy_grouped_replis %>% mutate(Month = case_when(
  Timepoint_cat == 't1' ~ "2.5",
  Timepoint_cat == 't2' ~ "5" ,
  Timepoint_cat == 't3'~ "7.5" ,
  Timepoint_cat == 't4'~ "10" ,
  Timepoint_cat == 't5' ~ "22" ,
)
)

head(tidy_grouped_replis)
colnames(tidy_grouped_replis)

write_csv2(tidy_grouped_replis, "../Analysis/tidy_data_OpenBio_Prokaryotes_202303_complete.csv")  

tidy_grouped_replis_filt <- tidy_grouped_replis %>%  select(OTU, Sample, Description, Abundance,
                                                       Polymer, Timepoint_cat, Month, Replicate,
                                                       Habitat, TixPol, TixHa, PolxHa, TixPolxHa, Desintegration,
                                                        degraded_pct, Kingdom, Phylum, Class, Order, Family, Genus, Species,  
                                                       Sample_rel_abund, Sample_st_dev, rep_rel_abund, Rep_st_dev, Kingdom_rel_abund_Sample, st_dev_Kingdom_abund, 
                                                       Kingdom_rep_rel_abund, Phylum_rel_abund_Sample, st_dev_Phylum_abund, Phylum_rep_rel_abund, Class_rel_abund_Sample,
                                                       st_dev_Class_abund, Class_rep_rel_abund, Order_rel_abund_Sample, st_dev_Order_abund, Order_rep_rel_abund, 
                                                       Family_rel_abund_Sample, st_dev_Family_abund, Family_rep_rel_abund,    
                                                       Genus_rel_abund_Sample, st_dev_Genus_abund, Genus_rep_rel_abund)

write_csv2(tidy_grouped_replis_filt, "../Analysis/tidy_data_OpenBio_Prokaryotes_202303.csv")  

####_______________________________________________________________________________________####
####                  Creating table per taxonomic level 
####_______________________________________________________________________________________####

# Create tables per taxonomic level, to export and for plotting 
#_KINGDOM_#
Kingdom <- tidy_grouped_replis  %>%  select(Timepoint_cat, Polymer, Habitat, Month,
                                            TixPol, TixHa, PolxHa, TixPolxHa, 
                                            Kingdom, Kingdom_rep_rel_abund, st_dev_Kingdom_abund)%>%
  distinct()
head(Kingdom)
write_csv(Kingdom, "../Analysis/Kingdom_Prok_202302.csv")

#_PHYLUM_#
Phyla <- tidy_grouped_replis  %>%  select(Timepoint_cat, Polymer, Habitat, Month,
                                          TixPol, TixHa, PolxHa, TixPolxHa, Phylum, Phylum_rep_rel_abund, st_dev_Phylum_abund)%>%
  distinct() 
head(Phyla)
write_csv(Phyla, "../Analysis/Phyla_Prok_202302.csv")

#_CLASS_#
Classes <- tidy_grouped_replis  %>%  select(Timepoint_cat, Polymer, Habitat,
                                            TixPol, TixHa, PolxHa, TixPolxHa, Month,
                                            Class, Class_rep_rel_abund, st_dev_Class_abund)%>%
  distinct() 
head(Classes)
write_csv(Classes, "../Analysis/Classes_Prok_202302.csv")

#_ORDER_#
Orders <- tidy_grouped_replis  %>%  select(Timepoint_cat, Polymer, Habitat, Month,
                                           TixPol, TixHa, PolxHa, TixPolxHa, 
                                           Order,Order_rep_rel_abund, st_dev_Order_abund)%>%
  distinct() 
head(Orders)
write_csv(Orders, "../Analysis/Orders_Prok_202302.csv")

#_FAMILY_#
Families <- tidy_grouped_replis  %>%  select(Timepoint_cat, Polymer, Habitat,  Month,
                                                   TixPol, TixHa, PolxHa, TixPolxHa,  Family, Family_rep_rel_abund, st_dev_Family_abund)%>%
  distinct() 
head(Families)
write_csv(Families, "../Analysis/Families_Prok_202302.csv")

#_GENUS_#
Genera <- tidy_grouped_replis  %>%   select(Timepoint_cat, Polymer, Habitat,  Month,
                                          TixPol, TixHa, PolxHa, TixPolxHa, 
                                         Genus, Genus_rep_rel_abund, st_dev_Genus_abund)%>%
  distinct() 
head(Genera)
write_csv(Genera, "../Analysis/Genera_Prok_202302.csv")

#_SPECIES_#
Specys <- tidy_grouped_replis  %>%  select(Timepoint_cat, Polymer, Habitat, Season, Month,
                                          TixPol, TixHa, PolxHa, TixPolxHa,
                                           Species, Species_rep_rel_abund, st_dev_Species_abund)%>%
  distinct() 
head(Specys)
write_csv(Specys, "../Analysis/Species_Prok_202302.csv")


#__ASVs__#
ASV <-  tidy_grouped_replis  %>%  select(TixPolxHa, Timepoint_cat, Polymer, Habitat,
                                         TixPol, TixHa, PolxHa, TixPolxHa, 
                                         Temp_avg_month, OTU, Sample_rel_abund, Sample_st_dev)%>% 
  distinct() 
head(ASV)
write_csv(ASV, "../Analysis/ASV_Prok_202302.csv")

#__ASVs_replis__#
ASV_reps<-  tidy_grouped_replis  %>%  select(TixPolxHa,Timepoint_cat, Polymer, Habitat,
                                             TixPol, TixHa, PolxHa, TixPolxHa, 
                                             Temp_avg_month, OTU, rep_rel_abund, Rep_st_dev)%>% 
  distinct() 
head(ASV_reps)
write_csv(ASV_reps, "../Analysis/ASV_reps_Prok_202302.csv")


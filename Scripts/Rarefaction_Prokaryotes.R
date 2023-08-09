##%######################################################%##
#                                                          #
####             Diversity Analysis of               ####
####            16S Amplicon Sequencing Data            ####
#                                                          #
##%######################################################%##

##%######################################################%##
#                                                          #
####           Project: OpenBio_ Prokaryotes            ####
####            Adapted by Maaike April 2022            ####
####              (this is the right one)               ####
####               Corrected Primers                    ####
#                                                          #
##%######################################################%##
############################################################
##  version.string R version 4.2.1 (2022-06-23 ucrt)
##  nickname       Funny-Looking Kid   

##%######################################################%##
#                                                          #
####              Rarefaction                            ####
#                                                          #
##%######################################################%##

###%#_______________________________________________________________________________________#%###
####                 Working Directory                 ####
###%#_______________________________________________________________________________________#%###
setwd('C:/Users/mgoudriaan/Documents/R-files/Projects/OpenBio/Scripts')

####_______________________________________________________________________________________#%###
####                   Load libraries                                                      ####
####_______________________________________________________________________________________#%###
library(tidyverse)
library(vegan)
library(knitr)
library("ggh4x")
library(ggpubr)
library("paletteer")
library("ggthemes")
library("Polychrome")
library("cowplot")
library(phyloseq)
library("ggrepel")
library(microbiome)


##%#########################################################################%##
##%#_______________________________________________________________________#%###
####                  Import Data                                           ####
###%#______________________________________________________________________#%###
# Transform data into 3 needed objects for physeq
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

###_______________________________________________________________________________________####
####           Create & check physeq object            
####_______________________________________________________________________________________####
#create physeq object
physeq_object = merge_phyloseq(otu, tax, map) 
summarize_phyloseq(physeq_object)

# Turn into a tidy tibble
source("tidy_tibble_maker.R")
tidy_physeq <- tidy_tibble_maker(physeq_object)

# Get rid of columns we no longer need for plotting
tidy_physeq$LinkerPrimerSequence <- NULL
tidy_physeq$ReversePrimerSequence <- NULL
tidy_physeq$InputFileName <- NULL
tidy_physeq$BarcodeSequence <- NULL
tidy_physeq$BarcodeSequence_1 <- NULL

colnames(tidy_physeq)

##############################################################################################
###%#______________________________________________________________________________________#%###
####                Rarefaction              
###%#______________________________________________________________________________________#%###
#Select colums we need for the rarefaction
ASV.ABS.count <- tidy_physeq %>%  select(TixPolxHa, Polymer, Habitat, Description, OTU, Abundance) %>%
  distinct() 
head(ASV.ABS.count)

#Adjust the number of significant digets to get from vegan in rarefactions
#back to old settings
old <- options(pillar.sigfig = 8)
options(old)

# Transform count table to get samples per row and ASV per column
# Not perse necessary, but this means we can use standard MARGIN in transfromation and ordination so makes life easier
#We only want 3 columns for rarefaction, because that's all the rarecurve function can take
ASV.ABS_df1 <- ASV.ABS.count %>% select(Description, OTU, Abundance)  %>% mutate(across(c(OTU),factor)) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance, values_fill = 0) %>% replace(is.na(.), 0)  %>% arrange(Description) %>% 
  column_to_rownames(var = "Description") %>% as.data.frame()
str(ASV.ABS_df1)
rownames(ASV.ABS_df1)
colnames(ASV.ABS_df1)

# Select the metadata we need for plotting
meta_df1 <- ASV.ABS.count %>% select(Description, Polymer, Habitat, OTU, Abundance)  %>% mutate(across(c(OTU),factor)) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance, values_fill = 0) %>% replace(is.na(.), 0)  %>% arrange(Description) %>% 
  column_to_rownames(var = "Description") %>% as.data.frame()%>% select(Polymer, Habitat)
rownames(meta_df1)
colnames(meta_df1)

# # Find the lowest amount of sequences in all samples with this
# # Can be used for rarefaction depth
# min_seqs_1 <- ASV.ABS.count%>%
#   group_by(Description) %>%
#   summarize(n_seqs = sum(Abundance))  %>% 
#   summarize(min = min(n_seqs)) %>%
#   pull(min)

# #Get a tibble of the number of sequences per sample and put in ascending order
# # To check amounts
# n_seqs_1 <- ASV.ABS.count%>%
#   group_by(Description) %>%
#   summarize(n_seqs = sum(Abundance))  %>% 
#   arrange(n_seqs)
# head(n_seqs_1)
# 
# Description n_seqs
# <chr>        <int>
#   1 5B1P      10916
# 2 5A5P         11684
# 3 4A1P         13992
# 4 4B2P         14547
# 5 5A5B         15094
# 6 5B5P         15662

# > min_seqs_1 
# [1] 10916

# # Get a tibble of rarefied data in case you want/need this
# rare_1 <- rarefy(ASV.ABS_df1, min_seqs_1) %>% 
#   as_tibble(rownames  = "Description") %>% 
#   select(Description, vegan=value)

# #Transform count table to get samples per row and ASV per column
# # In this case we make the table per sample w tripiclates pulled, instead individual samples
# ASV.ABS_df2 <- ASV.ABS.count %>% select(TixPolxHa,  OTU, Abundance)  %>%
#   group_by(TixPolxHa, OTU) %>% summarize(count = as.integer(mean(Abundance))) %>% 
#   mutate(across(c(OTU),factor)) %>%
#   pivot_wider(names_from = OTU, values_from = count, values_fill = 0)  %>% arrange(TixPolxHa) %>% 
#   column_to_rownames(var = "TixPolxHa") %>% as.data.frame()
# str(ASV.ABS_df2)
# 
# # Find the lowest amount of sequences per triplo of samples
# min_seqs_2 <- ASV.ABS.count%>%
#   group_by(TixPolxHa, Description) %>%
#   summarize(n_seqs = sum(Abundance)) %>% ungroup() %>% 
#   group_by(TixPolxHa) %>% 
#   summarize(av_seqs = mean(n_seqs)) %>% 
#   summarize(min = min(av_seqs)) %>%
#   pull(min) %>% as.integer()

# #Get a tibble of the number of sequences per sample and put in ascending order
# n_seqs_2 <- ASV.ABS.count%>%
#   group_by(TixPolxHa, Description) %>%
#   summarize(n_seqs = sum(Abundance)) %>% ungroup() %>% 
#   group_by(TixPolxHa) %>% 
#   summarize(av_seqs = mean(n_seqs)) %>% 
#   arrange(av_seqs)
# 
# head(n_seqs_2)

# The lowest number is very low, we use the 2nd one
# > head(n_seqs_2)
# # A tibble: 6 × 2
# TixPolxHa       av_seqs
# <chr>             <dbl>
#   1 5_PHB_Pelagic    13289 
# 2 5_PBSeT_Pelagic  14803 
# 3 3_PHB_Pelagic    19696 
# 4 3_LDPE_Pelagic   20370.
# 5 4_PBSeT_Pelagic  21422 
# 6 4_PHB_Pelagic    23576.

# # Get a tibble of rarefied data in case you want/need this
# rare_2 <- rarefy(ASV.ABS_df2, 12746) %>% 
#   as_tibble(rownames  = "TixPolxHa") %>% 
#   select(TixPolxHa, vegan=value)

#Remove the first column for rarefaction
ASV.ABS.df1 <- ASV.ABS_df1[,-1]
# #another way to get a tibble with rarefied data
# samples <- vegan::rarefy(ASV.ABS.df1, 10916) %>%
#   as_tibble(rownames = "Sample") %>%
#   select(Sample, samp.val = value)

ASV.ABS.df2 <- ASV.ABS_df2[,-1]
# reps<- rarefy(ASV.ABS_df2, 12746 ) %>%
#   as_tibble(rownames = "Description") %>%
#   select(Description, samp.val = value)

#Create the rarefaction curves for both individual samples, and the replicates combined
samples.rarecurve <- vegan::rarecurve(ASV.ABS.df1, step = 1000) 
reps.rarecurve <- vegan::rarecurve(ASV.ABS.df2, step = 100)
#check structure of the object
str(samples.rarecurve)

#We need to create a dataframe form the rarecurve data, to be able to use it in ggplot2
# map_dfr binds all individual list elements into one df 
# Bind the columns of metadata to the df, for the colors and facets of the plot in this
# change rownames for plotting 
samp.rare.data <- map_dfr(samples.rarecurve, bind_rows) %>% 
  bind_cols(meta_df1, Sample = rownames(ASV.ABS_df1)) %>% 
  pivot_longer(-c("Polymer", "Sample", "Habitat")) %>% 
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>% 
  select(-name)

head(samp.rare.data)
colnames(samp.rare.data)

# We add an extra column with samplenames to be used for the label
# Group by Sample and make sure there is only a label for the last value of the sample
samp.rare.data <- samp.rare.data %>% group_by(Sample) %>% mutate(label = if_else(value == max(value), as.character(Sample), NA_character_))

# This is our general color palette for polymers as we have been using
pal.pol.line <- c('#EE7733','#117733', '#66CCEE' )

#Now we can plot the data with ggplot2 
rareplot.1 <- ggplot(samp.rare.data, aes(x = n_seqs, y = value, group = Sample)) +
  geom_line(aes(color = Polymer)) +
  facet_grid(. ~ Habitat) +
  # geom_vline(xintercept = 10916, color="black", linewidth = 1) +    ## We do not rarefy in this case, so we do need a line for rarefaction depth
  xlim(0, 100000) +   ## We set an x-lim to make it easier to digest the plot
  geom_label_repel(aes(label = label, color = Polymer), nudge_x = 1, na.rm = T) +
  scale_colour_manual(values = pal.pol.line) +
  xlab("# of sequences") + ylab("# of ASVs") +
  theme_classic()

rareplot.1

rep.rare.data <- map_dfr(reps.rarecurve, bind_rows) %>% 
  bind_cols(TixPolxHa = rownames(ASV.ABS_df2), .) %>% 
  pivot_longer(-TixPolxHa) %>% 
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>% 
  select(-name)

head(rep.rare.data)

rareplot.2 <- ggplot(rep.rare.data, aes(x = n_seqs, y = value, color = TixPolxHa)) +
  geom_line( linewidth = 1.5) +
  geom_vline(xintercept = min_seqs_2, linewidth = 1.5) +
  xlim(0, 50000) +
  guides(color= guide_legend(ncol = 2)) +
  theme_classic()

rareplot.2

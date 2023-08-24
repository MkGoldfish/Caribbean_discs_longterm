##%#####################################################################################%##
# NIOZ164 Statia Discs - FA SIP Fingerprint Plots                                     #####                                            
# Author: Maaike Goudriaan, NIOZ, MMB                                                     #
#                                                                                         #
# Purpose:                                                                                # 
# NIOZ164 is an illumina sequencing lane consisting of different samplesets.              # 
#                                                                                         #
# Discs covered in polymer films that were incubated at 2 different                       #
# locations in 2 different waterdepths (seafloor and watercolumn) in                      #
# coastal Caribbean waters close to the island of St. EUstatius.                          #
# The tested polymers were PE, PP, PS, PEt and Nylon, PE-13C and PP-13C.                  #
# Of each polymere there was an UV pretreated and a non-treated version.                  #
# In addition to amplicon analysis, we performed SIP by extracting membrane FA's          #
# and analyzing the extracts with FID, MS and IRMS. This script is used to make plots of  #
# the resulting fingerprints                                                              #
##%#####################################################################################%##

# Date: 2023 - 08 - 24
# R-version: 4.3.1 

# Set working directory ------------------------------------------------------------------
setwd("C:/Users/mgoudriaan/Documents/GitHub/Caribbean_discs_longterm/Scripts")

# Load libraries -------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library("ggh4x")
library(ggpubr)
library("ggsci")
library("cowplot")

# Import data ----------------------------------------------------------------------------
Areas <- read.delim("../Data/Lipid_areas_RA.txt", sep = '\t', dec = ".")
d13C <- read.delim("../Data/Lipid_d13C_corrected_signals.txt", sep = '\t', na.strings = c(" "), dec = ".")
metadata  <- read.delim("../Data/Metadata_lipids_R_plots.txt", sep = '\t')

# Create one tibble for easy plotting ----------------------------------------------------
# Remove the row of the standard and pivot
Areas.long <- Areas %>% filter(Fatty.Acid != "C19:0 (standard!)") %>% pivot_longer( -Fatty.Acid, names_to = "Sample", values_to = "Rel.Abund")
d13C.long <- d13C %>% filter(Fatty.Acid != "C19:0 (standard!)") %>% pivot_longer( -Fatty.Acid, names_to = "Sample", values_to = "d13C")
# Combine tibbles and add metadata
df <- inner_join(Areas.long, d13C.long)
df.m <-inner_join(df, metadata, by = "Sample") 
colnames(df.m)

df.m %>% filter(Fatty.Acid == "C16:0")
str(df.m)

##%#####################################################################################%##
#NIOZ164 Statia Discs - Alpha Diversity                                              #####                                            
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

## Load libraries -------------------------------------------------------------------------
library(devtools)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(ggh4x)
library(cowplot)

## Colors for plotting --------------------------------------------------------------------
pal.loc <- c("#FF6DB6FF" , "#004949FF",  "#66A61E")
# CB, CC, Zeelandia
pal.habs<- c("#CC5800FF", "#51C3CCFF")
# Benthic, Pelagic
pal.loc.hab <- c("#FF6DB6FF","#FFB6DBFF", "#004949FF", "#009292FF", "#66A61E")
# CB_P, CB_B, CC_P, CB_B, Zeelandia

## Import data ----------------------------------------------------------------------------
CC_15m <- read.delim("../Data/Temperature_Data_CC_15m.txt", sep = '\t', dec = ".")
CC_10m <- read.delim("../Data/Temperature_Data_CC_10m.txt", sep = '\t', dec = ".")

# Add locations to create tibble
CC_15m <- CC_15m[-1] %>% mutate(Location = "Benthic")
CC_10m <- CC_10m[-1] %>% mutate(Location = "Pelagic")

#remove first column with only numbers, bind 2 df's into 1
Temps <- bind_rows(CC_15m,  CC_10m, id = NULL) 

# Splite Date and Time, and calculate temp in C instead of F, and round
Temps_d <- Temps %>% separate(Date_Time, sep = " ", into = c("Date", "Time"), remove = F)
Temps_C <- Temps_d %>%  mutate(Temp_C = (Temp - 32) * (5/9)) %>%  mutate(Temp_C = round(Temp_C, digits = 2))

CC_15m_C <- CC_15m  %>%  mutate(Temp_C = (Temp - 32) * (5/9)) %>%  mutate(Temp_C = round(Temp_C, digits = 2))
CC_10m_C <- CC_10m %>% mutate(Temp_C = (Temp - 32) * (5/9)) %>%  mutate(Temp_C = round(Temp_C, digits = 2))

# %>% group_by(Date) %>% summarise(Temp_avg = mean(Temp)) %>%  ungroup()

Plot_B <- ggplot(CC_15m_C, aes(x= Date_Time, y = Temp_C, group = 1)) +
  geom_line(colour = "#CC5800FF") +
  labs(colour = NULL, x = "Time", y = "Temperature (C)") 


Plot_P <- ggplot(CC_10m_C, aes(x= Date_Time, y = Temp_C, group = 1)) +
  geom_line(colour = "#51C3CCFF") +
  labs(colour = NULL, x = "Time", y = "Temperature (C)") 

plot_grid(Plot_P, Plot_B,
          ncol = 1,
          nrow = 2,
          labels = c('A', 'B'))


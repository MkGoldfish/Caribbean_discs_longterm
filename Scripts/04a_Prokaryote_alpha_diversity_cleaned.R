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
####              Alpha and Beta Diversity              ####
#                                                          #
##%######################################################%##
# BiocManager::install("MicrobiotaProcess")
install.packages("ggforce")
install.packages("ggdist")
devtools::install_github('erocoar/gghalves')

###_______________________________________________________________________________________####
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
library("microbiome")
# library("MicrobiotaProcess")

library(ggpubr)
library("plotrix")
library("FactoMineR")
library("factoextra")
# library(mia)
library(patchwork)
# library(tidySummarizedExperiment)
# library(ANCOMBC)
# library(ALDEx2)
# library(Maaslin2)
library(knitr)
library(usedist)
library("heatmaply")
# library("TreeSummarizedExperiment")
library(ggforce)
library(ggdist)
library(gghalves)
library("ggh4x")

library(dichromat)
library("wesanderson")
library("colorspace")
library("RColorBrewer")
library("ggsci")
library("paletteer")
library(ochRe)
library("IslamicArt")
library("ggthemes")
library(harrypotter)
library(trekcolors)
library("dutchmasters")
library(palettetown)
#library(popthemes)
#library(rockthemes)
library(lisa)
#library(LaCroixColoR)
library("Polychrome")
library("colorblindr")
library("ggsci")
library("paletteer")
library(rvg)
library(officer)

library(phyloseqCompanion)

####_______________________________________________________________________________________####
####                  Import Data &                    
####_______________________________________________________________________________________####
tax <- as.matrix(read.delim('../Data/representative_seq_set_tax_assignments.txt', row.names = 1, na.strings = c(" ")))
tax <- tax_table(tax)
otu <- read.delim('../Data/20220225_PE_asvTable_noSingletons_Decon_Filt.txt', row.names = 1)
otu <- otu %>%  select(-taxonomy) %>% as.matrix() 
otu <- otu_table(otu, taxa_are_rows = T)
map <- sample_data(read.delim('../Data/OpenBio_Sample_Metadata3.txt', row.names = 1, na.strings = c("", "NA")))

set.seed(42)
####_______________________________________________________________________________________####
####           Create & check physeq object            
####_______________________________________________________________________________________####
physeq_object = merge_phyloseq(otu, tax, map) 
summarize_phyloseq(physeq_object)

# check the physeq object
source("basic_info_physeq_object.R")
basic_info_physeq_object(physeq_object)


####_______________________________________________________________________________________####
####      Check & correct taxonomy assignments         
####_______________________________________________________________________________________####
source("wonky_tonky_taxonomy.R")
wonky_tonky_taxonomy(physeq_object)

# Pruning to get rid of singletons
physeq_object <- prune_taxa(taxa_sums(physeq_object) > 1, physeq_object) #no singletons
summarize_phyloseq(physeq_object)


# transform into dataframe
taxo <- as.data.frame(physeq_object@tax_table)

####_Loop to redefine weird taxonomy to a single common character string "unassigned"_____####
for (i in 1:nrow(taxo)) {
  for (y in 1:ncol(taxo)) {
    if 
    (any(str_detect(taxo[i,y], c("uncultured","Uncultured","metagenome","Metagenome","unknown","Unknown","NA")))) {taxo[i,y] <- "unassigned" }
  }
} 

# re-define as tax table object
taxo <- tax_table(as.matrix(taxo))

# merge updated taxonomy
physeq_object <- merge_phyloseq(physeq_object@otu_table, taxo, map) 


# check the physeq object
source("basic_info_physeq_object.R")
basic_info_physeq_object(physeq_object)
summarize_phyloseq(physeq_object)

# Different colorpallettes to choose from
colors_by_Fons <- c("#004C9AFF","#FF0000FF","#FFEE00FF","#8F7700FF","#00FF24FF"
                    ,"#FD534CFF","#22DDBBFF","#0024FFFF","#0092FFFF","#00FFFFFF"
                    ,"#FFAAAAFF", "#9999FFFF","#00A000FF","#DDDDDDFF","#999999FF"
                    ,"#A73030FF","#FF00DBFF","#FF6D00FF","#AA7789FF","#BBFFBBFF"
                    ,"#EFC000FF","#B600B6FF")

colors_by_Emna <- c("#199442", "#ED1F1F", "#F5EE2C", "#B636D6", "#3D68E0", "#EBA53D", 
                    "#00688B", "#CDCD00", "#EE3A8C", "#00EE76", "#CD9B9B", "#00BFFF", 
                    "#FFF68F", "#FF7F50", "#68228B", "#ADFF2F", "#CD0000", "#0000FF", 
                    "#CD9B1D", "#FF34B3", "#BBFFFF", "#191970", "#14A821", "#E6DB45", 
                    "#EB2C2C", "#4BEE8", "#C66EE6")

colors_by_Maaike <- c("#004e64", "#ecc8af", "#F2AF29", "#436436", "#00a5cf", 
                      "#c18c5d", "#5f0f40", "#DC602E", "#495867", "#A29F15", "#570000", 
                      "#FFF5B2", "#20221B", "#9fffcb", "#c08497", "#8D6346", "#FF4B3E", "#149911", 
                      "#472d30", "#ce796b", "#25a18e", "#BC412B", "#95D9DA", 
                      "#B10F2E", "#0E273C", "#E3FDD8", "#353535", "#e7ad99", "#0F8B8D", 
                      "#7ae582", "#F2AF29", "#606c38", "#3d405b", "#94d2bd",
                      "#772e25", "#344e41", "#0047E0", "#6c584c", "#5f0f40","#D7F171", 
                      "#c89f9c", "#339989", "#faf3dd", "#04724d", "#98B9AB",
                      "#b09e99", "#AD343E", "#F2AF29", "#362C28", "#5171A5",
                      "#F7FE72", "#F4978E", "#7A9B76", "#8A7E72", "#143642", "#662C91")



pal_isme <- c("#006d77", "#ffddd2", "#00C49A", "#e29578", "#83c5be")

pal.hab.fill <- c('#EE99AA', '#6699CC', '#EECC66')
pal.hab.line <- c('#994455', '#004488', '#997700')

pal.pol.fill <- c('#EE8866','#AAAA00', '#99DDFF' )
pal.pol.line <- c('#EE7733','#117733', '#66CCEE' )

pal.time <- c('#AA4499', '#009988', '#EE3377','#0077BB','#CC3311')



####_______________________________________________________________________________________####
####           Alpha Diversity          
####_______________________________________________________________________________________####
Prevalence <- read_pptx()
Prevalence <- read_pptx("../Reports/Prevalence.pptx")

summarize_phyloseq(physeq_object)
alpha_tab <-microbiome::alpha(physeq_object, index = "all")
write.csv(alpha_tab, file = "../alpha_div/alpha_div_indexes_OpenBio_Prok.csv")

physeq_Benthic <- subset_samples(physeq_object, Habitat == "Benthic")
physeq_Pelagic <- subset_samples(physeq_object, Habitat == "Pelagic")
physeq_Eulittoral <- subset_samples(physeq_object, Habitat == "Eulittoral")

physeq_LDPE <- subset_samples(physeq_object, Polymer == "_LDPE_")
physeq_PBSeT <- subset_samples(physeq_object, Polymer == "_PBSeT_")
physeq_PHB <- subset_samples(physeq_object, Polymer == "_PHB_")


Prevalence_plot <- microbiome::plot_taxa_prevalence(physeq_PBSeT, "Phylum") + theme(legend.position = "none") + ylab("Prevalence") + labs(title = "Prevalence PBSet Prokaryotic Phyla") #prevalence
Prevalence_plot

Prevalence_plot <- microbiome::plot_taxa_prevalence(physeq_PHB, "Phylum") + theme(legend.position = "none") + ylab("Prevalence") + labs(title = "Prevalence PHB Prokaryotic Phyla") #prevalence
Prevalence_plot


Prevalence_plot <- microbiome::plot_taxa_prevalence(physeq_object, "Order")
Prevalence_plot

editable_graph <- dml(ggobj = Prevalence_plot)
Prevalence <- add_slide(Prevalence) 
Prevalence <- ph_with(x = Prevalence, editable_graph,location = ph_location_type(type = "body") )
print(Prevalence, target = "../Reports/Prevalence.pptx")


####___Alpha diversity with phyloseq________________________________________________________________####
Alpha.OpenBio <- phyloseq::estimate_richness(physeq_object, measures=c("Observed", "Simpson", "Shannon", "Chao1") )
Alpha.OpenBio <- alpha_phylo %>% cbind(data.frame(physeq_object@sam_data))
head(Alpha.OpenBio)

write_csv(Alpha.OpenBio , file = "../alpha_div/alpha_div_no_prune.csv")
Alpha.OpenBio  <- Alpha.OpenBio  %>% cbind(data.frame(physeq_object@sam_data)) 
Alpha.OpenBio$InputFileName <- NULL
Alpha.OpenBio$BarcodeSequence <- NULL
Alpha.OpenBio$BarcodeSequence_1 <- NULL
Alpha.OpenBio$LinkerPrimerSequence <- NULL
Alpha.OpenBio$ReversePrimerSequence <- NULL
Alpha.OpenBio <- Alpha.OpenBio %>% rownames_to_column(var = "Sample") %>% column_to_rownames(var = "Description")
head(Alpha.OpenBio )

#Two options for plotting
#1, select the different columns to obtain the different measures in different columns
obs.ft <- Alpha.OpenBio %>%  select(Observed, TixPolxHa, Time_cat, Season, Polymer, Habitat, TixPol, TixHa, PolxHa, TixPolxHa)
chao1 <- Alpha.OpenBio %>%  select(Chao1, TixPolxHa, Time_cat,Season, Polymer, Habitat, TixPol, TixHa, PolxHa, TixPolxHa)
Shannon <- Alpha.OpenBio %>%  select(Shannon, TixPolxHa, Time_cat,Season, Polymer, Habitat, TixPol, TixHa, PolxHa, TixPolxHa)
Simpson <- Alpha.OpenBio %>%  select(Simpson, TixPolxHa, Time_cat,Season, Polymer, Habitat, TixPol, TixHa, PolxHa, TixPolxHa)

#2, transform the complete df in long format
Alpha.OpenBio.long <- Alpha.OpenBio %>% select(!se.chao1) %>%  pivot_longer(cols = c("Observed", "Chao1", "Shannon", "Simpson"), names_to = "Diversity_Index", 
                                                                        values_to = "Value")

# Summarizing function to calcuate mean and SD
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#1 Summarize per measure
#Define groups per which you want to have the mean and SD and you want to keep for plotting
obs.ft.c <-  summarySE(obs.ft, measurevar="Observed", groupvars=c("Time_cat", "Season", "Polymer", "Habitat", "TixPol", "TixHa", "PolxHa", "TixPolxHa"))
chao1.c <-  summarySE(chao1, measurevar="Chao1", groupvars=c("Time_cat", "Season", "Polymer", "Habitat", "TixPol", "TixHa", "PolxHa", "TixPolxHa"))
Shannon.c <- summarySE(Shannon, measurevar="Shannon", groupvars=c("Time_cat", "Season", "Polymer", "Habitat", "TixPol", "TixHa", "PolxHa", "TixPolxHa"))
Simpson.c <- summarySE(Simpson, measurevar="Simpson", groupvars=c("Time_cat", "Season", "Polymer", "Habitat", "TixPol", "TixHa", "PolxHa", "TixPolxHa"))

# And bind into one df
AlphaOpenBio.c <- bind_rows(obs.ft.c, chao1.c, Shannon.c, Simpson.c)

#Summarize long format
#Here you also need to include column Diversity_Index in calculation
Alpha.OpenBio.long.c <-  summarySE(Alpha.OpenBio.long, measurevar="Value", groupvars=c("Time_cat", "Season", "Polymer", "Habitat", "TixPol", "TixHa", "PolxHa", "TixPolxHa","Diversity_Index"))
#Divide into to df's, for plotting with similar y-axes
Richness.long.c <- Alpha.OpenBio.long.c %>% filter(Diversity_Index %in% c("Observed", "Chao1"))
Diversity.long.c <- Alpha.OpenBio.long.c %>% filter(Diversity_Index %in% c("Shannon", "Simpson"))

#### Simpson evenness
Alpha.Div.OpenBio.Prok <- read_pptx()
Alpha.Div.OpenBio.Prok  <- read_pptx("../Reports/Alpha.Div.OpenBio.Prok.pptx")

## Make plot with GGPLOT
P1 <- ggplot(obs.ft.c,         #Pick data to plot
             aes(x = Time_cat, y = Observed, fill = Time_cat, color = Time_cat)) + #Pick factors to use
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se, width=.1)) +
  geom_point(size = 3) +
  facet_nested(~ fct_relevel(Habitat, 'Pelagic', 'Benthic', 'Eulittoral'  ) + Polymer, drop = T, 
               scales = "free_y", space = "free_y", switch = "y",
               axes = 'margins', 
               nest_line = element_line(),
               strip = strip_nested(
                 background_x = elem_list_rect(fill = c("#6699CC", "#EE99AA","#EECC66", rep_len("white", 9))),
                 by_layer_x = F
               )) + 
  theme_pubclean()+
  theme(legend.position = "top",
        legend.key = element_rect(fill = "white", colour = "white"),legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        plot.title = element_text(size = 20, hjust = 0.5),
        panel.border = element_rect(color = "grey90", fill = NA),
        strip.background = element_rect(color = "darkgrey"),
        panel.background = element_rect(color = "darkgrey")) +
  scale_colour_manual(values = pal.time) +
  scale_fill_manual(values = pal.time) +
  labs( title = "Observed features")

P1

# plotP1 <- ggboxplot(Alpha.OpenBio, x = "read_numb_cat_up", y = "Simpson",
#                     fill = "read_numb_cat_up", palette ="Dark2", add = "jitter", 
#                      size = 1) +facet_grid(  ~ Timepoint_num)
# plotP1 <- ggpar(plotP1,xlab ="Readnumber", title = "Simpson eveness index")
# plotP1

editable_graph <- dml(ggobj = P1)
Alphadiv_Prok<- add_slide(Alpha.Div.OpenBio.Prok) 
Alphadiv_Prok <- ph_with(x = Alpha.Div.OpenBio.Prok, editable_graph,location = ph_location_type(type = "body") )
print(Alphadiv_Prok, target = "../Reports/Alpha.Div.OpenBio.Prok.pptx")

#### Shannon Evenness
Shannon_phylo <- read_pptx()
Shannon_phylo <- read_pptx("../Reports/Shannon_phylo.pptx")

plotP2 <- ggboxplot(alpha_phylo, x = "Timepoint_num", y = "Shannon",
                    fill = "Timepoint_num", palette ="Dark2",
                     size = 1) +facet_grid( read_numb_cat_up ~ Polymer)
plotP2 <- ggpar(plotP2,legend = "none",xlab ="", title = "Shannon eveness index")
plotP2

editable_graph <- dml(ggobj = plotP2)
Shannon_phylo <- add_slide(Shannon_phylo) 
Shannon_phylo <- ph_with(x = Shannon_phylo, editable_graph,location = ph_location_type(type = "body") )
print(Shannon_phylo, target = "../Reports/Shannon_phylo.pptx")

####________________________________________________________________________________________________####
####___Alpha diversity with microbiome______________________________________________________________####
Alphadiv_micro <- read_pptx()
Alphadiv_micro <- read_pptx("../Reports/Alphadiv_micro.pptx")


alpha_micro <-microbiome::alpha(physeq_object, index = "all")
write_csv(alpha_micro, file = "../alpha_div/alpha_div_microbiome.csv")
alpha_micro <- alpha_micro %>% cbind(data.frame(physeq_object@sam_data))
head(alpha_micro)


#### Simpson evenness
plotM1 <- ggboxplot(alpha_micro, x = "Material", y = "evenness_simpson",
                    fill = "Material", palette ="Dark2",
                    add = "jitter", size = 1) +facet_grid(treatment ~ timepoint)
plotM1 <- ggpar(plotM1,legend = "none",xlab ="Material", title = "Simpson eveness index")
plotM1

editable_graph <- dml(ggobj = plotM1)
Alphadiv_micro <- add_slide(Alphadiv_micro) 
Alphadiv_micro <- ph_with(x = Alphadiv_micro, editable_graph,location = ph_location_type(type = "body") )
print(Alphadiv_micro, target = "../Reports/Alphadiv_micro.pptx")

#### Shannon Diversity
plotM2 <- ggboxplot(alpha_micro, x = "Material", y = "diversity_shannon",
                    fill = "Material", palette ="Dark2",
                    add = "jitter", size = 1) +facet_grid(treatment ~ timepoint)
plotM2 <- ggpar(plotM2,legend = "none",xlab ="Material", title = "Shannon eveness index")
plotM2

#### Chao1 index
Chao1_micro <- read_pptx()
Chao1_micro <- read_pptx("../Reports/Chao1_micro.pptx")

plotM3 <- ggboxplot(alpha_micro, x = "read_numb_cat_up", y = "chao1",
                    fill = "read_numb_cat_up", palette ="Dark2",
                    add = "jitter", size = 1) +facet_grid( ~Timepoint_num)
plotM3 <- ggpar(plotM3,legend = "none",xlab ="Readnumber", title = "Chao1 features")
plotM3

editable_graph <- dml(ggobj = plotM3)
Chao1_micro <- add_slide(Chao1_micro) 
Chao1_micro <- ph_with(x = Chao1_micro, editable_graph,location = ph_location_type(type = "body") )
print(Chao1_micro, target = "../Reports/Chao1_micro.pptx")

#### Pielou evenness
Pilou_micro <- read_pptx()
Pilou_micro <- read_pptx("../Reports/Pilou_micro.pptx")

plotM4 <- ggboxplot(alpha_micro, x = "Polymer", y = "evenness_pielou",
                    fill = "Polymer", palette ="Dark2",
                    add = "jitter", size = 1) +facet_grid( Timepoint_num ~Habitat)
plotM4 <- ggpar(plotM4,legend = "none",xlab ="", title = "Pielou evenness index")
plotM4

editable_graph <- dml(ggobj = plotM4)
Pilou_micro <- add_slide(Pilou_micro) 
Pilou_micro <- ph_with(x = Pilou_micro, editable_graph,location = ph_location_type(type = "body") )
print(Pilou_micro, target = "../Reports/Pilou_micro.pptx")



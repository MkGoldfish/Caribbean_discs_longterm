#Load libraries.
library("dplyr")
library("microDecon")

#Set working directory.
setwd('C:/Users/fdevogel/OneDrive - NIOZ/OpenBio/Molecular analysis/R_analysis')

#Load sequencing results (make sure header is correct).
All_samples_NS <- read.table("20220225_PE_asvTable_noSingletons.txt",sep = '\t', header = TRUE)


#Split samples into the different PCR batches. 
##Note: Switched column positions of sample NIOZ197.135.130 and NIOZ197.137.132.
PCR_batch1_samples <- All_samples_NS %>% select(OTU.ID,NIOZ197.253.248,NIOZ197.001.252,NIOZ197.003.254,NIOZ197.005.256,NIOZ197.007.002,NIOZ197.009.004,NIOZ197.011.006,NIOZ197.013.008,NIOZ197.015.010,NIOZ197.017.012,NIOZ197.019.014,NIOZ197.021.016,NIOZ197.023.018,NIOZ197.025.020,NIOZ197.027.022,NIOZ197.029.024,NIOZ197.031.026,NIOZ197.033.028,NIOZ197.035.030,NIOZ197.037.032,NIOZ197.039.034,NIOZ197.041.036,NIOZ197.043.038,NIOZ197.045.040,NIOZ197.047.042,NIOZ197.049.044,NIOZ197.051.046,NIOZ197.053.048,NIOZ197.055.050,NIOZ197.057.052,NIOZ197.059.054,NIOZ197.061.056,NIOZ197.063.058,NIOZ197.065.060,NIOZ197.067.062,NIOZ197.069.064,NIOZ197.071.066,NIOZ197.073.068,NIOZ197.075.070,NIOZ197.077.072,NIOZ197.079.074,NIOZ197.081.076,NIOZ197.083.078,NIOZ197.085.080,NIOZ197.087.082,taxonomy)
PCR_batch2_samples <- All_samples_NS %>% select(OTU.ID,NIOZ197.251.246,NIOZ197.255.250,NIOZ197.089.084,NIOZ197.091.086,NIOZ197.093.088,NIOZ197.095.090,NIOZ197.097.092,NIOZ197.099.094,NIOZ197.101.096,NIOZ197.103.098,NIOZ197.105.100,NIOZ197.107.102,NIOZ197.109.104,NIOZ197.111.106,NIOZ197.113.108,NIOZ197.115.110,NIOZ197.117.112,NIOZ197.119.114,NIOZ197.121.116,NIOZ197.123.118,NIOZ197.125.120,NIOZ197.127.122,NIOZ197.129.124,NIOZ197.131.126,NIOZ197.133.128,NIOZ197.137.132,NIOZ197.135.130,NIOZ197.139.134,NIOZ197.141.136,NIOZ197.143.138,NIOZ197.145.140,NIOZ197.147.142,NIOZ197.149.144,NIOZ197.151.146,NIOZ197.153.148,NIOZ197.155.150,NIOZ197.157.152,NIOZ197.159.154,NIOZ197.161.156,NIOZ197.163.158,NIOZ197.165.160,NIOZ197.167.162,NIOZ197.169.164,NIOZ197.171.166,taxonomy)
PCR_batch3_samples <- All_samples_NS %>% select(OTU.ID,NIOZ197.249.244,NIOZ197.173.168,NIOZ197.175.170,NIOZ197.177.172,NIOZ197.179.174,NIOZ197.181.176,NIOZ197.183.178,NIOZ197.185.180,NIOZ197.187.182,NIOZ197.189.184,NIOZ197.191.186,NIOZ197.193.188,NIOZ197.195.190,NIOZ197.197.192,NIOZ197.199.194,NIOZ197.201.196,NIOZ197.203.198,NIOZ197.205.200,NIOZ197.207.202,NIOZ197.209.204,NIOZ197.211.206,NIOZ197.213.208,NIOZ197.215.210,NIOZ197.217.212,NIOZ197.219.214,NIOZ197.221.216,NIOZ197.223.218,NIOZ197.225.220,NIOZ197.227.222,NIOZ197.229.224,NIOZ197.231.226,NIOZ197.233.228,NIOZ197.235.230,NIOZ197.237.232,NIOZ197.239.234,NIOZ197.241.236,NIOZ197.243.238,NIOZ197.245.240,taxonomy)
##For checking the contamination in the water that was used to dilute samples, combine all blanks with the "NTC dil. Water" NIOZ197.247.242 in a single file, as this sample is ran in all the PCR batches. Reorder columns, to put PCR blanks in front. 
PCR_dilutionwater_samples <- All_samples_NS %>% select(OTU.ID,NIOZ197.249.244,NIOZ197.251.246,NIOZ197.253.248,NIOZ197.247.242,taxonomy)


#Decontaminate samples, making use of the appropriate blanks, using MicroDecon. 
##Function needs to be split to not lose rows in the different batches, which complicates merging again later. 
PCR_batch1_Decon1 <- remove.cont(data = PCR_batch1_samples,numb.blanks=1,taxa=T)
PCR_batch2_Decon1 <- remove.cont(data = PCR_batch2_samples,numb.blanks=1,taxa=T)
PCR_batch3_Decon1 <- remove.cont(data = PCR_batch3_samples,numb.blanks=1,taxa=T)
PCR_dilutionwater_Decon1 <- remove.cont(data = PCR_dilutionwater_samples,numb.blanks=3,taxa=T)

PCR_batch1_Decon2 <- remove.thresh(data=PCR_batch1_Decon1,numb.ind=c(3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2),taxa=T)
PCR_batch2_Decon2 <- remove.thresh(data=PCR_batch2_Decon1,numb.ind=c(1,1,3,3,3,3,3,3,3,2,1,3,3,3,2,2,2,2),taxa=T)
PCR_batch3_Decon2 <- remove.thresh(data=PCR_batch3_Decon1,numb.ind=c(1,3,3,3,3,3,3,3,3,3,3,3,3),taxa=T)
PCR_dilutionwater_Decon2 <- remove.thresh(data=PCR_dilutionwater_Decon1,numb.ind=c(1),taxa=T)


#Check MicroDecon results. This function requires the numb.ind argument, which groups samples, even though remove.cont treats each samples individually. 
PCR_batch1_Decon1_report <- decon.diff(data=PCR_batch1_samples,output=PCR_batch1_Decon1,numb.blanks=1,numb.ind=c(3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2),taxa=T)
PCR_batch1_Decon1.reads.removed <- PCR_batch1_Decon1_report$reads.removed
####PCR_batch1_Decon1.difference.sum <- PCR_batch1_Decon1_report$sum.per.group
####PCR_batch1_Decon1.difference.mean <- PCR_batch1_Decon1_report$mean.per.group
PCR_batch1_Decon1.ASVs.removed <- PCR_batch1_Decon1_report$OTUs.removed

PCR_batch1_Decon2_report <- decon.diff(data=PCR_batch1_samples,output=PCR_batch1_Decon2,numb.blanks=1,numb.ind=c(3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2),taxa=T)
PCR_batch1_Decon2.reads.removed <- PCR_batch1_Decon2_report$reads.removed
####PCR_batch1_Decon2.difference.sum <- PCR_batch1_Decon2_report$sum.per.group
####PCR_batch1_Decon2.difference.mean <- PCR_batch1_Decon2_report$mean.per.group
PCR_batch1_Decon2.ASVs.removed <- PCR_batch1_Decon2_report$OTUs.removed

PCR_batch2_Decon1_report <- decon.diff(data=PCR_batch2_samples,output=PCR_batch2_Decon1,numb.blanks=1,numb.ind=c(1,1,3,3,3,3,3,3,3,2,1,3,3,3,2,2,2,2),taxa=T)
PCR_batch2_Decon1.reads.removed <- PCR_batch2_Decon1_report$reads.removed
####PCR_batch2_Decon1.difference.sum <- PCR_batch2_Decon1_report$sum.per.group
####PCR_batch2_Decon1.difference.mean <- PCR_batch2_Decon1_report$mean.per.group
PCR_batch2_Decon1.ASVs.removed <- PCR_batch2_Decon1_report$OTUs.removed 

PCR_batch2_Decon2_report <- decon.diff(data=PCR_batch2_samples,output=PCR_batch2_Decon2,numb.blanks=1,numb.ind=c(1,1,3,3,3,3,3,3,3,2,1,3,3,3,2,2,2,2),taxa=T)
PCR_batch2_Decon2.reads.removed <- PCR_batch2_Decon2_report$reads.removed
####PCR_batch2_Decon2.difference.sum <- PCR_batch2_Decon2_report$sum.per.group
####PCR_batch2_Decon2.difference.mean <- PCR_batch2_Decon2_report$mean.per.group
PCR_batch2_Decon2.ASVs.removed <- PCR_batch2_Decon2_report$OTUs.removed 

PCR_batch3_Decon1_report <- decon.diff(data=PCR_batch3_samples,output=PCR_batch3_Decon1,numb.blanks=1,numb.ind=c(1,3,3,3,3,3,3,3,3,3,3,3,3),taxa=T)
PCR_batch3_Decon1.reads.removed <- PCR_batch3_Decon1_report$reads.removed
####PCR_batch3_Decon1.difference.sum <- PCR_batch3_Decon1_report$sum.per.group
####PCR_batch3_Decon1.difference.mean <- PCR_batch3_Decon1_report$mean.per.group
PCR_batch3_Decon1.ASVs.removed <- PCR_batch3_Decon1_report$OTUs.removed 

PCR_batch3_Decon2_report <- decon.diff(data=PCR_batch3_samples,output=PCR_batch3_Decon2,numb.blanks=1,numb.ind=c(1,3,3,3,3,3,3,3,3,3,3,3,3),taxa=T)
PCR_batch3_Decon2.reads.removed <- PCR_batch3_Decon2_report$reads.removed
####PCR_batch3_Decon2.difference.sum <- PCR_batch3_Decon2_report$sum.per.group
####PCR_batch3_Decon2.difference.mean <- PCR_batch3_Decon2_report$mean.per.group
PCR_batch3_Decon2.ASVs.removed <- PCR_batch3_Decon2_report$OTUs.removed 

PCR_dilutionwater_Decon1_report <- decon.diff(data=PCR_dilutionwater_samples,output=PCR_dilutionwater_Decon1,numb.blanks=3,numb.ind=c(1),taxa=T)
PCR_dilutionwater_Decon1.reads.removed <- PCR_dilutionwater_Decon1_report$reads.removed
####PCR_dilutionwater_Decon1.difference.sum  <- PCR_dilutionwater_Decon1_report$sum.per.group
####PCR_dilutionwater_Decon1.difference.mean  <- PCR_dilutionwater_Decon1_report$mean.per.group
PCR_dilutionwater_Decon1.ASVs.removed <- PCR_dilutionwater_Decon1_report$OTUs.removed 

PCR_dilutionwater_Decon2_report <- decon.diff(data=PCR_dilutionwater_samples,output=PCR_dilutionwater_Decon2,numb.blanks=3,numb.ind=c(1),taxa=T)
PCR_dilutionwater_Decon2.reads.removed <- PCR_dilutionwater_Decon2_report$reads.removed
####PCR_dilutionwater_Decon2.difference.sum  <- PCR_dilutionwater_Decon2_report$sum.per.group
####PCR_dilutionwater_Decon2.difference.mean  <- PCR_dilutionwater_Decon2_report$mean.per.group
PCR_dilutionwater_Decon2.ASVs.removed <- PCR_dilutionwater_Decon2_report$OTUs.removed 


#Export MicroDecon Report dataframes and save as tab-separated .txt file.
write.table(PCR_batch1_Decon1.reads.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch1_Decon1.reads.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_batch1_Decon2.reads.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch1_Decon2.reads.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_batch2_Decon1.reads.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch2_Decon1.reads.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_batch2_Decon2.reads.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch2_Decon2.reads.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_batch3_Decon1.reads.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch3_Decon1.reads.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_batch3_Decon2.reads.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch3_Decon2.reads.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_dilutionwater_Decon1.reads.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_dilutionwater_Decon1.reads.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_dilutionwater_Decon2.reads.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_dilutionwater_Decon2.reads.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)

write.table(PCR_batch1_Decon1.ASVs.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch1_Decon1.ASVs.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_batch1_Decon2.ASVs.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch1_Decon2.ASVs.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_batch2_Decon1.ASVs.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch2_Decon1.ASVs.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_batch2_Decon2.ASVs.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch2_Decon2.ASVs.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_batch3_Decon1.ASVs.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch3_Decon1.ASVs.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_batch3_Decon2.ASVs.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_batch3_Decon2.ASVs.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_dilutionwater_Decon1.ASVs.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_dilutionwater_Decon1.ASVs.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_dilutionwater_Decon2.ASVs.removed,"01_Decontamination_Filtering_Prok_OpenBio/20220225_PE_asvTable_noSingletons_PCR_dilutionwater_Decon2.ASVs.removed.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)


#Combine reads again into a single table.
##Create seperate dataframes for samples and PCR controls.
All_samples_NS_Dec1 <- cbind(PCR_batch1_Decon2,PCR_batch2_Decon2,PCR_batch3_Decon2)
All_samples_NS_Dec2 <- All_samples_NS_Dec1[c(-2,-47,-48,-49,-50,-93,-94,-95)]
##Option: Decontaminate the diluted samples again for "NTC water blank" NIOZ197.247.242, using MicroDecon? 
###This impacts only 3 samples fully (NIOZ197.067.062,	NIOZ197.069.064 &	NIOZ197.071.066), and 1 reaction of the PCR triplicates of 2 other samples (NIOZ197.235.230 & NIOZ197.135.130). 
###Although the 3 samples fully impacted were biological replicates , they do not necessarily show a trend compared to other samples having more reads in the contaminated ASVs at first sight. This is similar to the two other samples.
###Decided not to decontaminate these samples further.


#Export decontaminated dataframes and save as tab-separated .txt file.
write.table(All_samples_NS_Dec2,"20220225_PE_asvTable_noSingletons_Decon.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(PCR_dilutionwater_Decon2,"20220225_PE_asvTable_noSingletons_Decon_dilutionwater.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)


#Filter reads
##Remove ASVs that have 0 reads left in them after applying MicroDecon
All_samples_NS_Dec2_F1 <- All_samples_NS_Dec2[rowSums(All_samples_NS_Dec2[, 2:124])>0, ]

##Retain only ASVs classified as Archaea and Bacteria 
All_samples_NS_Dec2_F2 <- All_samples_NS_Dec2_F1[grep("Bacteria|Archaea", All_samples_NS_Dec2_F1$taxonomy), ]

##Remove ASVs containing taxonomic identities related to chloroplasts and mitochondria.
All_samples_NS_Dec2_F3 <- All_samples_NS_Dec2_F2[!grepl("Chloroplast|Mitochondria", All_samples_NS_Dec2_F2$taxonomy),]


#Export filtered dataframes and save as tab-separated .txt file.
write.table(All_samples_NS_Dec2_F3,"20220225_PE_asvTable_noSingletons_Decon_Filt.txt",sep="\t",dec = ".",col.names=TRUE,row.names=FALSE,quote=FALSE)

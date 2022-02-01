#CST210428
#Thepurposeofthisistofiltervariantsandthenexportfilestodesktopforclumping/analysis. TheIVisdifferentforeachbloodtraitsoneedtoclumpforeachone.
# DonotnecessarilyneedRforthis, thoughalreadyhavethisinRsogoingtouseit
# Openallfilesupfront, filteronce(likelyalwayswillbetoBASO), thensimply

### StartUp
#library(MVMR)
#library(TwoSampleMR)
#library(MRPRESSO)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MR.BMI.Blood/")

### Load data and filter. Have arbitrarily chosen MPV as the blood trait of choice here to initially filter (other blood traits done in 2nd part)

e1<- fread("/project/voight_datasets/GWAS/04_giant/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz", header=T)
e1$Phenotype<- "BMI"
e1$SNP<- gsub("[A-Z]|:", "", e1$SNP) #bcWHRadjBMIandWHRhaversID:Allele:Allele
e1$ChrPosRefAlt<- paste0(e1$CHR, ":", e1$POS, "_", e1$Tested_Allele, "_", e1$Other_Allele)
e1$ChrPosAltRef<- paste0(e1$CHR, ":", e1$POS, "_", e1$Other_Allele, "_", e1$Tested_Allele)

#e2<- fread("/project/voight_datasets/GWAS/04_giant/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T)
#e2$Phenotype<- "BMI"

e3<- fread("/project/voight_datasets/GWAS/59_2020BloodTraits/BCX2_MPV_EA_GWAMA.out.gz", header=T)
e3$Phenotype<- "MPV"

#e1 <- fread("/project/voight_datasets/GWAS/61_dbGAP_MVP/79613/MVP.Lipids/MVP.EUR.TC.gwas.dbGAP.txt.gz", header=T)
#e1$Phenotype<- "TC"
#e1$ChrPosRefAlt<- paste0(e1$Chromosome, ":", e1$Position, "_", e1$Allele1, "_", e1$Allele2)
#e1$ChrPosAltRef<- paste0(e1$Chromosome, ":", e1$Position, "_", e1$Allele2, "_", e1$Allele1)

#e1_comm<- e1[e1$SNP %in% e2$SNP,] 
#e1_comm<- e1_comm[e1_comm$SNP %in% out$SNP_ID,]
e1_a<- e1[e1$ChrPosRefAlt %in% e3$rs_number,]
e1_b<- e1[e1$ChrPosAltRef %in% e3$rs_number,]
e1_comm<- bind_rows(e1_a,e1_b)

#e2_comm<- e2[e2$SNP %in% e1_comm$SNP_ID,]
#out_comm<- out[out$SNP_ID %in% e1_comm$SNP,]

e3_joina<- left_join(e1_a, e3, by= c("ChrPosRefAlt" = "rs_number"))
e3_joinb<- left_join(e1_b, e3, by= c("ChrPosAltRef" = "rs_number"))
e3_comm<- bind_rows(e3_joina, e3_joinb)

write.table(e1_comm, "BMI.filtered.txt", sep="\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(e2_comm, "BMI.filtered.txt", sep="\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
write.table(e3_comm, "MPV.filtered.txt", sep="\t", quote= FALSE, row.names= FALSE, col.names= TRUE) # this should match what is later output as 'MPV.filtered'
#write.table(out_comm, "TC.filtered.txt", sep="\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

quit()
#don't need to do blood right now

### 2nd part... Load for each blood trait, merge with WHR.filtered.txt because this has ChrPosRefAlt and ChrPosAltRef in it (which match rs_number in blood traits)

files <- list.files(path="/project/voight_datasets/GWAS/59_2020BloodTraits/", pattern="*EA_GWAMA.out.gz", full.names=FALSE, recursive=FALSE)

lapply(files, function(x) {

  #load blood trait file and assign trait
  
  t <- fread(paste0("/project/voight_datasets/GWAS/59_2020BloodTraits/", x), header=TRUE)
  t1 <- inner_join(e1_comm, t, by=c("ChrPosRefAlt" = "rs_number"))
  t2 <- inner_join(e1_comm, t, by=c("ChrPosAltRef" = "rs_number"))
  t <- bind_rows(t1, t2)

  #split at _ and use 1st sapply for real
  spl <- strsplit(x, "_" ,fixed=TRUE)
  trait <- sapply(spl, "[", 2) #first element is BCX2, 2nd is the trait
  t$Phenotype <- trait

  write.table(t, paste0(trait, ".filtered"), sep="\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

})

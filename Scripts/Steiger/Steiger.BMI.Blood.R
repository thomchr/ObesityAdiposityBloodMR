library(TwoSampleMR)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
library(psych)

setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MR.Blood.BMI/")



### Load exposure data 
####do this once, then open each blood trait file

##open exposure
e1<- fread("/project/voight_datasets/GWAS/04_giant/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz", header=T)
e1$Phenotype<- "BMI"
e1$SNP<- gsub("[A-Z]|:", "", e1$SNP) #bcWHRadjBMIandWHRhaversID:Allele:Allele
e1$ChrPosRefAlt<- paste0(e1$CHR, ":", e1$POS, "_", e1$Tested_Allele, "_", e1$Other_Allele)
e1$ChrPosAltRef<- paste0(e1$CHR, ":", e1$POS, "_", e1$Other_Allele, "_", e1$Tested_Allele)

e1_sig <- e1[e1$P <= 5*10^-8,]



# find blood trait clumped files
files <- list.files(path="/project/voight_datasets/GWAS/59_2020BloodTraits/", pattern="*EA_GWAMA.out.gz", full.names=FALSE, recursive=FALSE)

lapply(files, function(x) {

  #load blood trait file and assign trait, need to filter to make this go faster
  t <- fread(paste0("/project/voight_datasets/GWAS/59_2020BloodTraits/", x), header=TRUE)
  
  #split at _ and use 1st sapply for real ; use   . and 2 for quick plt
  spl <- strsplit(x, "_" ,fixed=TRUE)
  trait <- sapply(spl, "[", 2) #first element is BCX2, 2nd is the trait
  t$Phenotype <- trait
    
  ### take the intersection of significant SNPs
  t_sig <- t[t$'p-value' <= 5*10^-8,]
  
  t1 <- inner_join(e1_sig, t_sig, by=c("ChrPosRefAlt" = "rs_number"))
  t2 <- inner_join(e1_sig, t_sig, by=c("ChrPosAltRef" = "rs_number"))
  t_index <- bind_rows(t1, t2)

  ##check format - likely this can be used for e1 and t formatting as outputs##           

  cat("BMI on ", trait, "\n")
  cat("Nrows t_index\n")
  print(nrow(t_index))
  
  #format exposure/lipid
  e1_int <- format_data(t_index, type="exposure", snps=t_index$SNP, header=T,
                        snp_col="SNP", phenotype_col = "Phenotype.x",
                        beta_col="BETA", se_col="SE",
                        pval_col="P", samplesize_col ="N",
                        eaf_col="Freq_Tested_Allele_in_HRS",
                        effect_allele_col="Tested_Allele",
                        other_allele_col="Other_Allele")

  e1_out <- format_data(t_index, type="outcome", snps=t_index$SNP, header=T,
  	 snp_col="SNP", phenotype_col = "Phenotype.x",
                        beta_col="BETA", se_col="SE",
                        pval_col="P", samplesize_col ="N",
                        eaf_col="Freq_Tested_Allele_in_HRS",
                        effect_allele_col="Tested_Allele",
                        other_allele_col="Other_Allele")


  # format blood trait
  t_int <- format_data(t_index, type="exposure", snps=t_index$SNP, header=T,
                        snp_col="SNP", phenotype_col = "Phenotype.y",
                        beta_col="beta", se_col="se",
                        pval_col="p-value", samplesize_col ="n_samples",
                        eaf_col="eaf",
                        effect_allele_col="reference_allele",
                        other_allele_col="other_allele")
  
  t_out <- format_data(t_index, type="outcome", snps=t_index$SNP, header=T,
  	snp_col="SNP", phenotype_col = "Phenotype.y",
                        beta_col="beta", se_col="se",
                        pval_col="p-value", samplesize_col ="n_samples",
                        eaf_col="eaf",
                        effect_allele_col="reference_allele",
                        other_allele_col="other_allele")


  ### MR Steiger test of directionality

  # harmonising
  e1_harm <- harmonise_data(exposure_dat = e1_int, outcome_dat = t_out)
  print(nrow(e1_harm))

  t_harm <- harmonise_data(exposure_dat = t_int, outcome_dat = e1_out)
  print(nrow(t_harm))

  write.table(e1_harm, paste0(trait, ".Steiger.Intx.txt"), quote=F, col.names=T, row.names=F, sep="\t")


  ### MR Steiger test of directionality
  ## Testing e1->t (for BMI on blood, this is BMI->blood)
  cat("\n\nTesting e1->t\n")
  
  print(directionality_test(e1_harm))
  print(mr_steiger(p_exp = e1_harm$pval.exposure, p_out = e1_harm$pval.outcome, 
           n_exp = e1_harm$samplesize.exposure, n_out = e1_harm$samplesize.outcome, 
           r_exp = e1_harm$beta.exposure, r_out = e1_harm$beta.outcome, 
           r_xxo = 1, r_yyo = 1))
	   
  ## Testing t->e1
  cat("\n\nTesting t->e1\n")
  
  print(directionality_test(t_harm))
  print(mr_steiger(p_exp = t_harm$pval.exposure, p_out = t_harm$pval.outcome, 
           n_exp = t_harm$samplesize.exposure, n_out = t_harm$samplesize.outcome, 
           r_exp = t_harm$beta.exposure, r_out = t_harm$beta.outcome, 
           r_xxo = 1, r_yyo = 1))
  
  #write.table("Testing e1->t", paste0("Steiger.TC.",trait,".txt"), quote=F, col.names=T, row.names=F, sep="\t")
  #write.table("directionality_e1", paste0("Steiger.TC.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table(directionality_e1, paste0("Steiger.TC.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table("steiger_e1", paste0("Steiger.TC.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table(steiger_e1, paste0("Steiger.TC.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table("\n\n\n\nTesting t->e1", paste0("Steiger.TC.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table("directionality_t", paste0("Steiger.TC.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table(directionality_t, paste0("Steiger.TC.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table("steiger_t", paste0("Steiger.TC.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table(steiger_t, paste0("Steiger.TC.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")

})

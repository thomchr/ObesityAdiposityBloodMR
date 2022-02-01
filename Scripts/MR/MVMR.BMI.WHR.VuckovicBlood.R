#This is designed to take a pre-clumped BMI IV (SNPs significant for BMI), condition on WHR, and report effect sizes on Vuckovic 2020 Blood traits

#Assumes BMI IV already produced from filtering SNPs common to all 3 data sets and then clumping on Desktop

### Start Up
library(MVMR)
library(TwoSampleMR)
#library(MRPRESSO)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.BMI.WHR.VuckovicBlood/")

### Load data
e2 <- fread("BMI.WHR.Blood.Clumped", header=T)
e2$Phenotype <- "BMI"
e2_data <- format_data(e2, type="exposure", snps=e2$SNP, header=T, phenotype_col = "Phenotype", 
  	       snp_col="SNP", beta_col="BETA", se_col="SE", eaf_col="Freq_Tested_Allele_in_HRS", 
	       effect_allele_col="Tested_Allele", other_allele_col="Other_Allele", pval_col="P")

e1 <- fread("WHR.Clumped", header=T)
e1$Phenotype <- "WHR"
e1_data <- format_data(e1, type="exposure", snps=e1$SNP, header=T, phenotype_col = "Phenotype", 
               snp_col="SNP", beta_col="beta.exposure", se_col="se.exposure", eaf_col="eaf.exposure", 
               effect_allele_col="effect_allele.exposure", other_allele_col="other_allele.exposure", pval_col="pval.exposure")


### Load and run MR for each blood trait
files <- list.files(path="/project/voight_datasets/GWAS/59_2020BloodTraits/", pattern="*EA_GWAMA.out.gz", full.names=FALSE, recursive=FALSE)

lapply(files, function(x) {

  #load blood trait file and assign trait, need to filter to make this go faster
  #real one
  t <- read.table(paste0("/project/voight_datasets/GWAS/59_2020BloodTraits/", x), header=TRUE)
  t1 <- inner_join(e2, t, by=c("ChrPosRefAlt" = "rs_number"))
  t2 <- inner_join(e2, t, by=c("ChrPosAltRef" = "rs_number"))
  t <- bind_rows(t1, t2)           
  
  #split at _ and use 1st sapply for real ; use   . and 2 for quick plt
  spl <- strsplit(x, "_" ,fixed=TRUE)
  trait <- sapply(spl, "[", 2) #first element is BCX2, 2nd is the trait
  t$Phenotype <- trait
    
  # format blood trait as outcome
  t_data <- format_data(t, type="outcome", snps=t$SNP, header=T,
                       snp_col="SNP", phenotype_col = "Phenotype",
                       beta_col="beta", se_col="se",   # when already formatted might be 'beta'
                       eaf_col="eaf", effect_allele_col="reference_allele", #or EAF, A1 or effect_allele, eaf
                       other_allele_col="other_allele", pval_col="p.value")  # A2, P or other_allele and pval
    
  ### MVMR
  e1.harm <- harmonise_data(exposure_dat = e1_data, outcome_dat = t_data)
  e2.harm <- harmonise_data(exposure_dat = e2_data, outcome_dat = t_data)

  combined.e1_e2 <- merge(select(e1.harm, SNP, beta.exposure, se.exposure, eaf.exposure, beta.outcome, se.outcome, eaf.outcome), 
                          select(e2.harm, SNP, beta.exposure, se.exposure, eaf.exposure, beta.outcome, se.outcome, eaf.outcome),
                          by="SNP")
                          
  # ensure consistent harmonisation
  inconsistent <- combined.e1_e2$beta.outcome.x != combined.e1_e2$beta.outcome.y
  combined.e1_e2[inconsistent,]$beta.exposure.y <- combined.e1_e2[inconsistent,]$beta.exposure.y * -1
    
  dat1 <- combined.e1_e2
  #exposure.x is WHR, exposure.y is BMI
  dat1$r2 <- 2*(dat1$beta.exposure.y)^2 * dat1$eaf.exposure.y * (1-dat1$eaf.exposure.y) /
  (2*(dat1$beta.exposure.y)^2 * dat1$eaf.exposure.y * (1-dat1$eaf.exposure.y) + 
     (dat1$se.exposure.y)^2*2*700000*dat1$eaf.exposure.y * (1-dat1$eaf.exposure.y))
     
  cat("r2:   ")
  cat("   ")
  cat(sum(dat1$r2))
  cat("    ")
  
  write.table(dat1, paste0("MVMR.BMI.WHR.", trait, ".IV.txt"), quote=F, col.names=T, row.names=F, sep="\t")

  #write.table(combined.e1_e2, paste0("MVMR.BMI.WHR.", trait, ".IV.txt"), row.names=F)
  
  mvmr_in <- format_mvmr(BXGs = select(combined.e1_e2, beta.exposure.x, beta.exposure.y), 
  	            BYG = combined.e1_e2$beta.outcome.x, 
          	      seBXGs = select(combined.e1_e2, se.exposure.x, se.exposure.y),
            	      seBYG =  combined.e1_e2$se.outcome.x,
            	      RSID = combined.e1_e2$SNP)
  rownames(mvmr_in) <- mvmr_in$SNP
  
  cat(paste0("BMI.WHR.", trait, ".MVMRIVW.txt"))
  mvmr_out <- mvmr(mvmr_in, gencov = 0, weights = 1)
  cat("\n\n\n")
  #out <- cat(mvmr_out)
  #write.table(out, paste0("BMI.WHR.", trait, ".MVMRIVW.txt"), quote=F, col.names=T, row.names=F, sep="\t")

})
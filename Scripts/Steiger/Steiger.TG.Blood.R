library(TwoSampleMR)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
library(psych)

setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.TG.BMI.Blood/")



### Load exposure data 
####do this once, then open each blood trait file

##open exposure
e1 <- fread("/project/voight_datasets/GWAS/61_dbGAP_MVP/79613/MVP.Lipids/MVP.EUR.TG.gwas.dbGAP.txt.gz", header=T)
e1$Phenotype <- "TG"
e1$ChrPosRefAlt<- paste0(e1$Chromosome, ":", e1$Position, "_", e1$Allele1, "_", e1$Allele2)
e1$ChrPosAltRef<- paste0(e1$Chromosome, ":", e1$Position, "_", e1$Allele2, "_", e1$Allele1)

e1_sig <- e1[e1$'PValue' <= 5*10^-8,]



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

  cat("TG on ", trait, "\n")
  cat("Nrows t_index\n")
  print(nrow(t_index))
  
  #format exposure/lipid
  e1_int <- format_data(t_index, type="exposure", snps=t_index$SNP_ID, header=T,
                        snp_col="SNP_ID", phenotype_col = "Phenotype.x",
                        beta_col="Effect", se_col="SE",
                        pval_col="PValue", samplesize_col ="SampleSize",
                        eaf_col="EAF",
                        effect_allele_col="Allele1",
                        other_allele_col="Allele2")

  e1_out <- format_data(t_index, type="outcome", snps=t_index$SNP_ID, header=T,
  	    		snp_col="SNP_ID", phenotype_col = "Phenotype.x",
                        beta_col="Effect", se_col="SE",
                        pval_col="PValue", samplesize_col ="SampleSize",
                        eaf_col="EAF",
                        effect_allele_col="Allele1",
                        other_allele_col="Allele2")


  # format blood trait
  t_int <- format_data(t_index, type="exposure", snps=t_index$SNP_ID, header=T,
                        snp_col="SNP_ID", phenotype_col = "Phenotype.y",
                        beta_col="beta", se_col="se",
                        pval_col="p-value", samplesize_col ="n_samples",
                        eaf_col="eaf",
                        effect_allele_col="reference_allele",
                        other_allele_col="other_allele")
  
  t_out <- format_data(t_index, type="outcome", snps=t_index$SNP_ID, header=T,
  	snp_col="SNP_ID", phenotype_col = "Phenotype.y",
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
  ## Testing e1->t (for TG on blood, this is TG->blood)
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
  
  #write.table("Testing e1->t", paste0("Steiger.TG.",trait,".txt"), quote=F, col.names=T, row.names=F, sep="\t")
  #write.table("directionality_e1", paste0("Steiger.TG.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table(directionality_e1, paste0("Steiger.TG.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table("steiger_e1", paste0("Steiger.TG.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table(steiger_e1, paste0("Steiger.TG.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table("\n\n\n\nTesting t->e1", paste0("Steiger.TG.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table("directionality_t", paste0("Steiger.TG.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table(directionality_t, paste0("Steiger.TG.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table("steiger_t", paste0("Steiger.TG.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  #write.table(steiger_t, paste0("Steiger.TG.",trait,".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")

})

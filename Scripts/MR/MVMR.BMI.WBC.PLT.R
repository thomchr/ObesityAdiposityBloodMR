library(MVMR)
library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.BMI.WBC.PLT/")

### Load data #do this once, then open each blood trait file
e1_clump <- fread("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.BMI.PLT.WBC/BMI.Clumped", header=T)
#head(e1_clump)

e1_comm <- fread("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.BMI.PLT.WBC/BMI.filtered.txt", header=T)
#head(e1_comm)

e2 <- fread("/project/voight_datasets/GWAS/59_2020BloodTraits/BCX2_WBC_EA_GWAMA.out.gz", header=T)
e2$Phenotype <- "WBC"
  t1 <- inner_join(e1_comm, e2, by=c("ChrPosRefAlt" = "rs_number"))
  t2 <- inner_join(e1_comm, e2, by=c("ChrPosAltRef" = "rs_number"))
  e2_comm <- bind_rows(t1, t2)
  #head(e2_comm)
  
out <- fread("/project/voight_datasets/GWAS/59_2020BloodTraits/BCX2_PLT_EA_GWAMA.out.gz", header=T)
out$Phenotype <- "PLT"
  t1 <- inner_join(e1_comm, out, by=c("ChrPosRefAlt" = "rs_number"))
  t2 <- inner_join(e1_comm, out, by=c("ChrPosAltRef" = "rs_number"))
  out_comm <- bind_rows(t1, t2)
  #head(out_comm)

e1_clump <- e1_comm %>% filter(SNP %in% e1_clump$SNP)
e1_clump <- format_data(e1_clump, type="exposure", snps=e1_clump$SNP, header=T, phenotype_col = "Phenotype", 
  	      snp_col="SNP", beta_col="BETA", se_col="SE", 
	      eaf_col="Freq_Tested_Allele", effect_allele_col="Tested_Allele", other_allele_col="Other_Allele", pval_col="P")
#head(e1_clump)

e2_clump <- e2_comm[e2_comm$SNP %in% e1_clump$SNP]
e2_clump <- format_data(e2_clump, type="exposure", snps=e2_clump$SNP, header=T,
                       snp_col="SNP", phenotype_col = "Phenotype.y",
                       beta_col="beta", se_col="se",   # when already formatted might be 'beta'
                       eaf_col="eaf", effect_allele_col="reference_allele", #or EAF, A1 or effect_allele, eaf
                       other_allele_col="other_allele", pval_col="p-value")  # A2, P or other_allele and pvalue
#head(e2_clump)

out_clump <- out_comm[out_comm$SNP %in% e1_clump$SNP]
out_clump <- format_data(out_clump, type="outcome", snps=out_clump$SNP, header=T,
                       snp_col="SNP", phenotype_col = "Phenotype.y",
                       beta_col="beta", se_col="se",   # when already formatted might be 'beta'
                       eaf_col="eaf", effect_allele_col="reference_allele", #or EAF, A1 or effect_allele, eaf
                       other_allele_col="other_allele", pval_col="p-value")  # A2, P or other_allele and pvalue
#head(out_clump)

### MVMR
e1.harm <- harmonise_data(exposure_dat = e1_clump, outcome_dat = out_clump)
e2.harm <- harmonise_data(exposure_dat = e2_clump, outcome_dat = out_clump)
#e1.harm <- harmonise_data(exposure_dat = e3_clump, outcome_dat = out_clump)

combined.e1_e2 <- merge(select(e1.harm, SNP, beta.exposure, se.exposure, beta.outcome, se.outcome), 
                          select(e2.harm, SNP, beta.exposure, se.exposure, beta.outcome, se.outcome),
                          by="SNP")
#combined.e1_e2_e3 <- merge(combined.e1_e2,
#                           select(e3.harm, SNP, beta.exposure, se.exposure, beta.outcome, se.outcome),
#                           by="SNP")

# ensure consistent harmonisation
inconsistent <- combined.e1_e2$beta.outcome.x != combined.e1_e2$beta.outcome.y
combined.e1_e2[inconsistent,]$beta.exposure.y <- combined.e1_e2[inconsistent,]$beta.exposure.y * -1

#inconsistent <- combined.e1_e2_e3$beta.outcome.x != combined.e1_e2_e3$beta.outcome
#combined.e1_e2_e3[inconsistent,]$beta.exposure <- combined.e1_e2_e3[inconsistent,]$beta.exposure * -1



write.table(combined.e1_e2, paste0("BMI.WBC.PLT.MVMR.IV.txt"), sep="\t", row.names=F, col.names=TRUE, quote=FALSE)

mvmr_in <- format_mvmr(BXGs = select(combined.e1_e2, beta.exposure.x, beta.exposure.y), 
            BYG = combined.e1_e2$beta.outcome.x, 
            seBXGs = select(combined.e1_e2, se.exposure.x, se.exposure.y),
            seBYG = combined.e1_e2$se.outcome.x,
            RSID = combined.e1_e2$SNP)
rownames(mvmr_in) <- mvmr_in$SNP

cat(paste0("BMI.WBC.PLT.MVMRIVW"))
mvmr_out <- mvmr(mvmr_in, gencov = 0, weights = 1)




library(TwoSampleMR)
library(dplyr)
library(data.table)

setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MR.BMI.Blood/")


# load Exposure data and an example Astle trait to filter and clump/prepare the instrument. This is already formatted
EX <- fread("BMI.Clumped", header=T)
filter <- fread("BMI.filtered.txt", header=T)

f <- filter[filter$SNP %in% EX$SNP]


#Load and run MR for SmkInit->Each blood trait

files <- list.files(path="/project/voight_datasets/GWAS/59_2020BloodTraits/", pattern="*EA_GWAMA.out.gz", full.names=FALSE, recursive=FALSE)

lapply(files, function(x) {

  #load blood trait file and assign trait, need to filter to make this go faster
  t <- read.table(paste0("/project/voight_datasets/GWAS/59_2020BloodTraits/", x), header=TRUE)
  t1 <- inner_join(f, t, by=c("ChrPosRefAlt" = "rs_number"))
  t2 <- inner_join(f, t, by=c("ChrPosAltRef" = "rs_number"))
  t <- bind_rows(t1, t2)	       
  
  #splt at _ and use 1st sapply for real ; use   . and 2 for quick plt
  spl <- strsplit(x, "_" ,fixed=TRUE)
  trait <- sapply(spl, "[", 2) #first element is BCX2, 2nd is the trait
  t$Phenotype <- trait
    
  # format blood trait as outcome
  t_data <- format_data(t, type="outcome", snps=t$SNP, header=T,
                       snp_col="SNP", phenotype_col = "Phenotype",
                       beta_col="beta", se_col="se",   # when already formatted might be 'beta'
                       eaf_col="eaf", effect_allele_col="reference_allele", 
                       other_allele_col="other_allele", pval_col="p.value")  # A2, P .... or other_allele and pval
    



  #harmonize
  dat1 <- harmonise_data(exposure_dat = EX, outcome_dat = t_data)

  print("this is R2 for the F-statistic based on Shim et al - plug in to mRnd")
  # sample size for CAD is 547261, SmkInit is 557337 per Dan Hui (paper says up to 1,232,091) --- Using Dan's estimate, Pulit BMI is 484680, klarin lipids = , vuckovic = 
  dat1$r2 <- 2*(dat1$beta.exposure)^2 * dat1$eaf.exposure * (1-dat1$eaf.exposure) /
  	  (2*(dat1$beta.exposure)^2 * dat1$eaf.exposure * (1-dat1$eaf.exposure) +
     	  (dat1$se.exposure)^2*2*484680*dat1$eaf.exposure * (1-dat1$eaf.exposure))
  #r2 <- sum(dat1$r2)
  cat("r2\t")
  cat(sum(dat1$r2))

  write.table(dat1, paste0("BMI.", trait, ".TwoSample.IV.txt"), quote=F,col.names=T,row.names=F,sep="\t")
    
  #run MR
  res <- mr(dat1)
  
  #generate HTML file
  #mr_report(dat1)
  
  #print output
  write.table(res, paste0("BMI.", trait, ".MR.txt"), quote=F,col.names=T,row.names=F,sep="\t")
  
  #heterogeneity test
  het <- mr_heterogeneity(dat1)
  write.table(het, paste0("BMI.", trait, ".het.txt"), quote=F,col.names=T,row.names=F,sep="\t")
  #horiz pleiotropy test
  hp <- mr_pleiotropy_test(dat1)
  write.table(hp, paste0("BMI.", trait, ".hp.txt"), quote=F,col.names=T,row.names=F,sep="\t")
  



  #Plots

  #Scatterplot
  pdf(paste0("BMI.", trait, ".Scatterplot.pdf"))
  p1 <- mr_scatter_plot(res, dat1)
  p1[[1]]
  dev.off()
  
  #Forest plot
  pdf(paste0("BMI.", trait, ".Scatterplot.pdf"))
  res_single <- mr_singlesnp(dat1)
  p2 <- mr_forest_plot(res_single)
  p2[[1]]
  dev.off()
  
  pdf(paste0("BMI.", trait, ".Scatterplot.pdf"))
  res_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_two_sample_ml"))
  p3 <- mr_forest_plot(res_single)
  p3[[1]]
  dev.off()
  
  #Leave One Out Plot
  pdf(paste0("BMI.", trait, ".Scatterplot.pdf"))
  res_loo <- mr_leaveoneout(dat1)
  p4 <- mr_leaveoneout_plot(res_loo)
  p4[[1]]
  dev.off()
  
    
})



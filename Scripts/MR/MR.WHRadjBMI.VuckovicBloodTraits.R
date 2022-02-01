#This version is geared to analyze Liu et al Smoking traits or Wootten et al Lifetime Smoking and VUCKOVIC 2020 CELL blood traits


###########################################
##  Below is the iterative run for all blood traits
###############################

library(TwoSampleMR)
library(dplyr)
setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MR.WHRadjBMI.Blood/")


# load Exposure data and an example Astle trait to filter and clump/prepare the instrument. This is already formatted
####if doing on a new trait, re-do the stuff above in order to format and clump!
EX <- read.table("WHRadjBMI.HGB.Clumped", header=T)
EX$Phenotype <- "WHRadjBMI"
EX_data <- format_data(EX, type="exposure", snps=EX$SNP, header=T,
                       snp_col="SNP", phenotype_col = "Phenotype.x",
                       beta_col="BETA", se_col="SE",
                       eaf_col="Freq_Tested_Allele", effect_allele_col="Tested_Allele",
                       other_allele_col="Other_Allele", pval_col="P")

#Load and run MR for SmkInit->Each blood trait

files <- list.files(path="/project/voight_datasets/GWAS/59_2020BloodTraits/", pattern="*EA_GWAMA.out.gz", full.names=FALSE, recursive=FALSE)

lapply(files, function(x) {

  #load blood trait file and assign trait, need to filter to make this go faster
  #real one
  t <- read.table(paste0("/project/voight_datasets/GWAS/59_2020BloodTraits/", x), header=TRUE)
  t1 <- inner_join(EX, t, by=c("ChrPosRefAlt" = "rs_number"))
  t2 <- inner_join(EX, t, by=c("ChrPosAltRef" = "rs_number"))
  t <- bind_rows(t1, t2)	       
  
  #quick one
  #t <- read.table(paste0("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/test/", x), header=TRUE)
  
  #splt at _ and use 1st sapply for real ; use   . and 2 for quick plt
  spl <- strsplit(x, "_" ,fixed=TRUE)
  trait <- sapply(spl, "[", 2) #first element is BCX2, 2nd is the trait
  t$Phenotype <- trait
    
  # format blood trait as outcome
  t_data <- format_data(t, type="outcome", snps=t$SNP, header=T,
                       snp_col="SNP", phenotype_col = "Phenotype",   #note these are all .y because colnames similar WHRadjBMI and VuckovicBlood
                       beta_col="beta.y", se_col="se.y",   # when already formatted might be 'beta'
                       eaf_col="eaf.y", effect_allele_col="reference_allele.y", #or EAF, A1  ... or effect_allele, eaf
                       other_allele_col="other_allele.y", pval_col="p-value")  # A2, P .... or other_allele and pval
    
  # clump and harmonize SmkInit_data to blood trait
  dat1 <- harmonise_data(exposure_dat = EX_data, outcome_dat = t_data)
    
  #print out instrument stats (for supp table) and F statistic (based on Shim et al PMID 25898129). 
  #print("this is R2 for the F-statistic based on Shim et al - plug in to mRnd") 
  # sample size for Astle is ~ 131564, Wootton et al LfSmk = 462690, Liu et al SmkInit = 1232091, YEngo BMI is 700000
  dat1$r2 <- 2*(dat1$beta.exposure)^2 * dat1$eaf.exposure * (1-dat1$eaf.exposure) /
    (2*(dat1$beta.exposure)^2 * dat1$eaf.exposure * (1-dat1$eaf.exposure) +
       (dat1$se.exposure)^2*2*484563*dat1$eaf.exposure * (1-dat1$eaf.exposure))
  r2 <- sum(dat1$r2)
  cat("\nr2\n")
  cat(sum(dat1$r2))
  cat("\n\n")

  write.table(dat1, paste0("WHRadjBMI.", trait, ".IV.txt"), quote=F,col.names=T,row.names=F,sep="\t")
  
  #run MR
  res <- mr(dat1)
  
  #generate HTML file
  #mr_report(dat1)
  
  #print output
  write.table(res, paste0("WHRadjBMI.", trait, ".MR.txt"), quote=F,col.names=T,row.names=F,sep="\t")
  
  #heterogeneity test
  het <- mr_heterogeneity(dat1)
  write.table(het, paste0("WHRadjBMI.", trait, ".het.txt"), quote=F,col.names=T,row.names=F,sep="\t")
  #horiz pleiotropy test
  hp <- mr_pleiotropy_test(dat1)
  write.table(hp, paste0("WHRadjBMI.", trait, ".hp.txt"), quote=F,col.names=T,row.names=F,sep="\t")
  



  #Plots

  #Scatterplot
  pdf(paste0("WHRadjBMI.", trait, ".Scatterplot.pdf"))
  p1 <- mr_scatter_plot(res, dat1)
  p1[[1]]
  dev.off()
  
  #Forest plot
  pdf(paste0("WHRadjBMI.", trait, ".Scatterplot.pdf"))
  res_single <- mr_singlesnp(dat1)
  p2 <- mr_forest_plot(res_single)
  p2[[1]]
  dev.off()
  
  pdf(paste0("WHRadjBMI.", trait, ".Scatterplot.pdf"))
  res_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_two_sample_ml"))
  p3 <- mr_forest_plot(res_single)
  p3[[1]]
  dev.off()
  
  #Leave One Out Plot
  pdf(paste0("WHRadjBMI.", trait, ".Scatterplot.pdf"))
  res_loo <- mr_leaveoneout(dat1)
  p4 <- mr_leaveoneout_plot(res_loo)
  p4[[1]]
  dev.off()
  
  #Good forest plot -- taking this out for now, as seems to fail. This worked well with dplyr/ggplot2 in past though
  #pdf(paste0("SmkInit.", trait, ".GoodForest.pdf"))
  #res_single <- mr_singlesnp(dat1, all_method = c("mr_ivw",
  #                                                "mr_egger_regression",
  #                                                "mr_weighted_median"))
  #singlesnp_results <- res_single
  #exponentiate <- FALSE
  # 
  #requireNamespace("ggplot2", quietly = TRUE)
  #requireNamespace("plyr", quietly = TRUE)
  #res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d) {
  #  d <- plyr::mutate(d)
  #  
  #  if (sum(!grepl("All", d$SNP)) < 2) {
  #    return(blank_plot("Insufficient number of SNPs"))
  #    }
  #  
  #  levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "Inverse variance weighted"
  #  levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "MR Egger"
  #  levels(d$SNP)[levels(d$SNP) == "All - Weighted median"] <- "Weighted median"
  #  #d[d$SNP == "All - Inverse variance weighted", "SNP"] <- "Inverse variance weighted"
  #  #d[d$SNP == "All - MR Egger", "SNP"] <- "MR Egger"
  #  #d[d$SNP == "All - Weighted median", "SNP"] <- "Inverse variance weighted"
  #  am <- grep("All", d$SNP, value = TRUE)
  #  d$up <- d$b + 1.96 * d$se #or 
  #  d$lo <- d$b - 1.96 * d$se
  #                   
  #  # change unit
  #  # continuous gets exp(beta), binary is exp(ln(2)*beta)
  #  d$b <- exp(d$b) # d$b * log(2) for binary exposure
  #  d$up <- exp(d$up) # d$b * log(2) for binary exposure
  #                   d$lo <- exp(d$lo) # d$b * log(2) for binary exposure
  #                    
  #                   d$tot <- 0.01
  #                   d$tot[d$SNP %in% am] <- 1
  #                   d$SNP <- as.character(d$SNP)
  #                   nom <- d$SNP[!d$SNP %in% am]
  #                   nom <- nom[order(d$b)]
  #                   d <- rbind(d, d[nrow(d), ])
  #                   #d$SNP[nrow(d) - 1] <- ""
  #                   #d$b[nrow(d) - 1] <- NA
  #                   #d$up[nrow(d) - 1] <- NA
  #                   #d$lo[nrow(d) - 1] <- NA
  #                   d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
  #                   xint <- 0
  #                   if (exponentiate) {
  #                     d$b <- exp(d$b)
  #                     d$up <- exp(d$up)
  #                     d$lo <- exp(d$lo)
  #                     xint <- 1
  #                   }
  #                   #print(tail(d, 4))
  #                   d <- tail(d, 4)
  #                   d <- head(d, 3)
  #                   d[d$SNP == "Inverse variance weighted", "samplesize"] <- 3
  #                   d[d$SNP == "Weighted median", "samplesize"] <- 2
  #                   d[d$SNP == "MR Egger", "samplesize"] <- 1
  #                   d$SNP2 <- reorder(d$SNP, d$samplesize)
  #                   print(d)
  #                   ggplot2::ggplot(d, aes(y = SNP2, x = b)) + 
  #                     ggplot2::geom_vline(xintercept = 1, linetype = "dotted") + 
  #                     ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = up, size=as.factor(tot),colour = as.factor(tot)), size=4,
  #                              height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot)), size=16)  + 
  #                     ggplot2::scale_colour_manual(values = c("black", "red")) + 
  #                     ggplot2::scale_size_manual(values = c(0.3, 1)) + 
  #                     ggplot2::theme(legend.position = "none", 
  #                      axis.text.y = ggplot2::element_text(size = 14), 
  #                      axis.ticks.y = ggplot2::element_line(size = 0), 
  #                      axis.title.x = ggplot2::element_text(size = 14)) + 
  #                     ggplot2::labs(y = "", x = paste0("MR odds ratio for\n'", 
  #                                                      d$exposure[1], "' on '", d$outcome[1], "'"))
  #                 })
  #  print(res) #this prints if within Rstudio
    
})






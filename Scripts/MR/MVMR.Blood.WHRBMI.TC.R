#CST 210428
#The purpose of this is to take files that have already been filtered in Filter.Blood.WHRBMI.PAD, and then clumped on Desktop using something like ClumpSNPs.R, and run MVMR for all blood traits -WHR-BMI--> T2D or other outcomes. 
#The IV is different for each blood trait so need to run separate analyses for each one.

### Start Up
library(MVMR)
library(TwoSampleMR)
#library(MRPRESSO)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.Blood.WHR.BMI.TC/")



### Load data #do this once, then open each 
e1_comm <- fread("WHR.filtered.txt", header=T)
e2_comm <- fread("BMI.filtered.txt", header=T)
out_comm <- fread("TC.filtered.txt", header=T)

# find blood trait clumped files
files <- list.files(path="/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.Blood.WHR.BMI.TC/", pattern="*Clumped", full.names=FALSE, recursive=FALSE)

lapply(files, function(x) {

  e3_clump <- fread(paste0("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.Blood.WHR.BMI.TC/", x), header=TRUE)

  #split at . and use 1st sapply 
  spl <- strsplit(x, "." ,fixed=TRUE)
  trait <- sapply(spl, "[", 1) #first element is trait, second is 'filtered, third 'txt'
  #e3_comm$Phenotype <- trait  #this file has already been clumped and is in 'exposure' format, with Phenotype as 'exposure' and Beta as 'beta.exposure' etc
  

e1_clump <- e1_comm %>% filter(SNP %in% e3_clump$SNP)
e1_clump <- format_data(e1_clump, type="exposure", snps=e1_clump$SNP, header=T, phenotype_col = "Phenotype", 
  	      snp_col="SNP", beta_col="BETA", se_col="SE", 
	      eaf_col="Freq_Tested_Allele", effect_allele_col="Tested_Allele", other_allele_col="Other_Allele", pval_col="P")
#write.table(e1_clump, file = "WHR.Clumped", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


e2_clump <- e2_comm[e2_comm$SNP %in% e3_clump$SNP]
e2_clump <- format_data(e2_clump, type="exposure", snps=e2_clump$SNP, header=T, phenotype_col = "Phenotype", 
  	      snp_col="SNP", beta_col="BETA", se_col="SE", 
				  eaf_col="Freq_Tested_Allele", effect_allele_col="Tested_Allele", other_allele_col="Other_Allele", pval_col="P")
#write.table(e2_clump, file = "BMI.Clumped", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

out_clump <- out_comm %>% filter(SNP_ID %in% e3_clump$SNP)
out_clump <- format_data(out_clump, type="outcome", snps=out_clump$SNP_ID, header=T, phenotype_col = "Phenotype",
                     snp_col="SNP_ID",
                     beta_col="Effect", se_col="SE",
                     eaf_col="EAF", effect_allele_col="Allele1",
                     other_allele_col="Allele2", pval_col="PValue")
#write.table(out_clump, file = "TC.Clumped", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

### MVMR
e1.harm <- harmonise_data(exposure_dat = e1_clump, outcome_dat = out_clump)
e2.harm <- harmonise_data(exposure_dat = e2_clump, outcome_dat = out_clump)
e3.harm <- harmonise_data(exposure_dat = e3_clump, outcome_dat = out_clump)

combined.e1_e2 <- merge(select(e1.harm, SNP, beta.exposure, se.exposure, beta.outcome, se.outcome), 
                          select(e2.harm, SNP, beta.exposure, se.exposure, beta.outcome, se.outcome),
                          by="SNP")
combined.e1_e2_e3 <- merge(combined.e1_e2,
                           select(e3.harm, SNP, beta.exposure, se.exposure, beta.outcome, se.outcome),
                           by="SNP")

# ensure consistent harmonisation
inconsistent <- combined.e1_e2_e3$beta.outcome.x != combined.e1_e2_e3$beta.outcome.y
combined.e1_e2_e3[inconsistent,]$beta.exposure.y <- combined.e1_e2_e3[inconsistent,]$beta.exposure.y * -1

inconsistent <- combined.e1_e2_e3$beta.outcome.x != combined.e1_e2_e3$beta.outcome
combined.e1_e2_e3[inconsistent,]$beta.exposure <- combined.e1_e2_e3[inconsistent,]$beta.exposure * -1



write.table(combined.e1_e2_e3, paste0(trait, ".WHRBMI.TC.MVMR.IV.txt"), quote=F,col.names=T,row.names=F,sep="\t")

mvmr_in <- format_mvmr(BXGs = select(combined.e1_e2_e3, beta.exposure.x, beta.exposure.y, beta.exposure), 
            BYG = combined.e1_e2_e3$beta.outcome.x, 
            seBXGs = select(combined.e1_e2_e3, se.exposure.x, se.exposure.y, se.exposure),
            seBYG = combined.e1_e2_e3$se.outcome.x,
            RSID = combined.e1_e2_e3$SNP)
rownames(mvmr_in) <- mvmr_in$SNP

cat(paste0(trait, ".TC.MVMRIVW"))
mvmr_out <- mvmr(mvmr_in, gencov = 0, weights = 1)

###TwoSample MR based on MVMR filtered SNPs

#harmonize
dat1 <- harmonise_data(exposure_dat = e3_clump, outcome_dat = out_clump)

#print out instrument stats (for supp table) and R2 for F statistic (based on Shim et al PMID 25898129). 
print("this is R2 for the F-statistic based on Shim et al - plug in to mRnd") 
# sample size for CAD is 547261, SmkInit is 557337 per Dan Hui (paper says up to 1,232,091) --- Using Dan's estimate, Yengo BMI is 700000, MPV is 750000ish
dat1$r2 <- 2*(dat1$beta.exposure)^2 * dat1$eaf.exposure * (1-dat1$eaf.exposure) /
  (2*(dat1$beta.exposure)^2 * dat1$eaf.exposure * (1-dat1$eaf.exposure) + 
     (dat1$se.exposure)^2*2*563085*dat1$eaf.exposure * (1-dat1$eaf.exposure))
#or could use generic, with fstat$samplesize.outcome
#dat$r2 <- 2*(dat$beta.outcome)^2 * dat$eaf.outcome * (1-dat$eaf.outcome) / (2*(dat$beta.outcome)^2 * dat$eaf.outcome * (1-dat$eaf.outcome) + (dat$se.outcome)^2*2*dat$samplesize.outcome*dat$eaf.outcome * (1-dat$eaf.outcome))
cat("r2")
cat(sum(dat1$r2))

write.table(dat1, paste0(trait, ".TC.TwoSample.IV.txt"), quote=F,col.names=T,row.names=F,sep="\t")


#run MR
res <- mr(dat1)
write.table(res,paste0(trait,".TC.TwoSample.MR.txt"),quote=F,col.names=T,row.names=F,sep="\t")
het <- mr_heterogeneity(dat1)
write.table(het,paste0(trait,".TC.TwoSample.heterogeneity.txt"),quote=F,col.names=T,row.names=F,sep="\t")
hp <- mr_pleiotropy_test(dat1)
write.table(hp,paste0(trait,".TC.TwoSample.pleiotropy.txt"),quote=F,col.names=T,row.names=F,sep="\t")


##Plots
#Scatterplot
pdf(paste0(trait,".TC.TwoSample.MR.Scatterplot.pdf"))
mr_scatter_plot(res, dat1)
dev.off()

#Forest plot
pdf(paste0(trait,".TC.TwoSample.MR.Forestplot.pdf"))
res_single <- mr_singlesnp(dat1)
mr_forest_plot(res_single)
dev.off()

#pdf("MR.OtherForestPlot.pdf")
#res_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_two_sample_ml"))
#mr_forest_plot(res_single)
#dev.off()

#Leave One Out Plot
pdf(paste0(trait,".TC.TwoSample.MR.LeaveOneOutplot.pdf"))
res_loo <- mr_leaveoneout(dat1)
mr_leaveoneout_plot(res_loo)
dev.off()

#Good forest plot with actual estimates
pdf(paste0(trait,".TC.TwoSample.Good.forestplot.pdf"))
res_single <- mr_singlesnp(dat1, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
singlesnp_results <- res_single
exponentiate <- FALSE

requireNamespace("ggplot2", quietly = TRUE)
requireNamespace("plyr", quietly = TRUE)
res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d) {
  d <- plyr::mutate(d)
  if (sum(!grepl("All", d$SNP)) < 2) {
    return(blank_plot("Insufficient number of SNPs"))
  }
  levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "Inverse variance weighted"
  levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "MR Egger"
  levels(d$SNP)[levels(d$SNP) == "All - Weighted median"] <- "Weighted median"
  #d[d$SNP == "All - Inverse variance weighted", "SNP"] <- "Inverse variance weighted"
  #d[d$SNP == "All - MR Egger", "SNP"] <- "MR Egger"
  #d[d$SNP == "All - Weighted median", "SNP"] <- "Inverse variance weighted"
                     am <- grep("All", d$SNP, value = TRUE)

                     d$up <- d$b + 1.96 * d$se
                     d$lo <- d$b - 1.96 * d$se
                     
                     ####### change unit depending on binary/continuous!!!
                     # binary: e.g. smoking initiation
                     #d$b <- exp(d$b * log(2))
                     #d$up <- exp(d$up * log(2))
                     #d$lo <- exp(d$lo * log(2))
                     
                     # continuous: e.g., LfSmk, BMI
                     d$b <- exp(d$b)
                     d$up <- exp(d$up)
                     d$lo <- exp(d$lo)
                     
                     d$tot <- 0.01
                     d$tot[d$SNP %in% am] <- 1
                     d$SNP <- as.character(d$SNP)
                     nom <- d$SNP[!d$SNP %in% am]
                     nom <- nom[order(d$b)]
                     d <- rbind(d, d[nrow(d), ])
                     #d$SNP[nrow(d) - 1] <- ""
                     #d$b[nrow(d) - 1] <- NA
                     #d$up[nrow(d) - 1] <- NA
                     #d$lo[nrow(d) - 1] <- NA
                     d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
                     xint <- 0
                     if (exponentiate) {
                       d$b <- exp(d$b)
                       d$up <- exp(d$up)
                       d$lo <- exp(d$lo)
                       xint <- 1
                     }
                     #print(tail(d, 4))
                     d <- tail(d, 4)
                     d <- head(d, 3)
                     d[d$SNP == "Inverse variance weighted", "samplesize"] <- 3
                     d[d$SNP == "Weighted median", "samplesize"] <- 2
                     d[d$SNP == "MR Egger", "samplesize"] <- 1
                     d$SNP2 <- reorder(d$SNP, d$samplesize)
                     print(d)
                     ggplot2::ggplot(d, aes(y = SNP2, x = b)) + 
                       ggplot2::geom_vline(xintercept = 1, linetype = "dotted") + 
                       ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = up, size=as.factor(tot),colour = as.factor(tot)), size=4,
                                height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot)), size=16)  + 
                       ggplot2::scale_colour_manual(values = c("black", "red")) + 
                       ggplot2::scale_size_manual(values = c(0.3, 1)) + 
                       ggplot2::theme(legend.position = "none", 
                        axis.text.y = ggplot2::element_text(size = 14), 
                        axis.ticks.y = ggplot2::element_line(size = 0), 
                        axis.title.x = ggplot2::element_text(size = 14)) + 
                       ggplot2::labs(y = "", x = paste0("MR odds ratio for\n'", 
                                                        d$exposure[1], "' on '", d$outcome[1], "'"))
                   })
p2 <- res

p2[[1]]
dev.off()


}) #this is the end of the lapply from the start of MVMR


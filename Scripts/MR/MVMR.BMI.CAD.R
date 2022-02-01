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

setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MR.BMI.CAD/")



### Load data #do this once, then open each 
#e1_comm <- fread("WHR.filtered.txt", header=T)
e2_clump <- fread("BMI.Clumped", header=T)
out_comm <- fread("CAD.filtered.txt", header=T)


out_clump <- out_comm %>% filter(oldID %in% e2_clump$SNP)
out_clump <- format_data(out_clump, type="outcome", snps=out_clump$oldID, header=T, phenotype_col = "Phenotype",
                     snp_col="oldID",
                     beta_col="Effect", se_col="StdErr",
                     eaf_col="Freq1", effect_allele_col="Allele1",
                     other_allele_col="Allele2", pval_col="P-value")
write.table(out_clump, file = "CAD.Clumped", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


###TwoSample MR based on MVMR filtered SNPs

#harmonize
dat1 <- harmonise_data(exposure_dat = e2_clump, outcome_dat = out_clump)

#print out instrument stats (for supp table) and R2 for F statistic (based on Shim et al PMID 25898129). 
print("this is R2 for the F-statistic based on Shim et al - plug in to mRnd") 
# sample size for CAD is 547261, SmkInit is 557337 per Dan Hui (paper says up to 1,232,091) --- Using Dan's estimate, Yengo BMI is 700000, Vuckovic blood is 563085, Pulit BMI is 484680
dat1$r2 <- 2*(dat1$beta.exposure)^2 * dat1$eaf.exposure * (1-dat1$eaf.exposure) /
  (2*(dat1$beta.exposure)^2 * dat1$eaf.exposure * (1-dat1$eaf.exposure) + 
     (dat1$se.exposure)^2*2*484680*dat1$eaf.exposure * (1-dat1$eaf.exposure))
#or could use generic, with fstat$samplesize.outcome
#dat$r2 <- 2*(dat$beta.outcome)^2 * dat$eaf.outcome * (1-dat$eaf.outcome) / (2*(dat$beta.outcome)^2 * dat$eaf.outcome * (1-dat$eaf.outcome) + (dat$se.outcome)^2*2*dat$samplesize.outcome*dat$eaf.outcome * (1-dat$eaf.outcome))
cat("r2")
cat(sum(dat1$r2))

write.table(dat1, "BMI.CAD.TwoSample.IV.txt", quote=F,col.names=T,row.names=F,sep="\t")


#run MR
res <- mr(dat1)
write.table(res,"BMI.CAD.TwoSample.MR.txt",quote=F,col.names=T,row.names=F,sep="\t")
het <- mr_heterogeneity(dat1)
write.table(het,"BMI.CAD.TwoSample.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
hp <- mr_pleiotropy_test(dat1)
write.table(hp,"BMI.CAD.TwoSample.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


##Plots
#Scatterplot
pdf("BMI.CAD.TwoSample.MR.Scatterplot.pdf")
mr_scatter_plot(res, dat1)
dev.off()

#Forest plot
pdf("BMI.CAD.TwoSample.MR.Forestplot.pdf")
res_single <- mr_singlesnp(dat1)
mr_forest_plot(res_single)
dev.off()

#pdf("MR.OtherForestPlot.pdf")
#res_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_two_sample_ml"))
#mr_forest_plot(res_single)
#dev.off()

#Leave One Out Plot
pdf("BMI.CAD.TwoSample.MR.LeaveOneOutplot.pdf")
res_loo <- mr_leaveoneout(dat1)
mr_leaveoneout_plot(res_loo)
dev.off()

#Good forest plot with actual estimates
pdf("BMI.CAD.TwoSample.Good.forestplot.pdf")
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





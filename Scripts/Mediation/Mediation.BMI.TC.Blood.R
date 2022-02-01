library(MVMR)
library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.BMI.TC.Blood/")

#Load data
files <- list.files(path="/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.BMI.TC.Blood/", pattern="*MVMR.IV.txt", full.names=FALSE, recursive=FALSE)

lapply(files, function(x) {

  #load blood trait file and assign trait
  spl <- strsplit(x, "." ,fixed=TRUE)
  trait <- sapply(spl, "[", 3) #the 3rd element is the blood trait
  #t$Phenotype <- trait
  
  
  mvmr <- fread(x, stringsAsFactors = FALSE)
  snpcount <- nrow(mvmr)

  # setups
  # assume idenpendent SNPs after clumping
  rho <- diag(snpcount)
  
  ## TC is x, BMI is y, lipids/blood are 3rd which is just beta.exposure or se.exposure
  # snp -> risk factor
  betaXG <- mvmr$beta.exposure.y
  sebetaXG <- mvmr$se.exposure.y
  
  # snp -> mediator
  betaMG <- mvmr$beta.exposure.x
  sebetaMG <- mvmr$se.exposure.x
  
  # snp -> outcome
  betaYG <- mvmr$beta.outcome.x
  sebetaYG <- mvmr$se.outcome.x
  
  # calculation
  Omega <- sebetaYG%o%sebetaYG*rho
  total.effect.correl <- solve(t(betaXG)%*%solve(Omega)%*%betaXG)*t(betaXG)%*%solve(Omega)%*%betaYG
  
  se.total.effect.fixed <- sqrt(solve(t(betaXG)%*%solve(Omega)%*%betaXG))
  
  resid.total <- betaYG-as.vector(total.effect.correl)*betaXG
  
  se.total.effect.random <- sqrt(solve(t(betaXG)%*%solve(Omega)%*%betaXG))*
  max(sqrt(t(resid.total)%*%solve(Omega)%*%resid.total/(length(betaXG)-1)),1)
  
  direct.effect.correl <- solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%
  cbind(betaXG, betaMG))%*%t(cbind(betaXG, betaMG))%*%solve(Omega)%*%betaYG
  
  se.direct.effect.fixed <- sqrt(solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%cbind(betaXG, betaMG))[1,1])
  
  resid.direct <- betaYG-direct.effect.correl[1]*betaXG-direct.effect.correl[2]*betaMG
  
  se.direct.effect.random <- sqrt(solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%cbind(betaXG, betaMG))[1,1])*
  max(sqrt(t(resid.direct)%*%solve(Omega)%*%resid.direct/(length(betaXG)-2)),1)
  
  write.table("total.effect.correl", paste0("Mediation.BMI.TC.", trait, ".txt"), quote=F, col.names=T, row.names=F, sep="\t")
  write.table(total.effect.correl, paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table("\n\nse.total.effect.fixed", paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table(se.total.effect.fixed, paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table("\n\nse.total.effect.random", paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table(se.total.effect.random, paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table("\n\ndirect.effect.correl", paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table(direct.effect.correl, paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table("\n\nse.direct.effect.fixed", paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table(se.direct.effect.fixed, paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table("\n\nse.direct.effect.random", paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table(se.direct.effect.random, paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table("\n\nresid.total", paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")
  write.table(resid.total, paste0("Mediation.BMI.TC.", trait, ".txt"), append=T, quote=F, col.names=T, row.names=F, sep="\t")


  #print(total.effect.correl)
  #print(se.total.effect.fixed)
  #print(resid.total)
  #print(se.total.effect.random)
  #print(direct.effect.correl)
  #print(se.direct.effect.fixed)
  #print(se.direct.effect.random)
  
})



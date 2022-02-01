library(MVMR)
library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

setwd("/project/voight_ML/thomc/thomc_results/MendelianRandomization/MR.Laptop.Backup.200502/MR.Smk.Blood.ASCVD/2020VuckovicBloodData/MVMR.BMI.WHR.VuckovicBlood/")

#Load data
mvmr <- fread("MVMR.BMI.WHR.HGB.IV.txt", stringsAsFactors = FALSE)
snpcount <- nrow(mvmr)

# setups
# assume idenpendent SNPs after clumping
rho <- diag(snpcount)

## WHR is x, BMI is y... this is BMI as exposure so using that as 'variable X' for betaXG here
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


print(total.effect.correl)
print(se.total.effect.fixed)
print(resid.total)
print(se.total.effect.random)
print(direct.effect.correl)
print(se.direct.effect.fixed)
print(se.direct.effect.random)


## Script for Generating plots when network is fixed at baseline or evaluated at the prior time-period ##
## Carly's backbone script with James model specification ##

---
title: "Publication Facing Script"
output: html_notebook
---
  
## Load our libraries
  
#install.packages("broom.mixed",repos=NULL,contriburl="file:/opt/software/cran/src/contrib")
library(lme4)
library(nlme)
library(Gmisc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(broom.mixed)

## Load in the Data

#set some new pathways
outdir<-"/drives/drive1/54054dua/idata/cbobak/"
outdir_spillover<-paste(outdir,"spillover_analysis/",sep="")
data<-read.csv(paste(outdir_spillover,"acp_data_with_spill_strict_prior_stochastic.csv",sep=""))
data$chosp_step <- as.factor(data$hosp_step)
data$ctrialmonth <- as.factor(data$trialmonth)
data$cbregion <- as.factor(data$bregion)
data$cageg <- as.factor(data$ageg)
data$ccovid_status <- as.factor(data$covid_status)
data$ccharl_sc <- as.factor(data$charl_sc)
data$csq <- as.factor(data$sq)
data$gameTime<-data$total_game_time
#So that effect-size corresponds to 10 percentage-points, multiply the following by 10
data$physWPost_strict <- (10/8)*data$physWPost_tinv #8 = max internal physicians encountered. Standardize this way to get similarly-sized estimates to those for phase. However, its max value is < 1 = max(physWPost_strict_ext_tinv) because we divide by the total number of intervened physicians and the internal peers were not intervened on until the same step as the focial physician.
data$physWPost_strict_ext <- 10*data$physWPost_strict_ext_tinv #Multiplying by 10 means unit change (0.1) approximates mean(physWPost_strict_ext_tinv)
#Instead of the above two lines, use the following two lines when allow network to vary over follow-up:
#data$physWPost_strict <- (10/8)*data$physWPost_strict #8 = max internal physicians encountered. Standardize this way to get similarly-sized estimates to those for phase. However, its max value is < 1 = max(physWPost_strict_ext_tinv) because we divide by the total number of intervened physicians and the internal peers were not intervened on until the same step as the focal physician.
#data$physWPost_strict_ext <- 10*data$physWPost_strict_ext

#Hard code an $post \time retention$ and $expose \times post \time retention$ variables:
  
data$phaseret <- data$phase*data$day_since_intervention/90 #Dividing by 90 makes coefficient correspond to a time-unit of 90 days (1 quarter)
data$phaseexpret <- data$phase*data$physWPost_strict_ext*data$day_since_intervention/90
data$n_intret <- data$n_intervened_mds*data$day_since_intervention/90
data$n_intexpret <- data$n_intervened_mds*data$physWPost_strict*data$day_since_intervention/90
hist(data$phaseret)

#Tabulate and plot the distribution of the peer-exposure variables (phase and n_intervened) by step.
tapply(data$phase,data$hosp_step,'mean')
tapply(data$phase,data$trialmonth,'mean')
table(data$hosp_step,data$trialmonth)

#Repeat the above for number intervened
tapply(data$n_intervened_mds,data$hosp_step,'mean')
tapply(data$n_intervened_mds,data$trialmonth,'mean') #n_intervened is non-constant over last 3 months as patients may see different numbers of intervened on physicians during their stays?

#Now lets look at the distribution of intervened-on peers
tapply(data$physWPost_strict_ext,data$hosp_step,'mean') #Expect external connections to be similar across steps
tapply(data$physWPost_strict_ext,data$trialmonth,'mean')
#Do this by step to check initial randomization with respect to spillover
tapply(data$physWPost_strict_ext[data$hosp_step==1],data$trialmonth[data$hosp_step==1],'mean')
tapply(data$physWPost_strict_ext[data$hosp_step==2],data$trialmonth[data$hosp_step==2],'mean')
tapply(data$physWPost_strict_ext[data$hosp_step==3],data$trialmonth[data$hosp_step==3],'mean')
tapply(data$physWPost_strict_ext[data$hosp_step==4],data$trialmonth[data$hosp_step==4],'mean')
tapply(data$physWPost_strict_ext[data$hosp_step==5],data$trialmonth[data$hosp_step==5],'mean')
#Seems like step 1 hospitals were way more interconnected?

#Plot diffusion of intervention and exposure to spillover across time
mnphase <- tapply(data$phase,data$trialmonth,'mean')
mnspill <- tapply(data$physWPost_strict_ext,data$trialmonth,'mean')
resetlag <- 1
if (resetlag==1) {
  mnspill <- c(mnspill[2:11],mean(mnspill[9:11]))
}
mnspill <- mnspill/mnspill[11]
mx <- max(mnspill)
mnth <- as.numeric(names(mnphase))
pdf(file = "TimePlotAllNoLag_tinv.pdf",width=6,height=6)
plot(mnth,mnphase,type='l',xlab="Trial month",ylab="Proportion exposed",main='Direct and indirect exposure',ylim=c(0,mx))
lines(mnth,mnspill,type='l',col='red')
legend(8,0.4,legend=c("Direct","Indirect"),col=c('black','red'),lty=c(1,1),bty="n")
dev.off()

pdf(file = "TimePlotByStepNoLag_tinv.pdf",width=6,height=6)
par(mfrow=c(2,1),mar=c(5,4,4,2),mai=c(0.36,0.8,0.2,0))
plot(c(1,11),c(0,1),type="n",xlab="Trial month",ylab="Proportion exposed",main='Direct exposure by trial month')
for (i in 1:5) {
  mnexp <- tapply(data$phase[data$hosp_step==i],data$trialmonth[data$hosp_step==i],'mean')
  lines(mnth,mnexp,type='l',col=i)
}
legend(1,0.95,legend=c("1","2","3","4","5"),col=seq(1,5),lty=rep(1,5),bty="n")
plot(c(1,11),c(0,4.5),type="n",xlab="Trial month",ylab="Num peer exposed",main='Indirect exposure by trial month')
for (i in 1:5) {
  mnspill <- tapply(data$physWPost_strict_ext[data$hosp_step==i],data$trialmonth[data$hosp_step==i],'mean')
  if (resetlag==1) {
    mnspill <- c(mnspill[2:11],mean(mnspill[9:11]))
  }
  lines(mnth,mnspill,type='l',col=i)
}
legend(1,4.5,legend=seq("1","5"),col=seq(1,5),lty=rep(1,5),bty="n")
dev.off()

tapply(data$day_since_intervention,data$trialmonth,'mean')

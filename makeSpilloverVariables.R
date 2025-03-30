rsource <- "/drives/54054-Linux/54054dua/programs/shared/Phys_Survey/"
survey <- "/drives/54054-Linux/54054dua/idata/idata/Phys_Survey"

#Load libraries
library(sas7bdat)
library(lme4)
library(nlme)
#library(optimx)

data <- read.sas7bdat(paste(survey,"trialdataset.sas7bdat",sep="/"))

#Define categorical predictors
data$chosp_step <- as.factor(data$hosp_step)
data$ctrialmonth <- as.factor(data$trialmonth)
data$cbregion <- as.factor(data$bregion)
data$cageg <- as.factor(data$ageg)
data$ccovid_status <- as.factor(data$covid_status)
data$ncovid_status <- as.numeric(as.vector(data$ccovid_status)=="Yes")
data$ccharl_sc <- as.factor(data$charl_sc)
data$sq[data$sq=='M'] <- 'N' #Megan did not recode this
data$csq <- as.factor(data$sq)
data$nsq <- as.numeric(as.vector(data$csq)=="Y")
data$ctrialday <- (data$trialday-mean(data$trialday))/30

#Descriptive statistics
sumstats <- summary(data[,c("acp","phase","hosp_step","trialmonth","bregion","ageg","trialday","ncovid_status","cad_cd","dem_cd","heart_cd","diab_cd","lung_cd","charl_sc","n_w2mds","covid_stays","Q2_Y20_Adj_ACP","Delta_ACP","nsq")])

#Run models
#mod0 <- glmer(acp~phase+chosp_step+cageg+ccovid_status+cad_cd+heart_cd+dem_cd+lung_cd+diab_cd+ccharl_sc+csq+cbregion+n_w2mds+covid_stays+Q2_Y20_Adj_ACP+Delta_ACP+(1|HospitalID),family=binomial(link="logit"),data=data)
mod1 <- glmer(acp~phase+chosp_step+ctrialmonth+ctrialday+cbregion+cageg+ccovid_status+cad_cd+dem_cd+heart_cd+diab_cd+lung_cd+ccharl_sc+n_w2mds+covid_stays+Q2_Y20_Adj_ACP+Delta_ACP+csq+(1|HospitalID),family=binomial(link="logit"),data=data)
#mod3 <- glmer(acp~phase+chosp_step+phase*chosp_step+ctrialmonth+trialday+cageg+ccovid_status+cad_cd+heart_cd+dem_cd+lung_cd+diab_cd+ccharl_sc+csq+cbregion+n_w2mds+covid_stays+Q2_Y20_Adj_ACP+Delta_ACP+(1|HospitalID),family=binomial(link="logit"),data=data)
summary(mod1)
#note convergence failure: are we concerned?

#link hospWPost and physWPost
data$hospWPost<-NA
data$physWPost<-NA

for(i in 1:nrow(data)){
  hospid<-data$HospitalID[i]
  step<-as.numeric(data$chosp_step[i])
  stay<-data$patient_stay_id[i]
  
  data$hospWPost[i]<-hosp_trial_effects[hosp_trial_effects$hosp_id==hospid,step+1]
  
  encountered_npi<-as.character(unlist(unique(trial_encounters[trial_encounters$patient_stay_id==stay,"npi"])))
  data$physWPost[i]<-sum(phys_trial_effects[phys_trial_effects$NPI%in%encountered_npi,step+1],na.rm =T)
}

write.csv(data,paste(outdir_spillover,"acp_data_with_spill.csv",sep=""))


hosp_trial_effects_rw <- read.csv(paste(outdir_hosp,"hospital_level_spillover_effects_rw.csv",sep="/"),row.names = 1)
phys_trial_effects_rw <- read.csv(paste(outdir_phys,"physician_level_spillover_effects_rw.csv",sep="/"),row.names = 1)

#link hospWPost and physWPost
data$hospWPost<-NA
data$physWPost<-NA

for(i in 1:nrow(data)){
  hospid<-data$HospitalID[i]
  step<-as.numeric(data$chosp_step[i])
  stay<-data$patient_stay_id[i]
  
  data$hospWPost[i]<-hosp_trial_effects_rw[hosp_trial_effects_rw$hosp_id==hospid,step+1]
  
  encountered_npi<-as.character(unlist(unique(trial_encounters[trial_encounters$patient_stay_id==stay,"npi"])))
  data$physWPost[i]<-sum(phys_trial_effects_rw[phys_trial_effects_rw$NPI%in%encountered_npi,step+1],na.rm =T)
}

write.csv(data,paste(outdir_spillover,"acp_data_with_spill_stochastic.csv",sep=""))

hosp_trial_effects_rw <- read.csv(paste(outdir_hosp,"hospital_level_spillover_effects_rw_strict_prior.csv",sep="/"),row.names = 1)
phys_trial_effects_rw <- read.csv(paste(outdir_phys,"physician_level_spillover_effects_rw_strict_prior.csv",sep="/"),row.names = 1)

#link hospWPost and physWPost
data$hospWPost<-NA
data$physWPost<-NA

for(i in 1:nrow(data)){
  hospid<-data$HospitalID[i]
  step<-as.numeric(data$chosp_step[i])
  stay<-data$patient_stay_id[i]
  
  data$hospWPost[i]<-hosp_trial_effects_rw[hosp_trial_effects_rw$hosp_id==hospid,step+1]
  
  encountered_npi<-as.character(unlist(unique(trial_encounters[trial_encounters$patient_stay_id==stay,"npi"])))
  data$physWPost[i]<-sum(phys_trial_effects_rw[phys_trial_effects_rw$NPI%in%encountered_npi,step+1],na.rm =T)
}

write.csv(data,paste(outdir_spillover,"acp_data_with_spill_strict_prior_stochastic.csv",sep=""))


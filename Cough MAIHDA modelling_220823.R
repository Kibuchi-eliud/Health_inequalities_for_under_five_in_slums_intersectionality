#====================================================================================
#MAIHDA using Stan *brms package)
#Author: Eliud Kibuchi 
#========================================================================================
rm(list=ls())

library(foreign)
library(dplyr)
library(naniar)
library(plyr)
library(gmodels)
library(lsr)
library(ProbBayes)
library(dplyr)
library(ggplot2)
library(tidyverse) # for data manipulation and plots
library(haven) #for reading sav data
library(sjstats) #for calculating intra-class correlation (ICC)
library(ROCR) #for calculating area under the curve (AUC) statistics
library(brms) #for Bayesian (multilevel) generalised linear modelling
library(modelr) #for data manipulation
library(tidybayes) #for analysis of posterior draws of a Bayesian model
library(crosstable)
library(gmodels) #for cross tabulation
library(ggpubr)

setwd("")

#Import data
ncss12_full <- read.csv("NCSSS2012_children cough_21092022.csv")
names(ncss12_full)
#====================================================================================================================
summary(ncss12_full$child_agegrp)
summary(ncss12_full$child_gender)
summary(ncss12_full$hhhead_agegrp)
summary(ncss12_full$hhead_ethnicgrp)
summary(ncss12_full$women_agegrp)
names (ncss12_full)


#=========================================================================================================
#Univariate analyses
ncss12_full$hhd_wealthquintile <- as.factor(ncss12_full$hhd_wealthquintile)
ncss12_full$hhd_wealthquintile <- relevel(ncss12_full$hhd_wealthquintile, ref = "rich")

#Food availability 
summary(ncss12_full$hhfood_availability); levels(ncss12_full$hhfood_availability)
ncss12_full$hhfood_availability<-  as.factor(ncss12_full$hhfood_availability)

#Highest education 
ncss12_full$hh_educgrp <- as.factor(ncss12_full$hh_educgrp)
ncss12_full$hh_educgrp <- relevel(ncss12_full$hh_educgrp, ref = "none")
summary(ncss12_full$hh_educgrp); levels(ncss12_full$hh_educgrp)

#disability 
ncss12_full$hh_disability <- as.factor(ncss12_full$hh_disability)
ncss12_full$hh_disability <- relevel(ncss12_full$hh_disability, ref = "no")
summary(ncss12_full$hh_disability); levels(ncss12_full$hh_disability)

#women age group
ncss12_full$woman_educgrp <- as.factor(ncss12_full$woman_educgrp); summary(ncss12_full$woman_educgrp)
ncss12_full$woman_educgrp<- relevel(ncss12_full$woman_educgrp, ref = "primary")


#-======================================================================
ncss12_full$child_cough<-as.factor(ncss12_full$child_cough)
ncss12_full$cough <- ncss12_full$child_cough
summary(ncss12_full$cough)
#========================================================================
CrossTable(ncss12_full$child_agegrp, ncss12_full$cough)
CrossTable(ncss12_full$child_gender, ncss12_full$cough)

CrossTable(ncss12_full$woman_agegrp, ncss12_full$cough)
CrossTable(ncss12_full$hoh_sex, ncss12_full$cough)
CrossTable(ncss12_full$hoh_ethnicgrp, ncss12_full$cough)
CrossTable(ncss12_full$hoh_ethnicgrp, ncss12_full$cough)
CrossTable(ncss12_full$hoh_agegrp, ncss12_full$cough)
CrossTable(ncss12_full$hh_educgrp, ncss12_full$cough)
CrossTable(ncss12_full$hhd_wealthquintile, ncss12_full$cough)
CrossTable(ncss12_full$hh_staygrp2, ncss12_full$cough)
CrossTable(ncss12_full$hh_religiongrp, ncss12_full$cough)
CrossTable(ncss12_full$hh_disability, ncss12_full$cough)
CrossTable(ncss12_full$hhtenure, ncss12_full$cough)
CrossTable(ncss12_full$hhfood_availability, ncss12_full$cough)
CrossTable(ncss12_full$Income_genActivitygrp, ncss12_full$cough)
CrossTable(ncss12_full$hh_insurance_grouped, ncss12_full$cough)
CrossTable(ncss12_full$hhCHE40, ncss12_full$cough)
CrossTable(ncss12_full$woman_educgrp, ncss12_full$cough)

#========================================================================
#child_agegrp

ncss12_full$child_agegrp<-as.factor(ncss12_full$child_agegrp)
child_agegrp_model1<- glm(cough ~ 1 + child_agegrp, family=binomial(link='logit'), data = ncss12_full) 
print(summary(child_agegrp_model1), digits = 2)

#child_gender
child_gender_model1<- glm(cough ~ 1 + child_gender, family=binomial(link='logit'), data = ncss12_full) 
print(summary(child_gender_model1), digits = 3)

# hhhead_sex
ncss12_full$hoh_sex <- as.factor(ncss12_full$hoh_sex )
hhhead_sex_model1<- glm(cough ~ 1 + hoh_sex, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hhhead_sex_model1), digits = 3)

#hhhead_agegrp
ncss12_full$hoh_agegrp <- as.factor(ncss12_full$hoh_agegrp)
hhhead_agegrp_model1<- glm(cough ~ 1 + hoh_agegrp, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hhhead_agegrp_model1), digits = 3)

#hhead_ethnicgrp
ncss12_full$hoh_ethnicgrp <- as.factor(ncss12_full$hoh_ethnicgrp)
hhead_ethnicgrp_model1<- glm(cough ~ 1 + hoh_ethnicgrp, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hhead_ethnicgrp_model1), digits = 3)

#hhd_wealthquintile
ncss12_full$hhd_wealthquintile <- as.factor(ncss12_full$hhd_wealthquintile)
hhd_wealthquintile_model1<- glm(cough ~ 1 + hhd_wealthquintile, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hhd_wealthquintile_model1), digits = 3)

#hh_staygrp2
hh_staygrp2_model1<- glm(cough ~ 1 + hh_staygrp2, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hh_staygrp2_model1), digits = 3)

#hh_insurance_grouped
hh_insurance_grouped_model1<- glm(cough ~ 1 + hh_insurance_grouped, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hh_insurance_grouped_model1), digits = 3)

#hhCHE40
hhCHE40_model1<- glm(cough ~ 1 + hhCHE40, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hhCHE40_model1), digits = 3)

#hhfood_availability
hhfood_availability_model1<- glm(cough ~ 1 + hhfood_availability, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hhfood_availability_model1), digits = 3)

#Income_genActivitygrp
Income_genActivitygrp_model1<- glm(cough ~ 1 + Income_genActivitygrp, family=binomial(link='logit'), data = ncss12_full) 
print(summary(Income_genActivitygrp_model1), digits = 3)

#Highest household education 
hh_highestEdu_model1<- glm(cough ~ 1 + hh_educgrp, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hh_highestEdu_model1), digits = 3)

#"hh_religiongrp" 
hh_religiongrp_model1<- glm(cough ~ 1 + hh_religiongrp, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hh_religiongrp_model1), digits = 3)

#hh_disability
ncss12_full$hh_disability <- as.factor(ncss12_full$hh_disability)
summary(ncss12_full$hh_disability)
hh_disability_model1<- glm(cough ~ 1 + hh_disability, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hh_disability_model1), digits = 3)

#Women age 
women_agegrp_model1<- glm(cough ~ 1 + woman_agegrp, family=binomial(link='logit'), data = ncss12_full) 
print(summary(women_agegrp_model1), digits = 3)

#women education 
woman_educgrp_model1<- glm(cough ~ 1 + woman_educgrp, family=binomial(link='logit'), data = ncss12_full) 
print(summary(woman_educgrp_model1), digits = 3)

#Tenure  
summary(ncss12_full$hhtenure)
hhtenure_model1<- glm(cough ~ 1 + hhtenure, family=binomial(link='logit'), data = ncss12_full) 
print(summary(hhtenure_model1), digits = 3)

#============================================================================================================================
#============================================================================================================================
#Create intersection strata for all significant univariate variables 
#hoh_ethnicgrp,hhd_wealthquintile,hh_staygrp2,hhCHE40

ncss12_full$sig_pid <- ncss12_full%>% group_indices(hoh_ethnicgrp,hhd_wealthquintile,hh_staygrp2,hhCHE40) 
summary(as.factor(ncss12_full$sig_pid))
summary(ncss12_full$cough)

#============================================================================================================================
#Intercept model
#============================================================================================================================
sig_model1 <- brm(cough ~ 1 + (1|sig_pid), data = ncss12_full,family = bernoulli,warmup = 2000, iter = 20000, control = list(adapt_delta = 0.90),
                  cores = 2, chains = 2, 
                  seed = 123)

print(summary(sig_model1), digits = 3)
model1_results <- exp(fixef(sig_model1)); model1_results
write.csv(model1_results, "model1_results.csv")

#=====================C======================
#AUC (area under the curve).
#The AUC measures discrimination, that is, the ability of the test to correctly classify those with and without the target response. 
#(i.e. ability of model to classify under-five with and without health condition as a function of under-five's predicted probabilities, thus measuring discriminatory accuracy (DA) )
#In the current data, the target response is having a cough. 
#0.5 = This suggests no discrimination, so we might as well flip a coin (fail).
#0.5-0.7 = We consider this poor discrimination, not much better than a coin toss (poor).
#0.7-0.8 = Acceptable discrimination (fair)
#0.8-0.9= Excellent discrimination (good)
#>0.9 = Outstanding discrimination (excellent)

# Compute AUC for predicting Class with the model
Prob <- predict(sig_model1, type="response")
Prob <- Prob[,1]
Pred <- prediction(Prob, as.vector(pull(ncss12_full, cough)))
AUC <- performance(Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC

#With an AUC score of  0.884, the model discrimination is excellent
#=====================================================================================================
#VPC Indicates the share of the total individual varinace  in the propensity for having cough that is accounted for at the intersectional strata level.
variance_group <- 0.329^2; variance_group 
VPC_model1 <- 0.329^2 / (0.329^2 +3.29)*100; VPC_model1
#model1_vpc =  6.5% - indicates low discriminatory accuracy of intersection strata 
#VPC by 100 and interpreted it as the percentage share of the individual variance which lies between strata
#===================================================================================================
#Plot
#extract posterior distributions of all the random effect terms
data_RandomEffect <- ranef(sig_model1)

#extract posterior distributions of `sd(Intercept)`
r_Intercept <- data_RandomEffect$sig_pid[, , 1] %>%
  as_tibble() %>%
  rownames_to_column(var = "StrataID") %>%
  mutate(Variable = "sd(Intercept)")


#arrange in ascending order 
names(r_Intercept)
r_Intercept2 <-r_Intercept[order(r_Intercept$Estimate),]
r_Intercept2$Strata_ID <- seq.int(nrow(r_Intercept2))
head(r_Intercept2)

#plot
plot1 <- r_Intercept2 %>%
  mutate(Contain_Zero = if_else(Q2.5*Q97.5 > 0, "no", "yes")) %>%
  ggplot(aes( x = Strata_ID, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5)) +
  scale_x_continuous(breaks = round(seq(min(r_Intercept2$Strata_ID), max(r_Intercept2$Strata_ID), by = 5),1)) +
  scale_y_continuous(breaks = round(seq(-2,2, by = 0.5),1))+
  coord_cartesian(ylim = c(-2, 2))+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  xlab("Stratum rank") +
  labs(title = "A")+
  ylab("Estimated intersectional effects") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position="none")

r_Intercept2$sig_pid<- r_Intercept2$StrataID

cough_group1_model1 <- merge(ncss12_full, r_Intercept2, by = "sig_pid") 

#Export data
write.csv(cough_group1_model1, file = "cough_sig_model1.csv")

#===================================================================================================================
#model 2 : includes all significant variables (main effects) used to construct the intersections strata 
#Disentangles tha main (additive) effects from interaction effects 
#main effects -are used to model averages and 2) random effects - used to model differences (e.g., intersectional variance) 
#====================================================================================

sig_model2 <- brm(cough ~ 1 + hoh_ethnicgrp + hhd_wealthquintile + hh_staygrp2 + hhCHE40 + (1|sig_pid),
                       data = ncss12_full,family = bernoulli,warmup = 2000, iter = 20000,control = list(adapt_delta = 0.95),cores = 2, chains = 2, 
                       seed = 123)

summary(sig_model2)
sigmodel2_exp <- exp(fixef(sig_model2)); sigmodel2_exp 
write.csv(sigmodel2_exp, "model2_results.csv"
          )
#===============================================================================================================
#Correct Classification Rate
#====================================================================================================================
model2_Prob <- predict(sig_model2, type="response")
model2_Prob <- model2_Prob[,1]
modell2_Pred <- prediction(model2_Prob, as.vector(pull(ncss12_full, cough)))
AUC <- performance(modell2_Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC

#===================================================================================
#Partially adjusted intersection model (PCV)
#PCV quantifies the degree to which the different dimensions used to construct the intersectional strata contributed to the 
#between-stratum variance observed in previous model 1
#PCV = (model1_variance - model2_variance)/ model1_variance 

#VPC
model2_variance <- 0.13^2; model2_variance 
model2_vpc <-   0.13^2 / (0.13^2 +3.29)*100; model2_vpc



#PCV
model1_variance <-  0.329^2
model2_variance <-  0.13^2
model2_PCV <- ( 0.329^2 -  0.13^2)/ 0.329^2 ; model2_PCV
model2_PCV <- model2_PCV*100; model2_PCV
#PCV=84.29%
#This indicates 15.71% of variance was not explained by adding fixed effects 
#=================================================================================================
#Plot
data2_RandomEffect <- ranef(sig_model2)
#extract posterior distributions of `sd(Intercept)`
r2_Intercept <- data2_RandomEffect$sig_pid[, , 1] %>%
  as_tibble() %>%
  rownames_to_column(var = "StrataID") %>%
  mutate(Variable = "sd(Intercept)")

#arrange in ascending order 
names(r2_Intercept)
r2_Intercept2 <-r2_Intercept[order(r2_Intercept$Estimate),]
r2_Intercept2$Strata_ID <- seq.int(nrow(r2_Intercept2))
head(r2_Intercept2)

#plot
plot2 <- r2_Intercept2 %>%
  mutate(Contain_Zero = if_else(Q2.5*Q97.5 > 0, "no", "yes")) %>%
  ggplot(aes( x = Strata_ID, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5)) +
  scale_x_continuous(breaks = round(seq(min(r2_Intercept2$Strata_ID), max(r2_Intercept2$Strata_ID), by = 5),1)) +
  scale_y_continuous(breaks = round(seq(-2,2, by = 0.5),1))+
  coord_cartesian(ylim = c(-2, 2))+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  xlab("Stratum rank") +
  labs(title = "B")+
  ylab("Estimated intersectional effects") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position="none")

r2_Intercept2$personal_group_id<- r2_Intercept2$StrataID
names(r2_Intercept2)
drops <- c("StrataID" ,"Variable","Est.Error")
r2_Intercept2<- r2_Intercept2[ , !(names(r2_Intercept2) %in% drops)]
names(r2_Intercept2)

colnames(r2_Intercept2) <- c("Estimate_s1m2","Q2.5_s1m2","Q97.5_s1m2","Strata_ID_g1m2","sig_pid")


cough_group1_model2 <- merge(ncss12_full, r2_Intercept2, by = "sig_pid") 


#Export data
write.csv(cough_group1_model2, file = "cough_sig_model2.csv")

#================

#===================================================================================================================
#model 3 : includes all variables 

#Disentangles tha main (additive) effects from interaction effects 
#main effects -are used to model averages and 2) random effects - used to model differences (e.g., intersectional variance) 
#====================================================================================
names(ncss12_full)
sig_model3 <- brm(cough~ 1 + child_agegrp + child_gender + woman_agegrp + woman_educgrp + hoh_sex + hoh_agegrp + hoh_ethnicgrp +
                    hh_educgrp + hhd_wealthquintile + hh_staygrp2 +  hh_religiongrp + hh_insurance_grouped 
                  + hh_disability + hhtenure + hhfood_availability + Income_genActivitygrp + hhCHE40 + (1|sig_pid),
                  data = ncss12_full,family = bernoulli,warmup = 2000, iter = 20000,control = list(adapt_delta = 0.95),cores = 2, chains = 2, 
                  seed = 123)
summary(sig_model3)
model3_results <- exp(fixef(sig_model3)); model3_results
write.csv(model3_results, "model3_results.csv")
#=====================================================================================================================
#Correct Classification Rate
#====================================================================================================================
model3_Prob <- predict(sig_model3, type="response")
model3_Prob <- model3_Prob[,1]
model3_Pred <- prediction(model3_Prob, as.vector(pull(ncss12_full, cough)))
AUC <- performance(model3_Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC
#===================================================================================
#Partially adjusted intersectional model (PCV)
#PCV quantifies the degree to which the diffeent dimensions used to construct the intersectional strata contributed to the 
#between-stratum variance observed in previous model 1
#PCV = (model1_variance - model2_variance)/ model1_variance 

model3_variance <-    0.17^2; model3_variance 
model3_vpc <-   0.17^2 / ( 0.17^2 +3.29)*100; model3_vpc
#modell2_vpc =9.86% 

#PCV
model1_variance <- 0.329^2
model3_PCV <- (0.329^2 -  0.17^2)/0.329^2 ; model3_PCV
model3_PCV <- model3_PCV*100; model3_PCV
#=================================================================================================
#Plot
data3_RandomEffect <- ranef(sig_model3)
#extract posterior distributions of `sd(Intercept)`
r3_Intercept <- data3_RandomEffect$sig_pid[, , 1] %>%
  as_tibble() %>%
  rownames_to_column(var = "StrataID") %>%
  mutate(Variable = "sd(Intercept)")

#arrange in ascending order 
names(r3_Intercept)
r3_Intercept3 <-r3_Intercept[order(r3_Intercept$Estimate),]
r3_Intercept3$Strata_ID <- seq.int(nrow(r3_Intercept3))
head(r3_Intercept3)

#plot
#tiff("plot2.tiff", units="in", width=5, height=5, res=300)
plot3 <- r3_Intercept3 %>%
  mutate(Contain_Zero = if_else(Q2.5*Q97.5 > 0, "no", "yes")) %>%
  ggplot(aes( x = Strata_ID, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5)) +
  scale_x_continuous(breaks = round(seq(min(r_Intercept2$Strata_ID), max(r_Intercept2$Strata_ID), by = 30),1)) +
  scale_y_continuous(breaks = round(seq(-2,2, by = 0.5),1))+
  coord_cartesian(ylim = c(-2, 2))+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  xlab("Stratum rank") +
  labs(title = "C")+
  ylab("Estimated intersectional effects") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position="none")

#dev.off()

r3_Intercept3$personal_group_id<- r3_Intercept3$StrataID
names(r3_Intercept3)
drops <- c("StrataID" ,"Variable","Est.Error")
r3_Intercept3<- r3_Intercept3[ , !(names(r3_Intercept3) %in% drops)]
names(r3_Intercept3)

colnames(r3_Intercept3) <- c("Estimate_s1m3","Q2.5_s1m3","Q97.5_s1m3","Strata_ID_g1m3","sig_pid")


cough_group1_model3 <- merge(ncss12_full, r3_Intercept3, by = "sig_pid") 


#Export data
write.csv(cough_group1_model2, file = "cough_sig_model3.csv")
#================
#Plots 
grid.arrange(plot1, plot2,plot3, nrow=3, ncol=1)

cough_plot <- ggpubr::ggarrange(plot1, plot2,plot3, ncol = 1, nrow = 3)

ggsave('cough.pdf',cough_plot, width =18,height=15,dpi=1200)
#=============================================================================
#================
#


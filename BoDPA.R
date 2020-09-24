setwd("O:/Helhed/BoD food allergy/data")

rm(list=ls())
dev.off()

memory.limit(200000)

####Packages####
library(mc2d)
library(readr)
library(EnvStats)
library(plyr)
library(tidyr)
#library(fitdistrplus)
#library(goftest)
library(ggplot2)
library(reshape2)

####helpers#### 
mean_median_ci <-
  function(x) {
    c(mean = mean(x),
      median = median(x),
      quantile(x, probs = c(0.025, 0.975)))
  }

####Settings####
nvar <- 1e5 #antal iterationer i "var"
nunc <- 1e6 #antal iterationer i "unc"

set.seed(123)

##Datasets
pa_prev <- read.csv2("prevalence.csv")
pop_agegr <- read.csv2("pop_q1_2016.csv")
mortality <- read.csv2("mortality.csv")
lifeexpect <- read.csv2("national life expectancies 2015_2016.csv")

##Population statistics 

#Danish population per Q1, 2016
pop <- 5707251 

#Mean population size of one-year agegroups from 0-5, per Q1 2016
pop_agegr$sum <- pop_agegr$male + pop_agegr$female
mean_pop <- mean(pop_agegr$sum) 
mean_pop <- pop_agegr[c(1:6),c(4)]
mean_pop <- mean(mean_pop) 

#National life expectancies for YLD calculations
#both 2015:2016
lifeexpect$mean <- (lifeexpect$male + lifeexpect$female)/2
LE_0_5 <- lifeexpect[c(1:6),c(2,1:4)]
LE_0_5 <- mean(LE_0_5$mean) #78.66917

#Standard Expected Years Of Life Lost 
SEYLL_01_04 <- 89.41
SEYLL_05_09 <- 84.52
SEYLL_10_14 <- 79.53
SEYLL_15_19 <- 74.54
SEYLL_20_24 <- 69.54
SEYLL_25_29 <- 64.60
SEYLL_30_34 <- 59.63
SEYLL_35_39 <- 54.67
SEYLL_40_44 <- 49.73
SEYLL_45_49 <- 44.81
SEYLL_50_54 <- 39.92
SEYLL_55_59 <- 35.07
SEYLL_60_64 <- 30.25
SEYLL_65_69 <- 25.49
SEYLL_70_74 <- 20.77
SEYLL_75_79 <- 16.43

#Age of death - uncertain 
aod <- rpert(nunc, min = 9, mode = 32, max = 78)
str(aod)
hist(aod)
mean_median_ci(aod)

SEYLL <- 0

for (i in 1:nunc) {
  
  if(aod[i] < 5)
  {SEYLL[i] <- 0}
  if(aod[i] > 04 & aod[i] < 10){SEYLL[i] <- SEYLL_05_09}
  if(aod[i] > 09 & aod[i] < 15){SEYLL[i] <- SEYLL_10_14}
  if(aod[i] > 14 & aod[i] < 20){SEYLL[i] <- SEYLL_15_19}
  if(aod[i] > 19 & aod[i] < 25){SEYLL[i] <- SEYLL_20_24}
  if(aod[i] > 24 & aod[i] < 30){SEYLL[i] <- SEYLL_25_29}
  if(aod[i] > 29 & aod[i] < 35){SEYLL[i] <- SEYLL_30_34}
  if(aod[i] > 34 & aod[i] < 40){SEYLL[i] <- SEYLL_35_39}
  if(aod[i] > 39 & aod[i] < 45){SEYLL[i] <- SEYLL_40_44}
  if(aod[i] > 44 & aod[i] < 50){SEYLL[i] <- SEYLL_45_49}
  if(aod[i] > 49 & aod[i] < 55){SEYLL[i] <- SEYLL_50_54}
  if(aod[i] > 54 & aod[i] < 60){SEYLL[i] <- SEYLL_55_59}
  if(aod[i] > 59 & aod[i] < 65){SEYLL[i] <- SEYLL_60_64}
  if(aod[i] > 64 & aod[i] < 70){SEYLL[i] <- SEYLL_65_69}
  if(aod[i] > 69 & aod[i] < 75){SEYLL[i] <- SEYLL_70_74}
  if(aod[i] > 74 & aod[i] < 80){SEYLL[i] <- SEYLL_75_79}
}

str(SEYLL)
hist(SEYLL)
mean_median_ci(SEYLL)

#### Prevalence ####

#Beta distribution to derive prevalence: rbeta(s+1, n-s+1)
n <- sum(pa_prev$n)
s <- sum(pa_prev$confirmed)
shape1 <- s+1
shape2 <- n-s+1

pa_prev <- rbeta(nunc, shape1, shape2) 
str(pa_prev)
hist(pa_prev)
mean_median_ci(pa_prev)

####Incidence ####

#total incidence, 2016 of PA in Denmark
pa_inc <- pa_prev * mean_pop 
str(pa_inc)
hist(pa_inc)
mean_median_ci(pa_inc)

#incidence per 100,000 individuals, 2016
pa_inc_100000 <- pa_inc/pop*100000
str(pa_inc_100000)
hist(pa_inc_100000)
mean_median_ci(pa_inc_100000)

pa_inc_100000_below5 <- pa_inc/mean_pop*100000
pa_inc_100000_below5 <- pa_inc/mean_pop
str(pa_inc_100000_below5)
hist(pa_inc_100000_below5)
mean_median_ci(pa_inc_100000_below5)

#### Disability weight of living with PA ####

#Defining severity distribution - MIRABEL study, Deschildre et al. (2015) 
Pr_non_severe <- 0.7
Pr_severe <- 0.3

#Define uncertainty around dw (salomon et al. 2014)
controlled <- rpert(nunc, min = 0.007, mode = 0.015, max = 0.026)
partially_controlled <- rpert(nunc, min = 0.022, mode = 0.036, max = 0.055)
uncontrolled <- rpert(nunc, min = 0.086, mode = 0.133, max = 0.192)

dw_asthma <- (controlled + partially_controlled + uncontrolled)/3
mean_median_ci(dw_asthma)

dw_non_severe <- rpert(nunc, min = 0.031, mode = 0.049, max = 0.072) #generic uncomplicated disease, worry and daily meds = uncomplicated diabetes (Salomon et al. 2013)
dw_severe <- dw_asthma #combined dw controlled, partially controlled and uncontrolled asthma, partially controlled (IHME 2017) 

dw <- Pr_non_severe * dw_non_severe + Pr_severe * dw_severe
str(dw)
mean_median_ci(dw)

####DALY per fatal case = SEYLL * dw####

DALY_fatal <- SEYLL * 1
mean_median_ci(DALY_fatal)

####DALY per resolving PA case = d * DW####

#mean age of onset

AoO <- 2.5 #assumption - but Sicherer et al. (1998) median age of first reaction: 24 mo (range: 6 - 108 mo)
AoR <- 10 #assumption - better data?

d_res <- AoR - AoO

DALY_resolve <- d_res*dw
mean_median_ci(DALY_resolve)

####DALY per permanent case =  d * dw ####

DALY_permanent <- (LE_0_5)*dw
mean_median_ci(DALY_permanent)

####Transition probabilities####

#Mortality due to peanut allergy in Denmark #
#Pert distribution to describe the uncertainty around the mortality; min = 0 mortality, mode = swedish mortality, max = maximum identified in lit.
minimum <- 1.09/1000000 #Umasunthar 2013
mode <- 2.13/1000000 #Umasunthar 2013
maximum <- 4.6/1000000 #Umasunthar 2013

p_fatal <- rpert(nunc, minimum, mode, maximum)

p_resolve <- rpert(nunc, min = 0.10, mode = 0.22, max = 0.27) ##Peters et al. (2014)

p_permanent <- (1-(p_fatal + p_resolve)) 


#### DALY ####

#mean DALY per case of PA
DALY_case = (p_fatal * DALY_fatal) + (p_resolve * DALY_resolve) + (p_permanent * DALY_permanent)
mean_median_ci(DALY_case)

#YLL
YLL = pa_inc * p_fatal * DALY_fatal
mean_median_ci(YLL)

#YLD contribution from permanent PA
YLD_permanent = pa_inc * p_permanent * DALY_permanent 
mean_median_ci(YLD_permanent)

#YLD contribution from resolving PA
YLD_resolve = pa_inc * p_resolve * DALY_resolve
mean_median_ci(YLD_resolve)

#Total YLD
YLD <- YLD_permanent + YLD_resolve
mean_median_ci(YLD)

#Total DALY
DALY <- YLD + YLL 
str(DALY)
hist(DALY)
mean_median_ci(DALY)

#DALY per 100,000 
DALY_100000 <- DALY/pop*100000
mean_median_ci(DALY_100000)

#YLL contribution to DALY
YLL_contribution <- YLL/DALY*100
mean_median_ci(YLL_contribution)

#YLD contribution to DALY
YLD_contribution <- YLD/DALY*100
mean_median_ci(YLD_contribution)

YLD_permanent_contribution <- YLD_permanent/DALY*100
mean_median_ci(YLD_permanent_contribution)

#number of cases with each health outcome
n_permanent <- pa_inc * p_permanent
mean_median_ci(n_permanent)

n_resolve <- pa_inc * p_resolve
mean_median_ci(n_resolve)

n_fatal <- pa_inc * p_fatal
mean_median_ci(n_fatal)

####SA1 - lower mean age of death####
#Age of death - uncertain 
#aod <- rnorm(nunc, mean = 15, sd = 4)
aod = 15
str(aod)
hist(aod)
mean_median_ci(aod)

SEYLL_SA1 <- SEYLL_15_19
#SEYLL <- 0

#for (i in 1:nunc) {
  
#  if(aod[i] < 5)
#  {SEYLL[i] <- 0}
#  if(aod[i] > 04 & aod[i] < 10){SEYLL[i] <- SEYLL_05_09}
#  if(aod[i] > 09 & aod[i] < 15){SEYLL[i] <- SEYLL_10_14}
#  if(aod[i] > 14 & aod[i] < 20){SEYLL[i] <- SEYLL_15_19}
#  if(aod[i] > 19 & aod[i] < 25){SEYLL[i] <- SEYLL_20_24}
#  if(aod[i] > 24 & aod[i] < 30){SEYLL[i] <- SEYLL_25_29}
#  if(aod[i] > 29 & aod[i] < 35){SEYLL[i] <- SEYLL_30_34}
#  if(aod[i] > 34 & aod[i] < 40){SEYLL[i] <- SEYLL_35_39}
#  if(aod[i] > 39 & aod[i] < 45){SEYLL[i] <- SEYLL_40_44}
#  if(aod[i] > 44 & aod[i] < 50){SEYLL[i] <- SEYLL_45_49}
#  if(aod[i] > 49 & aod[i] < 55){SEYLL[i] <- SEYLL_50_54}
#  if(aod[i] > 54 & aod[i] < 60){SEYLL[i] <- SEYLL_55_59}
#  if(aod[i] > 59 & aod[i] < 65){SEYLL[i] <- SEYLL_60_64}
#  if(aod[i] > 64 & aod[i] < 70){SEYLL[i] <- SEYLL_65_69}
#  if(aod[i] > 69 & aod[i] < 75){SEYLL[i] <- SEYLL_70_74}
#  if(aod[i] > 74 & aod[i] < 80){SEYLL[i] <- SEYLL_75_79}
#}

str(SEYLL_SA1)
hist(SEYLL_SA1)
mean_median_ci(SEYLL_SA1)

DALY_fatal <- SEYLL_SA1 * 1
mean_median_ci(DALY_fatal)

#SA1:mean DALY per case of PA
DALY_case_SA1 = (p_fatal * DALY_fatal) + (p_resolve * DALY_resolve) + (p_permanent * DALY_permanent)
mean_median_ci(DALY_case_SA1)

#SA1:YLL
YLL_SA1 = pa_inc * p_fatal * DALY_fatal
mean_median_ci(YLL_SA1)

#SA1:YLD contribution from permanent PA
YLD_permanent_SA1 = pa_inc * p_permanent * DALY_permanent 
mean_median_ci(YLD_permanent_SA1)

#SA1:YLD contribution from resolving PA
YLD_resolve_SA1 = pa_inc * p_resolve * DALY_resolve
mean_median_ci(YLD_resolve_SA1)

#SA1:Total YLD
YLD_SA1 <- YLD_permanent + YLD_resolve
mean_median_ci(YLD_SA1)

#SA1:Total DALY
DALY_SA1 <- YLD_SA1 + YLL_SA1 
str(DALY_SA1)
hist(DALY_SA1)
mean_median_ci(DALY_SA1)

#SA1:YLL contribution to DALY
YLL_contribution_SA1 <- YLL_SA1/DALY_SA1*100
mean_median_ci(YLL_contribution_SA1)

#SA1:YLD contribution to DALY
YLD_contribution_SA1 <- YLD_SA1/DALY_SA1*100
mean_median_ci(YLD_contribution_SA1)

YLD_permanent_contribution_SA1 <- YLD_permanent_SA1/DALY_SA1*100
mean_median_ci(YLD_permanent_contribution_SA1)


####SA2 - risk of death is zero####

p_fatal <- 0

p_resolve <- rpert(nunc, min = 0.10, mode = 0.22, max = 0.27) ##Peters et al. (2014)

p_permanent <- (1-(p_fatal + p_resolve)) 

#SA2:mean DALY per case of PA
DALY_case_SA2 = (p_fatal * DALY_fatal) + (p_resolve * DALY_resolve) + (p_permanent * DALY_permanent)
mean_median_ci(DALY_case_SA2)

#SA2:YLL
YLL_SA2 = pa_inc * p_fatal * DALY_fatal
mean_median_ci(YLL_SA2)

#SA2:YLD contribution from permanent PA
YLD_permanent_SA2 = pa_inc * p_permanent * DALY_permanent 
mean_median_ci(YLD_permanent_SA2)

#SA2:YLD contribution from resolving PA
YLD_resolve_SA2 = pa_inc * p_resolve * DALY_resolve
mean_median_ci(YLD_resolve_SA2)

#SA2:Total YLD
YLD_SA2 <- YLD_permanent_SA2 + YLD_resolve_SA2
mean_median_ci(YLD_SA2)

#SA2:Total DALY
DALY_SA2 <- YLD_SA2 + YLL_SA2 
str(DALY_SA2)
hist(DALY_SA2)
mean_median_ci(DALY_SA2)

#SA2:YLL contribution to DALY
YLL_contribution_SA2 <- YLL_SA2/DALY_SA2*100
mean_median_ci(YLL_contribution_SA2)

#SA2:YLD contribution to DALY
YLD_contribution_SA2 <- YLD_SA2/DALY_SA2*100
mean_median_ci(YLD_contribution_SA2)

YLD_permanent_contribution_SA2 <- YLD_permanent_SA2/DALY_SA2*100
mean_median_ci(YLD_permanent_contribution_SA2)

####SA3a - lower dw ####
#SA3.a - dw = DM = 0.049

#dw equal to mild reumathoid 
dw_SA3a <- dw_non_severe
mean_median_ci(dw_SA3a)

aod <- rpert(nunc, min = 9, mode = 32, max = 78)
str(aod)
hist(aod)
mean_median_ci(aod)

SEYLL <- 0

for (i in 1:nunc) {
  
  if(aod[i] < 5)
  {SEYLL[i] <- 0}
  if(aod[i] > 04 & aod[i] < 10){SEYLL[i] <- SEYLL_05_09}
  if(aod[i] > 09 & aod[i] < 15){SEYLL[i] <- SEYLL_10_14}
  if(aod[i] > 14 & aod[i] < 20){SEYLL[i] <- SEYLL_15_19}
  if(aod[i] > 19 & aod[i] < 25){SEYLL[i] <- SEYLL_20_24}
  if(aod[i] > 24 & aod[i] < 30){SEYLL[i] <- SEYLL_25_29}
  if(aod[i] > 29 & aod[i] < 35){SEYLL[i] <- SEYLL_30_34}
  if(aod[i] > 34 & aod[i] < 40){SEYLL[i] <- SEYLL_35_39}
  if(aod[i] > 39 & aod[i] < 45){SEYLL[i] <- SEYLL_40_44}
  if(aod[i] > 44 & aod[i] < 50){SEYLL[i] <- SEYLL_45_49}
  if(aod[i] > 49 & aod[i] < 55){SEYLL[i] <- SEYLL_50_54}
  if(aod[i] > 54 & aod[i] < 60){SEYLL[i] <- SEYLL_55_59}
  if(aod[i] > 59 & aod[i] < 65){SEYLL[i] <- SEYLL_60_64}
  if(aod[i] > 64 & aod[i] < 70){SEYLL[i] <- SEYLL_65_69}
  if(aod[i] > 69 & aod[i] < 75){SEYLL[i] <- SEYLL_70_74}
  if(aod[i] > 74 & aod[i] < 80){SEYLL[i] <- SEYLL_75_79}
}

str(SEYLL)
hist(SEYLL)
mean_median_ci(SEYLL)

#DALY per fatal PA case = SEYLL * 1#
DALY_fatal <- SEYLL * 1
mean_median_ci(DALY_fatal)

#DALY per resolving PA case = d * DW#
d_res <- AoR - AoO

DALY_resolve_SA3a <- d_res*dw_SA3a
mean_median_ci(DALY_resolve_SA3a)

#DALY per permanent case =  d * dw #

DALY_permanent_SA3a <- (LE_0_5)*dw_SA3a
mean_median_ci(DALY_permanent_SA3a)

#Transition probabilities#

#Mortality due to peanut allergy in Denmark #
#Pert distribution to describe the uncertainty around the mortality; min = 0 mortality, mode = swedish mortality, max = maximum identified in lit.
minimum <- 1.09/1000000 #Umasunthar 2013
mode <- 2.13/1000000 #Umasunthar 2013
maximum <- 4.6/1000000 #Umasunthar 2013

p_fatal <- rpert(nunc, minimum, mode, maximum)

p_resolve <- rpert(nunc, min = 0.10, mode = 0.22, max = 0.27) ##Peters et al. (2014)

p_permanent <- (1-(p_fatal + p_resolve)) 

#DALY #

#mean DALY per case of PA
DALY_case_SA3a = (p_fatal * DALY_fatal) + (p_resolve * DALY_resolve_SA3a) + (p_permanent * DALY_permanent_SA3a)
mean_median_ci(DALY_case_SA3a)

#YLL
YLL_SA3a = pa_inc * p_fatal * DALY_fatal
mean_median_ci(YLL_SA3a)

#YLD contribution from permanent PA
YLD_permanent_SA3a = pa_inc * p_permanent * DALY_permanent_SA3a 
mean_median_ci(YLD_permanent_SA3a)

#YLD contribution from resolving PA
YLD_resolve_SA3a = pa_inc * p_resolve * DALY_resolve_SA3a
mean_median_ci(YLD_resolve_SA3a)

#Total YLD
YLD_SA3a <- YLD_permanent_SA3a + YLD_resolve_SA3a
mean_median_ci(YLD_SA3a)

#Total DALY
DALY_SA3a <- YLD_SA3a + YLL_SA3a 
str(DALY_SA3a)
hist(DALY_SA3a)
mean_median_ci(DALY_SA3a)

#YLL contribution to DALY
YLL_contribution_SA3a <- YLL_SA3a/DALY_SA3a*100
mean_median_ci(YLL_contribution_SA3a)

#YLD contribution to DALY
YLD_contribution_SA3a <- YLD_SA3a/DALY_SA3a*100
mean_median_ci(YLD_contribution_SA3a)

YLD_permanent_contribution_SA3a <- YLD_permanent_SA3a/DALY_SA3a*100
mean_median_ci(YLD_permanent_contribution_SA3a)

####SA3b - higher dw ####
#SA3.a - dw = DM = 0.049

#dw equal to mild reumathoid 
dw_SA3b <- dw_severe
mean_median_ci(dw_SA3b)

aod <- rpert(nunc, min = 9, mode = 32, max = 78)
str(aod)
hist(aod)
mean_median_ci(aod)

#DALY per fatal PA case = SEYLL * 1#
DALY_fatal <- SEYLL * 1
mean_median_ci(DALY_fatal)

#DALY per resolving PA case = d * DW#
d_res <- AoR - AoO

DALY_resolve_SA3b <- d_res*dw_SA3b
mean_median_ci(DALY_resolve_SA3b)

#DALY per permanent case =  d * dw #

DALY_permanent_SA3b <- (LE_0_5)*dw_SA3b
mean_median_ci(DALY_permanent_SA3b)

#Transition probabilities#

#Mortality due to peanut allergy in Denmark #
#Pert distribution to describe the uncertainty around the mortality; min = 0 mortality, mode = swedish mortality, max = maximum identified in lit.
minimum <- 1.09/1000000 #Umasunthar 2013
mode <- 2.13/1000000 #Umasunthar 2013
maximum <- 4.6/1000000 #Umasunthar 2013

p_fatal <- rpert(nunc, minimum, mode, maximum)

p_resolve <- rpert(nunc, min = 0.10, mode = 0.22, max = 0.27) ##Peters et al. (2014)

p_permanent <- (1-(p_fatal + p_resolve)) 

#DALY #

#mean DALY per case of PA
DALY_case_SA3b = (p_fatal * DALY_fatal) + (p_resolve * DALY_resolve_SA3b) + (p_permanent * DALY_permanent_SA3b)
mean_median_ci(DALY_case_SA3b)

#YLL
YLL_SA3b = pa_inc * p_fatal * DALY_fatal
mean_median_ci(YLL_SA3b)

#YLD contribution from permanent PA
YLD_permanent_SA3b = pa_inc * p_permanent * DALY_permanent_SA3b 
mean_median_ci(YLD_permanent_SA3b)

#YLD contribution from resolving PA
YLD_resolve_SA3b = pa_inc * p_resolve * DALY_resolve_SA3b
mean_median_ci(YLD_resolve_SA3b)

#Total YLD
YLD_SA3b <- YLD_permanent_SA3b + YLD_resolve_SA3b
mean_median_ci(YLD_SA3b)

#Total DALY
DALY_SA3b <- YLD_SA3b + YLL_SA3b 
str(DALY_SA3b)
hist(DALY_SA3b)
mean_median_ci(DALY_SA3b)

#YLL contribution to DALY
YLL_contribution_SA3b <- YLL_SA3b/DALY_SA3b*100
mean_median_ci(YLL_contribution_SA3b)

#YLD contribution to DALY
YLD_contribution_SA3b <- YLD_SA3b/DALY_SA3b*100
mean_median_ci(YLD_contribution_SA3b)

YLD_permanent_contribution_SA3b <- YLD_permanent_SA3b/DALY_SA3b*100
mean_median_ci(YLD_permanent_contribution_SA3b)

####SA3c - extreme dw ####
#SA3.a - dw = DM = 0.049

#dw is the mean of mild, moderate, severe RA and Crohn's disease. 
RA_mild <- rpert(nunc, min = 0.08, mode = 0.117, max = 0.163)
RA_moderate <- rpert(nunc, min = 0.216, mode = 0.317, max = 0.440)
RA_severe <- rpert(nunc, min = 0.403, mode = 0.581, max = 0.739)
Crohns <- rpert(nunc, min = 0.156, mode = 0.213, max = 0.320)

dw_SA3c <- (RA_mild + RA_moderate + RA_severe + Crohns)/4
mean_median_ci(dw_SA3c)

aod <- rpert(nunc, min = 9, mode = 32, max = 78)
str(aod)
hist(aod)
mean_median_ci(aod)

#DALY per fatal PA case = SEYLL * 1#
DALY_fatal <- SEYLL * 1
mean_median_ci(DALY_fatal)

#DALY per resolving PA case = d * DW#
d_res <- AoR - AoO

DALY_resolve_SA3c <- d_res*dw_SA3c
mean_median_ci(DALY_resolve_SA3c)

#DALY per permanent case =  d * dw #

DALY_permanent_SA3c <- (LE_0_5)*dw_SA3c
mean_median_ci(DALY_permanent_SA3c)

#Transition probabilities#

#Mortality due to peanut allergy in Denmark #
#Pert distribution to describe the uncertainty around the mortality; min = 0 mortality, mode = swedish mortality, max = maximum identified in lit.
minimum <- 1.09/1000000 #Umasunthar 2013
mode <- 2.13/1000000 #Umasunthar 2013
maximum <- 4.6/1000000 #Umasunthar 2013

p_fatal <- rpert(nunc, minimum, mode, maximum)

p_resolve <- rpert(nunc, min = 0.10, mode = 0.22, max = 0.27) ##Peters et al. (2014)

p_permanent <- (1-(p_fatal + p_resolve)) 

#DALY #

#mean DALY per case of PA
DALY_case_SA3c = (p_fatal * DALY_fatal) + (p_resolve * DALY_resolve_SA3c) + (p_permanent * DALY_permanent_SA3c)
mean_median_ci(DALY_case_SA3c)

#YLL
YLL_SA3c = pa_inc * p_fatal * DALY_fatal
mean_median_ci(YLL_SA3c)

#YLD contribution from permanent PA
YLD_permanent_SA3c = pa_inc * p_permanent * DALY_permanent_SA3c 
mean_median_ci(YLD_permanent_SA3c)

#YLD contribution from resolving PA
YLD_resolve_SA3c = pa_inc * p_resolve * DALY_resolve_SA3c
mean_median_ci(YLD_resolve_SA3c)

#Total YLD
YLD_SA3c <- YLD_permanent_SA3c + YLD_resolve_SA3c
mean_median_ci(YLD_SA3c)

#Total DALY
DALY_SA3c <- YLD_SA3c + YLL_SA3c 
str(DALY_SA3c)
hist(DALY_SA3c)
mean_median_ci(DALY_SA3c)

#YLL contribution to DALY
YLL_contribution_SA3c <- YLL_SA3c/DALY_SA3c*100
mean_median_ci(YLL_contribution_SA3c)

#YLD contribution to DALY
YLD_contribution_SA3c <- YLD_SA3c/DALY_SA3c*100
mean_median_ci(YLD_contribution_SA3c)

YLD_permanent_contribution_SA3c <- YLD_permanent_SA3c/DALY_SA3c*100
mean_median_ci(YLD_permanent_contribution_SA3c)

##### Results plot #####
#DALY plot
DALY_mean <- mean(DALY)
DALY_LL <- quantile(DALY, probs = 0.025)
DALY_UL <- quantile(DALY, probs = 0.975)

DALY_SA1_mean <- mean(DALY_SA1)
DALY_SA1_LL <- quantile(DALY_SA1, probs = 0.025)
DALY_SA1_UL <- quantile(DALY_SA1, probs = 0.975)

DALY_SA2_mean <- mean(DALY_SA2)
DALY_SA2_LL <- quantile(DALY_SA2, probs = 0.025)
DALY_SA2_UL <- quantile(DALY_SA2, probs = 0.975)

DALY_SA3a_mean <- mean(DALY_SA3a)
DALY_SA3a_LL <- quantile(DALY_SA3a, probs = 0.025)
DALY_SA3a_UL <- quantile(DALY_SA3a, probs = 0.975)

DALY_SA3b_mean <- mean(DALY_SA3b)
DALY_SA3b_LL <- quantile(DALY_SA3b, probs = 0.025)
DALY_SA3b_UL <- quantile(DALY_SA3b, probs = 0.975)

DALY_SA3c_mean <- mean(DALY_SA3c)
DALY_SA3c_LL <- quantile(DALY_SA3c, probs = 0.025)
DALY_SA3c_UL <- quantile(DALY_SA3c, probs = 0.975)

DALY_scenarios <- rbind(DALY_mean, DALY_SA1_mean, DALY_SA2_mean, DALY_SA3a_mean, DALY_SA3b_mean, DALY_SA3c_mean)

scenarios <-  c(1:6)

DALY_scenarios <- as.data.frame(cbind(scenarios, DALY_scenarios))

DALY_scenarios_melt <- melt(DALY_scenarios, id.vars = 'scenarios')

DALY_scenarios_LL <- rbind(DALY_LL, DALY_SA1_LL, DALY_SA2_LL, DALY_SA3a_LL, DALY_SA3b_LL, DALY_SA3c_LL)

DALY_scenarios_LL <- as.data.frame(cbind(scenarios, DALY_scenarios_LL))
DALY_scenarios_LL_melt <- melt(DALY_scenarios_LL, id.vars = 'scenarios')

DALY_scenarios_UL <- rbind(DALY_UL, DALY_SA1_UL, DALY_SA2_UL, DALY_SA3a_UL, DALY_SA3b_UL, DALY_SA3c_UL)

DALY_scenarios_UL <- as.data.frame(cbind(scenarios, DALY_scenarios_UL))
DALY_scenarios_UL_melt <- melt(DALY_scenarios_UL, id.vars = 'scenarios')

DALY_scenarios_UI <- merge(DALY_scenarios_melt, DALY_scenarios_LL_melt, by = c('scenarios'))
names(DALY_scenarios_UI)[5] <- 'LL' 

DALY_scenarios_UI <- merge(DALY_scenarios_UI, DALY_scenarios_UL_melt, by = c('scenarios'))
names(DALY_scenarios_UI)[7] <- 'UL' 

names(DALY_scenarios_UI)[3] <- 'value'

DALY_scenarios_UI <- DALY_scenarios_UI[, c(1,3,5,7)]

DALY_scenarios_UI$scenarios <- as.factor(DALY_scenarios_UI$scenarios)

#tiff("Fig_DALY_PA_scenarios.tiff", width = 9, height = 6, units = "in", res = 300, compress = "lzw")
c <-  ggplot(DALY_scenarios_UI, aes(y = value, x=scenarios))

c +
  geom_bar(stat = "identity", color = 'gray30', fill='lightgrey') +
  geom_errorbar(aes(ymin = DALY_scenarios_UI$LL, ymax = DALY_scenarios_UI$UL), width=.2)+ 
  labs(x= '23', y = 'DALY') +
  scale_x_discrete(name = "", limits=c("1", "2", "3", "4", "5", "6"),
                   labels = c('Baseline', 'SA1', 'SA2', "SA3a", "SA3b", "SA3c")) +
  scale_fill_manual(name = '', 
                    values = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BF7D", "#00BFC4")) +
  # scale_y_continuous(labels = scales::comma +
  theme_bw()+
  theme(axis.ticks.x = element_blank())
#graphics.off()


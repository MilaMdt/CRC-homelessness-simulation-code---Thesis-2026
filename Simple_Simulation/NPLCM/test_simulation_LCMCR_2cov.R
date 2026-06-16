rm(list=ls())
gc()

set.seed(42)

library(LCMCR)
library(matrixStats)
library(ggplot2)
library(foreign)
library(cat)
library(haven)
library(combinat)
library(stringr)



## EM algorithm

# we create a function that will be useful later
# it inserts ":" in between the letters of the coefs vector, to include it in the formula updating
insert_colon = function(x){ paste(strsplit(x, "")[[1]], collapse=":") }



MC = 625
popSize = 20000
Nsamples = 3

results = matrix(NA, MC, 3)
colnames(results)=c("itt","n000","N_est_standard_cov")
results[,"itt"] = 1:MC


get_prob= function(age, famous,sex, sample){
  base_prob = 0.2
  
  if(sex==2 & sample==3){
    base_prob = base_prob + 0.15
  }
    
  base_prob = base_prob + 0.15*(famous==1)
  
  #if(age==2){base_prob = base_prob + 0.1}
  
  #if((age==1)|(age==3)){base_prob = base_prob - 0.1}
  
  pmax(pmin(base_prob, 1), 0)
}

for(mc in 1:MC){
  
  famous = sample(c(1,2), popSize, replace=TRUE, prob=c(0.5,0.5)) # 1 for famous
  sex    = sample(c(1,2), popSize, replace=TRUE, prob=c(0.4,0.6)) # 1 for woman, 2 for man
  age    = sample(c(1,2,3), popSize, replace=TRUE, prob=c(0.2,0.6,0.2))
  # created the population
  
  Samples = matrix(NA, popSize, Nsamples)
  inclusionProbs = matrix(NA, popSize, Nsamples)
  
  for (s in 1:Nsamples){
    inclusionProbs[,s] = mapply(get_prob,age=age, famous = famous, sex=sex, MoreArgs=list(sample=s))
    Samples[,s] = rbinom(popSize, 1, inclusionProbs[,s])
  }
  colnames(Samples) = c("G","W","Z")
  freq = 1
  dfile = cbind(Samples,age, famous, sex)
  colnames(dfile) = c("G","W","Z","A", "F","S")
  dfile = as.data.frame(dfile)
  tb = as.data.frame(table(dfile))
  for (col in c('G','W','Z')){
    tb[[col]] = as.factor(tb[[col]])
    dfile[[col]] = as.factor(dfile[[col]])
    levels(tb[[col]]) = c("2","1") # change to 1-2 factor variables, with 2 'not in list'
    levels(dfile[[col]]) = c(2,1)
  }
  
  struc.zero <- ifelse(tb[,1] == 2 & tb[,2] == 2 & tb[,3] == 2, 0, 1)
  tb	     <- cbind(tb, as.data.frame(struc.zero)) 
  colnames(tb) = c(colnames(tb)[-ncol(tb)],"struc.zero")
  
  filtered_dfile = subset(dfile, G !=2 | Z!=2 | W!=2)
  for (col in c('G','W','Z','A','F','S')){
    filtered_dfile[[col]] = as.numeric(filtered_dfile[[col]])
  }
  dmat.list    	 <- as.matrix(filtered_dfile)
  
  null_obs = sum(tb[tb[,1] == 2 & tb[,2] == 2 & tb[,3] == 2,]$Freq)
  tb[tb[,1] == 2 & tb[,2] == 2 & tb[,3] == 2,]$Freq = 0
  
  # first with sex==1 and famous == 1
  tb.s1f1 = as.data.frame(table(dfile[(dfile$`F`==1) & (dfile$S==1),1:3]))
  null_obs_s1f1 = sum(tb.s1f1[(tb.s1f1[,1] == 2) & (tb.s1f1[,2] == 2) & (tb.s1f1[,3] == 2),]$Freq)
  filtered.tb.s1f1 = subset(tb.s1f1, G !=2 | Z!=2 | W!=2)
  
  sampler.s1f1 <- lcmCR(captures = filtered.tb.s1f1,tabular=TRUE, in_list_label = '1',
                   not_in_list_label = '2', K = 10, a_alpha = 0.25, b_alpha = 0.25,
                   seed = 'auto', buffer_size = 10000, thinning = 100)
  N.s1f1 <- lcmCR_PostSampl(sampler.s1f1, burnin = 10000, samples = 1000, thinning = 100, output = FALSE)
  
  # with sex==1 and famous==2
  
  tb.s1f2 = as.data.frame(table(dfile[(dfile$`F`==2) & (dfile$S==1),1:3]))
  null_obs_s1f2 = sum(tb.s1f2[(tb.s1f2[,1] == 2) & (tb.s1f2[,2] == 2) & (tb.s1f2[,3] == 2),]$Freq)
  filtered.tb.s1f2 = subset(tb.s1f2, G !=2 | Z!=2 | W!=2)
  
  sampler.s1f2 <- lcmCR(captures = filtered.tb.s1f2,tabular=TRUE, in_list_label = '1',
                        not_in_list_label = '2', K = 10, a_alpha = 0.25, b_alpha = 0.25,
                        seed = 'auto', buffer_size = 10000, thinning = 100)
  N.s1f2 <- lcmCR_PostSampl(sampler.s1f2, burnin = 10000, samples = 1000, thinning = 100, output = FALSE)
  
  # with sex==2 and famous==1
  
  tb.s2f1 = as.data.frame(table(dfile[(dfile$`F`==1) & (dfile$S==2),1:3]))
  null_obs_s2f1 = sum(tb.s2f1[(tb.s2f1[,1] == 2) & (tb.s2f1[,2] == 2) & (tb.s2f1[,3] == 2),]$Freq)
  filtered.tb.s2f1 = subset(tb.s2f1, G !=2 | Z!=2 | W!=2)
  
  sampler.s2f1 <- lcmCR(captures = filtered.tb.s2f1,tabular=TRUE, in_list_label = '1',
                        not_in_list_label = '2', K = 10, a_alpha = 0.25, b_alpha = 0.25,
                        seed = 'auto', buffer_size = 10000, thinning = 100)
  N.s2f1 <- lcmCR_PostSampl(sampler.s2f1, burnin = 10000, samples = 1000, thinning = 100, output = FALSE)
  
  # with sex==2 and famous==2
  
  tb.s2f2 = as.data.frame(table(dfile[(dfile$`F`==2) & (dfile$S==2),1:3]))
  null_obs_s2f2 = sum(tb.s1f1[(tb.s2f2[,1] == 2) & (tb.s2f2[,2] == 2) & (tb.s2f2[,3] == 2),]$Freq)
  filtered.tb.s2f2 = subset(tb.s2f2, G !=2 | Z!=2 | W!=2)
  
  sampler.s2f2 <- lcmCR(captures = filtered.tb.s2f2,tabular=TRUE, in_list_label = '1',
                        not_in_list_label = '2', K = 10, a_alpha = 0.25, b_alpha = 0.25,
                        seed = 'auto', buffer_size = 10000, thinning = 100)
  N.s2f2 <- lcmCR_PostSampl(sampler.s2f2, burnin = 10000, samples = 1000, thinning = 100, output = FALSE)
  
  # final estimate of population size:
  
  N.LCMCR = median(N.s1f1) + median(N.s1f2) + median(N.s2f1) + median(N.s2f2)
  
  
  results[mc,2] = null_obs_s1f1 + null_obs_s1f2 + null_obs_s2f1 + null_obs_s2f2
  results[mc,3] = N.LCMCR
  
  
}

# from the first estimate we observe that N.s1f1 and N.s2f1 are good estimates
# but for s1f2 and s2f2 we get estimates that off by a few thousands : 
# for s2f2 the difference is more than 2000
# for s1f2 the difference is very similar with also more than 2000
# this is because we do not have many people who are not famous and still captured

colMeans(results[,c("n000", "N_est_standard_cov")])
colSds(  results[,c("n000", "N_est_standard_cov")])/(MC^0.5)

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

values = c(a = 0,
           g = 1,
           w = 2,
           z = 3,
           f = 4,
           s = 5)


# we create a function that will be useful later
# it inserts ":" in between the letters of the coefs vector, to include it in the formula updating
insert_colon = function(x){ paste(strsplit(x, "")[[1]], collapse=":") }

# another function for the specs based on the names
names2specs = function(x, values){
  low = str_to_lower(x[-1]) # take out the intercept
  x = Reduce(function(x,y) c(x,"a",y), low)
  letters = unlist(strsplit(x, split=""))
  specs = values[letters]
  names(specs)=NULL
  specs
}



MC = 625
popSize = 20000
Nsamples = 3

results = matrix(NA, MC, 3)
colnames(results)=c("itt","n000","N_est_standard_cov")
results[,"itt"] = 1:MC


get_prob= function(age,famous,sex, sample){
  base_prob = 0.2
  if(sex==2 & sample==3){
    base_prob = base_prob + 0.15
  }
  
  base_prob = base_prob + 0.15*(famous==1)
  
  #if(age==2){
  #  base_prob = base_prob + 0.1
  #}
  
  if((age==1)|(age==3)){
    base_prob = base_prob - 0.1
  }
  
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
  
  # first without any covariate
  tb.nocov = as.data.frame(table(dfile[,1:3]))
  null_obs_nocov = tb.nocov[(tb.nocov[,1] == 2) & (tb.nocov[,2] == 2) & (tb.nocov[,3] == 2),]$Freq
  tb.nocov[(tb.nocov[,1] == 2) & (tb.nocov[,2] == 2) & (tb.nocov[,3] == 2),]$Freq = 0
  
  
  sampler <- lcmCR(captures = tb.nocov,tabular=TRUE, in_list_label = '1',
                   not_in_list_label = '2', K = 10, a_alpha = 0.25, b_alpha = 0.25,
                   seed = 'auto', buffer_size = 10000, thinning = 100)
  # obtain 1000 samples from the posterior distribution of N#
  N_LCMCR <- lcmCR_PostSampl(sampler, burnin = 10000, samples = 1000, thinning = 100, output = FALSE)
  
  results[mc,2]=null_obs_nocov
  results[mc,3]= median(N_LCMCR)
  

}

colMeans(results[,c("n000", "N_est_standard_cov")]) 
colSds(  results[,c("n000", "N_est_standard_cov")])/(MC^0.5)
rm(list=ls())
gc()

# in this code I would like to try to delete all parameters with non-significant pvalues at once, per k-way

set.seed(345)

library(nimble)
library(MCMCvis)
library(nimbleMacros)
library(combinat)


values = c(b = 0,
           g = 1,
           w = 2,
           z = 3,
           a = 4,
           f = 5,
           s = 6)

## EM algorithm



# we create a function that will be useful later
# it inserts ":" in between the letters of the coefs vector, to include it in the formula updating
insert_colon = function(x){ paste(strsplit(x, "")[[1]], collapse=":") }

# another function for the specs based on the names
names2specs = function(x, values){
  low = str_to_lower(x[-1]) # take out the intercept
  x = Reduce(function(x,y) c(x,"b",y), low)
  letters = unlist(strsplit(x, split=""))
  specs = values[letters]
  names(specs)=NULL
  specs
}

popSize = 20000
Nsamples = 3


MC = 100
results = matrix(NA, MC, 2)
colnames(results) = c("itt","est.pop.size")
results[,1] = 1:MC

get_prob= function(age, famous,sex, sample){
  base_prob = 0.2
  
  if(sex==2 & sample==3) base_prob = base_prob + 0.15
  
  base_prob = base_prob + 0.15 *(famous==1)
  
  #base_prob = base_prob + 0.1*(age==2)
  
  #base_prob = base_prob - 0.1*(age == 1 | age==3)
  
  pmax(pmin(base_prob, 1), 0) # keep base_prob between 0 and 1
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

get_logit_prob <- function(age, famous, sex, sample) {
  lp <- -1.5  # corresponds to a basic probability 1/(1+exp(-(-1.5))) = 0.18
  
  # Covariates and sample effects
  lp <- lp + 0.3 * (sex == 1)
  lp <- lp + 0.6 * (famous == 1)
  
  lp <- lp + 0.4 * (age == 2)
  lp <- lp - 0.3 * (age == 1 | age == 3)
  
  if (sample == 3 & sex == 1) lp <- lp + 0.3
  
  return(lp)
}

#for(mc in 1:MC){

famous = sample(c(1,2), popSize, replace=TRUE, prob=c(0.4,0.6)) # 1 for famous
sex    = sample(c(1,2), popSize, replace=TRUE, prob=c(0.5,0.5)) # 1 for woman, 2 for man
age    = sample(c(1,2,3), popSize, replace=TRUE, prob=c(0.2,0.6,0.2))
# created the population

Samples = matrix(NA, popSize, Nsamples)
inclusionProbs = matrix(NA, popSize, Nsamples)

for (s in 1:Nsamples) {
  #all_ind_prob <- mapply(get_logit_prob, sex, famous, age, MoreArgs=list(sample=s))
  #inclusionProbs[, s] <- inv_logit(all_ind_prob)
  inclusionProbs[,s] = mapply(get_prob, age=age, famous=famous, sex=sex, MoreArgs=list(sample=s))
  Samples[, s] <- rbinom(popSize, 1, inclusionProbs[, s])
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

A2 = model.matrix(~A, tb)[,2]
A3 = model.matrix(~A, tb)[,3]
names(A2) = NULL
names(A3) = NULL
tb$A2 = as.factor(A2)
tb$A3 = as.factor(A3)

# PRIORS
# prior for abundance, capture/inclusion probability,
# in this code not need for a prior for model selection as estimation of only one model
# scale prior for abundance N
# beta prior for capture probabilities p_j with hyperparameters a and b
# jeffrey's prior for unobserved population sizes
# vague gamma prior for sigma square
# independent normal priors with 0 mean and sigma² for parameters in the model

# maybe we could create a distribution to get N based on the formula exp(beta0) + exp(beta0+beta1) + ... = \hat{N}

# model building
model <- nimbleCode({
  # prior for coefficients
  betaG ~ dnorm(0, tau)
  betaW ~ dnorm(0, tau)
  betaZ ~ dnorm(0, tau)
  betaF ~ dnorm(0, tau)
  betaS ~ dnorm(0, tau)
  betaA2 ~ dnorm(0, tau)
  betaA3 ~ dnorm(0, tau)
  betaZS ~ dnorm(0, tau)
  betaGF ~ dnorm(0, tau)
  betaZF ~ dnorm(0, tau)
  betaWF ~ dnorm(0, tau)
  betaFS ~ dnorm(0, tau)
  
  # we need a different one for beta0, as it is at the heart of the unobserved population estimation
  # we want for now a non-informative prior, but positive and of mean large (around 10000), as we already observe around 8000
  beta0 ~ dnorm(0, tau=0.01)
  
  # non-informative prior for tau
  tau ~ dgamma(0.01,0.01) # maybe a t-half or uniform prior is better
  
  
  n000 <- exp(beta0) + exp(beta0+betaF) + exp(beta0+betaS) + exp(beta0+betaA2) +exp(beta0+betaA3) +
                      exp(beta0+betaF+betaS+betaFS) + exp(beta0+betaF+betaA2) + exp(beta0+betaF+betaA3) + exp(beta0+betaA2+betaS) +exp(beta0+betaA3+betaS) +
                      exp(beta0+betaF+betaS+betaA2 + betaFS) + exp(beta0+betaF+betaS+betaA3+betaFS)
  
  # model
  for(i in 1:n) {
    freq[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + betaG * G[i] + betaZ * Z[i] + betaW * W[i] +
      betaF * f[i] + betaS * S[i] + betaA2 * A2[i] + betaA3 * A3[i] +
      betaZS * Z[i] * S[i] + betaGF * G[i] * f[i] + betaZF * Z[i] * f[i] + betaWF * W[i] * f[i] + betaFS * S[i] * f[i]# true model 
    # "A" had to be dummy coded because it has 3 possible values
    # because value 1 is the reference, we only need 2 dummy variables A2 and A3
    
  }
  
})
# you can either write the code manually 
# or use LM() and setPrior() from the nimbleMacros package

# specify constants and read data
tb.struc1 = tb[tb$struc.zero==1,]
my.data <- list(freq = tb.struc1$Freq) # should we include the struc.zero==0??
n = nrow(tb.struc1)
glmConstants <- list(n = n[1], G = tb.struc1[,1], W = tb.struc1[,2], Z = tb.struc1[,3],
                     f = tb.struc1[,5], S = tb.struc1[,6], A2 = tb.struc1[,9], A3 = tb.struc1[,10])
# pick initial values
initial.values <- list(beta0 = 2, betaG = 0, betaZ = 0, betaW = 0,
              betaF = 0,betaS = 0,betaA2 = 0,betaA3 = 0,
              betaZS = 0, betaGF = 0, betaZF = 0, betaWF = 0,betaFS = 0,tau = rgamma(1,0.01,0.01))

parameters.to.save = c("beta0","betaG","lambda", "n000")
# specify MCMC details
n.iter <- 13000
n.burnin <- 7000
n.chains <- 3
# run NIMBLE
mcmc.output <- nimbleMCMC(code = model,
                          data = my.data,
                          monitors=parameters.to.save,
                          inits = initial.values,
                          constants = glmConstants,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
n000 = as.numeric(MCMCsummary(object=mcmc.output, params="n000")[4])
#results[mc,2] = n000 + sum(tb$Freq)

#}
# calculate numerical summaries
MCMCsummary(object = mcmc.output, round = 2)
# visualize parameter posterior distribution
MCMCplot(object = mcmc.output, 
         params = 'n000')
MCMCplot(object = mcmc.output, 
         params = 'beta0')
# check convergence
MCMCtrace(object = mcmc.output,
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "n000")
MCMCtrace(object = mcmc.output,
          pdf = TRUE,
          ind = TRUE,
          params = "beta0")

####################
## longer version

unobserved <- nimbleModel(code = model,
                        data = my.data,
                        constants = glmConstants,
                        inits = initial.values)
Cunobserved <- compileNimble(unobserved)
# create a MCMC configuration
unobservedConf <- configureMCMC(unobserved)
# add lifespan to list of parameters to monitor
unobservedConf$addMonitors(c("n000"))
# create a MCMC function and compile it
unobservedMCMC <- buildMCMC(unobservedConf)
CunobservedMCMC <- compileNimble(unobservedMCMC, project=unobserved)
# specify MCMC details
n.iter <- 5000; n.burnin <- 1000; n.chains <- 3
# run NIMBLE
samples <- runMCMC(mcmc = CunobservedMCMC, 
                   niter = n.iter,
                   nburnin = n.burnin,
                   nchain = n.chains)

niter_ad <- 6000
CunobservedMCMC$run(niter_ad, reset = FALSE)

more_samples <- as.matrix(CunobservedMCMC$mvSamples)
samplesSummary(more_samples)

MCMCplot(object = more_samples, 
         params = 'n000')
# check convergence
MCMCtrace(object = more_samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "n000")



rm(list=ls())
gc()

set.seed(42)

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

## EM algorithm

dcat.fit <- function(MODEL, MARGINS){ 
  dcat.em         <- em.cat (dcat.list, start = struc.zero, showits=FALSE) 
  dcat.ecm        <- ecm.cat(dcat.list, start = struc.zero, margins = MARGINS, showits=FALSE, eps = 1e-7) 
  # it seems like the start argument could be taken out and the results would be the same
  # this might be because dcat.list is based on dfile and not tb, and the 222 combinations are taken out
  
  
  tb.ecm          <- tb # grouped data
  tb.ecm$Freq     <- as.numeric(dcat.ecm*dcat.list$n) # results with log-linear EM : new frequencies
  
  # this two lines are not in the current Rscript:
  #   tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==1,"Freq"] = tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==1,"Freq"] + 1/6
  #   tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==2,"Freq"] = tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==2,"Freq"] + 1/3  #new

    tb.fit.glm      <- glm(MODEL, family = poisson, subset = (struc.zero == 1), data = tb.ecm)  #original
    tb.summary      <- summary(tb.fit.glm) # summary of log-linear model based on EM log-linear imputation
    tb.model.matrix <- model.frame(MODEL, data = tb.ecm)
    tb.fitted   <- predict(tb.fit.glm, tb.model.matrix, type = "response")
    tb.ecm$est = tb.fitted
    #  tb.ecm[tb.ecm$struc.zero==1,"est"] = tb.ecm[tb.ecm$struc.zero==1,"Freq"] # not in current code
    norm.deviance   <- -2*(logpost.cat(dcat.list,dcat.ecm)-logpost.cat(dcat.list,dcat.em))
    npar            <- tb.fit.glm$df.null- tb.fit.glm$df.residual + 1
    norm.aic        <- norm.deviance + 2*npar
    d.o.f           <- tb.fit.glm$df.residual
    tbf             <- tb
    tbf$Freq        <- tb.fitted 
    tb.missed       <- xtabs(Freq ~ struc.zero, data = tbf )[1]
  
  
  return(list(
    norm.deviance = norm.deviance,
    d.o.f         = d.o.f,
    norm.aic      = norm.aic,
    tb.missed     = tb.missed, 
    tb.summary    = tb.summary,
    tb.fit.glm    = tb.fit.glm,
    tb.fitted = tb.fitted,
    tb.ecm = tb.ecm)) # table with information on the model
} 

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
colnames(results)=c("itt","n00","N_est_standard_cov")
results[,"itt"] = 1:MC


get_prob= function(sex, famous,age, sample){
  base_prob = 0.1
  
  if(sex==1 & sample==3){
    base_prob = base_prob - 0.07
  }
  
  if(famous == 1){
    base_prob = base_prob + 0.2
  }
  
  if(age==2){
    base_prob = base_prob + 0.1
  }
  
  if((age==1)|(age==3)){
    base_prob = base_prob - 0.1
  }
  
  pmax(pmin(base_prob, 1), 0) # keep base_prob between 0 and 1
}

for(mc in 1:MC){
  
  famous = sample(c(1,2), popSize, replace=TRUE, prob=c(0.2,0.8)) # 1 for famous
  sex    = sample(c(1,2), popSize, replace=TRUE, prob=c(0.4,0.6)) # 1 for woman, 2 for man
  age    = sample(c(1,2,3), popSize, replace=TRUE, prob=c(0.2,0.6,0.2))
  # created the population
  
  Samples = matrix(NA, popSize, Nsamples)
  inclusionProbs = matrix(NA, popSize, Nsamples)
  
  for (s in 1:Nsamples){
    for (i in 1:popSize){
      inclusionProbs[i,s] = get_prob(sex=sex[i], famous = famous[i], age=age[i], sample=s)
      Samples[i,s] = rbinom(1, 1, inclusionProbs[i,s])
    }
  }
  colnames(Samples) = c("G","Z","W")
  freq = 1
  dfile = cbind(Samples, famous, sex, age)
  colnames(dfile) = c("G","Z","W", "F","S","A")
  dfile = as.data.frame(dfile)
  tb = as.data.frame(table(dfile))
  for (col in c('G','Z','W')){
    tb[[col]] = as.factor(tb[[col]])
    dfile[[col]] = as.factor(dfile[[col]])
    levels(tb[[col]]) = c("2","1") # change to 1-2 factor variables, with 2 'not in list'
    levels(dfile[[col]]) = c(2,1)
  }
  
  struc.zero <- ifelse(tb[,1] == 2 & tb[,2] == 2 & tb[,3] == 2, 0, 1)
  tb	     <- cbind(tb, as.data.frame(struc.zero)) 
  colnames(tb) = c(colnames(tb)[-ncol(tb)],"struc.zero")
  
  filtered_dfile = subset(dfile, G !=2 | Z!=2 | W!=2)
  for (col in c('G','Z','W','F','S','A')){
    filtered_dfile[[col]] = as.numeric(filtered_dfile[[col]])
  }
  dmat.list    	 <- as.matrix(filtered_dfile)
  
  null_obs = sum(tb[tb[,1] == 2 & tb[,2] == 2 & tb[,3] == 2,]$Freq)
  tb[tb[,1] == 2 & tb[,2] == 2 & tb[,3] == 2,]$Freq = 0
  
  dcat.list	        <- prelim.cat(dmat.list) 
  struc.zero	    	<- array(struc.zero, dim = dcat.list$d)/sum(struc.zero) 
  
  test.model <- list()
  test.specs <- list()
  test.names <- list()
  
  
  # need to change the start model to the one with all terms up to 4-way
  # should we take out 0 cell counts?
  
  # main effects
  me = c("G","Z","W","F","S","A")
  # interactions in the start model
  twoe   = unlist(lapply(combn(me,2,simplify=FALSE),paste0, collapse = ""))
  threee = unlist(lapply(combn(me,3,simplify=FALSE),paste0, collapse = ""))
  foure  = unlist(lapply(combn(me,4,simplify=FALSE),paste0, collapse = ""))
  
  all_coef    = c(me,twoe,threee,foure)
  all_formula = sapply(all_coef, insert_colon)
  all_model   = paste(all_formula, collapse="+")
  all_specs   = names2specs(all_coef, values=values)
  
  start.model = as.formula(paste("Freq~",all_model))
  start.specs = all_specs
  start.names = paste(all_coef, collapse=" ")
  
  start.FIT <- dcat.fit(start.model, start.specs)
  start.pars 		<- round(start.FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)
  
  start.pars
  pars.names = rownames(start.pars)
  # first let's transform this list into a list of interactions
  start.coef = unique(gsub("[^a-zA-Z]","",pars.names)) # weirdly here the coefficient SLH is not in the vector but SLBH is
  
  
  test.model[[1]] = start.model
  test.specs[[1]] = names2specs(start.coef[-1], values=values) # unfortunately we do not get quite the same models but I do not know why
  test.names[[1]] = paste(start.coef[-1], collapse=" ")
  
  # now we want to focus on the 4-ways parameters
  coef4 = start.coef[nchar(start.coef)==4]
  formula4 = sapply(coef4, insert_colon)
  if(length(coef4)==0){
    best4.model = test.model[[1]]
    best4.specs = test.specs[[1]]
    best4.names = test.names[[1]]
  }
  if(length(coef4)!=0){
    test.pars4 = list()
    test.coef4 = list()
    i=1
    while (length(coef4)!=0){
      test.pars4[[i]] = round(dcat.fit(test.model[[i]], test.specs[[i]])$tb.summary$coefficients[,c(1,2,4)],
                              digits=2)
      coefs = unique(gsub("[^a-zA-Z]","",rownames(test.pars4[[i]])))
      test.coef4[[i]] = coefs[nchar(coefs)==4]
      coef4 = test.coef4[[i]]
      no4 = length(coefs)-length(test.coef4[[i]])
      pvalues4 = test.pars4[[i]][which(coefs %in% coef4),3]
      out = which.max(pvalues4)
      if(pvalues4[out]>0.05){
        new.coef = coefs[-(out+no4)]
        # create new model
        test.model[[i+1]] = as.formula(paste("Freq~",paste(sapply(new.coef[-1], insert_colon), collapse="+")))
        test.specs[[i+1]] = names2specs(new.coef, values=values)
        test.names[[i+1]] = paste(new.coef[-1], collapse=" ")
      } else{
        new.coef = coefs
        best4.model = test.model[[i]]
        best4.specs = test.specs[[i]]
        best4.names = test.names[[i]]
        break
      }
      coef4 = new.coef[nchar(new.coef)==4]
      best4.model = test.model[[i+1]]
      best4.specs = test.specs[[i+1]]
      best4.names = test.names[[i+1]]
      i=i+1
    }
    
    }
  
  
  # we want to keep all interactions that are included in a higher term kept in the model
  
  # what are the 4-way terms still in the model?
  if(length(coef4)==0) {
    coef3 = new.coef[nchar(new.coef)==3]
    out_test = 0
  } # well then we're good
  
  # if not, we can make a list of interactions not to take out
  # if GZBL is in the model for example, we want all 3-way interactions from this term to stay in the model
  # such as GZL, GZB, but also GLZ and ZGL, because we do not know in what order they are in the model names
  allcomb = function(coefs, n=3){
    letters = unlist(strsplit(coefs, split=""))
    combs = combn(letters, n, simplify=TRUE) # all combinations per 4-way term
    perms = unlist(apply(combs,2, function(x){
      lapply(permn(x),paste0, collapse = "") # all permutations from each combination
    }))
    return(perms)
  }
  
  if(length(coef4)>=1){
    exclude = unlist(lapply(coef4, allcomb))
    coef3 = new.coef[nchar(new.coef)==3]
    out_test = which(coef3 %in% exclude)
    coef3 = coef3[-out]
  }
  
  # now we can apply the same algorithm to the 3-way interactions
  # you still need to adapt it to out_test
  
  test3.model <- list()
  test3.specs <- list()
  test3.names <- list()
  
  formula3 = sapply(coef3, insert_colon)
  
  best4.FIT <- dcat.fit(best4.model, best4.specs)
  test3.model[[1]]= best4.model
  test3.specs[[1]]= best4.specs
  test3.names[[1]]= best4.names
  
  best4.coef = unlist(strsplit(best4.names, split=" "))
  
  # while the preferred one is a new one, should stop if it is the last one created
  if(length(coef3)==0){
    best3.model = test3.model[[1]]
    best3.specs = test3.specs[[1]]
    best3.names = test3.names[[1]]
  }
  if(length(coef3)!=0){
    test.pars3 = list()
    test.coef3 = list()
    i=1
    while (length(coef3)!=0){
      test.pars3[[i]] = round(dcat.fit(test3.model[[i]], test3.specs[[i]])$tb.summary$coefficients[,c(1,2,4)],
                              digits=2)
      coefs = unique(gsub("[^a-zA-Z]","",rownames(test.pars3[[i]])))
      test.coef3[[i]] = coefs[nchar(coefs)==3]
      #coef3 = test.coef3[[i]][which(test.coef3[[i]] %in% coef3)]
      pvalues3 = test.pars3[[i]][which(coefs %in% coef3),3]
      out = which.max(pvalues3)
      if(pvalues3[out]>0.05){
        no3 = which(test.coef3[[i]] == coef3[out])
        new.coef = c(coefs[nchar(coefs)<=2], test.coef3[[i]][-no3],best4.coef[nchar(best4.coef)>3])
        # create new model
        test3.model[[i+1]] = as.formula(paste("Freq~",paste(sapply(new.coef, insert_colon), collapse="+")))
        test3.specs[[i+1]] = names2specs(new.coef, values=values)
        test3.names[[i+1]] = paste(new.coef, collapse=" ")
      } else{
        new.coef = coefs
        best3.model = test3.model[[i]]
        best3.specs = test3.specs[[i]]
        best3.names = test3.names[[i]]
        break
      }
      coef3 = coef3[-out]
      best3.model = test3.model[[i+1]]
      best3.specs = test3.specs[[i+1]]
      best3.names = test3.names[[i+1]]
      i=i+1
    }
    
  }
  
  
  if(length(coef3)==0) {
    coef2 = new.coef[nchar(new.coef)==2]
    out_test = 0
  } # keeping all 2-way interactions that were in the last version of newcoef
  # they should also be the 2-way interactions that were in start.coef and in best4.coef
  
  
  if(length(coef3)>=1){
    exclude = unlist(lapply(coef3, allcomb, n=2))
    coef2 = new.coef[nchar(new.coef)==2]
    out_test = which(coef2 %in% exclude)
    coef2 = coef2[-out_test]
  }
  
  # 2-way interactions
  
  
  test2.model <- list()
  test2.specs <- list()
  test2.names <- list()
  
  best3.FIT <- dcat.fit(best3.model, best3.specs)
  test2.model[[1]]= best3.model
  test2.specs[[1]]= best3.specs
  test2.names[[1]]= best3.names
  
  best3.coef = unlist(strsplit(best3.names, split=" "))
  
  # while the preferred one is a new one, should stop if it is the last one created
  if(length(coef2)==0){
    best2.model = test2.model[[1]]
    best2.specs = test2.specs[[1]]
    best2.names = test2.names[[1]]
  }
  if(length(coef2)!=0){
    test.pars2 = list()
    test.coef2 = list()
    i=1
    while (length(coef2)!=0){
      test.pars2[[i]] = round(dcat.fit(test2.model[[i]], test2.specs[[i]])$tb.summary$coefficients[,c(1,2,4)],
                              digits=2)
      coefs = unique(gsub("[^a-zA-Z]","",rownames(test.pars2[[i]])))
      test.coef2[[i]] = coefs[nchar(coefs)==2]
      pvalues2 = test.pars2[[i]][which(coefs %in% coef2),3]
      out = which.max(pvalues2)
      if(pvalues2[out]>0.05){
        no2 = which(test.coef2[[i]] == coef2[out])
        new.coef = c(coefs[nchar(coefs)<=1], test.coef2[[i]][-no2],best3.coef[nchar(best3.coef)>=3])
        # create new model
        test2.model[[i+1]] = as.formula(paste("Freq~",paste(sapply(new.coef, insert_colon), collapse="+")))
        test2.specs[[i+1]] = names2specs(new.coef, values=values)
        test2.names[[i+1]] = paste(new.coef, collapse=" ")
      } else{
        new.coef = coefs
        best2.model = test2.model[[i]]
        best2.specs = test2.specs[[i]]
        best2.names = test2.names[[i]]
        break
      }
      coef2 = coef2[-out]
      best2.model = test2.model[[i+1]]
      best2.specs = test2.specs[[i+1]]
      best2.names = test2.names[[i+1]]
      i=i+1
    }
    
  }
  
  # best model should have all main effects of G,W,Z and S,L,B,H.
  # the algorithm ends here and the best model is :
  best2.model
  
  best2.FIT <- dcat.fit(best2.model, best2.specs)
  best2.pars = round(best2.FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)
  best2.freqs = cbind(tb[,1:length(tb)-1],round(best2.FIT$tb.fitted,2))
  coefs.origineel = best2.FIT$tb.summary$coefficients
  tb.ecm.origineel = best2.FIT$tb.ecm
  
  fitted = best2.freqs
  colnames(fitted)	= c("G","W","Z","F","S","A","obs","est")
  
  est.pop.size = sum(fitted$est)
  obs.pop.size = sum(fitted$obs)
  
  dark.figure = est.pop.size-obs.pop.size 
  results[mc,2]=dark.figure
  results[mc,3]=est.pop.size
  
}

colMeans(results[,c("n00", "N_est_standard_cov")])
colSds(  results[,c("n00", "N_est_standard_cov")])/(MC^0.5)

# we get 10026.843, with only famous and sex, and base_prob=0.2.
# with also age and base_prob = 0.1, we get 9663.891, which is quite further (by more than 300)
# with sd 22.21862, with only famous and sex.
# with also age, sd is 44.39


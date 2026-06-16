# Backward BIC 

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
           b = 4,
           f = 5,
           s = 6)

a = 0
g = 1
w = 2
z = 3
b = 4
f = 5
s = 6

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
  n               <- nrow(tb.ecm)
  norm.bic        <- norm.deviance + (log(n))*npar
  d.o.f           <- tb.fit.glm$df.residual
  tbf             <- tb
  tbf$Freq        <- tb.fitted 
  tb.missed       <- xtabs(Freq ~ struc.zero, data = tbf )[1]
  
  
  return(list(
    norm.deviance = norm.deviance,
    d.o.f         = d.o.f,
    norm.bic      = norm.bic,
    tb.missed     = tb.missed, 
    tb.summary    = tb.summary,
    tb.fit.glm    = tb.fit.glm,
    tb.fitted = tb.fitted,
    tb.ecm = tb.ecm)) # table with information on the model
} 

# we create a function that will be useful later
# inserts ":" in between the letters of the coefs vector, to include it in the formula updating
insert_colon = function(x){ paste(strsplit(x, "")[[1]], collapse=":") }
insert_star  = function(x){ paste(strsplit(x, "")[[1]], collapse="*") }

# get the specs based on the names
names2specs = function(x, values){
  if(length(x)>1){
    low = str_to_lower(x) 
    vec = Reduce(function(x,y) c(x,"a",y), low)
    letters = unlist(strsplit(vec, split=""))
    specs = values[letters]
    names(specs)=NULL
  }
  if(length(x)==1){
    low = str_to_lower(x)
    letters = unlist(strsplit(low, split=""))
    specs = values[letters]
    names(specs)=NULL
  }
  return(specs)
}

model2specs = function(formula, values){
  # right side of the formula
  sum = sub(".*~", "", formula)[3]
  # ":" taken out of interactions
  sum2 = gsub(":", "", sum)
  # + replaced by a, whose value is 0
  vec = gsub("\\+", "a", sum2)
  low = str_to_lower(vec)
  # transform into vector of letters
  low2 = unlist(strsplit(low, " "))
  letters = unlist(strsplit(low2, ""))
  specs = values[letters]
  names(specs) = NULL
  specs
}

model2names = function(formula){
  # right side of the formula
  sum = sub(".*~", "", formula)[3]
  # remove all white space
  sum = gsub(" ", "", sum)
  # ":" taken out of interactions
  sum2 = gsub(":", "", sum)
  names_long = gsub("\\+", " ", sum2)
  # remove the eventual \n if the formula way too long...
  names = gsub(" \n", " ", names_long)
  names
}

allcomb = function(coefs, n=3){
  letters = unlist(strsplit(coefs, split=""))
  combs = combn(letters, n, simplify=TRUE) # all combinations per 4-way term
  perms = unlist(apply(combs,2, function(x){
    lapply(permn(x),paste0, collapse = "") # all permutations from each combination
  }))
  return(perms)
}


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
  
  famous = sample(c(1,2), popSize, replace=TRUE, prob=c(0.4,0.6)) # 1 for famous
  sex    = sample(c(1,2), popSize, replace=TRUE, prob=c(0.5,0.5)) # 1 for woman, 2 for man
  age    = sample(c(1,2,3), popSize, replace=TRUE, prob=c(0.2,0.6,0.2))
  # created the population
  
  Samples = matrix(NA, popSize, Nsamples)
  inclusionProbs = matrix(NA, popSize, Nsamples)
  
  for (s in 1:Nsamples){
      inclusionProbs[,s] = mapply(get_prob,age=age, famous = famous,sex= sex, MoreArgs=list(sample=s))
      Samples[,s] = rbinom(popSize, 1, inclusionProbs[,s])
  }
  colnames(Samples) = c("G","W","Z")
  freq = 1
  dfile = cbind(Samples,age, famous, sex)
  colnames(dfile) = c("G","W","Z","B", "F","S")
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
  for (col in c('G','W','Z','B','F','S')){
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
  me = c("G","W","Z","B","F","S")
  # interactions in the start model
  twoe   = unlist(lapply(combn(me,2,simplify=FALSE),paste0, collapse = ""))
  threee = unlist(lapply(combn(me,3,simplify=FALSE),paste0, collapse = ""))
  foure  = unlist(lapply(combn(me,4,simplify=FALSE),paste0, collapse = ""))
  
  # we do not want GZW to be in the 3rd and 4th orders
  threee = threee[-which(threee == "GWZ")]
  foure = foure[-which(foure %in% c("GWZB", "GWZS", "GWZF"))]
  
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
    BIC = numeric(length(coef4)+1)
    BIC[1]=start.FIT$norm.bic
    min.BIC = 1
    # while the preferred one is a new one, should stop if it is the last one created
    for (i in 1:length(coef4)){
      test.model[[i+1]]= update.formula(test.model[[1]],paste("~.-",formula4[i]))
      test.specs[[i+1]]= names2specs(start.coef[-c(1,which(start.coef==coef4[i]))], values=values)
      test.names[[i+1]] = model2names(test.model[[i+1]])
      BIC[i+1]=dcat.fit(test.model[[i+1]], test.specs[[i+1]])$norm.bic
    }
    if(which.min(BIC)==1){
      best4.model = test.model[[1]]
      best4.specs = test.specs[[1]]
      best4.names = test.names[[1]]
    }else{
      min.BIC = which.min(BIC)
      best4.model = test.model[[min.BIC]]
      best4.specs = test.specs[[min.BIC]]
      best4.names = test.names[[min.BIC]]
      newcoef = unlist(strsplit(test.names[[min.BIC]], split=" "))
      coef4 = newcoef[nchar(newcoef)==4]
      formula4 = sapply(coef4, insert_colon)
      while(length(coef4)>0){
        j=length(test.model)
        min.BIC = which.min(BIC)
        for (i in 1:length(coef4)){
          j=j+1
          test.model[[j]]= update.formula(test.model[[min.BIC]],paste("~.-",formula4[i]))
          test.specs[[j]]= names2specs(newcoef[-which(newcoef==coef4[i])], values=values)
          test.names[[j]] = model2names(test.model[[j]])
          BIC[j]=dcat.fit(test.model[[j]], test.specs[[j]])$norm.bic
        }
        if (min.BIC == which.min(BIC)){break}
        best4.model = test.model[[which.min(BIC)]]
        best4.specs = test.specs[[which.min(BIC)]]
        best4.names = test.names[[which.min(BIC)]]
        newcoef = unlist(strsplit(best4.names, split=" "))
        coef4 = newcoef[nchar(newcoef)==4]
        formula4 = sapply(coef4, insert_colon)
      }
    }
  }
  
  # we want to keep all interactions that are included in a higher term kept in the model
  
  # what are the 4-way terms still in the model?
  if(length(coef4)==0) {
    coef3 = start.coef[nchar(start.coef)==3]
    out_test = 0
  } # well then we're good
  
  # if not, we can make a list of interactions not to take out
  # if GZBL is in the model for example, we want all 3-way interactions from this term to stay in the model
  # such as GZL, GZB, but also GLZ and ZGL, because we do not know in what order they are in the model names
  
  
  if(length(coef4)>=1){
    exclude = unlist(lapply(coef4, allcomb))
    coef3 = start.coef[nchar(start.coef)==3]
    out_test = which(coef3 %in% exclude)
    coef3 = coef3[-out_test]
  }
  
  # now we can apply the same algorithm to the 3-way interactions
  
  test3.model <- list()
  test3.specs <- list()
  test3.names <- list()
  
  formula3 = sapply(coef3, insert_colon)
  BIC3 = numeric(length(coef3)+1)
  min.BIC3 = 1
  
  best4.FIT <- dcat.fit(best4.model, best4.specs)
  test3.model[[1]]= best4.model
  test3.specs[[1]]= best4.specs
  test3.names[[1]]= best4.names
  
  best4.coef = unlist(strsplit(best4.names, split=" "))
  
  BIC3[1]=best4.FIT$norm.bic
  # while the preferred one is a new one, should stop if it is the last one created
  if(length(coef3)==0){
    best3.model = test3.model[[1]]
    best3.specs = test3.specs[[1]]
    best3.names = test3.names[[1]]
  }
  if(length(coef3)!=0){
    for (i in 1:length(coef3)){
      test3.model[[i+1]]= update.formula(test3.model[[1]],paste("~.-",formula3[i]))
      test3.specs[[i+1]]= names2specs(best4.coef[-which(best4.coef==coef3[i])], values=values)
      test3.names[[i+1]]= model2names(test3.model[[i+1]])
      BIC3[i+1]=dcat.fit(test3.model[[i+1]], test3.specs[[i+1]])$norm.bic
    }
    if(which.min(BIC3)==1){
      best3.model = test3.model[[1]]
      best3.specs = test3.specs[[1]]
      best3.names = test3.names[[1]]
      newcoef = unlist(strsplit(test3.names[[1]], split=" "))
    }else{
      min.BIC3 = which.min(BIC3)
      best3.model = test3.model[[min.BIC3]]
      best3.specs = test3.specs[[min.BIC3]]
      best3.names = test3.names[[min.BIC3]]
      newcoef = unlist(strsplit(test3.names[[min.BIC3]], split=" "))
      coef3 = newcoef[nchar(newcoef)==3]
      if(sum(out_test)!=0){coef3 = coef3[-which(coef3 %in% exclude)]}
      formula3 = sapply(coef3, insert_colon)
      last.min = 0
      while(length(coef3)>0){
        j=length(test3.model)
        min.BIC3 = which.min(BIC3)
        for (i in 1:length(coef3)){
          j=j+1
          test3.model[[j]]= update.formula(best3.model,paste("~.-",formula3[i]))
          test3.specs[[j]]= model2specs(test3.model[[j]], values=values)
          test3.names[[j]] = model2names(test3.model[[j]])
          BIC3[j]=dcat.fit(test3.model[[j]], test3.specs[[j]])$norm.bic
        }
        BIC3 = round(BIC3, 2)
        if ((min.BIC3 == which.min(BIC3))&(length(which(BIC3== BIC3[min.BIC3]))==1)){break} # if after the loop the best model is still the one we had before the loop,
        # then it means that we want to keep this model
        if ((min.BIC3 == which.min(BIC3))&(length(which(BIC3== BIC3[min.BIC3]))>1)){
          if(last.min == which(BIC3== BIC3[min.BIC3])[length(which(BIC3 == BIC3[min.BIC3]))]){break}
          best3.model = test3.model[[which(BIC3== BIC3[min.BIC3])[length(which(BIC3 == BIC3[min.BIC3]))]]]
          best3.specs = test3.specs[[which(BIC3== BIC3[min.BIC3])[length(which(BIC3 == BIC3[min.BIC3]))]]]
          best3.names = test3.names[[which(BIC3== BIC3[min.BIC3])[length(which(BIC3 == BIC3[min.BIC3]))]]]
          newcoef = unlist(strsplit(best3.names, split=" "))
          coef3 = newcoef[nchar(newcoef)==3]
          if(sum(out_test)!=0){coef3 = coef3[-which(coef3 %in% exclude)]}
          formula3 = sapply(coef3, insert_colon)
          last.min = which(BIC3== BIC3[min.BIC3])[length(which(BIC3 == BIC3[min.BIC3]))]
        } else {
          best3.model = test3.model[[which.min(BIC3)]]
          best3.specs = test3.specs[[which.min(BIC3)]]
          best3.names = test3.names[[which.min(BIC3)]]
          newcoef = unlist(strsplit(best3.names, split=" "))
          coef3 = newcoef[nchar(newcoef)==3]
          if(sum(out_test)!=0){coef3 = coef3[-which(coef3 %in% exclude)]}
          formula3 = sapply(coef3, insert_colon)
        }
      }
    }
  }
  coef3 = newcoef[nchar(newcoef)==3]
  if(length(coef3)==0) {
    coef2 = newcoef[nchar(newcoef)==2]
    out_test = 0
  } # keeping all 2-way interactions that were in the last version of newcoef
  # they should also be the 2-way interactions that were in start.coef and in best4.coef
  
  
  if(length(coef3)>=1){
    exclude = unlist(lapply(coef3, allcomb, n=2))
    coef2 = newcoef[nchar(newcoef)==2]
    out_test = which(coef2 %in% exclude)
    coef2 = coef2[-out_test]
  }
  
  # 2-way interactions
  
  
  test2.model <- list()
  test2.specs <- list()
  test2.names <- list()
  
  formula2 = sapply(coef2, insert_colon)
  BIC2 = numeric(length(coef2)+1)
  
  best3.FIT <- dcat.fit(best3.model, best3.specs)
  test2.model[[1]]= best3.model
  test2.specs[[1]]= best3.specs
  test2.names[[1]]= best3.names
  
  best3.coef = unlist(strsplit(best3.names, split=" "))
  
  BIC2[1]=best3.FIT$norm.bic
  min.BIC2 = 1
  if(length(coef2)==0){
    best2.model = test2.model[[1]]
    best2.specs = test2.specs[[1]]
    best2.names = test2.names[[1]]
  }
  if(length(coef2)!=0){
    # while the preferred one is a new one, should stop if it is the last one created
    for (i in 1:length(coef2)){
      test2.model[[i+1]]= update.formula(test2.model[[1]],paste("~.-",formula2[i]))
      test2.specs[[i+1]]= names2specs(best3.coef[-which(best3.coef==coef2[i])], values=values)
      test2.names[[i+1]]= model2names(test2.model[[i+1]])
      BIC2[i+1]=dcat.fit(test2.model[[i+1]], test2.specs[[i+1]])$norm.bic
    }
    if(which.min(BIC2)==1){
      best2.model = test2.model[[1]]
      best2.specs = test2.specs[[1]]
      best2.names = test2.names[[1]]
    }else{
      min.BIC2 = which.min(BIC2)
      best2.model = test2.model[[min.BIC2]]
      best2.specs = test2.specs[[min.BIC2]]
      best2.names = test2.names[[min.BIC2]]
      newcoef = unlist(strsplit(test2.names[[min.BIC2]], split=" "))
      coef2 = newcoef[nchar(newcoef)==2]
      if(sum(out_test)!=0){coef2 = coef2[-which(coef2 %in% exclude)]} # keeping only interactions that we want out
      formula2 = sapply(coef2, insert_colon)
      last.min = 0
      while(length(coef2)>0){
        j=length(test2.model)
        min.BIC2 = which.min(BIC2)
        for (i in 1:length(coef2)){
          j=j+1
          test2.model[[j]]= update.formula(best2.model,paste("~.-",formula2[i]))
          test2.specs[[j]]= model2specs(test2.model[[j]], values = values)
          test2.names[[j]] = model2names(test2.model[[j]])
          BIC2[j]=dcat.fit(test2.model[[j]], test2.specs[[j]])$norm.bic
        }
        BIC2 = round(BIC2,2)
        if ((min.BIC2 == which.min(BIC2))&(length(which(BIC2== BIC2[min.BIC2]))==1)){break} # if after the loop the best model is still the one we had before the loop,
        # then it means that we want to keep this model
        if ((min.BIC2 == which.min(BIC2))&(length(which(BIC2== BIC2[min.BIC2]))>1)){
          if(last.min == which(BIC2== BIC2[min.BIC2])[length(which(BIC2 == BIC2[min.BIC2]))])
          best2.model = test2.model[[which(BIC2== BIC2[min.BIC2])[length(which(BIC2 == BIC2[min.BIC2]))]]]
          best2.specs = test2.specs[[which(BIC2== BIC2[min.BIC2])[length(which(BIC2 == BIC2[min.BIC2]))]]]
          best2.names = test2.names[[which(BIC2== BIC2[min.BIC2])[length(which(BIC2 == BIC2[min.BIC2]))]]]
          newcoef = unlist(strsplit(best2.names, split=" "))
          coef2 = newcoef[nchar(newcoef)==2]
          if(sum(out_test)!=0){coef2 = coef2[-which(coef2 %in% exclude)]}
          formula2 = sapply(coef2, insert_colon)
          last.min = which(BIC2== BIC2[min.BIC2])[length(which(BIC2 == BIC2[min.BIC2]))]
        } else {
          best2.model = test2.model[[which.min(BIC2)]]
          best2.specs = test2.specs[[which.min(BIC2)]]
          best2.names = test2.names[[which.min(BIC2)]]
          newcoef = unlist(strsplit(best2.names, split=" "))
          coef2 = newcoef[nchar(newcoef)==2]
          if(sum(out_test)!=0){coef2 = coef2[-which(coef2 %in% exclude)]}
          formula2 = sapply(coef2, insert_colon)
        }
      }
    }
  }
  # best model should have all main effects of G,W,Z and S,L,B,H.
  # the algorithm ends here and the best model is :
  best2.model
  BIC2[min.BIC2] 
  
  best2.FIT <- dcat.fit(best2.model, best2.specs)
  best2.pars = round(best2.FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)
  best2.freqs = cbind(tb[,1:length(tb)-1],round(best2.FIT$tb.fitted,2))
  coefs.origineel = best2.FIT$tb.summary$coefficients
  tb.ecm.origineel = best2.FIT$tb.ecm
  
  fitted = best2.freqs
  colnames(fitted)	= c("G","W","Z","B","F","S","obs","est")
  
  est.pop.size = sum(fitted$est)
  obs.pop.size = sum(fitted$obs)
  
  dark.figure = est.pop.size-obs.pop.size 
  results[mc,2]=null_obs
  results[mc,3]=est.pop.size
  
}

colMeans(results[,c("n000", "N_est_standard_cov")]) 
colSds(  results[,c("n000", "N_est_standard_cov")])/(MC^0.5)

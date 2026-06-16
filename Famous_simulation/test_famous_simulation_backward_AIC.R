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
library(matlib)

values = c(a = 0,
           g = 1,
           w = 2,
           z = 3,
           b = 4,
           f = 5,
           s = 6)

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

results = matrix(NA, MC, 5)
colnames(results)=c("itt","n000","est.pop.size", "est.chap.size","truen000")
results[,"itt"] = 1:MC


get_prob= function(age, famous,sex, sample){
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
  
  pmax(pmin(base_prob, 1), 0)
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
      inclusionProbs[i,s] = get_prob(age=age[i], famous = famous[i],sex= sex[i], sample=s)
      Samples[i,s] = rbinom(1, 1, inclusionProbs[i,s])
    }
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
  AIC = numeric(length(coef4)+1)
  AIC[1]=start.FIT$norm.aic
  min.AIC = which.min(AIC)
  # while the preferred one is a new one, should stop if it is the last one created
  for (i in 1:length(coef4)){
    test.model[[i+1]]= update.formula(test.model[[1]],paste("~.-",formula4[i]))
    test.specs[[i+1]]= names2specs(start.coef[-c(1,which(start.coef==coef4[i]))], values=values)
    test.names[[i+1]] = paste(start.coef[-c(1,which(start.coef == coef4[i]))], collapse=" ")
    AIC[i+1]=dcat.fit(test.model[[i+1]], test.specs[[i+1]])$norm.aic
  }
  if(which.min(AIC)==1){
    best4.model = test.model[[1]]
    best4.specs = test.specs[[1]]
    best4.names = test.names[[1]]
  }else{
    min.AIC = which.min(AIC)
    best4.model = test.model[[min.AIC]]
    best4.specs = test.specs[[min.AIC]]
    best4.names = test.names[[min.AIC]]
    newcoef = unlist(strsplit(test.names[[min.AIC]], split=" "))
    coef4 = newcoef[nchar(newcoef)==4]
    formula4 = sapply(coef4, insert_colon)
    while(length(coef4)>0){
      j=length(test.model)
      min.AIC = which.min(AIC)
      for (i in 1:length(coef4)){
        j=j+1
        test.model[[j]]= update.formula(test.model[[min.AIC]],paste("~.-",formula4[i]))
        test.specs[[j]]= names2specs(newcoef[-which(newcoef==coef4[i])], values=values)
        test.names[[j]] = paste(newcoef[-which(newcoef == coef4[i])], collapse=" ")
        AIC[j]=dcat.fit(test.model[[j]], test.specs[[j]])$norm.aic
      }
      if (min.AIC == which.min(AIC)){break}
      best4.model = test.model[[which.min(AIC)]]
      best4.specs = test.specs[[which.min(AIC)]]
      best4.names = test.names[[which.min(AIC)]]
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
    out = 0
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
    coef3 = start.coef[nchar(start.coef)==3]
    out = which(coef3 %in% exclude)
    coef3 = coef3[-out]
  }
  
  # now we can apply the same algorithm to the 3-way interactions
  
  test3.model <- list()
  test3.specs <- list()
  test3.names <- list()
  
  formula3 = sapply(coef3, insert_colon)
  AIC3 = numeric(length(coef3)+1)
  min.AIC = which.min(AIC3)
  
  best4.FIT <- dcat.fit(best4.model, best4.specs)
  test3.model[[1]]= best4.model
  test3.specs[[1]]= best4.specs
  test3.names[[1]]= best4.names
  
  best4.coef = unlist(strsplit(best4.names, split=" "))
  
  AIC3[1]=best4.FIT$norm.aic
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
    test3.names[[i+1]]= paste(best4.coef[-which(best4.coef == coef3[i])], collapse=" ")
    AIC3[i+1]=dcat.fit(test3.model[[i+1]], test3.specs[[i+1]])$norm.aic
  }
  if(which.min(AIC3)==1){
    best3.model = test3.model[[1]]
    best3.specs = test3.specs[[1]]
    best3.names = test3.names[[1]]
    newcoef = unlist(strsplit(test3.names[[min.AIC]], split=" "))
  }else{
    min.AIC = which.min(AIC3)
    best3.model = test3.model[[min.AIC]]
    best3.specs = test3.specs[[min.AIC]]
    best3.names = test3.names[[min.AIC]]
    newcoef = unlist(strsplit(test3.names[[min.AIC]], split=" "))
    coef3 = newcoef[nchar(newcoef)==3]
    if(sum(out)!=0){coef3 = coef3[-out]}
    formula3 = sapply(coef3, insert_colon)
    while(length(coef3)>0){
      j=length(test3.model)
      min.AIC = which.min(AIC3)
      for (i in 1:length(coef3)){
        j=j+1
        test3.model[[j]]= update.formula(test3.model[[min.AIC]],paste("~.-",formula3[i]))
        test3.specs[[j]]= names2specs(newcoef[-which(newcoef==coef3[i])], values=values)
        test3.names[[j]] = paste(newcoef[-which(newcoef == coef3[i])], collapse=" ")
        AIC3[j]=dcat.fit(test3.model[[j]], test3.specs[[j]])$norm.aic
      }
      if (min.AIC == which.min(AIC3)){break} # if after the loop the best model is still the one we had before the loop,
      # then it means that we want to keep this model
      best3.model = test3.model[[which.min(AIC3)]]
      best3.specs = test3.specs[[which.min(AIC3)]]
      best3.names = test3.names[[which.min(AIC3)]]
      newcoef = unlist(strsplit(best3.names, split=" "))
      coef3 = newcoef[nchar(newcoef)==3]
      if(sum(out)!=0){coef3 = coef3[-out]}
      formula3 = sapply(coef3, insert_colon)
    }
  }
  }
  
  if(length(coef3)==0) {
    coef2 = newcoef[nchar(newcoef)==2]
    out = 0
  } # keeping all 2-way interactions that were in the last version of newcoef
  # they should also be the 2-way interactions that were in start.coef and in best4.coef
  
  
  if(length(coef3)>=1){
    exclude = unlist(lapply(coef3, allcomb, n=2))
    coef2 = newcoef[nchar(newcoef)==2]
    out = which(coef2 %in% exclude)
    coef2 = coef2[-out]
  }
  
  # 2-way interactions
  
  
  test2.model <- list()
  test2.specs <- list()
  test2.names <- list()
  
  formula2 = sapply(coef2, insert_colon)
  AIC2 = numeric(length(coef2)+1)
  
  best3.FIT <- dcat.fit(best3.model, best3.specs)
  test2.model[[1]]= best3.model
  test2.specs[[1]]= best3.specs
  test2.names[[1]]= best3.names
  
  best3.coef = unlist(strsplit(best3.names, split=" "))
  
  AIC2[1]=best3.FIT$norm.aic
  min.AIC = which.min(AIC2)
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
    test2.names[[i+1]]= paste(best3.coef[-which(best3.coef == coef2[i])], collapse=" ")
    AIC2[i+1]=dcat.fit(test2.model[[i+1]], test2.specs[[i+1]])$norm.aic
  }
  if(which.min(AIC2)==1){
    best2.model = test2.model[[1]]
    best2.specs = test2.specs[[1]]
    best2.names = test2.names[[1]]
  }else{
    min.AIC = which.min(AIC2)
    best2.model = test2.model[[min.AIC]]
    best2.specs = test2.specs[[min.AIC]]
    best2.names = test2.names[[min.AIC]]
    newcoef = unlist(strsplit(test2.names[[min.AIC]], split=" "))
    coef2 = newcoef[nchar(newcoef)==2]
    if(sum(out)!=0){coef2 = coef2[-out]} # keeping only interactions that we want out
    formula2 = sapply(coef2, insert_colon)
    while(length(coef2)>0){
      j=length(test2.model)
      min.AIC = which.min(AIC2)
      for (i in 1:length(coef2)){
        j=j+1
        test2.model[[j]]= update.formula(test2.model[[min.AIC]],paste("~.-",formula2[i]))
        test2.specs[[j]]= names2specs(newcoef[-which(newcoef==coef2[i])], values=values)
        test2.names[[j]] = paste(newcoef[-which(newcoef == coef2[i])], collapse=" ")
        AIC2[j]=dcat.fit(test2.model[[j]], test2.specs[[j]])$norm.aic
      }
      if (min.AIC == which.min(AIC2)){break} # if after the loop the best model is still the one we had before the loop,
      # then it means that we want to keep this model
      best2.model = test2.model[[which.min(AIC2)]]
      best2.specs = test2.specs[[which.min(AIC2)]]
      best2.names = test2.names[[which.min(AIC2)]]
      newcoef = unlist(strsplit(best2.names, split=" "))
      coef2 = newcoef[nchar(newcoef)==2]
      if(sum(out)!=0){coef2 = coef2[-out]}
      formula2 = sapply(coef2, insert_colon)
    }
  }
  }
  # best model should have all main effects of G,W,Z and S,L,B,H.
  # the algorithm ends here and the best model is :
  best2.model
  AIC2[min.AIC] # in 2023 we get a model with an AIC of around 565.98
  
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
  results[mc,2]=dark.figure
  results[mc,3]=est.pop.size
  
  formula = as.character(best2.model[3])
  formula_X = strsplit(formula, " ")[[1]]
  formula_X = formula_X[formula_X != "" & formula_X != "\n"]
  formula_X = formula_X[!grepl("S", formula_X)] # take out terms with S
  formula_X = formula_X[!grepl("B", formula_X)] # take out terms with B
  formula_X = formula_X[!grepl("F", formula_X)] # F
  formula_X = formula_X[!grepl("+", formula_X, fixed = TRUE)] 
  formula_X = formula_X[!grepl("*", formula_X, fixed = TRUE)] 
  # only things left: G,W,Z (but could be more if we had interactions of G,W,Z)
  formula_X = gsub("G:W", "I(G*W)", formula_X)
  formula_X = gsub("G:Z", "I(G*Z)", formula_X)
  formula_X = gsub("Z:W", "I(Z*W)", formula_X)
  formula_X = gsub("W:G", "I(W*G)", formula_X)
  formula_X = gsub("Z:G", "I(Z*G)", formula_X)
  formula_X = gsub("W:Z", "I(W*Z)", formula_X)
  formula_X = paste0("Freq~", paste0(c(formula_X), collapse="+"), collapse="")
  
  l <- rep(list(0:1), 3)
  X = as.data.frame(cbind(rep(1, 7),as.matrix(expand.grid(l))[-1,]))
  colnames(X) = c("Freq","G","W","Z")
  X # all combinations of G,Z,W except for 0,0,0, and all Freq set to 1
  X = as.matrix(model.frame(as.formula(formula_X), data = X))
  
  A = Inverse(t(X)%*%X)%*%t(X)
  a_ijk_Chapman = A[1,]
  a_ijk_Chapman[a_ijk_Chapman>0] = 0
  
  tb_Chapman = tb
  rows = ifelse(X[which(a_ijk_Chapman<0),] == 1, 1, 2)
  if (length(which(a_ijk_Chapman<0))==1){
    row = rows
    tb_Chapman[(tb_Chapman$G == row["G"])& (tb_Chapman$W == row["W"]) & (tb_Chapman$Z == row["Z"]),"Freq"] = tb_Chapman[(tb_Chapman$G == row["G"])& (tb_Chapman$W == row["W"]) & (tb_Chapman$Z == row["Z"]), "Freq"]- a_ijk_Chapman[which(a_ijk_Chapman<0)]
  } else{
    for (num in 1:dim(rows)[1]){
      row = rows[num,]
      tb_Chapman[(tb_Chapman$G == row["G"])& (tb_Chapman$W == row["W"]) & (tb_Chapman$Z == row["Z"]),"Freq"] = tb_Chapman[(tb_Chapman$G == row["G"])& (tb_Chapman$W == row["W"]) & (tb_Chapman$Z == row["Z"]), "Freq"]- a_ijk_Chapman[which(a_ijk_Chapman<0)][num] # no need to specify struc.zero!=0 because we took it out from X
    }
  }
  
  
  chap.fit.glm      <- glm(best2.model, family = poisson, subset = (struc.zero == 1), data = tb_Chapman)  #original
  chap.summary      <- summary(chap.fit.glm) # summary of log-linear model based on EM log-linear imputation
  chap.model.matrix <- model.frame(best2.model, data = tb_Chapman)
  chap.fitted   <- predict(chap.fit.glm, chap.model.matrix, type = "response") 
  best.chap.pars = round(chap.summary$coefficients[,c(1,2,4)],digits=2)
  best.chap.freqs = cbind(tb_Chapman[,1:length(tb)-1],round(chap.fitted,2))
  
  chap.fitted = best.chap.freqs
  colnames(chap.fitted)	= c("G","W","Z","B","F","S","obs","est")
  
  est.chap.size = sum(chap.fitted$est)
  results[mc,4] = est.chap.size
  results[mc,5] = null_obs
  
}

mean(results[,3])
sd(results[,3])/(MC^0.5)
mean(results[,4])
sd(results[,4])/(MC^0.5)

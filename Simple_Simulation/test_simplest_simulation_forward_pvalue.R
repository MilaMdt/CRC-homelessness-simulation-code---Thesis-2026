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

model2names = function(formula, values){
  # right side of the formula
  sum = sub(".*~", "", formula)[3]
  # remove all white space
  sum = gsub(" ", "", sum)
  # ":" taken out of interactions
  sum2 = gsub(":", "", sum)
  names = gsub("\\+", " ", sum2)
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

smallest_p = function(pvalues){
  if(length(pvalues)>1){
    names(pvalues) = gsub("[^a-zA-Z]", "", names(pvalues))
    pv = numeric()
    for (i in 1:length(unique(names(pvalues)))){
      pv[i] = pvalues[which(names(pvalues)==unique(names(pvalues))[i])][which.min(pvalues[which(names(pvalues)==unique(names(pvalues))[i])])]
      names(pv)[i] = insert_colon(unique(names(pvalues))[i])
    }
    return(pv)
  } else {
    pv = pvalues
    return(pv)
  }
}

three_way = function(coef, best2 = best2.coef){
  let = unlist(strsplit(coef, split=""))
  coef.1 = sapply(permn(c(let[1], let[2])), paste0, collapse="")
  coef.2 = sapply(permn(c(let[1], let[3])), paste0, collapse="")
  coef.3 = sapply(permn(c(let[2], let[3])), paste0, collapse="")
  if(sum(c(coef.1,coef.2,coef.3) %in% best2)==3) return(coef)
}

four_way = function(coef, best3 = best3.coef){
  let = unlist(strsplit(coef, split=""))
  coef.1 = sapply(permn(c(let[1], let[2], let[3])), paste0, collapse="")
  coef.2 = sapply(permn(c(let[1], let[3], let[4])), paste0, collapse="")
  coef.3 = sapply(permn(c(let[2], let[3], let[4])), paste0, collapse="")
  coef.4 = sapply(permn(c(let[1], let[2], let[4])), paste0, collapse="")
  if(sum(c(coef.1,coef.2,coef.3, coef.4) %in% best3)==4) return(coef)
}

MC = 625
popSize = 20000
Nsamples = 3

results = matrix(NA, MC, 5)
colnames(results)=c("itt","n000","est.pop.size","model", "chapman")
results[,"itt"] = 1:MC


get_prob= function(age, famous,sex, sample){
  base_prob = 0.2
  
  if(sex==1 & sample==3){base_prob = base_prob + 0.15}
  
  base_prob = base_prob + 0.15*(famous==1)
  
  #if(age==2){base_prob = base_prob + 0.1}
  
  #if((age==1)|(age==3)){base_prob = base_prob - 0.1}
  
  pmax(pmin(base_prob, 1), 0)
}

for(mc in 1:MC) {
  
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
  
  # main effects
  me = c("G","W","Z","B","F","S")
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
  
  start.model = as.formula("Freq~G+W+Z+B+F+S")
  start.specs = c(g,a,w,a,z,a,b,a,f,a,s)
  start.names = c("G W Z B F S")
  
  start.FIT <- dcat.fit(start.model, start.specs)
  start.pars 		<- round(start.FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)
  
  start.pars
  pars.names = rownames(start.pars)
  # first let's transform this list into a list of interactions
  start.coef = unique(gsub("[^a-zA-Z]","",pars.names)) 
  
  
  test.model[[1]] = start.model
  test.specs[[1]] = names2specs(start.coef[-1], values=values) 
  test.names[[1]] = paste(start.coef[-1], collapse=" ")
  
  # now we want to focus on the 2-ways parameters
  coef2 = twoe
  formula2 = sapply(coef2, insert_colon)
  pvalues = numeric()
  pvalues[1] = 1 # just a number for the start model, so that we can give it the same indices as 
  test.pars2 = list()
  test.pars2[[1]] = start.pars
  test.coef2 = list()
  test.coef2[[1]] = "None"
  # while the preferred one is a new one, should stop if it is the last one created
  for (i in 1:length(coef2)){
    test.model[[i+1]]= update.formula(test.model[[1]],paste("~.+",formula2[i]))
    test.specs[[i+1]]= c(test.specs[[1]],0,names2specs(coef2[i], values=values))
    test.names[[i+1]] = paste(c(start.coef[-1],coef2[i]), collapse=" ") 
    test.pars2[[i+1]] = round(dcat.fit(test.model[[i+1]], test.specs[[i+1]])$tb.summary$coefficients[,c(1,2,4)],
                            digits=2)
    coefs = unique(gsub("[^a-zA-Z]","", rownames(test.pars2[[i+1]])))
    test.coef2[[i+1]] = coefs[nchar(coefs)==2] # there should be either one or two (two if B is included as B2, B3)
    pv = test.pars2[[i+1]][which(gsub("[^a-zA-Z]", "", rownames(test.pars2[[i+1]])) %in% test.coef2[[i+1]]),3]
    pvalues[i+1]=smallest_p(pvalues=pv)
  }
  names(pvalues) = c("None",coef2)
  if(pvalues[which.min(pvalues)]>=0.001){ # no pvalues are significant
    best2.model = test.model[[1]]
    best2.specs = test.specs[[1]]
    best2.names = test.names[[1]]
    newcoef = unlist(strsplit(best2.names, split=" "))
  }else{
    min.pvalue = which.min(pvalues) # we only include one coefficient: the one with the lowest pvalue
    in.coef2 = rownames(test.pars2[[min.pvalue]])[which(nchar(rownames(test.pars2[[min.pvalue]]))==5)]
    best2.model = test.model[[min.pvalue]]
    best2.specs = test.specs[[min.pvalue]]
    best2.names = test.names[[min.pvalue]]
    newcoef = unlist(strsplit(test.names[[min.pvalue]], split=" "))
    coef2 = coef2[-which(coef2 %in% newcoef[nchar(newcoef)==2])]
    formula2 = sapply(coef2, insert_colon)
    while(length(coef2)>0){
      j=length(test.model) # to add new models
      l=length(test.model) # to keep track of link between pvalue indices and model indices
      pvalues = numeric()
      for (i in 1:length(coef2)){
        j=j+1
        test.model[[j]]= update.formula(test.model[[min.pvalue]],paste("~.+",formula2[i]))
        test.specs[[j]]= c(test.specs[[min.pvalue]],0,names2specs(coef2[i], values=values))
        test.names[[j]] = paste(c(newcoef,coef2[i]), collapse=" ")
        test.pars2[[j]] = round(dcat.fit(test.model[[j]], test.specs[[j]])$tb.summary$coefficients[,c(1,2,4)],
                                digits=2)
        coefs = unique(gsub("[^a-zA-Z]","", rownames(test.pars2[[j]])))
        test.coef2[[j]] =coefs[nchar(coefs)==2] 
        # but we do not want the ones that have been added before
        test.coef2[[j]] = test.coef2[[j]][-which(test.coef2[[j]] %in% gsub("[^a-zA-Z]", "", in.coef2))]
        pv = test.pars2[[j]][which(gsub("[^a-zA-Z]", "", rownames(test.pars2[[j]])) %in% test.coef2[[j]]),3]
        pvalues[i]=smallest_p(pvalues=pv)
      }
      names(pvalues) = coef2
      min.pvalue = which.min(pvalues) +l
      if (pvalues[which.min(pvalues)] >= 0.001){break}
      best2.model = test.model[[which.min(pvalues) + l]]
      best2.specs = test.specs[[which.min(pvalues) + l]]
      best2.names = test.names[[which.min(pvalues) + l]]
      newcoef = unlist(strsplit(best2.names, split=" "))
      coef2 = coef2[-which(coef2 %in% newcoef[nchar(newcoef)==2])]
      formula2 = sapply(coef2, insert_colon)
      in.coef2 = rownames(test.pars2[[min.pvalue]])[which(nchar(rownames(test.pars2[[min.pvalue]]))==5)]
    }
  }
  # we want to keep all interactions that are included in a higher term kept in the model
  
  # if 2 or less 2-way interaction included
  # no need to look into 3-way interactions because of hierarchy principle : no 3-way interaction can be formed
  if(length(coef2)>=(length(twoe)-2)){
    best.model = best2.model
    best.specs = best2.specs
    best.names = best2.names
  } 
  
  
  if(length(coef2)<(length(twoe)-2)){
    
    test3.model <- list()
    test3.specs <- list()
    test3.names <- list()
    
    best2.FIT <- dcat.fit(best2.model, best2.specs)
    test3.model[[1]]= best2.model
    test3.specs[[1]]= best2.specs
    test3.names[[1]]= best2.names
    
    best2.coef = unlist(strsplit(best2.names, split=" "))
    
    # only consider the 3-way terms derived from newcoef[nchar(newcoef)==2]
    coef3 = unlist(sapply(threee, three_way, best2=best2.coef))
    if(length(coef3)==0){ 
      best.model = best2.model
      best.specs = best2.specs
      best.names = best2.names
    } else {
    demon = which(coef3 %in% c("GZW", "GWZ", "WZG","WGZ","ZGW","ZWG"))
    if (length(demon) !=0) { coef3 = coef3[-demon]}
    possible_coef3 = coef3
    if(length(coef3)==0){ 
      best.model = best2.model
      best.specs = best2.specs
      best.names = best2.names
    }
    if(length(coef3)>0){
    formula3 = sapply(coef3, insert_colon)
    
    pvalues = numeric()
    pvalues[1] = 1 # just a number for the start model, so that we can give it the same indices as 
    test.pars3 = list()
    test.pars3[[1]] = round(best2.FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)
    test.coef3 = list()
    test.coef3[[1]] = "None"
    # while the preferred one is a new one, should stop if it is the last one created
    for (i in 1:length(coef3)){
      test3.model[[i+1]]= update.formula(test3.model[[1]],paste("~.+",formula3[i]))
      test3.specs[[i+1]]= c(test3.specs[[1]],0,names2specs(coef3[i], values=values))
      test3.names[[i+1]] = paste(c(best2.coef,coef3[i]), collapse=" ") # to test, could be paste(c(start.coef[-1],coef2[i]), collapse=" ") instead
      test.pars3[[i+1]] = round(dcat.fit(test3.model[[i+1]], test3.specs[[i+1]])$tb.summary$coefficients[,c(1,2,4)],
                                digits=2)
      coefs = unique(gsub("[^a-zA-Z]","", rownames(test.pars3[[i+1]])))
      test.coef3[[i+1]] = coefs[nchar(coefs)==3] # there should be either one or two (two if B is included as B2, B3)
      pv = test.pars3[[i+1]][which(gsub("[^a-zA-Z]", "", rownames(test.pars3[[i+1]])) %in% test.coef3[[i+1]]),3]
      pvalues[i+1]=smallest_p(pvalues=pv)
    }
    names(pvalues) = c("None",coef3)
    if(pvalues[which.min(pvalues)]>=0.001){ # no pvalues are significant
      best3.model = test3.model[[1]]
      best3.specs = test3.specs[[1]]
      best3.names = test3.names[[1]]
      newcoef = unlist(strsplit(best3.names, split=" "))
    }else{
      min.pvalue3 = which.min(pvalues) # we only include one coefficient: the one with the lowest pvalue
      in.coef3 = rownames(test.pars3[[min.pvalue3]])[which(nchar(rownames(test.pars3[[min.pvalue3]]))==8)]
      best3.model = test3.model[[min.pvalue3]]
      best3.specs = test3.specs[[min.pvalue3]]
      best3.names = test3.names[[min.pvalue3]]
      newcoef = unlist(strsplit(test3.names[[min.pvalue3]], split=" "))
      coef3 = coef3[-which(coef3 %in% newcoef[nchar(newcoef)==3])]
      formula3 = sapply(coef3, insert_colon)
      while(length(coef3)>0){
        j=length(test3.model) # to add new models
        l=length(test3.model) # to keep track of link between pvalue indices and model indices
        pvalues = numeric()
        for (i in 1:length(coef3)){
          j=j+1
          test3.model[[j]]= update.formula(test3.model[[min.pvalue3]],paste("~.+",formula3[i]))
          test3.specs[[j]]= c(test3.specs[[min.pvalue3]],0,names2specs(coef3[i], values=values))
          test3.names[[j]] = paste(c(newcoef,coef3[i]), collapse=" ")
          test.pars3[[j]] = round(dcat.fit(test3.model[[j]], test3.specs[[j]])$tb.summary$coefficients[,c(1,2,4)],
                                  digits=2)
          coefs = unique(gsub("[^a-zA-Z]","", rownames(test.pars3[[j]])))
          test.coef3[[j]] =coefs[nchar(coefs)==3] 
          # but we do not want the ones that have been added before
          test.coef3[[j]] = test.coef3[[j]][-which(test.coef3[[j]] %in% gsub("[^a-zA-Z]", "", in.coef3))]
          pv = test.pars3[[j]][which(gsub("[^a-zA-Z]", "", rownames(test.pars3[[j]])) %in% test.coef3[[j]]),3]
          pvalues[i]=smallest_p(pvalues=pv)
        }
        names(pvalues) = coef3
        min.pvalue3 = which.min(pvalues) +l
        if (pvalues[which.min(pvalues)] >= 0.001){break}
        best3.model = test3.model[[which.min(pvalues) + l]]
        best3.specs = test3.specs[[which.min(pvalues) + l]]
        best3.names = test3.names[[which.min(pvalues) + l]]
        newcoef = unlist(strsplit(best3.names, split=" "))
        coef3 = coef3[-which(coef3 %in% newcoef[nchar(newcoef)==3])]
        formula3 = sapply(coef3, insert_colon)
        in.coef3 = rownames(test.pars3[[min.pvalue3]])[which(nchar(rownames(test.pars3[[min.pvalue3]]))==8)]
      }
      }
    }
    
    if(length(coef3)>=(length(possible_coef3)-3)){ # need at least 4 three way to get a four way
      best.model = best3.model
      best.specs = best3.specs
      best.names = best3.names
    }
    
    
    if(length(coef3)<(length(possible_coef3)-3)){
      
      test4.model <- list()
      test4.specs <- list()
      test4.names <- list()
      
      best3.FIT <- dcat.fit(best3.model, best3.specs)
      test4.model[[1]]= best3.model
      test4.specs[[1]]= best3.specs
      test4.names[[1]]= best3.names
      
      best3.coef = unlist(strsplit(best3.names, split=" "))
      
      coef4 = unlist(sapply(foure, four_way))
      if(length(coef4) == 0){
        best.model = test4.model[[1]]
        best.specs = test4.specs[[1]]
        best.names = test4.names[[1]]
      } else{
      bigdemon = sapply(coef4, function(x) all(c("G", "Z", "W") %in% strsplit(x, "")[[1]]))
      # we try all combinations just in case: we don't want to miss any demon
      if (sum(bigdemon) !=0) { coef4 = coef4[-which(bigdemon)]}
      if(length(coef4) == 0){
        best.model = test4.model[[1]]
        best.specs = test4.specs[[1]]
        best.names = test4.names[[1]]
      }
      if(length(coef4)>0){
      formula4 = sapply(coef4, insert_colon)
      
      pvalues = numeric()
      pvalues[1] = 1 # just a number for the start model, so that we can give it the same indices as 
      test.pars4 = list()
      test.pars4[[1]] = round(best3.FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)
      test.coef4 = list()
      test.coef4[[1]] = "None"
      # while the preferred one is a new one, should stop if it is the last one created
      for (i in 1:length(coef4)){
        test4.model[[i+1]]= update.formula(test4.model[[1]],paste("~.+",formula4[i]))
        test4.specs[[i+1]]= c(test4.specs[[1]],0,names2specs(coef4[i], values=values))
        test4.names[[i+1]] = paste(c(best3.coef,coef4[i]), collapse=" ") # to test, could be paste(c(start.coef[-1],coef2[i]), collapse=" ") instead
        test.pars4[[i+1]] = round(dcat.fit(test4.model[[i+1]], test4.specs[[i+1]])$tb.summary$coefficients[,c(1,2,4)],
                                  digits=2)
        coefs = unique(gsub("[^a-zA-Z]","", rownames(test.pars4[[i+1]])))
        test.coef4[[i+1]] = coefs[nchar(coefs)==4] 
        pv = test.pars4[[i+1]][which(gsub("[^a-zA-Z]", "", rownames(test.pars4[[i+1]])) %in% test.coef4[[i+1]]),3]
        pvalues[i+1]=smallest_p(pvalues=pv)
      }
      names(pvalues) = c("None",coef4)
      if(pvalues[which.min(pvalues)]>=0.001){ # no pvalues are significant
        best4.model = test4.model[[1]]
        best4.specs = test4.specs[[1]]
        best4.names = test4.names[[1]]
        newcoef = unlist(strsplit(best4.names, split=" "))
      }else{
        min.pvalue4 = which.min(pvalues) # we only include one coefficient: the one with the lowest pvalue
        in.coef4 = rownames(test.pars4[[min.pvalue4]])[which(nchar(rownames(test.pars4[[min.pvalue4]]))==11)]
        best4.model = test4.model[[min.pvalue4]]
        best4.specs = test4.specs[[min.pvalue4]]
        best4.names = test4.names[[min.pvalue4]]
        newcoef = unlist(strsplit(test4.names[[min.pvalue4]], split=" "))
        coef4 = coef4[-which(coef4 %in% newcoef[nchar(newcoef)==4])]
        formula4 = sapply(coef4, insert_colon)
        while(length(coef4)>0){
          j=length(test4.model) # to add new models
          l=length(test4.model) # to keep track of link between pvalue indices and model indices
          pvalues = numeric()
          for (i in 1:length(coef4)){
            j=j+1
            test4.model[[j]]= update.formula(test4.model[[min.pvalue4]],paste("~.+",formula4[i]))
            test4.specs[[j]]= c(test4.specs[[min.pvalue4]],0,names2specs(coef4[i], values=values))
            test4.names[[j]] = paste(c(newcoef,coef4[i]), collapse=" ")
            test.pars4[[j]] = round(dcat.fit(test4.model[[j]], test4.specs[[j]])$tb.summary$coefficients[,c(1,2,4)],
                                    digits=2)
            coefs = unique(gsub("[^a-zA-Z]","", rownames(test.pars4[[j]])))
            test.coef4[[j]] =coefs[nchar(coefs)==4] 
            # but we do not want the ones that have been added before
            test.coef4[[j]] = test.coef4[[j]][-which(test.coef4[[j]] %in% gsub("[^a-zA-Z]", "", in.coef4))]
            pv = test.pars4[[j]][which(gsub("[^a-zA-Z]", "", rownames(test.pars4[[j]])) %in% test.coef4[[j]]),3]
            pvalues[i]=smallest_p(pvalues=pv)
          }
          names(pvalues) = coef4
          min.pvalue4 = which.min(pvalues) +l
          if (pvalues[which.min(pvalues)] >= 0.001){break}
          best4.model = test4.model[[which.min(pvalues) + l]]
          best4.specs = test4.specs[[which.min(pvalues) + l]]
          best4.names = test4.names[[which.min(pvalues) + l]]
          newcoef = unlist(strsplit(best4.names, split=" "))
          coef4 = coef4[-which(coef4 %in% newcoef[nchar(newcoef)==4])]
          formula4 = sapply(coef4, insert_colon)
          in.coef4 = rownames(test.pars4[[min.pvalue4]])[which(nchar(rownames(test.pars4[[min.pvalue4]]))==11)]
        } 
        }
      }
      
      best.model = best4.model
      best.specs = best4.specs
      best.names = best4.names 
      }
      }
    }
    }
  
  # best model should have all main effects of G,W,Z and S,L,B,H.
  # the algorithm ends here and the best model is :
  best.model 
  
  best.FIT <- dcat.fit(best.model, best.specs)
  best.pars = round(best.FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)
  best.freqs = cbind(tb[,1:length(tb)-1],round(best.FIT$tb.fitted,2))
  coefs.origineel = best.FIT$tb.summary$coefficients
  tb.ecm.origineel = best.FIT$tb.ecm
  
  fitted = best.freqs
  colnames(fitted)	= c("G","W","Z","B","F","S","obs","est")
  
  est.pop.size = sum(fitted$est)
  obs.pop.size = sum(fitted$obs)
  
  # chapman version of this model
  dcat.em         <- em.cat (dcat.list, start = struc.zero, showits=FALSE) 
  dcat.ecm        <- ecm.cat(dcat.list, start = struc.zero, margins = best.specs, showits=FALSE, eps = 1e-7) 
  # it seems like the start argument could be taken out and the results would be the same
  # this might be because dcat.list is based on dfile and not tb, and the 222 combinations are taken out
  
  
  tb.ecm          <- tb # grouped data
  tb.ecm$Freq     <- as.numeric(dcat.ecm*dcat.list$n) # results with log-linear EM : new frequencies
  
  
  formula_X = strsplit(as.character(best.model[3]), " ")[[1]]
  formula_X = formula_X[!grepl("S", formula_X)] # take out terms with S
  formula_X = formula_X[!grepl("F", formula_X)] # take out terms with F
  formula_X = formula_X[!grepl("B", formula_X)] # take out terms with B
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
  # X = as.matrix(X) : is this really useful??
  
  A = Inverse(t(X)%*%X)%*%t(X)
  a_ijk_Chapman = A[1,]
  a_ijk_Chapman[a_ijk_Chapman>0] = 0
  
  tb.ecm_Chapman = tb.ecm
  rows = ifelse(X[which(a_ijk_Chapman<0),] == 1, 1, 2)
  if (length(which(a_ijk_Chapman<0))==1){
    row = rows
    tb.ecm_Chapman[(tb.ecm$G == row["G"])& (tb.ecm$W == row["W"]) & (tb.ecm$Z == row["Z"]),"Freq"] = tb.ecm_Chapman[(tb.ecm$G == row["G"])& (tb.ecm$W == row["W"]) & (tb.ecm$Z == row["Z"]), "Freq"]- a_ijk_Chapman[which(a_ijk_Chapman<0)]
  } else{
    for (num in 1:dim(rows)[1]){
      row = rows[num,]
      tb.ecm_Chapman[(tb.ecm$G == row["G"])& (tb.ecm$W == row["W"]) & (tb.ecm$Z == row["Z"]),"Freq"] = tb.ecm_Chapman[(tb.ecm$G == row["G"])& (tb.ecm$W == row["W"]) & (tb.ecm$Z == row["Z"]), "Freq"]- a_ijk_Chapman[which(a_ijk_Chapman<0)][num] # no need to specify struc.zero!=0 because we took it out from X
    }
  }
  
  tb.fit.glm      <- glm(best.model, family=poisson(link="log"), subset=(struc.zero==1), data=tb.ecm_Chapman) 
  tb.summary      <- summary(tb.fit.glm) # summary of log-linear model based on EM log-linear imputation
  tb.model.matrix <- model.frame(best.model, data = tb.ecm)
  tb.fitted   <- predict(tb.fit.glm, tb.model.matrix, type = "response")
  tb.ecm$est = tb.fitted
  
  best.pars.chapman = round(tb.summary$coefficients[,c(1,2,4)],digits=2)
  best.freqs.chapman = cbind(tb[,1:length(tb)-1],round(tb.fitted,2))
  fitted.chapman = best.freqs.chapman
  colnames(fitted.chapman)	= c("G","W","Z","B","F","S","obs","est")
  
  est.pop.size.chapman = sum(fitted.chapman$est)
  obs.pop.size.chapman = sum(fitted.chapman$obs)
  
  dark.figure = est.pop.size-obs.pop.size 
  results[mc,2] = null_obs
  results[mc,3] = est.pop.size
  results[mc,4] = best.names
  results[mc,5] = est.pop.size.chapman
  
}

rownames(results) = c("it","n000", "est.pop.size","model")
results = as.data.frame(t(results))

colMeans(results[,c("n000", "est.pop.size")]) # 19,238.836
colSds(  as.matrix(results[,c("n000", "est.pop.size")]))/(MC^0.5) # 10.657742

# clearly biased 
table(as.vector(results[,4]))
mean(as.numeric(results[,3]))
sd(as.numeric(results[,3]))/(MC^0.5)
length(which(as.vector(results[,4])=="G W Z B F S GF WF ZF ZS FS ZFS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS WF GF ZF FS ZSF"))
length(which(as.vector(results[,4])=="G W Z B F S ZS ZF GF WF FS ZSF"))
length(which(as.vector(results[,4])=="G W Z B F S ZS GF WF ZF FS ZSF"))
length(which(as.vector(results[,4])=="G W Z B F S ZS ZF GF WF FS ZSF"))
length(which(as.vector(results[,4])=="G W Z B F S ZS ZF WF GF FS ZSF"))

length(which(as.vector(results[,4])=="G W Z B F S WF GF ZF ZS FS"))
length(which(as.vector(results[,4])=="G W Z B F S ZF GF WF ZS FS"))
length(which(as.vector(results[,4])=="G W Z B F S GF GS WF WS ZF ZS FS"))
length(which(as.vector(results[,4])=="G W Z B F S GF ZF WF ZS FS"))
length(which(as.vector(results[,4])=="G W Z B F S ZF GF WF ZS FS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS ZF WF GF FS"))

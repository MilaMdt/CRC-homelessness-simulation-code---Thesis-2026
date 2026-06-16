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


MC = 625
popSize = 20000
Nsamples = 3

results = matrix(NA, MC, 4)
colnames(results)=c("itt","n000","est.pop.size","model")
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

for(mc in 1:MC) {
  
  famous = sample(c(1,2), popSize, replace=TRUE, prob=c(0.2,0.8)) # 1 for famous
  sex    = sample(c(1,2), popSize, replace=TRUE, prob=c(0.4,0.6)) # 1 for woman, 2 for man
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
  BIC = numeric(length(coef2)+1)
  BIC[1]=start.FIT$norm.bic
  min.BIC = 1
  # while the preferred one is a new one, should stop if it is the last one created
  for (i in 1:length(coef2)){
    test.model[[i+1]]= update.formula(test.model[[1]],paste("~.+",formula2[i]))
    test.specs[[i+1]]= c(test.specs[[1]],0,names2specs(coef2[i], values=values))
    test.names[[i+1]] = paste(c(start.coef[-1],coef2[i]), collapse=" ") 
    BIC[i+1]=dcat.fit(test.model[[i+1]], test.specs[[i+1]])$norm.bic
  }
  if(which.min(BIC)==1){
    best2.model = test.model[[1]]
    best2.specs = test.specs[[1]]
    best2.names = test.names[[1]]
    newcoef = unlist(strsplit(best2.names, split=" "))
  }else{
    min.BIC = which.min(BIC)
    best2.model = test.model[[min.BIC]]
    best2.specs = test.specs[[min.BIC]]
    best2.names = test.names[[min.BIC]]
    newcoef = unlist(strsplit(test.names[[min.BIC]], split=" "))
    coef2 = coef2[-which(coef2 %in% newcoef[nchar(newcoef)==2])]
    formula2 = sapply(coef2, insert_colon)
    while(length(coef2)>0){
      j=length(test.model)
      min.BIC = which.min(BIC)
      for (i in 1:length(coef2)){
        j=j+1
        test.model[[j]]= update.formula(test.model[[min.BIC]],paste("~.+",formula2[i]))
        test.specs[[j]]= c(test.specs[[min.BIC]],0,names2specs(coef2[i], values=values))
        test.names[[j]] = paste(c(newcoef,coef2[i]), collapse=" ")
        BIC[j]=dcat.fit(test.model[[j]], test.specs[[j]])$norm.bic
      }
      if (min.BIC == which.min(BIC)){break}
      best2.model = test.model[[which.min(BIC)]]
      best2.specs = test.specs[[which.min(BIC)]]
      best2.names = test.names[[which.min(BIC)]]
      newcoef = unlist(strsplit(best2.names, split=" "))
      coef2 = coef2[-which(coef2 %in% newcoef[nchar(newcoef)==2])]
      formula2 = sapply(coef2, insert_colon)
    }
  }
  # we want to keep all interactions that are included in a higher term kept in the model
  
  # if 1 or less 2-way interaction included
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
        BIC3 = numeric(length(coef3)+1)
        BIC3[1]=best2.FIT$norm.bic
        min.BIC3 = 1
        # while the preferred one is a new one, should stop if it is the last one created
        for (i in 1:length(coef3)){
          test3.model[[i+1]]= update.formula(test3.model[[1]],paste("~.+",formula3[i]))
          test3.specs[[i+1]]= c(test3.specs[[1]],0,names2specs(coef3[i], values=values))
          test3.names[[i+1]] = paste(c(best2.coef,coef3[i]), collapse=" ") # to test, could be paste(c(start.coef[-1],coef2[i]), collapse=" ") instead
          BIC3[i+1]=dcat.fit(test3.model[[i+1]], test3.specs[[i+1]])$norm.bic
        }
        if(which.min(BIC3)==1){
          best3.model = test3.model[[1]]
          best3.specs = test3.specs[[1]]
          best3.names = test3.names[[1]]
          newcoef = unlist(strsplit(best3.names, split=" "))
        }else{
          min.BIC3 = which.min(BIC3)
          best3.model = test3.model[[min.BIC3]]
          best3.specs = test3.specs[[min.BIC3]]
          best3.names = test3.names[[min.BIC3]]
          newcoef = unlist(strsplit(test3.names[[min.BIC3]], split=" "))
          coef3 = coef3[-which(coef3 %in% newcoef[nchar(newcoef)==3])]
          formula3 = sapply(coef3, insert_colon)
          while(length(coef3)>0){
            j=length(test3.model)
            min.BIC3 = which.min(BIC3)
            for (i in 1:length(coef3)){
              j=j+1
              test3.model[[j]]= update.formula(test3.model[[min.BIC3]],paste("~.+",formula3[i]))
              test3.specs[[j]]= c(test3.specs[[min.BIC3]],0,names2specs(coef3[i], values=values))
              test3.names[[j]] = paste(c(newcoef,coef3[i]), collapse=" ")
              BIC3[j]=dcat.fit(test3.model[[j]], test3.specs[[j]])$norm.bic
            }
            if (min.BIC3 == which.min(BIC3)){break}
            best3.model = test3.model[[which.min(BIC3)]]
            best3.specs = test3.specs[[which.min(BIC3)]]
            best3.names = test3.names[[which.min(BIC3)]]
            newcoef = unlist(strsplit(best3.names, split=" "))
            coef3 = coef3[-which(coef3 %in% newcoef[nchar(newcoef)==3])]
            formula3 = sapply(coef3, insert_colon)
          }
        }
      }
      
      if(length(coef3)>=(length(possible_coef3)-3)){
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
            BIC4 = numeric(length(coef4)+1)
            BIC4[1]=best3.FIT$norm.bic
            min.BIC4 = 1
            # while the preferred one is a new one, should stop if it is the last one created
            for (i in 1:length(coef4)){
              test4.model[[i+1]]= update.formula(test4.model[[1]],paste("~.+",formula4[i]))
              test4.specs[[i+1]] = model2specs(test4.model[[i+1]], values = values)
              test4.names[[i+1]] = model2names(test4.model[[i+1]], values = values)
              BIC4[i+1]=dcat.fit(test4.model[[i+1]], test4.specs[[i+1]])$norm.bic
            }
            if(which.min(BIC4)==1){
              best4.model = test4.model[[1]]
              best4.specs = test4.specs[[1]]
              best4.names = test4.names[[1]]
            }else{
              min.BIC4 = which.min(BIC4)
              best4.model = test4.model[[min.BIC4]]
              best4.specs = test4.specs[[min.BIC4]]
              best4.names = test4.names[[min.BIC4]]
              newcoef = unlist(strsplit(test4.names[[min.BIC4]], split=" "))
              coef4 = coef4[-which(coef4 %in% newcoef[nchar(newcoef)==4])]
              formula4 = sapply(coef4, insert_star)
              while(length(coef4)>0){
                j=length(test4.model)
                min.BIC4 = which.min(BIC4)
                for (i in 1:length(coef4)){
                  j=j+1
                  test4.model[[j]] = update.formula(test4.model[[min.BIC4]],paste("~.+",formula4[i]))
                  test4.specs[[j]] = model2specs(test4.model[[j]], values = values)
                  test4.names[[j]] = model2names(test4.model[[j]], values = values)
                  BIC4[j]=dcat.fit(test4.model[[j]], test4.specs[[j]])$norm.bic
                }
                if (min.BIC4 == which.min(BIC4)){break}
                best4.model = test4.model[[which.min(BIC4)]]
                best4.specs = test4.specs[[which.min(BIC4)]]
                best4.names = test4.names[[which.min(BIC4)]]
                newcoef = unlist(strsplit(best4.names, split=" "))
                coef4 = coef4[-which(coef4 %in% newcoef[nchar(newcoef)==4])]
                formula4 = sapply(coef4, insert_star)
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
  
  dark.figure = est.pop.size-obs.pop.size 
  results[mc,2] = null_obs
  results[mc,3] = est.pop.size
  results[mc,4] = best.names
  
}

rownames(results) = c("n000", "est.pop.size")
results = as.data.frame(t(results))

colMeans(results[,c("n000", "est.pop.size")]) # 19,238.836
colSds(  as.matrix(results[,c("n000", "est.pop.size")]))/(MC^0.5) # 10.657742

# clearly biased 

table(as.vector(results[,4])) # length of 100 so never twice the same
length(which(as.vector(results[,4])=="G W Z B F S ZS WF GF ZF FS ZFS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS GF WF ZF FS ZFS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS ZF GF WF FS ZFS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS ZF WF GF FS ZFS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS WF ZF GF FS ZFS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS GF ZF WF FS ZFS"))

length(which(as.vector(results[,4])=="G W Z B F S ZS WF GF ZF FS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS GF WF ZF FS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS ZF GF WF FS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS ZF WF GF FS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS WF ZF GF FS"))
length(which(as.vector(results[,4])=="G W Z B F S ZS GF ZF WF FS"))
mean(as.numeric(results[,3]))
sd(as.numeric(results[,3]))/(MC^0.5)

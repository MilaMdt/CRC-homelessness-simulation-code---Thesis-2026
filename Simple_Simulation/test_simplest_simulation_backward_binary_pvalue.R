# Note : Here we do not use ECM because the start.specs have more than 400 values 
# and the R session aborts when it ecm.cat is run with it

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
           b2 = 4,
           b3 = 5,
           b4 = 6,
           b5 = 7,
           f = 8,
           s = 9)

a = 0
g = 1
w = 2
z = 3
b2 = 4
b3 = 5
b4 = 6
b5 = 7
f = 8
s = 9

## EM algorithm

dcat.fit <- function(MODEL, MARGINS){ 
  #dcat.em         <- em.cat (dcat.list, start = struc.zero, showits=FALSE) 
  #dcat.ecm        <- ecm.cat(dcat.list, start = struc.zero, margins = MARGINS, showits=FALSE, eps = 1e-7) 
  # it seems like the start argument could be taken out and the results would be the same
  # this might be because dcat.list is based on dfile and not tb, and the 222 combinations are taken out
  
  
  #tb.ecm          <- tb # grouped data
  #tb.ecm$Freq     <- as.numeric(dcat.ecm*dcat.list$n) # results with log-linear EM : new frequencies
  
  # this two lines are not in the current Rscript:
  #   tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==1,"Freq"] = tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==1,"Freq"] + 1/6
  #   tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==2,"Freq"] = tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==2,"Freq"] + 1/3  #new

  tb.fit.glm      <- glm(MODEL, family = poisson, subset = (struc.zero == 1), data = tb)  #original
  tb.summary      <- summary(tb.fit.glm) # summary of log-linear model based on EM log-linear imputation
  tb.model.matrix <- model.frame(MODEL, data = tb)
  tb.fitted   <- predict(tb.fit.glm, tb.model.matrix, type = "response") 
  # warning : In predict.lm(object, newdata, se.fit, scale = 1, type = if (type ==  :
  # prediction from rank-deficient fit; attr(*, "non-estim") has doubtful cases
  tb$est = tb.fitted
  #  tb.ecm[tb.ecm$struc.zero==1,"est"] = tb.ecm[tb.ecm$struc.zero==1,"Freq"] # not in current code
  norm.deviance   <- -2*logLik(tb.fit.glm)
  npar            <- tb.fit.glm$df.null- tb.fit.glm$df.residual + 1
  n               <- nrow(tb)
  norm.aic        <- as.numeric(norm.deviance + 2*npar)
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
    tb = tb)) # table with information on the model
} 

# we create a function that will be useful later
# it inserts ":" in between the letters of the coefs vector, to include it in the formula updating
insert_colon = function(x){
  var = regmatches(x, gregexpr("[A-Za-z][0-9]*", x))[[1]]
  paste0(var, collapse = ":")
}

allcomb = function(coef, n=3){
  pattern = "[A-Za-z][2-9]*"
  letters = unlist(regmatches(coef, gregexpr(pattern, coef)))
  combs = combn(letters, n, simplify=TRUE) # all combinations per 4-way term
  perms = unlist(apply(combs,2, function(x){
    lapply(permn(x),paste0, collapse = "") # all permutations from each combination
  }))
  return(perms)
}

# another function for the specs based on the names
names2specs = function(x, values){
  if(length(x)>1){
    low = str_to_lower(x) 
    vec = Reduce(function(x,y) c(x,"a",y), low)
    pattern = "[A-Za-z][2-9]*"
    letters = unlist(regmatches(vec, gregexpr(pattern, vec)))
    specs = values[letters]
    names(specs)=NULL
  }
  if(length(x)==1){
    low = str_to_lower(x)
    low = gsub(" ", "a", low)
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
smallest_p = function(pvalues){
  if(length(pvalues)>1){
    names(pvalues) = gsub("[^a-zA-Z][0-9]*","",names(pvalues))
    pv = numeric()
    for (i in 1:length(unique(names(pvalues)))){
      pv[i] = pvalues[which(names(pvalues)== unique(names(pvalues))[i])][which.min(pvalues[which(names(pvalues)== unique(names(pvalues))[i])])] 
      names(pv)[i]= insert_colon(unique(names(pvalues))[i])
    }
    return(pv)
  } else {
    pv = pvalues
    return(pv)
  }
}





MC = 625
popSize = 20000
Nsamples = 3

results = matrix(NA, MC, 5)
colnames(results)=c("itt","n000","est.pop.size", "model","chap.est.pop")
results[,"itt"] = 1:MC


get_prob= function(age, famous,sex, sample){
  base_prob = 0.2
  
  if(sex==2 & sample==3){
    base_prob = base_prob +0.15
  }
  
  base_prob = base_prob + 0.15*(famous==1)
  
  #if(age==2){base_prob = base_prob + 0.1}
  
  #if((age==1)|(age==3)){base_prob = base_prob - 0.1}
  
  pmax(pmin(base_prob, 1), 0) # keep base_prob between 0 and 1
}

for(mc in 1:MC){
  
  famous = sample(c(1,2), popSize, replace=TRUE, prob=c(0.4,0.6)) # 1 for famous
  sex    = sample(c(1,2), popSize, replace=TRUE, prob=c(0.5,0.5)) # 1 for woman, 2 for man
  age    = sample(c(1,2,3,4,5), popSize, replace=TRUE, prob=c(0.2,0.4,0.2,0.1,0.1))
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
  colnames(dfile) = c("G","W","Z","B","F","S")
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
  
  
  
  B2 = model.matrix(~B, tb)[,2]
  B3 = model.matrix(~B, tb)[,3]
  B4 = model.matrix(~B, tb)[,4]
  B5 = model.matrix(~B, tb)[,5]
  names(B2) = NULL
  names(B3) = NULL
  names(B4) = NULL
  names(B5) = NULL
  tb$B2 = as.factor(B2)
  tb$B3 = as.factor(B3)
  tb$B4 = as.factor(B4)
  tb$B5 = as.factor(B5)
  tb$B = NULL
  
  for (col in c('B2','B3','B4','B5','F','S')){
    tb[[col]] = as.numeric(tb[[col]])
  }
  
  for (col in c('B2','B3','B3','B4','S','F')){
    for (row in 1:nrow(tb)){
      if(tb[[col]][row] == 1) tb[[col]][row] = 0 else tb[[col]][row] = 1
    }
  }
  
  tb = tb[, c("G","W","Z","F","S","B2","B3","B4","B5","Freq","struc.zero")]
  
  # need to change the start model to the one with all terms up to 4-way
  # should we take out 0 cell counts?
  
  # main effects
  me = c("G","W","Z","B2","B3","B4","B5","F","S")
  # interactions in the start model
  twoe   = unlist(lapply(combn(me,2,simplify=FALSE),paste0, collapse = ""))
  threee = unlist(lapply(combn(me,3,simplify=FALSE),paste0, collapse = ""))
  foure  = unlist(lapply(combn(me,4,simplify=FALSE),paste0, collapse = ""))
  
  # we do not want GZW to be in the 3rd and 4th orders
  cov2_demon = grepl("(\\d.*\\d)", twoe)
  twoe = twoe[-which(cov2_demon)]
  threee = threee[-which(threee == "GWZ")]
  cov3_demon = grepl("(\\d.*\\d)", threee)
  threee = threee[-which(cov3_demon)]
  bigdemon = sapply(foure, function(x) all(c("G", "Z", "W") %in% strsplit(x, "")[[1]]))
  foure = foure[-which(bigdemon)]
  cov4_demon = grepl("(\\d.*\\d)", foure) # if more than one number we take it out, because two Bs cannot be estimated together
  foure = foure[-which(cov4_demon)]
  
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
  start.coef = unique(gsub(":","",pars.names)) 
  start.coef = gsub("1","",start.coef)
  
  
  test.model[[1]] = start.model
  test.specs[[1]] = names2specs(start.coef[-1], values=values) 
  test.names[[1]] = paste(start.coef[-1], collapse=" ")
  
  # now we want to focus on the 4-ways parameters
  coef4 = start.coef[nchar(gsub("[0-9]","",start.coef))==4]
  formula4 = sapply(coef4, insert_colon)
  if(length(coef4)==0){
    best4.model = test.model[[1]]
    best4.specs = test.specs[[1]]
    best4.names = test.names[[1]]
    new.coef = start.coef
  }
  if(length(coef4)!=0){
    test.pars4 = list()
    test.coef4 = list()
    i=1
    while (length(coef4)!=0){
      test.pars4[[i]] = round(dcat.fit(test.model[[i]], test.specs[[i]])$tb.summary$coefficients[,c(1,2,4)],
                              digits=2)
      # here one interaction can have two coefficients if one variable has more than 2 possible values
      # we want to take out an interaction only if none of the coefficients is significant
      coefs = unique(gsub("[:1]","",rownames(test.pars4[[i]])))
      test.coef4[[i]] = coefs[nchar(gsub("[0-9]","",coefs))==4]
      coef4 = test.coef4[[i]]
      #no4 = length(coefs)-length(test.coef4[[i]])
      pvalues4 = test.pars4[[i]][which(gsub("[:1]","",rownames(test.pars4[[i]])) %in% coef4),3]
      # small_pvalues4 = smallest_p(pvalues4) not needed anymore
      if(length(pvalues4)>1){
      out = which.max(pvalues4)
      if(pvalues4[out]>0.01){
        new.coef = coefs[-which(coefs %in% gsub("[:1]*","",names(out)))]
        # create new model
        test.model[[i+1]] = as.formula(paste("Freq~",paste(sapply(new.coef[-1], insert_colon), collapse="+")))
        test.specs[[i+1]] = names2specs(new.coef[-1], values=values)
        test.names[[i+1]] = paste(new.coef[-1], collapse=" ")
      } else{
        new.coef = coefs
        best4.model = test.model[[i]]
        best4.specs = test.specs[[i]]
        best4.names = test.names[[i]]
        break
      }
      coef4 = new.coef[nchar(gsub("[0-9]","",new.coef))==4]
      best4.model = test.model[[i+1]]
      best4.specs = test.specs[[i+1]]
      best4.names = test.names[[i+1]]
      i=i+1
      } else {
        if(pvalues4>0.01){
          new.coef = coefs[-which(coefs %in% coef4)] # you take out the only 4-way interaction
          # create new model
          test.model[[i+1]] = as.formula(paste("Freq~",paste(sapply(new.coef[-1], insert_colon), collapse="+")))
          test.specs[[i+1]] = names2specs(new.coef[-1], values=values)
          test.names[[i+1]] = paste(new.coef[-1], collapse=" ")
        } else{
          new.coef = coefs
          best4.model = test.model[[i]]
          best4.specs = test.specs[[i]]
          best4.names = test.names[[i]]
          break
        }
        coef4 = new.coef[nchar(gsub("[0-9]","",new.coef))==4]
        best4.model = test.model[[i+1]]
        best4.specs = test.specs[[i+1]]
        best4.names = test.names[[i+1]]
        break
      }
    }
    
    }
  
  
  # we want to keep all interactions that are included in a higher term kept in the model
  
  # what are the 4-way terms still in the model?
  if(length(coef4)==0) {
    coef3 = new.coef[nchar(gsub("[0-9]","",new.coef))==3]
    out_test = 0
  } # well then we're good
  
  # if not, we can make a list of interactions not to take out
  # if GZBL is in the model for example, we want all 3-way interactions from this term to stay in the model
  # such as GZL, GZB, but also GLZ and ZGL, because we do not know in what order they are in the model names

  
  if(length(coef4)>=1){
    exclude = unlist(lapply(coef4, allcomb))
    coef3 = new.coef[nchar(gsub("[0-9]","",new.coef))==3]
    out_test = which(coef3 %in% exclude)
    coef3 = coef3[-out_test]
  }
  
  # now we can apply the same algorithm to the 3-way interactions
  
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
    new.coef = best4.coef
  }
  if(length(coef3)!=0){
    test.pars3 = list()
    test.coef3 = list()
    i=1
    while (length(coef3)!=0){
      test.pars3[[i]] = round(dcat.fit(test3.model[[i]], test3.specs[[i]])$tb.summary$coefficients[,c(1,2,4)],
                              digits=2)
      coefs = unique(gsub("[:1]","",rownames(test.pars3[[i]])))
      test.coef3[[i]] = coefs[nchar(gsub("[0-9]","",coefs))==3]
      pvalues3 = test.pars3[[i]][which(gsub("[:1]","",rownames(test.pars3[[i]])) %in% coef3),3]
      # small_pvalues3 = smallest_p(pvalues3)
      if(length(pvalues3)>1){
      out = which.max(pvalues3)
      if(pvalues3[out]>0.01){
        new.coef = coefs[-which(coefs %in% gsub("[:1]","",names(out)))]
        # create new model
        test3.model[[i+1]] = as.formula(paste("Freq~",paste(sapply(new.coef[-1], insert_colon), collapse="+")))
        test3.specs[[i+1]] = names2specs(new.coef[-1], values=values)
        test3.names[[i+1]] = paste(new.coef[-1], collapse=" ")
      } else{
        new.coef = coefs
        best3.model = test3.model[[i]]
        best3.specs = test3.specs[[i]]
        best3.names = test3.names[[i]]
        break
      }
      coef3 = coef3[-which(coef3 %in% gsub("[:1]","",names(out)))]
      best3.model = test3.model[[i+1]]
      best3.specs = test3.specs[[i+1]]
      best3.names = test3.names[[i+1]]
      i=i+1
      }
      if(length(pvalues3)==1){
        if(pvalues3<=0.01){
          new.coef = coefs
          best3.model = test3.model[[i]]
          best3.specs = test3.specs[[i]]
          best3.names = test3.names[[i]]
          break
        }
        if(pvalues3>0.01){
          new.coef = coefs[-which(coefs %in% coef3)]
          # create new model
          test3.model[[i+1]] = as.formula(paste("Freq~",paste(sapply(new.coef[-1], insert_colon), collapse="+")))
          test3.specs[[i+1]] = names2specs(new.coef[-1], values=values)
          test3.names[[i+1]] = paste(new.coef[-1], collapse=" ")
          best3.model = test3.model[[i+1]]
          best3.specs = test3.specs[[i+1]]
          best3.names = test3.names[[i+1]]
          break
        }
      }
    }
    
  }
  coef3 = new.coef[nchar(gsub("[0-9]","",new.coef))==3]
  
  
  if(length(coef3)==0) {
    coef2 = new.coef[nchar(gsub("[0-9]","",new.coef))==2]
    out_test = 0
  } # keeping all 2-way interactions that were in the last version of newcoef
  # they should also be the 2-way interactions that were in start.coef and in best4.coef
  
  
  if(length(coef3)>=1){
    exclude = unlist(lapply(coef3, allcomb, n=2))
    coef2 = new.coef[nchar(gsub("[0-9]","",new.coef))==2]
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
  
  # after that we still have a problem : coef2 doesn't include all it should, which misleads pvalue2
  # why?
  
  # while the preferred one is a new one, should stop if it is the last one created
  if(length(coef2)==0){
    best2.model = test2.model[[1]]
    best2.specs = test2.specs[[1]]
    best2.names = test2.names[[1]]
    new.coef = best3.coef
  }
  if(length(coef2)!=0){
    test.pars2 = list()
    test.coef2 = list()
    i=1
    while (length(coef2)!=0){
      test.pars2[[i]] = round(dcat.fit(test2.model[[i]], test2.specs[[i]])$tb.summary$coefficients[,c(1,2,4)],
                              digits=2)
      coefs = unique(gsub("[:1]","",rownames(test.pars2[[i]])))
      test.coef2[[i]] = coefs[nchar(gsub("[0-9]","",coefs))==2]
      pvalues2 = test.pars2[[i]][which(gsub("[:1]","",rownames(test.pars2[[i]])) %in% coef2),3]
      # small_pvalues2 = smallest_p(pvalues2)
      if(length(pvalues2)>1){
      out = which.max(pvalues2)
      if(pvalues2[out]>0.01){
        new.coef = coefs[-which(coefs %in% gsub("[:1]","",names(out)))]
        # create new model
        test2.model[[i+1]] = as.formula(paste("Freq~",paste(sapply(new.coef[-1], insert_colon), collapse="+")))
        test2.specs[[i+1]] = names2specs(new.coef[-1], values=values)
        test2.names[[i+1]] = paste(new.coef[-1], collapse=" ")
      } else{
        new.coef = coefs
        best2.model = test2.model[[i]]
        best2.specs = test2.specs[[i]]
        best2.names = test2.names[[i]]
        break
      }
      coef2 = coef2[-which(coef2 %in% gsub("[:1]","",names(out)))]
      best2.model = test2.model[[i+1]]
      best2.specs = test2.specs[[i+1]]
      best2.names = test2.names[[i+1]]
      i=i+1
      }
      if(length(pvalues2)==1){
        if(pvalues2<=0.01){
          new.coef = coefs
          best2.model = test2.model[[i]]
          best2.specs = test2.specs[[i]]
          best2.names = test2.names[[i]]
          break
        }
        if(pvalues2>0.01){
          new.coef = coefs[-which(coefs %in% coef2)]
          # create new model
          test2.model[[i+1]] = as.formula(paste("Freq~",paste(sapply(new.coef[-1], insert_colon), collapse="+")))
          test2.specs[[i+1]] = names2specs(new.coef[-1], values=values)
          test2.names[[i+1]] = paste(new.coef[-1], collapse=" ")
          best2.model = test2.model[[i+1]]
          best2.specs = test2.specs[[i+1]]
          best2.names = test2.names[[i+1]]
          break
        }
      }
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
  colnames(fitted)	= c("G","W","Z","F","S","B2","B3","B4","B5","obs","est")
  
  est.pop.size = sum(fitted$est)
  obs.pop.size = sum(fitted$obs)
  
  
  formula = as.character(best2.model[3])
  formula_X = strsplit(formula, " ")[[1]]
  formula_X = formula_X[formula_X != "" & formula_X != "\n"]
  formula_X = formula_X[!grepl("S", formula_X)] # take out terms with S
  formula_X = formula_X[!grepl("F", formula_X)] # take out terms with F
  formula_X = formula_X[!grepl("B2", formula_X)] # take out terms with B2
  formula_X = formula_X[!grepl("B3", formula_X)] # take out terms with B3
  formula_X = formula_X[!grepl("B4", formula_X)] # take out terms with B4
  formula_X = formula_X[!grepl("B5", formula_X)] # take out terms with B5
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
  colnames(chap.fitted)	= c("G","W","Z","F","S","B2","B3","B4","B5","obs","est")
  
  est.chap.size = sum(chap.fitted$est)
  obs.chap.size = sum(chap.fitted$obs)
  
  dark.figure = est.pop.size-obs.pop.size 
  results[mc,2]=dark.figure
  results[mc,3]=est.pop.size
  results[mc,4]=best2.names
  results[mc,5]= est.chap.size
  
}

table(as.vector(results[,4]))
mean(as.numeric(results[,3]))
sd(as.numeric(results[,3]))/(MC^0.5)
mean(as.numeric(results[,5]))
sd(as.numeric(results[,5]))/(MC^0.5)

length(which(as.vector(results[,4])=="G W Z B2 B3 B4 B5 F S GF WF ZF ZS FS"))
length(which(as.vector(results[,4])=="G W Z B2 B3 B4 B5 F S GF WF ZF ZS FS ZFS"))

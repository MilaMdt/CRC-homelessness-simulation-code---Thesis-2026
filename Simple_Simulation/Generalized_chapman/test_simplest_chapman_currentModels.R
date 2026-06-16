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
           h = 5,
           l = 6,
           s = 7)

a=0
g=1
w=2
z=3
b=4
h=5
l=6
s=7



get_prob= function(age, famous,sex,town, sample){
  base_prob = 0.2
  
  if(sex==1 & sample==3){
    base_prob = base_prob + 0.15
  }
  
  base_prob = base_prob + 0.15*(famous==1)
  
  #if(age==2){base_prob = base_prob + 0.1}
  
  #if((age==1)|(age==3)){base_prob = base_prob - 0.1}
  
  #if(town == 2){base_prob = base_prob + 0.05}
  
  #if(town == 1){base_prob = base_prob - 0.05}
  
  pmax(pmin(base_prob, 1), 0) # keep base_prob between 0 and 1
}


## EM algorithm

dcat.fit <- function(MODEL, MARGINS){ 
  dcat.em         <- em.cat (dcat.list, start = struc.zero, showits=FALSE) 
  dcat.ecm        <- ecm.cat(dcat.list, start = struc.zero, margins = MARGINS, showits=FALSE, eps = 1e-7) 
  
  tb.ecm          <- tb # grouped data
  tb.ecm$Freq     <- as.numeric(dcat.ecm*dcat.list$n) # results with log-linear EM : new frequencies
  
  formula_X = strsplit(as.character(MODEL[3]), " ")[[1]]
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
  
  tb.fit.glm      <- glm(MODEL, family = poisson, subset = (struc.zero == 1), data = tb.ecm_Chapman)  #original
  tb.summary      <- summary(tb.fit.glm) # summary of log-linear model based on EM log-linear imputation
  tb.model.matrix <- model.frame(MODEL, data = tb.ecm)
  tb.ecm.fitted   <- predict(tb.fit.glm, tb.model.matrix, type = "response")
  tb.ecm$est = tb.ecm.fitted
  #  tb.ecm[tb.ecm$struc.zero==1,"est"] = tb.ecm[tb.ecm$struc.zero==1,"Freq"] # not in current code
  norm.deviance   <- -2*(logpost.cat(dcat.list,dcat.ecm)-logpost.cat(dcat.list,dcat.em))
  npar            <- tb.fit.glm$df.null- tb.fit.glm$df.residual + 1
  norm.aic        <- norm.deviance + 2*npar
  d.o.f           <- tb.fit.glm$df.residual
  tbf             <- tb
  tbf$Freq        <- tb.ecm.fitted 
  tb.missed       <- xtabs(Freq ~ struc.zero, data = tbf )[1]
  
  
  return(list(
    norm.deviance = norm.deviance,
    d.o.f         = d.o.f,
    norm.aic      = norm.aic,
    tb.missed     = tb.missed, 
    tb.summary    = tb.summary,
    tb.fit.glm    = tb.fit.glm,
    tb.ecm.fitted = tb.ecm.fitted,
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
colnames(results)=c("itt","n000","est.pop.size")
results[,"itt"] = 1:MC




for(mc in 1:MC){
  
  famous = sample(c(1,2), popSize, replace=TRUE, prob=c(0.4,0.6)) # 1 for famous
  sex    = sample(c(1,2), popSize, replace=TRUE, prob=c(0.5,0.5)) # 1 for woman, 2 for man
  age    = sample(c(1,2,3), popSize, replace=TRUE, prob=c(0.2,0.6,0.2))
  town = sample(c(1,2), popSize, replace=TRUE, prob=c(0.5,0.5)) 
  # created the population
  
  Samples = matrix(NA, popSize, Nsamples)
  inclusionProbs = matrix(NA, popSize, Nsamples)
  
  for (s in 1:Nsamples){
      inclusionProbs[,s] = mapply(get_prob,age=age, famous = famous, sex=sex, town = town, MoreArgs=list(sample=s))
      Samples[,s] = rbinom(popSize, 1, inclusionProbs[,s])
  }
  colnames(Samples) = c("G","W","Z")
  freq = 1
  dfile = cbind(Samples, age, famous, sex, town)
  colnames(dfile) = c("G","W","Z", "B","H","L","S")
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
  for (col in c('G','W','Z','B','H','L','S')){
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

  # basismodel 

test.model[[1]] <-  formula(Freq ~G*W*(S*B)+G*(Z*H)+G*L+W*(S*L*B + S*H*B)+Z*(S*L)+Z*(L*H)+B*(Z+L)+S*L*H)
test.specs[[1]] <- c(g,w,s,b,0,g,z,h,0,g,l,0,w,s,l,b,0,w,s,h,b,0,z,s,l,0,z,l,h,0,z,b,0,l,b,0,s,l,h)
test.names[[1]] <- c("GW*(S*B)+G(Z*H)+GL+W(S*L*B)+W(S*H*B)+Z(S*L)+Z(L*H)+ZB+LB+SLH")

# startmodel
test.model[[2]] <- formula(Freq ~ G*W*L*B+G*W*H+G*Z*L*B+G*Z*H+S*L*H+S*L*B+S*H*B)
test.specs[[2]] <- c(g,w,l,b,0,g,w,h,0,g,z,l,b,0,g,z,h,0,s,l,h,0,s,l,b,0,s,h,b)
test.names[[2]] <- c("GWLB GWH GZLB GZH SLH SLB SHB")

test.model[[3]] <- update.formula(test.model[[2]],~.-G:W:L:B)
test.specs[[3]] <- c(g,w,l,0,g,w,b,0,g,l,b,0,g,w,h,0,g,z,l,b,0,g,z,h,0,w,l,b,0,s,l,h,0,s,l,b,0,s,h,b)
test.names[[3]] <- c("GWL GWB GLB GWH GZLB GZH WLB SLH SLB SHB")

test.model[[4]] <- update.formula(test.model[[3]],~.-G:Z:L:B)
test.specs[[4]] <- c(g,w,l,0,g,w,b,0,g,l,b,0,g,w,h,0,g,z,l,0,g,z,b,0,g,l,b,0,g,z,h,0,w,l,b,0,z,l,b,0,s,l,h,0,s,l,b,0,s,h,b)
test.names[[4]] <- c("GWL GWB GLB GWH GZL GZB GLB GZH WLB ZLB SLH SLB SHB")

test.model[[5]] <- update.formula(test.model[[4]],~.-G:W:L)
test.specs[[5]] <- c(g,w,b,0,g,l,b,0,g,w,h,0,g,z,l,0,g,z,b,0,g,l,b,0,g,z,h,0,w,l,b,0,z,l,b,0,s,l,h,0,s,l,b,0,s,h,b)
test.names[[5]] <- c("GWB GLB GWH GZL GZB GLB GZH WLB ZLB SLH SLB SHB")

test.model[[6]] <- update.formula(test.model[[5]],~.-G:L:B)
test.specs[[6]] <- c(g,w,b,0,g,w,h,0,g,z,l,0,g,z,b,0,g,l,b,0,g,z,h,0,w,l,b,0,z,l,b,0,s,l,h,0,s,l,b,0,s,h,b)
test.names[[6]] <- c("GWB GWH GZL GZB GLB GZH WLB ZLB SLH SLB SHB")

test.model[[7]] <- update.formula(test.model[[6]],~.-G:W:H)
test.specs[[7]] <- c(g,w,b,0,g,z,l,0,g,z,b,0,g,l,b,0,g,z,h,0,w,h,0,w,l,b,0,z,l,b,0,s,l,h,0,s,l,b,0,s,h,b)
test.names[[7]] <- c("GWB GZL GZB GLB GZH WH WLB ZLB SLH SLB SHB")

test.model[[8]] <- update.formula(test.model[[7]],~.-G:B:Z)
test.specs[[8]] <- c(g,w,b,0,g,z,l,0,g,l,b,0,g,z,h,0,w,h,0,w,l,b,0,z,l,b,0,s,l,h,0,s,l,b,0,s,h,b)
test.names[[8]] <- c("GWB GZL GLB GZH WH WLB ZLB SLH SLB SHB")

test.model[[9]] <- update.formula(test.model[[8]],~.-G:L:Z)
test.specs[[9]] <- c(g,w,b,0,g,l,b,0,g,z,h,0,w,h,0,w,l,b,0,z,l,b,0,s,l,h,0,s,l,b,0,s,h,b)
test.names[[9]] <- c("GWB GLB GZH WH WLB ZLB SLH SLB SHB")

test.model[[10]] <- update.formula(test.model[[9]],~.-Z:B:L)
test.specs[[10]] <- c(g,w,b,0,g,l,b,0,g,z,h,0,w,h,0,w,l,b,0,z,l,0,z,b,0,s,l,h,0,s,l,b,0,s,h,b)
test.names[[10]] <- c("GWB GLB GZH WH WLB ZL ZB SLH SLB SHB")

test.model[[11]] <- update.formula(test.model[[10]],~.-S:B:H)
test.specs[[11]] <- c(g,w,b,0,g,l,b,0,g,z,h,0,w,h,0,w,l,b,0,z,l,0,z,b,0,s,l,h,0,s,l,b,0,h,b)
test.names[[11]] <- c("GWB GLB GZH WH WLB ZL ZB SLH SLB HB")

test.model[[12]] <- update.formula(test.model[[11]],~.-G:Z:H-G:H)
test.specs[[12]] <- c(g,w,b,0,g,l,b,0,g,z,0,w,h,0,w,l,b,0,z,l,0,z,h,0,z,b,0,s,l,h,0,s,l,b,0,h,b)
test.names[[12]] <- c("GWB GLB GZ WH WLB ZL ZH ZB SLH SLB HB")

test.model[[13]] <- update.formula(test.model[[12]],~.+G*W*S)
test.specs[[13]] <- c(g,w,b,0,g,l,b,0,g,w,s,0,g,z,0,w,h,0,w,l,b,0,z,l,0,z,h,0,z,b,0,s,l,h,0,s,l,b,0,h,b)
test.names[[13]] <- c("GWB GLB GWS GZ WH WLB ZL ZH ZB SLH SLB HB")

test.model[[14]] <- update.formula(test.model[[13]],~.+G*W*S)
test.specs[[14]] <- c(g,w,b,0,g,w,s,0,g,l,b,0,g,z,0,w,h,0,w,l,b,0,z,l,0,z,h,0,z,b,0,s,l,h,0,s,l,b,0,h,b)
test.names[[14]] <- c("GWB GWS GLB GZ WH WLB ZL ZH ZB SLH SLB HB")

test.model[[15]] <- update.formula(test.model[[14]],~.+G*W*S*B)
test.specs[[15]] <- c(g,w,s,b,0,g,l,b,0,g,z,0,w,h,0,w,l,b,0,z,l,0,z,h,0,z,b,0,s,l,h,0,s,l,b,0,h,b)
test.names[[15]] <- c("GWSB GLB GZ WH WLB ZL ZH ZB SLH SLB HB")

test.model[[16]] <- update.formula(test.model[[15]],~.+G*B*H)
test.specs[[16]] <- c(g,w,s,b,0,g,z,0,g,l,b,0,g,h,b,0,w,h,0,w,l,b,0,z,l,0,z,h,0,z,b,0,s,l,h,0,s,l,b,0,h,b)
test.names[[16]] <- c("GWSB GZ GLB GHB WH WLB ZL ZH ZB SLH SLB HB")

test.model[[17]] <- update.formula(test.model[[16]],~.+Z*S*L)
test.specs[[17]] <- c(g,w,s,b,0,g,z,0,g,l,b,0,g,h,b,0,w,h,0,w,l,b,0,z,s,l,0,z,h,0,z,b,0,s,l,h,0,s,l,b,0,h,b)
test.names[[17]] <- c("GWSB GZ GLB GHB WH WLB ZSL ZH ZB SLH SLB HB")

test.model[[18]] <- update.formula(test.model[[17]],~.+Z*L*H)
test.specs[[18]] <- c(g,w,s,b,0,g,z,0,g,l,b,0,g,h,b,0,w,h,0,w,l,b,0,z,s,l,0,z,l,h,0,z,b,0,s,l,h,0,s,l,b,0,h,b)
test.names[[18]] <- c("GWSB GZ GLB GHB WH WLB ZSL ZLH ZB SLH SLB HB")

test.model[[19]] <- update.formula(test.model[[18]],~.+Z*H*B)
test.specs[[19]] <- c(g,w,s,b,0,g,z,0,g,l,b,0,g,h,b,0,w,h,0,w,l,b,0,z,s,l,0,z,l,h,0,z,h,b,0,s,l,h,0,s,l,b)
test.names[[19]] <- c("GWSB GZ GLB GHB WH WLB ZSL ZLH ZHB SLH SLB")

test.model[[20]] <- update.formula(test.model[[19]],~.+G*W*H*B) 
test.specs[[20]] <- c(g,w,s,b,0,g,w,h,b,0,g,z,0,g,l,b,0,w,l,b,0,z,s,l,0,z,l,h,0,z,h,b,0,s,l,h,0,s,l,b)
test.names[[20]] <- c("GWSB GWHB GZ GLB WLB ZSL ZLH ZHB SLH SLB")

test.model[[21]] <- update.formula(test.model[[20]],~.-G:W:H:B) 
test.specs[[21]] <- c(g,w,s,b,0,g,w,h,0,g,w,s,0,g,z,0,g,l,b,0,g,h,b,0,w,l,b,0,w,h,b,0,z,s,l,0,z,l,h,0,z,h,b,0,s,l,h,0,s,l,b)
test.names[[21]] <- c("GWSB GWH GWS GZ GLB GHB WLB WHB ZSL ZLH ZHB SLH SLB")

test.model[[22]] <- update.formula(test.model[[21]],~.-B:H:Z) 
test.specs[[22]] <- c(g,w,s,b,0,g,w,h,0,g,w,s,0,g,z,0,g,l,b,0,g,h,b,0,w,l,b,0,w,h,b,0,z,s,l,0,z,l,h,0,z,h,0,z,b,0,s,l,h,0,s,l,b)
test.names[[22]] <- c("GWSB GWH GWS GZ GLB GHB WLB WHB ZSL ZLH ZH ZB SLH SLB")

test.model[[23]] <- update.formula(test.model[[22]],~.-G:B:H) 
test.specs[[23]] <- c(g,w,s,b,0,g,w,h,0,g,w,s,0,g,z,0,g,l,b,0,w,l,b,0,w,h,b,0,z,s,l,0,z,l,h,0,z,h,0,z,b,0,s,l,h,0,s,l,b)
test.names[[23]] <- c("GWSB GWH GWS GZ GLB WLB WHB ZSL ZLH ZH ZB SLH SLB")

test.model[[24]] <- update.formula(test.model[[23]],~.-G:W:H) 
test.specs[[24]] <- c(g,w,s,b,0,g,l,b,0,g,h,0,g,z,0,w,l,b,0,w,h,b,0,z,s,l,0,z,l,h,0,z,b,0,s,l,h,0,s,l,b)
test.names[[24]] <- c("GWSB GLB GH GZ WLB WHB ZSL ZLH ZB SLH SLB")

test.model[[25]] <- update.formula(test.model[[24]],~.+W*L*B*S) 
test.specs[[25]] <- c(g,w,s,b,0,g,l,b,0,g,h,0,g,z,0,w,s,l,b,0,w,h,b,0,z,s,l,0,z,l,h,0,z,b,0,s,l,h,0,s,l,b)
test.names[[25]] <- c("GWSB GLB GH GZ WSLB WHB ZSL ZLH ZB SLH SLB")

test.model[[26]] <- update.formula(test.model[[25]],~.+W*H*B*S) 
test.specs[[26]] <- c(g,w,s,b,0,g,l,b,0,g,h,0,g,z,0,w,s,l,b,0,w,h,b,s,0,z,s,l,0,z,l,h,0,z,b,0,s,l,h,0,s,l,b)
test.names[[26]] <- c("GWSB GLB GH GZ WSLB WHBS ZSL ZLH ZB SLH SLB")

test.model[[27]] <- update.formula(test.model[[26]],~.-B:S) 
test.specs[[27]] <- c(g,w,s,b,0,g,l,0,g,h,0,g,z,0,w,s,l,b,0,w,h,b,s,0,z,s,l,0,z,l,h,0,z,b,0,s,l,h)
test.names[[27]] <- c("GWSB GL GH GZ WSLB WHBS ZSL ZLH ZB SLH")

test.model[[28]] <- update.formula(test.model[[27]],~.+L:B) 
test.specs[[28]] <- c(g,w,s,b,0,g,l,0,g,h,0,g,z,0,w,s,l,b,0,w,h,b,s,0,z,s,l,0,z,l,h,0,z,b,0,l,b,0,s,l,h)
test.names[[28]] <- c("GWSB GL GH GZ WSLB WHBS ZSL ZLH ZB LB SLH")

test.model[[29]] <- update.formula(test.model[[28]],~.-G:H:Z) 
test.specs[[29]] <- c(g,w,s,b,0,g,l,0,w,s,l,b,0,w,h,b,s,0,z,s,l,0,z,l,h,0,z,b,0,l,b,0,s,l,h)
test.names[[29]] <- c("GWSB GL WSLB WHBS ZSL ZLH ZB LB SLH")

test.model[[30]] <- update.formula(test.model[[29]],~.+G*Z*H) 
test.specs[[30]] <- c(g,w,s,b,0,g,z,h,0,g,l,0,w,s,l,b,0,w,s,h,b,0,z,s,l,0,z,l,h,0,z,b,0,l,b,0,s,l,h)
test.names[[30]] <- c("GWSB GZH GL WSLB WSHB ZSL ZLH ZB LB SLH")

n.test <- length(test.model)

test.norm          <-matrix(nrow=n.test,ncol=5)
colnames(test.norm)<-c("nr","DEV","df","AIC","nHat0")
rownames(test.norm)<-test.names

test.pars <-list(n.test)
test.freqs<-list(n.test)

for (j in 1:n.test){
  FIT <- dcat.fit(test.model[[j]], test.specs[[j]])
  test.norm[j,]<-cbind(j,round(FIT$norm.deviance),FIT$d.o.f,round(FIT$norm.aic,3),  
                       round(as.vector(FIT$tb.missed),1))
  
  test.pars[[j]] <-  round(FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)
  test.freqs[[j]] <- cbind(tb[,1:length(tb)-1],round(FIT$tb.ecm.fitted,2))
}
test.norm

# Door -1 hieronder te kiezen worden altijd alle parameters afgedrukt, veranderd aan oorspronkelijke code.

#print(subset(test.pars[[n.test]],test.pars[[n.test]][,3]>-1))

# Deel 3, soms werkt bootstrap niet, een ander model kiezen hieronder kan helpen.

#code aangepast n.a.v. onderzoek Daan Zult.

#best.model.origineel.list[[y]] <- 18
  best.model.origineel <- test.norm[test.norm[, "AIC"] == min(test.norm[, "AIC"]), "nr"]

  FIT <- dcat.fit(test.model[[best.model.origineel]], test.specs[[best.model.origineel]])
  coefs.origineel = FIT$tb.summary$coefficients
  tb.ecm.origineel = FIT$tb.ecm

  fitted <- test.freqs[best.model.origineel][[1]]
  colnames(fitted)	<- c("G","W","Z","B","H","L","S","obs","est")

  est.pop.size = sum(fitted$est)
  obs.pop.size = sum(fitted$obs)

  dark.figure = est.pop.size-obs.pop.size 
  results[mc,2]=dark.figure
  results[mc,3]=est.pop.size
}

colMeans(results[,c("n00", "N_est_standard_cov")]) # n000 6928.61, est.pop 18886.58
colSds(  results[,c("n00", "N_est_standard_cov")])/(MC^0.5) # n000 13.04 est.pop 13.7936
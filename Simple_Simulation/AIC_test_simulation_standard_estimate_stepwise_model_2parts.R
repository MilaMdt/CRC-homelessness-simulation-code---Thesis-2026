library(foreign)
library(cat)
library(haven)
# install.packages(combinat)
library(combinat)
library(stringr)

### betekenis van de kolommen:
# G = GBA (mensen in laagdrempelig opvangadres) 1=ja, 2=nee
# W = WBB (daklozen met uitkering)
# Z = Zwervend (bestand van de reclassering)
# S = Sekse 1=man, 2=vrouw
# L = Leeftijd  1=18-29, 2= 30-49, 3= 50-64
# H = Herkomst 3=niet-westers allochtoon, 1= autochtoon, 2= westers-allochton
# B = Dummy voor G4 1=ja, 2=nee

var.names <- c("G","W","Z","S","L","H","B")
g <- 1
w <- 2
z <- 3
l <- 4
h <- 5
s <- 6
b <- 7

values = c(a = 0,
           g = 1,
           w = 2,
           z = 3,
           l = 4,
           h = 5,
           s = 6,
           b = 7)

allcomb = function(coefs, n=3){
  letters = unlist(strsplit(coefs, split=""))
  combs = combn(letters, n, simplify=TRUE) # all combinations per 4-way term
  perms = unlist(apply(combs,2, function(x){
    lapply(permn(x),paste0, collapse = "") # all permutations from each combination
  }))
  return(perms)
}

# we create a function that will be useful later
# it inserts ":" in between the letters of the coefs vector, to include it in the formula updating
insert_colon = function(x){ paste(strsplit(x, "")[[1]], collapse=":") }

# another function for the specs based on the names
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
  names = gsub("\\+", " ", sum2)
  names
}

MC = 3
popSize = 20000
Nsamples = 3

results = matrix(NA, MC, 3)
colnames(results)=c("itt","n000","N_est_standard_cov")
results[,"itt"] = 1:MC


get_prob= function(age,famous,sex,town, sample){
  base_prob = 0.2
  
  if(sex==1 & sample==3){
    base_prob = base_prob - 0.07
  }
  
  if(famous == 1){
    base_prob = base_prob + 0.15
  }
  
  pmax(pmin(base_prob, 1), 0) # keep base_prob between 0 and 1
}

## EM algorithm

dcat.fit <- function(MODEL, MARGINS){ 
  dcat.em         <- em.cat (dcat.list, start = struc.zero, showits=FALSE) 
  dcat.ecm        <- ecm.cat(dcat.list, start = struc.zero, margins = MARGINS, showits=FALSE, eps = 1e-7) 
  
  tb.ecm          <- tb # grouped data
  tb.ecm$Freq     <- as.numeric(dcat.ecm*dcat.list$n) # results with log-linear EM : new frequencies
  
  # this two lines are not in the current Rscript:
  #   tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==1,"Freq"] = tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==1,"Freq"] + 1/6
  #   tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==2,"Freq"] = tb.ecm[rowSums(tb.ecm[,c("G","W","Z")]==1)==2,"Freq"] + 1/3  #new
  
  tb.fit.glm      <- glm(MODEL, family = poisson, subset = (struc.zero == 1), data = tb.ecm)  #original
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

for(mc in 1:MC){
  
  famous = sample(c(1,2), popSize, replace=TRUE, prob=c(0.5,0.5)) # 1 for famous
  sex    = sample(c(1,2), popSize, replace=TRUE, prob=c(0.4,0.6)) # 1 for woman, 2 for man
  age    = sample(c(1,2,3), popSize, replace=TRUE, prob=c(0.2,0.6,0.2))
  town = sample(c(1,2), popSize, replace=TRUE, prob=c(0.5,0.5)) 
  # created the population
  
  Samples = matrix(NA, popSize, Nsamples)
  inclusionProbs = matrix(NA, popSize, Nsamples)
  
  for (s in 1:Nsamples){
    inclusionProbs[,s] = mapply(get_prob,age=age, famous = famous,sex= sex,town=town, MoreArgs=list(sample=s))
    Samples[,s] = rbinom(popSize, 1, inclusionProbs[,s])
  }
  colnames(Samples) = c("G","W","Z")
  freq = 1
  dfile = cbind(Samples, age, famous, sex, town)
  colnames(dfile) = c("G","W","Z", "L","H","S","B")
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
  for (col in c('G','W','Z','L','H','S','B')){
    filtered_dfile[[col]] = as.numeric(filtered_dfile[[col]])
  }
  dmat.list    	 <- as.matrix(filtered_dfile)
  
  null_obs = sum(tb[tb[,1] == 2 & tb[,2] == 2 & tb[,3] == 2,]$Freq)
  tb[tb[,1] == 2 & tb[,2] == 2 & tb[,3] == 2,]$Freq = 0
  
  dcat.list	        <- prelim.cat(dmat.list) 
  struc.zero	    	<- array(struc.zero, dim = dcat.list$d)/sum(struc.zero) 
  
final.models <- data.frame(nr=integer(),
                           DEV=numeric(),
                           df=numeric(),
                           AIC=numeric(),
                           nHat0=numeric(),
                           model=character(),
                           specs=character())
model.norm          	<- data.frame(nr=integer(),
                                   DEV=numeric(),
                                   df=numeric(),
                                   AIC=numeric(),
                                   nHat0=numeric(),
                                   model=character(),
                                   specs=character())
colnames(model.norm)	<- c("nr","DEV","df","AIC","nHat0","model","specs")



## Models


GWZ.model<-list()
GWZ.specs<-list()
GWZ.names<-list()

# the models are based on the grouped data, not on the original data

# only sources
GWZ.model[[1]] <- Freq~G+W+Z
GWZ.specs[[1]] <- c(g,0,w,0,z)
GWZ.names[[1]] <- c("(G+W+Z)")

GWZ.model[[2]] <- Freq~G+W*Z
GWZ.specs[[2]] <- c(g,0,w,z)
GWZ.names[[2]] <- c("(G+WZ)")

GWZ.model[[3]] <- Freq~W+G*Z
GWZ.specs[[3]] <- c(w,0,g,z)
GWZ.names[[3]] <- c("(GZ+W)")

GWZ.model[[4]] <- Freq~Z+G*W
GWZ.specs[[4]] <- c(z,0,g,w)
GWZ.names[[4]] <- c("(GW+Z)")

GWZ.model[[5]] <- Freq~G*W+Z*W
GWZ.specs[[5]] <- c(g,w,0,w,z)
GWZ.names[[5]] <- c("(GW+WZ)")

GWZ.model[[6]] <- Freq~G*Z+W*Z
GWZ.specs[[6]] <- c(g,z,0,w,z)
GWZ.names[[6]] <- c("(GZ+WZ)")

# why is there not one (GW+GZ)? Funny enough, this is what is chosen in the test.models!!
# coincidence? I don't think so

GWZ.model[[7]] <- Freq~G*Z+W*Z+G*W
GWZ.specs[[7]] <- c(g,z,0,w,z,0,g,w)
GWZ.names[[7]] <- c("(GZ+WZ+GW)")

GWZ.model[[8]] <- Freq~G*Z+G*W
GWZ.specs[[8]] <- c(g,z,0,g,w)
GWZ.names[[8]] <- c("(GZ+GW)")

GWZ.model[[9]] <- Freq~G*Z+W*Z+G*W
GWZ.specs[[9]] <- c(g,z,0,w,z,0,g,w)
GWZ.names[[9]] <- c("(GZ+WZ+GW)")

GWZ.model[[10]] <- Freq~G*W*Z
GWZ.specs[[10]] <- c(g,w,z)
GWZ.names[[10]] <- c("(GWZ)")



n.GWZ.models <- length(GWZ.model)

for(zz in 1:n.GWZ.models){
  
  m0        <- GWZ.names[[zz]] # name of the model
  
  GWZ       <- GWZ.model[[zz]]
  gwz       <- GWZ.specs[[zz]]
  
  
  cov.model <- list()
  cov.specs <- list()
  cov.names <- list()
  
  cov.model[[1]] <- update.formula(GWZ, ~.) # from model in GWZ to all variables (G+W+Z+ all interactions already included in GWZ)
  cov.specs[[1]] <- gwz
  cov.names[[1]] <- c("")
  
  # creating models with covariates
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s)
  cov.names[[length(cov.names)+1]] <- c("S")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+L)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,l)
  cov.names[[length(cov.names)+1]] <- c("L")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,b)
  cov.names[[length(cov.names)+1]] <- c("B")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,h)
  cov.names[[length(cov.names)+1]] <- c("H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S+L)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,0,l)
  cov.names[[length(cov.names)+1]] <- c("S+L")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S+B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,0,b)
  cov.names[[length(cov.names)+1]] <- c("S+B")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,0,h)
  cov.names[[length(cov.names)+1]] <- c("S+H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+L+B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,l,0,b)
  cov.names[[length(cov.names)+1]] <- c("L+B")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+L+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,l,0,h)
  cov.names[[length(cov.names)+1]] <- c("L+H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+B+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,b,0,h)
  cov.names[[length(cov.names)+1]] <- c("B+H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S+L+B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,0,l,0,b)
  cov.names[[length(cov.names)+1]] <- c("S+L+B")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S+L+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,0,l,0,h)
  cov.names[[length(cov.names)+1]] <- c("S+L+H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S+B+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,0,b,0,h)
  cov.names[[length(cov.names)+1]] <- c("S+B+H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+L+B+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,l,0,b,0,h)
  cov.names[[length(cov.names)+1]] <- c("L+B+H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S+L+B+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,0,l,0,b,0,h)
  cov.names[[length(cov.names)+1]] <- c("S+L+B+H")
  
  # now also introducing interactions of covariates
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l)
  cov.names[[length(cov.names)+1]] <- c("SL")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L+B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,0,b)
  cov.names[[length(cov.names)+1]] <- c("SL+B")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,0,h)
  cov.names[[length(cov.names)+1]] <- c("SL+H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,b)
  cov.names[[length(cov.names)+1]] <- c("SB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*B+L)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,b,0,l)
  cov.names[[length(cov.names)+1]] <- c("SB+L")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*B+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,b,0,h)
  cov.names[[length(cov.names)+1]] <- c("SB+H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*H+L)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,h,0,l)
  cov.names[[length(cov.names)+1]] <- c("SH+L")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*H+B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,h,0,b)
  cov.names[[length(cov.names)+1]] <- c("SH+B")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+L*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,l,b)
  cov.names[[length(cov.names)+1]] <- c("LB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+L*B+S)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,l,b,0,s)
  cov.names[[length(cov.names)+1]] <- c("LB+S")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+L*B+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,l,b,0,h)
  cov.names[[length(cov.names)+1]] <- c("LB+H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L+S*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,0,s,b)
  cov.names[[length(cov.names)+1]] <- c("SL+SB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L+S*H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,0,s,h)
  cov.names[[length(cov.names)+1]] <- c("SL+SH")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L+L*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,0,l,b)
  cov.names[[length(cov.names)+1]] <- c("SL+LB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L+L*H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,0,l,h)
  cov.names[[length(cov.names)+1]] <- c("SL+LH")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L+B*H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,0,b,h)
  cov.names[[length(cov.names)+1]] <- c("SL+BH")
  
  # now with 3 terms interactions of covariates
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h)
  cov.names[[length(cov.names)+1]] <- c("SLH")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,0,b)
  cov.names[[length(cov.names)+1]] <- c("SLH+B")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+B*S)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h,0,b,s)
  cov.names[[length(cov.names)+1]] <- c("SHL+BS")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+B*L)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h,0,b,l)
  cov.names[[length(cov.names)+1]] <- c("SHL+BL")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+S*H*B+B*L)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h,0,s,h,b,0,b,l)
  cov.names[[length(cov.names)+1]] <- c("SHL+SHB+BL")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+B*H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h,0,b,h)
  cov.names[[length(cov.names)+1]] <- c("SLH+BH")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,b)
  cov.names[[length(cov.names)+1]] <- c("SLB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*B+H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,b,0,h)
  cov.names[[length(cov.names)+1]] <- c("SLB+H")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*B+H*S)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,b,0,h,s)
  cov.names[[length(cov.names)+1]] <- c("SLB+HS")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*B+H*L)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,b,0,h,l)
  cov.names[[length(cov.names)+1]] <- c("SLB+HL")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*B+H*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,b,0,h,b)
  cov.names[[length(cov.names)+1]] <- c("SLB+HB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*H*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,h,b)
  cov.names[[length(cov.names)+1]] <- c("SHB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*H*B+L)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,h,b,0,l)
  cov.names[[length(cov.names)+1]] <- c("SHB+L")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*H*B+L*S)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,h,b,0,l,s)
  cov.names[[length(cov.names)+1]] <- c("SHB+LS")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*H*B+L*H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,h,b,0,l,h)
  cov.names[[length(cov.names)+1]] <- c("SHB+LH")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*H*B+L*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,h,b,0,l,b)
  cov.names[[length(cov.names)+1]] <- c("SHB+LB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+S*L*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h,0,s,l,b)
  cov.names[[length(cov.names)+1]] <- c("SLH+SLB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+S*H*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h,0,s,h,b)
  cov.names[[length(cov.names)+1]] <- c("SLH+SHB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*H*B+L*H*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,h,b,0,l,h,b)
  cov.names[[length(cov.names)+1]] <- c("SHB+LHB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+S*L*B+S*H*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h,0,s,l,b,0,s,h,b)
  cov.names[[length(cov.names)+1]] <- c("SLH+SLB+SHB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+S*L*B+L*H*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h,0,s,l,b,0,l,h,b)
  cov.names[[length(cov.names)+1]] <- c("SLH+SLB+LHB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+S*H*B+L*H*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h,0,s,h,b,0,l,h,b)
  cov.names[[length(cov.names)+1]] <- c("SLH+SHB+LHB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*B+S*H*B+L*H*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,b,0,s,h,b,0,l,h,b)
  cov.names[[length(cov.names)+1]] <- c("SLB+SHB+LHB")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*B+S*H*B+S*H*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,b,0,s,h,b,0,s,h,b)
  cov.names[[length(cov.names)+1]] <- c("SLB+SHB+SLH")
  
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*H+S*L*B+S*H*B+L*H*B)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,h,0,s,l,b,0,s,h,b,0,l,h,b)
  cov.names[[length(cov.names)+1]] <- c("SLH+SLB+SHB+LHB")
  
  # 4-way interaction
  cov.model[[length(cov.model)+1]] <- update.formula(GWZ, ~.+S*L*B*H)
  cov.specs[[length(cov.specs)+1]] <- c(gwz,0,s,l,b,h)
  cov.names[[length(cov.names)+1]] <- c("SLBH")
  
  n.cov.models <- length(cov.model)
  
  # Fit statistieken toevoegen aan modellen.
  
  cov.model.norm          	<- matrix(nrow=n.cov.models,ncol=5)
  colnames(cov.model.norm)	<- c("nr","DEV","df","AIC","nHat0")
  rownames(cov.model.norm)	<- cov.names
  
  cov.pars <-list(n.cov.models)
  
  for (j in 1:n.cov.models){
    
    FIT <- dcat.fit(cov.model[[j]], cov.specs[[j]]) # applying the function created 
    # so everytime we do a new model, we also redo the EM imputation
    cov.model.norm[j,]	<-cbind(j,
                               round(FIT$norm.deviance),
                               FIT$d.o.f,
                               round(FIT$norm.aic,3),  
                               round(as.vector(FIT$tb.missed),1)) # performance of the model
    
    cov.pars[[j]] 	<-  round(FIT$tb.summary$coefficients[,c(1,2,4)],digits=2) 
    # we get the columns Estimate, Std.Error and Pr(>|z|)
  }
  # all models fitted
  
  
  xpar			<- matrix(0,n.cov.models,1)
  colnames(xpar)	<- c("xpar")
  
  for (j in 1:n.cov.models){
    sel	<- as.data.frame(cov.pars[[j]])
    xpar[j]<- nrow(sel[sel[,2]>5,]) # number of coefficients estimated that are larger than 5
  }
  
  identified.cov.models 	<- subset(cov.model.norm,xpar==0) # performance of models with no coeff larger than 5
  rang   			<- round(rank(identified.cov.models[,4]),0) # gives the position from lower value to larger for AIC
  cov.model.nr		 	<- subset(identified.cov.models[,1],rang==1) # number of the best model regarding AIC
  
  # zz is een index om verschillende modellen te printen, m1 is andere naam voor GWZ model.
  # So zz refers to the interactions that are included between the sources (there are 8 of these basis)
  # depending on the zz, we have different effects to retain?
  m1				<- cov.names[[cov.model.nr]]
  
  if(zz==1){
    gwz.effects <- 3
    gwz1<-gwz[1]
    gwz2<-gwz[3]
    gwz3<-gwz[5]}
  if(zz==2){
    gwz.effects <- 2
    gwz1<-gwz[1]
    gwz2<-gwz[3:4]}
  if(zz==3){
    gwz.effects <- 2
    gwz1<-gwz[1]
    gwz2<-gwz[3:4]}
  if(zz==4){
    gwz.effects<-2
    gwz1<-gwz[1]
    gwz2<-gwz[3:4]}
  if(zz==5){
    gwz.effects<-2
    gwz1<-gwz[1:2]
    gwz2<-gwz[4:5]}
  if(zz==6){
    gwz.effects<-2
    gwz1<-gwz[1:2]
    gwz2<-gwz[4:5]}
  if(zz==7){
    gwz.effects<-3
    gwz1<-gwz[1:2]
    gwz2<-gwz[4:5]
    gwz3<-gwz[7:8]}
  if(zz==8){
    gwz.effects<-2
    gwz1<-gwz[1:2]
    gwz2<-gwz[4:5]}
  if(zz==9){
    gwz.effects<-3
    gwz1<-gwz[1:2]
    gwz2<-gwz[4:5]
    gwz3<-gwz[7:8]}
  if(zz==10){
    gwz.effects<-1
    gwz1<-gwz[1:3]
  }
  

  
  COVx	<- cov.model[[cov.model.nr]] # best model
  COV    <- paste(update(COVx,~.-G*W*Z))[3] # we take out the three way interaction of the sources
  Mcov	<- paste(COV," + ",sep="")
  Mx	<- paste(update(COVx,~.-S*L*H))[3] # I'm guessing these are unlikely interactions
  # it took out S, H, L, S:H,H:L (there was initially no S:H:L nor S:L in the model, at least in 2022)
  Mind	<- paste("(",Mx,")",sep="") # same as Mx but in brackets
  lgwz 	<- 2+length(gwz) 
  cov    <- cov.specs[[cov.model.nr]][lgwz:length(cov.specs[[cov.model.nr]])]
  
  model <- list()
  specs <- list()
  names <- list()
  
  model[[1]]   		<- as.formula(paste("Freq~",Mcov,Mind,sep=""))
  specs[[1]] 			<- if(gwz.effects==2){c(gwz1,0,gwz2,0,cov)}else if(gwz.effects==3){
    c(gwz1,0,gwz2,0,gwz3,0,cov)}else{c(gwz1,0,cov)} 
  names[[1]] 			<- paste(m0,"+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","S",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,0,gwz2,s,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,0,gwz2,s,0,gwz3,s,0,cov)}else{c(gwz1,s,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(S)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","L",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,l,0,gwz2,l,0,cov)}else if(gwz.effects==3){
    c(gwz1,l,0,gwz2,l,0,gwz3,l,0,cov)}else{c(gwz1,l,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(L)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","B",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,b,0,gwz2,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,b,0,gwz2,b,0,gwz3,b,0,cov)}else{c(gwz1,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(B)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","H",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{c(gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S+L)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,0,gwz2,s,0,gwz1,l,0,gwz2,l,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,0,gwz2,s,0,gwz3,s,0,gwz1,l,0,gwz2,l,0,gwz3,l,0,cov)}else{c(gwz1,s,0,gwz1,l,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(S+L)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S+B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,0,gwz2,s,0,gwz1,b,0,gwz2,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,0,gwz2,s,0,gwz3,s,0,gwz1,b,0,gwz2,b,0,gwz3,b,0,cov)}else{c(gwz1,s,0,gwz1,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(S+B)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,0,gwz2,s,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,0,gwz2,s,0,gwz3,s,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{c(gwz1,s,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(S+H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(L+B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,l,0,gwz2,l,0,gwz1,b,0,gwz2,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,l,0,gwz2,l,0,gwz3,l,0,gwz1,b,0,gwz2,b,0,gwz3,b,0,cov)}else{c(gwz1,l,0,gwz1,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(L+B)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(L+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,l,0,gwz2,l,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,l,0,gwz2,l,0,gwz3,l,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{c(gwz1,l,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(L+H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(B+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,b,0,gwz2,b,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,b,0,gwz2,b,0,gwz3,b,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{c(gwz1,b,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(B+H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S+L+B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,0,gwz2,s,0,gwz1,l,0,gwz2,l,0,gwz1,b,0,gwz2,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,0,gwz2,s,0,gwz3,s,0,gwz1,l,0,gwz2,l,0,gwz3,l,0,gwz1,b,0,gwz2,b,0,gwz3,b,0,cov)}else{
      c(gwz1,s,0,gwz1,l,0,gwz1,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(S+L+B)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S+L+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,0,gwz2,s,0,gwz1,l,0,gwz2,l,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,0,gwz2,s,0,gwz3,s,0,gwz1,l,0,gwz2,l,0,gwz3,l,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{
      c(gwz1,s,0,gwz1,l,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(S+L+H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S+B+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,0,gwz2,s,0,gwz1,b,0,gwz2,b,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,0,gwz2,s,0,gwz3,s,0,gwz1,b,0,gwz2,b,0,gwz3,b,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{
      c(gwz1,s,0,gwz1,b,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(S+B+H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(L+B+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,l,0,gwz2,l,0,gwz1,b,0,gwz2,b,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,l,0,gwz2,l,0,gwz3,l,0,gwz1,b,0,gwz2,b,0,gwz3,b,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{
      c(gwz1,l,0,gwz1,b,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(L+B+H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S+L+B+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,0,gwz2,s,0,gwz1,l,0,gwz2,l,0,gwz1,b,0,gwz2,b,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,0,gwz2,s,0,gwz3,s,0,gwz1,l,0,gwz2,l,0,gwz3,l,0,gwz1,b,0,gwz2,b,0,gwz3,b,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{
      c(gwz1,s,0,gwz1,l,0,gwz1,b,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(S+L+B+H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,0,gwz2,s,l,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,0,gwz2,s,l,0,gwz3,s,l,0,cov)}else{c(gwz1,s,l,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SL)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L+B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,0,gwz2,s,l,0,gwz1,b,0,gwz2,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,0,gwz2,s,l,0,gwz3,s,l,0,gwz1,b,0,gwz2,b,0,gwz3,b,0,cov)}else{c(gwz1,s,l,0,gwz1,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SL+B)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,0,gwz2,s,l,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,0,gwz2,s,l,0,gwz3,s,l,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{c(gwz1,s,l,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SL+H)","+",m1,sep="")
  
  ##bij model hieronder gaat het mis, checken!3 next models were taken out of current code
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,b,0,gwz2,s,b,0,cov)}else if(gwz.effects==3){
			c(gwz1,s,b,0,gwz2,s,b,0,gwz3,s,b,0,cov)}else{c(gwz1,s,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*B+L)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,b,0,gwz2,s,b,0,gwz1,l,0,gwz2,l,0,cov)}else if(gwz.effects==3){
  					c(gwz1,s,b,0,gwz2,s,b,0,gwz3,s,b,0,gwz1,l,0,gwz2,l,0,gwz3,l,0,cov)}else{c(gwz1,s,b,0,gwz1,l,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SB+L)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*B+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,b,0,gwz2,s,b,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
  					c(gwz1,s,b,0,gwz2,s,b,0,gwz3,s,b,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{c(gwz1,s,b,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SB+H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*H+L)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,h,0,gwz2,s,h,0,gwz1,l,0,gwz2,l,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,h,0,gwz2,s,h,0,gwz3,s,h,0,gwz1,l,0,gwz2,l,0,gwz3,l,0,cov)}else{c(gwz1,s,h,0,gwz1,l,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SH+L)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*H+B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,h,0,gwz2,s,h,0,gwz1,b,0,gwz2,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,h,0,gwz2,s,h,0,gwz3,s,h,0,gwz1,b,0,gwz2,b,0,gwz3,b,0,cov)}else{c(gwz1,s,h,0,gwz1,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SH+B)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(L*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,l,b,0,gwz2,l,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,l,b,0,gwz2,l,b,0,gwz3,l,b,0,cov)}else{c(gwz1,l,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(LB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(L*B+S)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,l,b,0,gwz2,l,b,0,gwz1,s,0,gwz2,s,0,cov)}else if(gwz.effects==3){
    c(gwz1,l,b,0,gwz2,l,b,0,gwz3,l,b,0,gwz1,s,0,gwz2,s,0,gwz3,s,0,cov)}else{
      c(gwz1,l,b,0,gwz1,s,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(LB+S)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(L*B+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,l,b,0,gwz2,l,b,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,l,b,0,gwz2,l,b,0,gwz3,l,b,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{
      c(gwz1,l,b,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(LB+H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(L*H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,l,h,0,gwz2,l,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,l,h,0,gwz2,l,h,0,gwz3,l,h,0,cov)}else{c(gwz1,l,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(LH)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(L*H+S)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,l,h,0,gwz2,l,h,0,0,gwz1,s,0,gwz2,s,0,cov)}else if(gwz.effects==3){
    c(gwz1,l,h,0,gwz2,l,h,0,gwz3,l,h,0,gwz1,s,0,gwz2,s,0,gwz3,s,0,cov)}else{
      c(gwz1,l,h,0,gwz1,s,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(LH+S)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(L*H+B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,l,h,0,gwz2,l,h,0,0,gwz1,b,0,gwz2,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,l,h,0,gwz2,l,h,0,gwz3,l,h,0,gwz1,b,0,gwz2,b,0,gwz3,b,0,cov)}else{
      c(gwz1,l,h,0,gwz1,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(LH+B)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L+S*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,0,gwz2,s,l,0,gwz1,s,b,0,gwz2,s,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,0,gwz2,s,l,0,gwz3,s,l,0,gwz1,s,b,0,gwz2,s,b,0,gwz3,s,b,0,cov)}else{
      c(gwz1,s,l,0,gwz1,s,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SL+SB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L+S*H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,0,gwz2,s,l,0,gwz1,s,h,0,gwz2,s,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,0,gwz2,s,l,0,gwz3,s,l,0,gwz1,s,h,0,gwz2,s,h,0,gwz3,s,h,0,cov)}else{
      c(gwz1,s,l,0,gwz1,s,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SL+SH)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L+L*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,0,gwz2,s,l,0,gwz1,l,b,0,gwz2,l,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,0,gwz2,s,l,0,gwz3,s,l,0,gwz1,l,b,0,gwz2,l,b,0,gwz3,l,b,0,cov)}else{
      c(gwz1,s,l,0,gwz1,l,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SL+LB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L+L*H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,0,gwz2,s,l,0,gwz1,l,h,0,gwz2,l,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,0,gwz2,s,l,0,gwz3,s,l,0,gwz1,l,h,0,gwz2,l,h,0,gwz3,l,h,0,cov)}else{
      c(gwz1,s,l,0,gwz1,gwz1,l,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SL+LH)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L+B*H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,0,gwz2,s,l,0,gwz1,b,h,0,gwz2,b,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,0,gwz2,s,l,0,gwz3,s,l,0,gwz1,b,h,0,gwz2,b,h,0,gwz3,b,h,0,cov)}else{
      c(gwz1,s,l,0,gwz1,b,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SL+BH)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,h,0,gwz2,s,l,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,cov)}else{c(gwz1,s,l,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLH)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,b,0,gwz2,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,b,0,gwz2,b,0,gwz3,b,0,cov)}else{
      c(gwz1,s,l,h,0,gwz1,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLH+B)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+B*S)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,b,s,0,gwz2,b,s,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,b,s,0,gwz2,b,s,0,gwz3,b,s,0,cov)}else{
      c(gwz1,s,h,l,0,gwz1,b,s,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SHL+BS)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+B*L)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,b,l,0,gwz2,b,l,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,b,l,0,gwz2,b,l,0,gwz3,b,l,0,cov)}else{
      c(gwz1,s,h,l,0,gwz1,b,l,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SHL+BL)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+B*H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,b,h,0,gwz2,b,h,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,b,h,0,gwz2,b,h,0,gwz3,b,h,0,cov)}else{
      c(gwz1,s,l,h,0,gwz1,b,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLH+BH)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,b,0,gwz2,s,l,b,0,cov)}else if(gwz.effects==3){
  						c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz3,s,l,b,0,cov)}else{c(gwz1,s,l,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*B+H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz1,h,0,gwz2,h,0,cov)}else if(gwz.effects==3){
  						c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz3,s,l,b,0,gwz1,h,0,gwz2,h,0,gwz3,h,0,cov)}else{
  						  c(gwz1,s,l,b,0,gwz1,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLB+H)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*B+H*S)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz1,h,s,0,gwz2,h,s,0,cov)}else if(gwz.effects==3){
  					c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz3,s,l,b,0,gwz1,h,s,0,gwz2,h,s,0,gwz3,h,s,0,cov)}else{
  					  c(gwz1,s,l,b,0,gwz1,h,s,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLB+HS)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*B+H*L)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz1,h,l,0,gwz2,h,l,0,cov)}else if(gwz.effects==3){
  				c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz3,s,l,b,0,gwz1,h,l,0,gwz2,h,l,0,gwz3,h,l,0,cov)}else{
  				  c(gwz1,s,l,h,0,gwz1,h,l,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLB+HL)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*B+H*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz1,h,b,0,gwz2,h,b,0,cov)}else if(gwz.effects==3){
  					c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz3,s,l,b,0,gwz1,h,b,0,gwz2,h,b,0,gwz3,h,b,0,cov)}else{
  					  c(gwz1,s,l,b,0,gwz1,h,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLB+HB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*H*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,h,b,0,gwz2,s,h,b,0,cov)}else if(gwz.effects==3){
  					 c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,cov)}else{c(gwz1,s,h,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SHB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*H*B+L)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz1,l,0,gwz2,l,0,cov)}else if(gwz.effects==3){
  					c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,gwz1,l,0,gwz2,l,0,gwz3,l,0,cov)}else{
  					  c(gwz1,s,h,b,0,gwz1,l,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SHB+L)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*H*B+L*S)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz1,l,s,0,gwz2,l,s,0,cov)}else if(gwz.effects==3){
  					c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,gwz1,l,s,0,gwz2,l,s,0,gwz3,l,s,0,cov)}else{
  					  c(gwz1,s,h,b,0,gwz1,l,s,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SHB+LS)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*H*B+L*H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz1,l,h,0,gwz2,l,h,0,cov)}else if(gwz.effects==3){
  					c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,gwz1,l,h,0,gwz2,l,h,0,gwz3,l,h,0,cov)}else{
  					  c(gwz1,s,h,b,0,gwz1,l,h,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SHB+LH)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*H*B+L*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz1,l,b,0,gwz2,l,b,0,cov)}else if(gwz.effects==3){
  				 c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,gwz1,l,b,0,gwz2,l,b,0,gwz3,l,b,0,cov)}else{
  				   c(gwz1,s,h,b,0,gwz1,l,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SHB+LB)","+",m1,sep="")
  
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+L*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,l,b,0,gwz2,l,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,l,b,0,gwz2,l,b,0,gwz3,l,b,0,cov)}else{
      c(gwz1,s,l,h,0,gwz1,l,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLH+LB)","+",m1,sep="")
  
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+S*L*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,s,l,b,0,gwz2,s,l,b,0,cov)}else if(gwz.effects==3){
  						c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,s,l,b,0,gwz2,s,l,b,0,gwz3,s,l,b,0,cov)}else{
  						  c(gwz1,s,l,h,0,gwz1,s,l,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLH+SLB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+S*H*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,s,h,b,0,gwz2,s,h,b,0,cov)}else if(gwz.effects==3){
  				c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,cov)}else{
  				  c(gwz1,s,l,h,0,gwz1,s,h,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLH+SHB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*H*B+L*H*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz1,l,h,b,0,gwz2,l,h,b,0,cov)}else if(gwz.effects==3){
  				c(gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,gwz1,l,h,b,0,gwz2,l,h,b,0,gwz3,l,h,b,0,cov)}else{
  				  c(gwz1,s,h,b,0,gwz1,l,h,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SHB+LHB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+S*L*B+L*H*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){
    c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,s,l,b,0,gwz2,s,l,b,0,gwz1,l,h,b,0,gwz2,l,h,b,0,cov)}else if(gwz.effects==3){
    c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,s,l,b,0,gwz2,s,l,b,0,gwz3,s,l,b,0,gwz1,l,h,b,0,gwz2,l,h,b,0,gwz3,l,h,b,0,cov)}else{
      c(gwz1,s,l,h,0,gwz1,s,l,b,0,gwz1,l,h,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLH+SLB+LHB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+S*L*B+S*H*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){
  	c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,s,l,b,0,gwz2,s,l,b,0,gwz1,s,h,b,0,gwz2,s,h,b,0,cov)}else if(gwz.effects==3){
  	c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,s,l,b,0,gwz2,s,l,b,0,gwz3,s,l,b,0,gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,cov)}else{
  	  c(gwz1,s,l,h,0,gwz1,s,l,b,0,gwz1,s,h,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLH+SLB+SHB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+S*H*B+L*H*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){
  	c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,s,h,b,0,gwz2,s,h,b,0,gwz1,l,h,b,0,gwz2,l,h,b,0,cov)}else if(gwz.effects==3){
  	c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,gwz1,l,h,b,0,gwz2,l,h,b,0,gwz3,l,h,b,0,cov)}else{
  	  c(gwz1,s,l,h,0,gwz1,s,h,b,0,gwz1,l,h,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLH+SHB+LHB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*B+S*H*B+L*H*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){
  	c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz1,s,h,b,0,gwz2,s,h,b,0,gwz1,l,h,b,0,gwz2,l,h,b,0,cov)}else if(gwz.effects==3){
  	c(gwz1,s,l,b,0,gwz2,s,l,b,0,gwz3,s,l,b,0,gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,gwz1,l,h,b,0,gwz2,l,h,b,0,gwz3,l,h,b,0,cov)}else{
  	  c(gwz1,s,l,b,0,gwz1,s,h,b,0,gwz1,l,h,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLB+SHB+LHB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*H+S*L*B+S*H*B+L*H*B)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){
  	c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz1,s,l,b,0,gwz2,s,l,b,0,gwz1,s,h,b,0,gwz2,s,h,b,0,gwz1,l,h,b,0,gwz2,l,h,b,0,cov)}else if(gwz.effects==3){
  	c(gwz1,s,l,h,0,gwz2,s,l,h,0,gwz3,s,l,h,0,gwz1,s,l,b,0,gwz2,s,l,b,0,gwz3,s,l,b,0,gwz1,s,h,b,0,gwz2,s,h,b,0,gwz3,s,h,b,0,gwz1,l,h,b,0,gwz2,l,h,b,0,gwz3,l,h,b,0,cov)}else{
  	  c(gwz1,s,l,h,0,gwz1,s,l,b,0,gwz1,s,h,b,0,gwz1,l,h,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLH+SLB+SHB+LHB)","+",m1,sep="")
  
  model[[length(model)+1]] <- as.formula((paste("Freq~",Mcov,Mind,"*","(S*L*B*H)",sep="")))
  specs[[length(specs)+1]] <- if(gwz.effects==2){
    c(gwz1,s,l,h,b,0,gwz2,s,l,h,b,0,cov)}else if(gwz.effects==3){
      c(gwz1,s,l,h,b,0,gwz2,s,l,h,b,0,gwz3,s,l,h,b,0,cov)}else{
        c(gwz1,s,l,h,b,0,cov)}
  names[[length(names)+1]] <- paste(m0,"(SLBH)","+",m1,sep="")
  
  n.models <- length(model)
  
  pars 			<- list(n.models)
  freqs			<- list(n.models)
  
  for (j in 1:n.models){
    
    FIT <- dcat.fit(model[[j]], specs[[j]])
    model.norm[j,]$nr=j
    model.norm[j,]$DEV=round(FIT$norm.deviance)
    model.norm[j,]$df=FIT$d.o.f
    model.norm[j,]$AIC=round(FIT$norm.aic,3)
    model.norm[j,]$nHat0=round(as.vector(FIT$tb.missed),1)
    model.norm[j,]$model=as.character(model[[j]])[3]
    model.norm[j,]$specs=paste(as.character(specs[[j]]),collapse=" ")
    
    pars[[j]] 		<- round(FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)
    freqs[[j]] 		<- cbind(tb[,1:length(tb)-1],round(FIT$tb.ecm.fitted,2))
  }
  rownames(model.norm)	<- names
  
  
  xpar			<- matrix(0,n.models,1)
  colnames(xpar)	<- c("xpar")
  
  for (j in 1:n.models){
    sel	<- as.data.frame(pars[[j]])
    xpar[j]<- nrow(sel[sel[,2]>5,])
  }

  if(sum(xpar==0)==length(xpar)){identified.models = model.norm}else{identified.models 	<- subset(model.norm,xpar==0)}
  rang   		<- round(rank(identified.models[,4]),0)
  rg3 = which(rang<4)
  
  final.models  	<- rbind(final.models,identified.models[rg3,])
  model.norm= model.norm[0,] # empty model.norm
}

final.models[,1] 	<- rank(final.models[,4])
best10 		<- data.frame(subset(final.models,final.models[,1]<11))
best10[order(best10$nr),-c(6,7)]

# best model first step:
bm = which(final.models[,1]==1)
bestmodel = final.models[bm,]

# best model specs
best.specs = as.numeric(unlist(regmatches(bestmodel$specs, gregexpr("[[:digit:]]+",bestmodel$specs))))



test.model <- list()
test.specs <- list()
test.names <- list()


# loop to apply stepwise selection
# ideally, this would include first deleting the 4 way interactions, then 3 way, etc
# ideally, we would also not exclude interactions when their terms are included in a higher interaction term

start.model = as.formula(paste("Freq~",bestmodel$model))
start.specs = best.specs
start.names = rownames(bestmodel)



start.FIT <- dcat.fit(start.model, start.specs)
start.pars 		<- round(start.FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)

start.pars
pars.names = rownames(start.pars)
# first let's transform this list into a list of interactions
start.coef = unique(gsub("[^a-zA-Z]","",pars.names)) [-1]

remove = start.coef[nchar(start.coef)>=5]
formula_remove = sapply(remove, insert_colon)
form_remove = paste(formula_remove, collapse= " -")
start2.model = update.formula(start.model, paste("~.-",form_remove ))
start2.names = model2names(start2.model)
start2.specs = names2specs(start2.names, values=values)

start2.FIT <- dcat.fit(start2.model, start2.specs)
start2.pars 		<- round(start2.FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)

start2.pars
pars2.names = rownames(start2.pars)
# first let's transform this list into a list of interactions
start2.coef = unique(gsub("[^a-zA-Z]","",pars2.names))

test.model[[1]] = start2.model
test.specs[[1]] = names2specs(start2.coef[-1], values=values) # unfortunately we do not get quite the same models but I do not know why
test.names[[1]] = paste(start2.coef[-1], collapse=" ")

# now we want to focus on the 4-ways parameters
coef4 = start2.coef[nchar(start2.coef)==4]
formula4 = sapply(coef4, insert_colon)
if(length(coef4)==0){
  best4.model = test.model[[1]]
  best4.specs = test.specs[[1]]
  best4.names = test.names[[1]]
}
if(length(coef4)!=0){
  AIC = numeric(length(coef4)+1)
  AIC[1]=start2.FIT$norm.aic
  # while the preferred one is a new one, should stop if it is the last one created
  for (i in 1:length(coef4)){
    test.model[[i+1]]= update.formula(test.model[[1]],paste("~.-",formula4[i]))
    test.specs[[i+1]]= names2specs(start2.coef[-c(1,which(start2.coef==coef4[i]))], values=values)
    test.names[[i+1]] = model2names(test.model[[i+1]])
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
        test.names[[j]] = model2names(test.model[[j]])
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
AIC3 = numeric(length(coef3)+1)

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
    test3.names[[i+1]]= model2names(test3.model[[i+1]])
    AIC3[i+1]=dcat.fit(test3.model[[i+1]], test3.specs[[i+1]])$norm.aic
  }
  if(which.min(AIC3)==1){
    best3.model = test3.model[[1]]
    best3.specs = test3.specs[[1]]
    best3.names = test3.names[[1]]
    newcoef = unlist(strsplit(test3.names[[1]], split=" "))
  }else{
    min.AIC3 = which.min(AIC3)
    best3.model = test3.model[[min.AIC3]]
    best3.specs = test3.specs[[min.AIC3]]
    best3.names = test3.names[[min.AIC3]]
    newcoef = unlist(strsplit(test3.names[[min.AIC3]], split=" "))
    coef3 = newcoef[nchar(newcoef)==3]
    if(sum(out_test)!=0){coef3 = coef3[-which(coef3 %in% exclude)]}
    formula3 = sapply(coef3, insert_colon)
    last.min = 0
    while(length(coef3)>0){
      j=length(test3.model)
      min.AIC3 = which.min(AIC3)
      for (i in 1:length(coef3)){
        j=j+1
        test3.model[[j]]= update.formula(best3.model,paste("~.-",formula3[i]))
        test3.specs[[j]]= model2specs(test3.model[[j]], values=values)
        test3.names[[j]] = model2names(test3.model[[j]])
        AIC3[j]=dcat.fit(test3.model[[j]], test3.specs[[j]])$norm.aic
      }
      AIC3 = round(AIC3, 2)
      if ((min.AIC3 == which.min(AIC3))&(length(which(AIC3== AIC3[min.AIC3]))==1)){break} # if after the loop the best model is still the one we had before the loop,
      # then it means that we want to keep this model
      if ((min.AIC3 == which.min(AIC3))&(length(which(AIC3== AIC3[min.AIC3]))>1)){
        if(last.min == which(AIC3== AIC3[min.AIC3])[length(which(AIC3 == AIC3[min.AIC3]))]){break}
        best3.model = test3.model[[which(AIC3== AIC3[min.AIC3])[length(which(AIC3 == AIC3[min.AIC3]))]]]
        best3.specs = test3.specs[[which(AIC3== AIC3[min.AIC3])[length(which(AIC3 == AIC3[min.AIC3]))]]]
        best3.names = test3.names[[which(AIC3== AIC3[min.AIC3])[length(which(AIC3 == AIC3[min.AIC3]))]]]
        newcoef = unlist(strsplit(best3.names, split=" "))
        coef3 = newcoef[nchar(newcoef)==3]
        if(sum(out_test)!=0){coef3 = coef3[-which(coef3 %in% exclude)]}
        formula3 = sapply(coef3, insert_colon)
        last.min = which(AIC3== AIC3[min.AIC3])[length(which(AIC3 == AIC3[min.AIC3]))]
      } else {
        best3.model = test3.model[[which.min(AIC3)]]
        best3.specs = test3.specs[[which.min(AIC3)]]
        best3.names = test3.names[[which.min(AIC3)]]
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
AIC2 = numeric(length(coef2)+1)

best3.FIT <- dcat.fit(best3.model, best3.specs)
test2.model[[1]]= best3.model
test2.specs[[1]]= best3.specs
test2.names[[1]]= best3.names

best3.coef = unlist(strsplit(best3.names, split=" "))

AIC2[1]=best3.FIT$norm.aic
# while the preferred one is a new one, should stop if it is the last one created
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
    AIC2[i+1]=dcat.fit(test2.model[[i+1]], test2.specs[[i+1]])$norm.aic
  }
  if(which.min(AIC2)==1){
    best2.model = test2.model[[1]]
    best2.specs = test2.specs[[1]]
    best2.names = test2.names[[1]]
  }else{
    min.AIC2 = which.min(AIC2)
    best2.model = test2.model[[min.AIC2]]
    best2.specs = test2.specs[[min.AIC2]]
    best2.names = test2.names[[min.AIC2]]
    newcoef = unlist(strsplit(test2.names[[min.AIC2]], split=" "))
    coef2 = newcoef[nchar(newcoef)==2]
    if(sum(out_test)!=0){coef2 = coef2[-which(coef2 %in% exclude)]} # keeping only interactions that we want out
    formula2 = sapply(coef2, insert_colon)
    last.min = 0
    while(length(coef2)>0){
      j=length(test2.model)
      min.AIC2 = which.min(AIC2)
      for (i in 1:length(coef2)){
        j=j+1
        test2.model[[j]]= update.formula(best2.model,paste("~.-",formula2[i]))
        test2.specs[[j]]= model2specs(test2.model[[j]], values = values)
        test2.names[[j]] = model2names(test2.model[[j]])
        AIC2[j]=dcat.fit(test2.model[[j]], test2.specs[[j]])$norm.aic
      }
      AIC2 = round(AIC2,2)
      if ((min.AIC2 == which.min(AIC2))&(length(which(AIC2== AIC2[min.AIC2]))==1)){break} # if after the loop the best model is still the one we had before the loop,
      # then it means that we want to keep this model
      if ((min.AIC2 == which.min(AIC2))&(length(which(AIC2== AIC2[min.AIC2]))>1)){
        if(last.min == which(AIC2== AIC2[min.AIC2])[length(which(AIC2 == AIC2[min.AIC2]))]){break}
        best2.model = test2.model[[which(AIC2== AIC2[min.AIC2])[length(which(AIC2 == AIC2[min.AIC2]))]]]
        best2.specs = test2.specs[[which(AIC2== AIC2[min.AIC2])[length(which(AIC2 == AIC2[min.AIC2]))]]]
        best2.names = test2.names[[which(AIC2== AIC2[min.AIC2])[length(which(AIC2 == AIC2[min.AIC2]))]]]
        newcoef = unlist(strsplit(best2.names, split=" "))
        coef2 = newcoef[nchar(newcoef)==2]
        if(sum(out_test)!=0){coef2 = coef2[-which(coef2 %in% exclude)]}
        formula2 = sapply(coef2, insert_colon)
        last.min = which(AIC2== AIC2[min.AIC2])[length(which(AIC2 == AIC2[min.AIC2]))]
      } else {
        best2.model = test2.model[[which.min(AIC2)]]
        best2.specs = test2.specs[[which.min(AIC2)]]
        best2.names = test2.names[[which.min(AIC2)]]
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
AIC2[min.AIC2] # in 2023 we get a model with an AIC of around 565.98

best2.FIT <- dcat.fit(best2.model, best2.specs)
best2.pars = round(best2.FIT$tb.summary$coefficients[,c(1,2,4)],digits=2)
best2.freqs = cbind(tb[,1:length(tb)-1],round(best2.FIT$tb.ecm.fitted,2))
coefs.origineel = best2.FIT$tb.summary$coefficients
tb.ecm.origineel = best2.FIT$tb.ecm

fitted = best2.freqs
colnames(fitted)	= c("G","W","Z","L","H","S","B","obs","est")

est.pop.size = sum(fitted$est)
obs.pop.size = sum(fitted$obs)

dark.figure = est.pop.size-obs.pop.size
results[mc,2]=dark.figure
results[mc,3]=est.pop.size

}

colMeans(results[,c("n000", "N_est_standard_cov")]) 
colSds(  results[,c("n000", "N_est_standard_cov")])/(MC^0.5) 

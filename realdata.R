rm(list = ls())

library(MASS)
library(INLA)
library(dplyr)
options(warn = -1)
set.seed(1)
load("./RData/fulldata.RData")
load("./RData/map.RData")


################################################################################
# Table D.1 in Supplementary



resfinal1 <- data.frame(matrix(0, nrow = 11+328*2, ncol = 100))
resfinal2 <- data.frame(matrix(0, nrow = 11+328, ncol = 100))
resfinal3 <- data.frame(matrix(0, nrow = 11, ncol = 100))
resfinal4 <- data.frame(matrix(0, nrow = 11+328*2, ncol = 100))

confint<-function(x,sigma=-1,alpha=0.05)
{
  n<-length(x)
  xb<-mean(x)
  if(sigma>=0)
  {
    tmp<-sigma*qnorm(1-alpha/2);df<-n
  }
  else{
    tmp<-sd(x)*qt(1-alpha/2,n-1);df<- n-1
  }
  data.frame(mean=xb,sd=sd(x),df=df,a=xb-tmp,b=xb+tmp)
}

# varying
i=1
while (i < 101) {
  print(paste0(1, ":", i))
  formula <- respiratory ~ smoking + fulldata[, paste0("temperature", i)] + fulldata[, paste0("rh", i)] + week + 
    hosdense + age + summer + autumn + winter + fulldata[, paste0("PM2.5", i)] + 
    f(idx1, fulldata[, paste0("PM2.5", i)], model="besag", graph=W.join, hyper=list(prec = list(prior = "loggamma", param = c(1, 10)))) +
    f(idx2, week, model="besag", graph=W.join, hyper=list(prec = list(prior = "loggamma", param = c(1, 1))))
  res <- try(inla(formula, data = fulldata, family="poisson", E=total),silent = T)
  if('try-error'%in%class(res)){
    next
  }else{
    resfinal1[(1:11), i] <- res$summary.fixed$mean
    resfinal1[(12:(11+328)), i] <- res$summary.random$idx1$mean
    resfinal1[((12+328):(11+328*2)), i] <- res$summary.random$idx2$mean
    i=i+1
  }
}

result1 <- as.matrix(resfinal1)
tableD <- rbind(confint(result1[1,]),confint(result1[2,]),confint(result1[3,]),confint(result1[4,]),
                confint(result1[5,]),confint(result1[6,]),confint(result1[7,]),confint(result1[8,]),
                confint(result1[9,]),confint(result1[10,]),confint(result1[11,]))


################################################################################
# Table C in Supplementary

# fixed coefficient
i=1
while (i < 101) {
  print(paste0(2, ":", i))
  formula <- respiratory ~ smoking + fulldata[, paste0("temperature", i)] + fulldata[, paste0("rh", i)] + week + 
    hosdense + age + summer + autumn + winter + fulldata[, paste0("PM2.5", i)] + 
    f(idx2, week, model="besag", graph=W.join, hyper=list(prec = list(prior = "loggamma", param = c(1, 1))))
  res <- try(inla(formula, data = fulldata, family="poisson", E=total),silent = T)
  if('try-error'%in%class(res)){
    next
  }else{
    resfinal2[(1:11), i] <- res$summary.fixed$mean
    resfinal2[(12:(11+328)), i] <- res$summary.random$idx2$mean
    i=i+1
  }
}

result2 <- as.matrix(resfinal2)
tableC3 <- rbind(confint(result2[1,]),confint(result2[2,]),confint(result2[3,]),confint(result2[4,]),
                 confint(result2[5,]),confint(result2[6,]),confint(result2[7,]),confint(result2[8,]),
                 confint(result2[9,]),confint(result2[10,]),confint(result2[11,]))

# bakar's model

### Spatial structure of Bakar's model
coordinate <- region_data[ids,2:3]
dists <- as.matrix(dist(coordinate))
phi_bakar <- 0.5
covar <- exp(-dists / phi_bakar)

i=1
while (i < 101) {
  print(paste0(4, ":", i))
  formula <- eta ~ smoking + fulldata[, paste0("temperature", i)] + fulldata[, paste0("rh", i)] + week + 
    hosdense + age + summer + autumn + winter + fulldata[, paste0("PM2.5", i)] + 
    f(idx1, fulldata[, paste0("PM2.5", i)], model="generic0", Cmatrix=covar, constr=TRUE,
      hyper=list(prec=list(prior="loggamma", param = c(1, 10)))) +
    f(idx2, week, model="generic0", Cmatrix=covar, constr=TRUE,
      hyper=list(prec = list(prior = "loggamma", param = c(1, 1))))
  res <- try(inla(formula, data = fulldata, family="gaussian", E=total),silent = T)
  if('try-error'%in%class(res)){
    next
  }else{
    resfinal4[(1:11), i] <- res$summary.fixed$mean
    resfinal4[(12:(11+328)), i] <- res$summary.random$idx1$mean
    resfinal4[((12+328):(11+328*2)), i] <- res$summary.random$idx2$mean
    i=i+1
  }
}

result4 <- as.matrix(resfinal4)
tableC4 <- rbind(confint(result4[1,]),confint(result4[2,]),confint(result4[3,]),confint(result4[4,]),
                 confint(result4[5,]),confint(result4[6,]),confint(result4[7,]),confint(result4[8,]),
                 confint(result4[9,]),confint(result4[10,]),confint(result4[11,]))


finalresult <- list(Varying = resfinal1,
                    Nonvarying = resfinal2,
                    NonSpatialEffect = resfinal3,
                    Bakar = resfinal4)


#save(finalresult,file="./paper.RData")
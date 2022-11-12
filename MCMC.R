rm(list = ls())

library(MASS)
library(INLA)
library(CARBayesST)
library(sf)
library(dplyr)
#options(warn = -1)
set.seed(0)
load("./RData/fulldata.RData")
load("./RData/map.RData")

################################################################################
# Table C.2 in Supplementary



fulldata$PM2.5_mean <- rowMeans(fulldata[,35:134])
fulldata$temperature_mean <- rowMeans(fulldata[,135:234])
fulldata$rh_mean <- rowMeans(fulldata[,235:334])

fulldata <- fulldata[order(fulldata$week,fulldata$idx1),]

formula <- respiratory ~ offset(log(total)) + smoking + temperature_mean + rh_mean +
  hosdense + age + summer + autumn + winter + PM2.5_mean

start <- Sys.time()
chain1 <- ST.CARlinear(formula, family = "poisson", data=fulldata, W=W[ids,ids],
                      burnin = 20000,
                      n.sample=220000, thin=100,
                      prior.tau2=c(1, 1), rho.slo=1, rho.int=1)
end <- Sys.time()
runtime1 = end - start

start <- Sys.time()
chain2 <- ST.CARanova(formula, family = "poisson", data=fulldata, W=W[ids,ids],
                      interaction=FALSE, burnin = 20000,
                      n.sample=220000, thin=100,
                      prior.tau2=c(1, 1), rho.S=1, rho.T=1)
end <- Sys.time()
runtime2 = end - start



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

tableC2 <- rbind(confint(chain1$samples$beta[,1]),confint(chain1$samples$beta[,2]),confint(chain1$samples$beta[,3]),confint(chain1$samples$beta[,4]),
                confint(chain1$samples$beta[,5]),confint(chain1$samples$beta[,6]),confint(chain1$samples$beta[,7]),confint(chain1$samples$beta[,8]),
                confint(chain1$samples$beta[,9]),confint(chain1$samples$beta[,10]),confint(chain1$samples$beta[,11]))

#save(chain1,runtime1,chain2,runtime2,file="./MCMC_result.RData")


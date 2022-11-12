rm(list = ls())

# load the region data
load("./RData/map.RData")

# packages
library(MASS)
library(INLA)

################################################################################
# Table B.2 in Supplementary



################################################################################
### fixed β_1 and varying β_2
################################################################################
one_varying <- function(M, rho, E_mean1, E_sd1, adjmatrix, coordinate) {
  # Week numbers
  times <- 21
  
  # Graph
  n <- nrow(adjmatrix)
  W <- as.matrix(adjmatrix)
  g <- inla.read.graph(W)
  
  # Result
  model_v <- model_nv <- model_ns <- model_bakar <- matrix(0, nrow=M*n, ncol=4)
  colnames(model_v) <- colnames(model_nv) <- colnames(model_ns) <- 
    colnames(model_bakar) <-  c("eta_gen","eta_data", "eta_pred", "Times")
  
  # Spatial structure of Bakar's model
  dists <- as.matrix(dist(coordinate))
  phi_bakar <- 0.5
  covar <- exp(-dists / phi_bakar)
  
  # Fixed coefficients 
  B0 <- -3
  B1 <- 0.5
  alpha <- 0.01
  
  # Generate beta_2
  d1 <- 0.5
  tau1 <- 10
  Q1 <- matrix(0, g$n, g$n)
  diag(Q1) <- tau1 * (d1 + g$nnbs)
  for(i in 1:g$n) {
    if (g$nnbs[i] > 0) {
      Q1[i, g$nbs[[i]]] <- -tau1
      Q1[g$nbs[[i]], i] <- -tau1
    }
  }
  R1 <- chol(Q1)
  B2 <- backsolve(R1, rnorm(g$n))
  
  # Generate phi_i
  d2 <- 1
  tau2 <- 1000
  Q2 <- matrix(0, g$n, g$n)
  diag(Q2) <- tau2 * (d2 + g$nnbs)
  for(i in 1:g$n) {
    if (g$nnbs[i] > 0) {
      Q2[i, g$nbs[[i]]] <- -tau2
      Q2[g$nbs[[i]], i] <- -tau2
    }
  }
  R2 <- chol(Q2)
  phi <- backsolve(R2, rnorm(g$n))
  
  # Index
  id1 <- rep(1:g$n, times=times)
  week <- rep(1:times, each=g$n)
  
  # Repeat M times
  for(t in 1:M) {
    print(paste("One varying","*********",t,"*********"))
    set.seed(t)
    
    ### Generate X_1 & X_2
    S <- matrix(c(1, rho, rho, 1), nrow=2)
    x <- mvrnorm(n * times, mu=rep(0, 2), Sigma=S)
    x1 <- x[, 1]
    x2 <- x[, 2]
    
    ### Generate the true disease rates
    lin <- B0 + B1 * x1 + B2 * x2 + (phi + alpha) * week
    eta <- exp(lin)
    
    ### Generate n_{it}
    E1 <- rnorm(n, mean=E_mean1, sd=E_sd1)  
    E1 <- round(E1, 0)
    lambda1 <- E1 * eta
    
    ### Generate Y_{it}
    y1 <- rpois(n=n * times, lambda=lambda1)
    
    ### Dataset
    dt1 <- data.frame(y1=y1, s.index1=id1, s.index2=id1,s.index3=id1,
                      week=week, E1=E1, x1=x1, x2=x2)
    dt1 <- within(dt1, eta <- log(y1 / E1))
    
    x1_pred <- dt1[dt1$week == 21, ]$x1
    x2_pred <- dt1[dt1$week == 21, ]$x2
    
    ############################################################################
    # Varying-coefficients + random effects
    formula <- y1 ~ x1 + week +
      f(s.index1, x2, model="besag", graph=g,
        hyper=list(prec=list(prior="loggamma", param=c(0.5, 2)))) + 
      f(s.index2, week, model="besag", graph=g,
        hyper=list(prec=list(prior="loggamma", param=c(1, 1))))
    res <- inla(formula, data=dt1[dt1$week != 21,], family="poisson", E=E1)
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.fixed$mean[2]
    b2_est <- res$summary.random$s.index1$mean
    b3_est <- res$summary.fixed$mean[3]
    b4_est <- res$summary.random$s.index2$mean
    eta_pred <- exp(b0_est + b1_est * x1_pred + b2_est * x2_pred + (b3_est + b4_est) * times)
    model_v[((t-1)*n+1):(t*n), "Times"] <- t
    model_v[((t-1)*n+1):(t*n), "eta_gen"] <- tail(eta,n)
    model_v[((t-1)*n+1):(t*n), "eta_data"] <- exp(dt1[dt1$week == 21,]$eta)
    model_v[((t-1)*n+1):(t*n), "eta_pred"] <- eta_pred
    ############################################################################
    
    ############################################################################
    # Nonvarying-coefficients + random effects
    formula <- y1 ~ x1 + x2 + week +
      f(s.index1, week, model="besag", graph=g,
        hyper=list(prec=list(prior="loggamma", param=c(1, 1))))
    res <- inla(formula, data=dt1[dt1$week != 21,], family="poisson", E=E1)
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.fixed$mean[2]
    b2_est <- res$summary.fixed$mean[3]
    b3_est <- res$summary.fixed$mean[4]
    b4_est <- res$summary.random$s.index1$mean
    eta_pred <- exp(b0_est + b1_est * x1_pred + b2_est * x2_pred + (b3_est + b4_est) * times)
    model_nv[((t-1)*n+1):(t*n), "Times"] <- t
    model_nv[((t-1)*n+1):(t*n), "eta_gen"] <- tail(eta,n)
    model_nv[((t-1)*n+1):(t*n), "eta_data"] <- exp(dt1[dt1$week == 21,]$eta)
    model_nv[((t-1)*n+1):(t*n), "eta_pred"] <- eta_pred
    ############################################################################    
    
    ############################################################################
    # Nonvarying-coefficients + no random effects
    formula <- y1 ~ x1 + x2 + week
    res <- inla(formula, data=dt1[dt1$week != 21,], family="poisson", E=E1)
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.fixed$mean[2]
    b2_est <- res$summary.fixed$mean[3]
    b3_est <- res$summary.fixed$mean[4]
    eta_pred <- exp(b0_est + b1_est * x1_pred + b2_est * x2_pred + b3_est * times)
    model_ns[((t-1)*n+1):(t*n), "Times"] <- t
    model_ns[((t-1)*n+1):(t*n), "eta_gen"] <- tail(eta,n)
    model_ns[((t-1)*n+1):(t*n), "eta_data"] <- exp(dt1[dt1$week == 21,]$eta)
    model_ns[((t-1)*n+1):(t*n), "eta_pred"] <- eta_pred
    ############################################################################     
    
    ############################################################################
    # Bakar's model
    formula <- eta ~ x1 + week +
      f(s.index1, x2, model="generic0", Cmatrix=covar, constr=TRUE,
        hyper=list(prec=list(prior="loggamma", param=c(0.5, 2)))) +
      f(s.index2, model="generic0", Cmatrix=covar, constr=TRUE,
        hyper=list(prec=list(prior="loggamma", param=c(1, 1))))
    res <- inla(formula, data=dt1[dt1$week != 21,], family="gaussian")
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.fixed$mean[2]
    b2_est <- res$summary.random$s.index1$mean
    b3_est <- res$summary.fixed$mean[3]
    b4_est <- res$summary.random$s.index2$mean
    eta_pred <- exp(b0_est + b1_est * x1_pred + b2_est * x2_pred + b3_est * times + b4_est)
    model_bakar[((t-1)*n+1):(t*n), "Times"] <- t
    model_bakar[((t-1)*n+1):(t*n), "eta_gen"] <- tail(eta,n)
    model_bakar[((t-1)*n+1):(t*n), "eta_data"] <- exp(dt1[dt1$week == 21,]$eta)
    model_bakar[((t-1)*n+1):(t*n), "eta_pred"] <- eta_pred
    ############################################################################ 
  }
  # Result 
  result <- list(model_Varying = model_v, model_Nonvarying = model_nv,
                 model_NoRandomEffect = model_ns, model_Bakar = model_bakar)
  
  MSE_1 <- list(MSE_Varying = mean((model_v[,"eta_data"] - model_v[,"eta_pred"])^2),
                MSE_Nonvarying = mean((model_nv[,"eta_data"] - model_nv[,"eta_pred"])^2),
                MSE_NoRandomEffect = mean((model_ns[,"eta_data"] - model_ns[,"eta_pred"])^2),
                MSE_Bakar = mean((model_bakar[,"eta_data"] - model_bakar[,"eta_pred"])^2))
  
  MSE_2 <- list(MSE_Varying = mean((model_v[,"eta_gen"] - model_v[,"eta_pred"])^2),
                MSE_Nonvarying = mean((model_nv[,"eta_gen"] - model_nv[,"eta_pred"])^2),
                MSE_NoRandomEffect = mean((model_ns[,"eta_gen"] - model_ns[,"eta_pred"])^2),
                MSE_Bakar = mean((model_bakar[,"eta_gen"] - model_bakar[,"eta_pred"])^2))
  
  # Output
  return(list(Result = result, MSE_data = MSE_1, MSE_gen = MSE_2))
}

################################################################################
### Varying β_1 and β_2 
################################################################################
two_varying <- function(M, rho, E_mean1, E_sd1, adjmatrix, coordinate) {
  # Week numbers
  times <- 21
  
  # Graph
  n <- nrow(adjmatrix)
  W <- as.matrix(adjmatrix)
  g <- inla.read.graph(W)
  
  # Result
  model_v <- model_nv <- model_ns <- model_bakar <- matrix(0, nrow=M*n, ncol=4)
  colnames(model_v) <- colnames(model_nv) <- colnames(model_ns) <- 
    colnames(model_bakar) <-  c("eta_gen","eta_data", "eta_pred", "Times")
  
  # Spatial structure of Bakar's model
  dists <- as.matrix(dist(coordinate))
  phi_bakar <- 0.5
  covar <- exp(-dists / phi_bakar)
  
  # Fixed coefficients 
  B0 <- -3
  alpha <- 0.01
  
  # Generate beta_1
  d1 <- 0.5
  tau1 <- 10
  Q1 <- matrix(0, g$n, g$n)
  diag(Q1) <- tau1 * (d1 + g$nnbs)
  for(i in 1:g$n) {
    if (g$nnbs[i] > 0) {
      Q1[i, g$nbs[[i]]] <- -tau1
      Q1[g$nbs[[i]], i] <- -tau1
    }
  }
  R1 <- chol(Q1)
  B1 <- backsolve(R1, rnorm(g$n))
  
  # Generate beta_2
  d2 <- 2
  tau2 <- 5
  Q2 <- matrix(0, g$n, g$n)
  diag(Q2) <- tau2 * (d2 + g$nnbs)
  for(i in 1:g$n) {
    if (g$nnbs[i] > 0) {
      Q2[i, g$nbs[[i]]] <- -tau2
      Q2[g$nbs[[i]], i] <- -tau2
    }
  }
  R2 <- chol(Q2)
  B2 <- backsolve(R2, rnorm(g$n))
  
  # Generate phi_i
  d3 <- 1
  tau3 <- 1000
  Q3 <- matrix(0, g$n, g$n)
  diag(Q3) <- tau3 * (d3 + g$nnbs)
  for(i in 1:g$n) {
    if (g$nnbs[i] > 0) {
      Q3[i, g$nbs[[i]]] <- -tau3
      Q3[g$nbs[[i]], i] <- -tau3
    }
  }
  R3 <- chol(Q3)
  phi <- backsolve(R3, rnorm(g$n))
  
  # Index
  id1 <- rep(1:g$n, times=times)
  week <- rep(1:times, each=g$n)
  
  # Repeat M times
  for(t in 1:M) {
    print(paste("Two varying","*********",t,"*********"))
    set.seed(t)
    
    ### Generate X_1 & X_2
    S <- matrix(c(1, rho, rho, 1), nrow=2)
    x <- mvrnorm(n * times, mu=rep(0, 2), Sigma=S)
    x1 <- x[, 1]
    x2 <- x[, 2]
    
    ### Generate the true disease rates
    lin <- B0 + B1 * x1 + B2 * x2 + (phi + alpha) * week
    eta <- exp(lin)
    
    ### Generate n_{it}
    E1 <- rnorm(n, mean=E_mean1, sd=E_sd1)  
    E1 <- round(E1, 0)
    lambda1 <- E1 * eta
    
    ### Generate Y_{it}
    y1 <- rpois(n=n * times, lambda=lambda1)
    
    ### Dataset
    dt1 <- data.frame(y1=y1, s.index1=id1, s.index2=id1,s.index3=id1,
                      week=week, E1=E1, x1=x1, x2=x2)
    dt1 <- within(dt1, eta <- log(y1 / E1))
    
    x1_pred <- dt1[dt1$week == 21, ]$x1
    x2_pred <- dt1[dt1$week == 21, ]$x2
    
    ############################################################################
    # Varying-coefficients + random effects
    formula <- y1 ~ week +
      f(s.index1, x1, model="besag", graph=g,
        hyper=list(prec=list(prior="loggamma", param=c(0.5, 1)))) +
      f(s.index2, x2, model="besag", graph=g,
        hyper=list(prec=list(prior="loggamma", param=c(2, 1)))) +
      f(s.index3, week, model="besag", graph=g,
        hyper=list(prec=list(prior="loggamma", param=c(1, 1))))
    res <- inla(formula, data=dt1[dt1$week != 21,], family="poisson", E=E1)
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.random$s.index1$mean
    b2_est <- res$summary.random$s.index2$mean
    b3_est <- res$summary.fixed$mean[2]
    b4_est <- res$summary.random$s.index3$mean
    eta_pred <- exp(b0_est + b1_est * x1_pred + b2_est * x2_pred + (b3_est + b4_est) * times)
    model_v[((t-1)*n+1):(t*n), "Times"] <- t
    model_v[((t-1)*n+1):(t*n), "eta_gen"] <- tail(eta,n)
    model_v[((t-1)*n+1):(t*n), "eta_data"] <- exp(dt1[dt1$week == 21,]$eta)
    model_v[((t-1)*n+1):(t*n), "eta_pred"] <- eta_pred
    ############################################################################
    
    ############################################################################
    # Nonvarying-coefficients + random effects
    formula <- y1 ~ x1 + x2 + week +
      f(s.index1, week, model="besag", graph=g,
        hyper=list(prec=list(prior="loggamma", param=c(1, 1))))
    res <- inla(formula, data=dt1[dt1$week != 21,], family="poisson", E=E1)
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.fixed$mean[2]
    b2_est <- res$summary.fixed$mean[3]
    b3_est <- res$summary.fixed$mean[4]
    b4_est <- res$summary.random$s.index1$mean
    eta_pred <- exp(b0_est + b1_est * x1_pred + b2_est * x2_pred + (b3_est + b4_est) * times)
    model_nv[((t-1)*n+1):(t*n), "Times"] <- t
    model_nv[((t-1)*n+1):(t*n), "eta_gen"] <- tail(eta,n)
    model_nv[((t-1)*n+1):(t*n), "eta_data"] <- exp(dt1[dt1$week == 21,]$eta)
    model_nv[((t-1)*n+1):(t*n), "eta_pred"] <- eta_pred
    ############################################################################    
    
    ############################################################################
    # Nonvarying-coefficients + no random effects
    formula <- y1 ~ x1 + x2 + week
    res <- inla(formula, data=dt1[dt1$week != 21,], family="poisson", E=E1)
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.fixed$mean[2]
    b2_est <- res$summary.fixed$mean[3]
    b3_est <- res$summary.fixed$mean[4]
    eta_pred <- exp(b0_est + b1_est * x1_pred + b2_est * x2_pred + b3_est * times)
    model_ns[((t-1)*n+1):(t*n), "Times"] <- t
    model_ns[((t-1)*n+1):(t*n), "eta_gen"] <- tail(eta,n)
    model_ns[((t-1)*n+1):(t*n), "eta_data"] <- exp(dt1[dt1$week == 21,]$eta)
    model_ns[((t-1)*n+1):(t*n), "eta_pred"] <- eta_pred
    ############################################################################     
    
    ############################################################################
    # Bakar's model
    formula <- eta ~ week +
      f(s.index1, x1, model="generic0", Cmatrix=covar, constr=TRUE,
        hyper=list(prec=list(prior="loggamma", param=c(0.5, 1)))) +
      f(s.index2, x2, model="generic0", Cmatrix=covar, constr=TRUE,
        hyper=list(prec=list(prior="loggamma", param=c(1, 1)))) +
      f(s.index3, model="generic0", Cmatrix=covar, constr=TRUE,
        hyper=list(prec=list(prior="loggamma", param=c(2, 1))))
    res <- inla(formula, data=dt1[dt1$week != 21,], family="gaussian")
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.random$s.index1$mean
    b2_est <- res$summary.random$s.index2$mean
    b3_est <- res$summary.fixed$mean[2]
    b4_est <- res$summary.random$s.index3$mean
    eta_pred <- exp(b0_est + b1_est * x1_pred + b2_est * x2_pred + b3_est * times + b4_est)
    model_bakar[((t-1)*n+1):(t*n), "Times"] <- t
    model_bakar[((t-1)*n+1):(t*n), "eta_gen"] <- tail(eta,n)
    model_bakar[((t-1)*n+1):(t*n), "eta_data"] <- exp(dt1[dt1$week == 21,]$eta)
    model_bakar[((t-1)*n+1):(t*n), "eta_pred"] <- eta_pred
    ############################################################################ 
  }
  # Result 
  result <- list(model_Varying = model_v, model_Nonvarying = model_nv,
                 model_NoRandomEffect = model_ns, model_Bakar = model_bakar)
  
  MSE_1 <- list(MSE_Varying = mean((model_v[,"eta_data"] - model_v[,"eta_pred"])^2),
                MSE_Nonvarying = mean((model_nv[,"eta_data"] - model_nv[,"eta_pred"])^2),
                MSE_NoRandomEffect = mean((model_ns[,"eta_data"] - model_ns[,"eta_pred"])^2),
                MSE_Bakar = mean((model_bakar[,"eta_data"] - model_bakar[,"eta_pred"])^2))
  
  MSE_2 <- list(MSE_Varying = mean((model_v[,"eta_gen"] - model_v[,"eta_pred"])^2),
                MSE_Nonvarying = mean((model_nv[,"eta_gen"] - model_nv[,"eta_pred"])^2),
                MSE_NoRandomEffect = mean((model_ns[,"eta_gen"] - model_ns[,"eta_pred"])^2),
                MSE_Bakar = mean((model_bakar[,"eta_gen"] - model_bakar[,"eta_pred"])^2))
  
  # Output
  return(list(Result = result, MSE_data = MSE_1, MSE_gen = MSE_2))
}

################################################################################
# One varying

# rho = 0, E_mean1 = 5000
finalresult_111 <- one_varying(M = 100, rho = 0, E_mean1 = 5000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# rho = 0.3, E_mean1 = 5000
finalresult_112 <- one_varying(M = 100, rho = 0.3, E_mean1 = 5000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# rho = 0.6, E_mean1 = 5000
finalresult_113 <- one_varying(M = 100, rho = 0.6, E_mean1 = 5000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# rho = 0, E_mean1 = 50000
finalresult_121 <- one_varying(M = 100, rho = 0, E_mean1 = 50000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# rho = 0.3, E_mean1 = 50000
finalresult_122 <- one_varying(M = 100, rho = 0.3, E_mean1 = 50000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# rho = 0.6, E_mean1 = 50000
finalresult_123 <- one_varying(M = 100, rho = 0.6, E_mean1 = 50000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

################################################################################
# Two varying

# rho = 0, E_mean1 = 5000
start <- Sys.time()
finalresult_211 <- two_varying(M = 100, rho = 0, E_mean1 = 5000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])
end <- Sys.time()
runtime_211 = end - start

# rho = 0.3, E_mean1 = 5000
finalresult_212 <- two_varying(M = 100, rho = 0.3, E_mean1 = 5000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# rho = 0.6, E_mean1 = 5000
finalresult_213 <- two_varying(M = 100, rho = 0.6, E_mean1 = 5000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# rho = 0, E_mean1 = 50000
finalresult_221 <- two_varying(M = 100, rho = 0, E_mean1 = 50000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# rho = 0.3, E_mean1 = 50000
finalresult_222 <- two_varying(M = 100, rho = 0.3, E_mean1 = 50000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# rho = 0.6, E_mean1 = 50000
finalresult_223 <- two_varying(M = 100, rho = 0.6, E_mean1 = 50000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])


#save(finalresult_111,finalresult_112,finalresult_113,finalresult_121,finalresult_122,finalresult_123,
#     finalresult_211,finalresult_212,finalresult_213,finalresult_221,finalresult_222,finalresult_223,file="./prediction_result.RData")
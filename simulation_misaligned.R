rm(list = ls())

# load the region data
load("./RData/map.RData")

# packages
library(MASS)
library(INLA)

################################################################################
# Table B.1 in Supplementary



################################################################################
### fixed β_1 and varying β_2
################################################################################
one_varying <- function(M, E_mean1, E_sd1, adjmatrix, coordinate) {
  # Week numbers
  times <- 20
  
  # Graph
  adjmatrix <- adjmatrix[ids,ids]
  n <- nrow(adjmatrix)
  W <- as.matrix(adjmatrix)
  g <- inla.read.graph(W)
  
  # Result
  model_v <- model_nv <- model_ns <- model_bakar <- matrix(0, nrow=M, ncol=2+n)
  
  # Generate data
  data_x1 <- matrix(0, nrow=M, ncol=n*times)

  # p_it
  eta_v <- eta_nv <- eta_ns <- eta_bakar <- matrix(0, nrow=M, ncol=n*times)
  eta_gen <- matrix(0, nrow=M, ncol=n*times)
  eta_data <- matrix(0, nrow=M, ncol=n*times)
  
  # Run time
  runtime_v <- runtime_nv <- runtime_ns <- runtime_bakar <- c()
  
  # Spatial structure of Bakar's model
  coordinate <- coordinate[ids,]
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
    x1 <- rnorm(n * times, 0, 1)
    x2 <- t1[t1$week<21,"PM2.5"]

    data_x1[t,] <- x1
    
    ### Generate the true disease rates
    lin <- B0 + B1 * x1 + B2 * x2 + (phi + alpha) * week
    eta <- exp(lin)
    
    eta_gen[t,] <- eta
    
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
    
    eta_data[t,] <- exp(dt1$eta)
    
    ############################################################################
    # Varying-coefficients + random effects
    formula <- y1 ~ x1 + week +
      f(s.index1, x2, model="besag", graph=g,
        hyper=list(prec=list(prior="loggamma", param=c(0.5, 2)))) + 
      f(s.index2, week, model="besag", graph=g,
        hyper=list(prec=list(prior="loggamma", param=c(1, 1))))
    start <- Sys.time()
    res <- inla(formula, data=dt1, family="poisson", E=E1)
    end <- Sys.time()
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.fixed$mean[2]
    b2_est <- res$summary.random$s.index1$mean
    b3_est <- res$summary.fixed$mean[3]
    b4_est <- res$summary.random$s.index2$mean
    eta_v[t,] <- exp(b0_est + b1_est * x1 + b2_est * x2 + (b3_est + b4_est) * week)
    model_v[t,] <- c(b0_est,b1_est,b2_est)
    runtime_v <- c(runtime_v,end-start)
    ############################################################################
    
    ############################################################################
    # Nonvarying-coefficients + random effects
    formula <- y1 ~ x1 + x2 + week +
      f(s.index1, week, model="besag", graph=g,
        hyper=list(prec=list(prior="loggamma", param=c(1, 1))))
    start <- Sys.time()
    res <- inla(formula, data=dt1, family="poisson", E=E1)
    end <- Sys.time()
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.fixed$mean[2]
    b2_est <- res$summary.fixed$mean[3]
    b3_est <- res$summary.fixed$mean[4]
    b4_est <- res$summary.random$s.index1$mean
    eta_nv[t,] <- exp(b0_est + b1_est * x1 + b2_est * x2 + (b3_est + b4_est) * week)
    model_nv[t,] <- c(b0_est,b1_est,rep(b2_est,n))
    runtime_nv <- c(runtime_nv,end-start)
    ############################################################################    
    
    ############################################################################
    # Nonvarying-coefficients + no random effects
    formula <- y1 ~ x1 + x2 + week
    start <- Sys.time()
    res <- inla(formula, data=dt1, family="poisson", E=E1)
    end <- Sys.time()
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.fixed$mean[2]
    b2_est <- res$summary.fixed$mean[3]
    b3_est <- res$summary.fixed$mean[4]
    eta_ns[t,] <- exp(b0_est + b1_est * x1 + b2_est * x2 + b3_est * week )
    model_ns[t,] <- c(b0_est,b1_est,rep(b2_est,n))
    runtime_ns <- c(runtime_ns,end-start)
    ############################################################################     
    
    ############################################################################
    # Bakar's model
    formula <- eta ~ x1 + week +
      f(s.index1, x2, model="generic0", Cmatrix=covar, constr=TRUE,
        hyper=list(prec=list(prior="loggamma", param=c(0.5, 2)))) +
      f(s.index2, model="generic0", Cmatrix=covar, constr=TRUE,
        hyper=list(prec=list(prior="loggamma", param=c(1, 1))))
    start <- Sys.time()
    res <- inla(formula, data=dt1, family="gaussian")
    end <- Sys.time()
    ### Get the result
    b0_est <- res$summary.fixed$mean[1]
    b1_est <- res$summary.fixed$mean[2]
    b2_est <- res$summary.random$s.index1$mean
    b3_est <- res$summary.fixed$mean[3]
    b4_est <- res$summary.random$s.index2$mean
    eta_bakar[t,] <- exp(b0_est + b1_est * x1 + b2_est * x2 + b3_est * week + b4_est)
    model_bakar[t,] <- c(b0_est,b1_est,b2_est)
    runtime_bakar <- c(runtime_bakar,end-start)
    ############################################################################ 
  }
  # Result 
  est <- list(model_Varying = model_v, model_Nonvarying = model_nv,
              model_NoRandomEffect = model_ns, model_Bakar = model_bakar)
  
  MSE_v <- c( mean((model_v[,1]-B0)^2), mean((model_v[,2]-B1)^2),
              mean((sweep(model_v[,3:(2+n)], 2, B2, '-'))^2),
              mean((eta_v - eta_gen)^2))
  MSE_nv <- c( mean((model_nv[,1]-B0)^2), mean((model_nv[,2]-B1)^2),
               mean((sweep(model_nv[,3:(2+n)], 2, B2, '-'))^2),
               mean((eta_nv - eta_gen)^2))
  MSE_ns <- c( mean((model_ns[,1]-B0)^2), mean((model_ns[,2]-B1)^2),
               mean((sweep(model_ns[,3:(2+n)], 2, B2, '-'))^2),
               mean((eta_ns - eta_gen)^2))
  MSE_bakar <- c( mean((model_bakar[,1]-B0)^2), mean((model_bakar[,2]-B1)^2),
                  mean((sweep(model_bakar[,3:(2+n)], 2, B2, '-'))^2),
                  mean((eta_bakar - eta_gen)^2))
  MSE <- list(MSE_Varying = MSE_v, MSE_Nonvarying = MSE_nv,
              MSE_NoRandomEffect = MSE_ns, MSE_Bakar = MSE_bakar)
  
  Var_v <- c( apply(model_v, 2, var)[1:2], mean(apply(model_v, 2, var)[3:(2+n)]),
              mean(apply(eta_v, 2, var)))
  Var_nv <- c( apply(model_nv, 2, var)[1:2], mean(apply(model_nv, 2, var)[3:(2+n)]),
               mean(apply(eta_nv, 2, var)))
  Var_ns <- c( apply(model_ns, 2, var)[1:2], mean(apply(model_ns, 2, var)[3:(2+n)]),
               mean(apply(eta_ns, 2, var)))
  Var_bakar <- c( apply(model_bakar, 2, var)[1:2], mean(apply(model_bakar, 2, var)[3:(2+n)]),
                  mean(apply(eta_bakar, 2, var)))
  Var <- list(Var_Varying = Var_v, Var_Nonvarying = Var_nv,
              Var_NoRandomEffect = Var_ns, Var_Bakar = Var_bakar)
  
  
  Runtime <- list(runtime_Varying = runtime_v,runtime_Nonvarying = runtime_nv,
                  runtime_NoRandomEffect = runtime_ns,runtime_Bakar = runtime_bakar)
  
  eta <- list(eta_Varying = eta_v,eta_Nonvarying = eta_nv,
              eta_NoRandomEffect = eta_ns,eta_Bakar = eta_bakar,
              eta_Generate = eta_gen,eta_Data = eta_data)
  
  data <- list(Beta_0 = B0, Beta_1 = B1, Beta_2 = B2, Phi = phi, Alpha = alpha, X1 = data_x1)

  # Output
  return(list(estimate = est, MSE = MSE, Variance = Var, Runtime = Runtime, eta = eta, Data = data))
}



################################################################################
# One varying

# E_mean1 = 5000
finalresult_111 <- one_varying(M = 100, E_mean1 = 5000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# E_mean1 = 50000
finalresult_112 <- one_varying(M = 100, E_mean1 = 50000, E_sd1 = 1000,
                               adjmatrix = W, coordinate = region_data[,2:3])

# E_mean1 = 500000
#finalresult_113 <- one_varying(M = 100, E_mean1 = 500000, E_sd1 = 1000,
#                               adjmatrix = W, coordinate = region_data[,2:3])


#save(finalresult_111,finalresult_112,file="./simulation_misaligned.RData")
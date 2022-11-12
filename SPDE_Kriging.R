rm(list = ls())

# load the data
load("./RData/data.RData")
load("./RData/SPDE.RData")

# packages
### INLA
library(MASS)
library(INLA)
library(sf)
### Kriging
library(raster)
library(sp)
library(gstat)

################################################################################
# Appendix B.1 in Supplementary






# function
interpolation <- function(times,coop,d){
  
  set.seed(12345)
  
  # 69 point-level observations for PM2.5 concentrations (52 time points) 
  polution_data <- pollution[,c(2:4,6:(5+times))]
  
  # Result
  result_SPDE <- matrix(NA, nrow=69, ncol=times)
  result_Kriging <- matrix(NA, nrow=69, ncol=times)
  
  for (t in 1:times) {
    
    print(t)
    
    for (index_leaveout in which(complete.cases(polution_data[,(3+t)]))) {
      
      # drop NA
      polution_data_temp <- polution_data[complete.cases(polution_data[,(3+t)]),c(1:3,(3+t))]
      
      # Leave-one-out
      polution_training <- polution_data_temp[-which(which(complete.cases(polution_data[,(3+t)]))==index_leaveout),]
      polution_test <- polution_data_temp[which(which(complete.cases(polution_data[,(3+t)]))==index_leaveout),]
      
      polution_test <- polution_test[rep(1,5),]
      polution_test[2,2] <- polution_test[2,2] + d
      polution_test[3,2] <- polution_test[3,2] - d
      polution_test[4,3] <- polution_test[4,3] + d
      polution_test[5,3] <- polution_test[5,3] - d
      
      # SPDE
      ##########################################################################
      coo_p <- data.frame(polution_training$LONG, polution_training$LAT)
      
      ### Constrained Refined Delaunay Triangulation (CRDT)
      mesh_p <- inla.mesh.2d(
        loc = coo_p, offset = c(0.35, 0.65),
        cutoff = 0.05, max.edge = c(0.2, 0.5)
      )
      
      ### Building the SPDE model on the mesh
      spde_p <- inla.spde2.matern(mesh = mesh_p, alpha = 2, constr = TRUE)
      
      ### Index set
      indexs_p <- inla.spde.make.index("s", spde_p$n.spde)
      
      ### Projection matrix
      A_p <- inla.spde.make.A(mesh = mesh_p, loc = as.matrix(coo_p))
      Ap_p <- inla.spde.make.A(mesh = mesh_p, loc = as.matrix(coop))
      
      # Estimation 
      ### stack for estimation stk.e
      stk.e <- inla.stack(
        tag = "est",
        data = list(y = polution_training[,4]), 
        A = list(1, A_p),
        effects = list(data.frame(b0 = rep(1, nrow(coo_p))), s = indexs_p)
      )
      ### stack for prediction stk.p
      stk.p <- inla.stack(
        tag = "pred",
        data = list(y = NA),
        A = list(1, Ap_p),
        effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs_p)
      )
      ### stk.full has stk.e and stk.p
      stk.full <- inla.stack(stk.e, stk.p)
      ### Model formula
      formula <- y ~ 0 + b0 + f(s, model = spde_p)
      ### inla call
      res <- inla(formula,
                  data = inla.stack.data(stk.full),
                  control.predictor = list(
                    compute = TRUE,
                    A = inla.stack.A(stk.full)
                  ), control.compute=list(config = TRUE)
      )
      ### Result
      sample <- res$summary.fitted.values$mean[70:(69+nrow(coop))]
      distcoop <- function(a, b) {
        return(dist(rbind(a, polution_test[b, c("LONG", "LAT")]))[1])
      }
      result_temp <- c()
      for (i in 1:nrow(polution_test)) {
        distance <- apply(coop, 1, distcoop, b = i)
        result_temp <- c(result_temp,sample[which(distance == min(distance))])
      }
      result_SPDE[index_leaveout,t] <- mean(result_temp)
      
      # Kriging
      ##########################################################################
      #newdata = polution_test[1,]
      #newdata[1,2] = newdata[1,2] + rnorm(1,0,d)
      #newdata[1,3] = newdata[1,3] + rnorm(1,0,d)
      newdata = polution_test
      v <- variogram(object = polution_training[,4]~1,data = polution_training ,locations =~LONG+LAT)
      v.fit = fit.variogram(v, vgm(model = "Mat"))
      TEM_krg <- krige(formula = polution_training[,4] ~ 1,
                       model = v.fit,
                       data = polution_training,
                       loc = ~LONG+LAT, 
                       newdata = newdata, 
                       nmax = 15, nmin = 5 
                       )
      result_Kriging[index_leaveout,t] <- mean(TEM_krg$var1.pred)
    }
    
  }
  return(list(Kriging = result_Kriging, SPDE = result_SPDE))
}

result <- interpolation(times = 10,coop = coop,d=0.05)


#save(result,file="./Kriging_SPDE_result.RData")
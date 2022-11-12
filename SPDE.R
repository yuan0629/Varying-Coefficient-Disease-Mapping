rm(list = ls())

# load the data
load("./RData/data.RData")
load("./RData/map.RData")

# packages
library(MASS)
library(INLA)
library(sf)

################################################################################
# Preparation for SPDE_Kriging



######################################################
################ PM2.5 concentrations ################
######################################################

coo_p <- data.frame(pollution$LONG, pollution$LAT)

# Constrained Refined Delaunay Triangulation (CRDT)
mesh_p <- inla.mesh.2d(
  loc = coo_p, offset = c(0.35, 0.65),
  cutoff = 0.05, max.edge = c(0.2, 0.5)
)

# Building the SPDE model on the mesh
spde_p <- inla.spde2.matern(mesh = mesh_p, alpha = 2, constr = TRUE)

# Index set
indexs_p <- inla.spde.make.index("s", spde_p$n.spde)

# Projection matrix
A_p <- inla.spde.make.A(mesh = mesh_p, loc = as.matrix(coo_p))

######################################################
#################### Temperature #####################
######################################################

coo_t <- data.frame(temperature$LONG, temperature$LAT)

mesh_t <- inla.mesh.2d(
  loc = coo_t, offset = c(0.35, 0.65),
  cutoff = 0.05, max.edge = c(0.2, 0.5)
)

spde_t <- inla.spde2.matern(mesh = mesh_t, alpha = 2, constr = TRUE)

indexs_t <- inla.spde.make.index("s", spde_t$n.spde)

A_t <- inla.spde.make.A(mesh = mesh_t, loc = as.matrix(coo_t))

######################################################
################# Relative humidity ##################
######################################################

coo_rh <- data.frame(rh$LONG, rh$LAT)

mesh_rh <- inla.mesh.2d(
  loc = coo_rh, offset = c(0.35, 0.65),
  cutoff = 0.05, max.edge = c(0.2, 0.5)
)

spde_rh <- inla.spde2.matern(mesh = mesh_rh, alpha = 2, constr = TRUE)

indexs_rh <- inla.spde.make.index("s", spde_rh$n.spde)

A_rh <- inla.spde.make.A(mesh = mesh_rh, loc = as.matrix(coo_rh))

########
# Grid #
########

# Coordinate range
bb <- st_bbox(map[-c(2,13),])

# Construct a regular grid with 50 × 50 grid points to cover the whole target region
x <- seq(bb[1] - 0.1, bb[3] + 0.1, length.out = 50)
y <- seq(bb[2] - 0.1, bb[4] + 0.1, length.out = 50)
coop <- as.matrix(expand.grid(x, y))

# Only the grid points in the Taiwan administrative region are kept,
# which leads to a total of m = 854 points left
inter <- function(a) {
  return(lengths(st_intersects(st_point(a), map)))
}
ind <- apply(coop, 1, inter)
coop <- coop[which(ind == 1), ]
interloc <- function(a) {
  return(as.numeric(st_intersects(st_point(a), map2)))# Administrative region
}
gridregion <- apply(coop, 1, interloc)
coop <- coop[!is.na(gridregion),]

# Region index
gridregion <- gridregion[!is.na(gridregion)]
gridregionname <- map2$oldLOC[(gridregion)]

# Projection matrix
Ap_p <- inla.spde.make.A(mesh = mesh_p, loc = as.matrix(coop))# PM2.5
Ap_t <- inla.spde.make.A(mesh = mesh_t, loc = as.matrix(coop))# Temperature
Ap_rh <- inla.spde.make.A(mesh = mesh_rh, loc = as.matrix(coop))# Relative humidity

##############
# Estimation #
##############

# To account for the uncertainty in predicting the reginal values, 
# multiple Monte Carlo integrations (100 times) are conducted

result_p <- data.frame(matrix(0, ncol = 100, nrow = 52 * 328))# PM2.5

for (t in 1:52) {
  # stack for estimation stk.e
  stk.e <- inla.stack(
    tag = "est",
    data = list(y = pollution[,(t+5)]), 
    A = list(1, A_p),
    effects = list(data.frame(b0 = rep(1, nrow(coo_p))), s = indexs_p)
  )
  # stack for prediction stk.p
  stk.p <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_p),
    effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs_p)
  )
  # stk.full has stk.e and stk.p
  stk.full <- inla.stack(stk.e, stk.p)
  # Model formula
  formula <- y ~ 0 + b0 + f(s, model = spde_p)
  # inla call
  res <- inla(formula,
              data = inla.stack.data(stk.full),
              control.predictor = list(
                compute = TRUE,
                A = inla.stack.A(stk.full)
              ), control.compute=list(config = TRUE)
  )
  # result
  samplefun <- function(a) {
    return(inla.rmarginal(100, res$marginals.fitted.values[[a+69]]))
  }
  #mean(samplefun(1))
  sample <- apply(matrix(c(1:nrow(coop)), ncol = 1), 1, samplefun) 
  for (j in 1:100){
    print(paste0(t, ":", j))
    new <- location
    new$PM2.5 <- 0
    distcoop <- function(a, b) {
      return(dist(rbind(a, new[b, c("long", "lat")]))[1])
    }
    for (i in 1:nrow(new)) {
      if (! is.na(mean(sample[j, gridregionname == new$oldLOC[i]]))) { 
        new$PM2.5[i] <- mean(sample[j, gridregionname == new$oldLOC[i]]) 
      } else {
        distance <- apply(coop, 1, distcoop, b = i)
        new$PM2.5[i] <- sample[j, which(distance == min(distance))]  
      }
    }
    result_p[(((t-1)*328 + 1):(t*328)),j] <- new$PM2.5
  }
}



result_t <- data.frame(matrix(0, ncol = 100, nrow = 52 * 328))# Temperature

for (t in 1:52) {
  # stack for estimation stk.e
  stk.e <- inla.stack(
    tag = "est",
    data = list(y = temperature[,(t+2)]), 
    A = list(1, A_t),
    effects = list(data.frame(b0 = rep(1, nrow(coo_t))), s = indexs_t)
  )
  # stack for prediction stk.p
  stk.p <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_t),
    effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs_t)
  )
  # stk.full has stk.e and stk.p
  stk.full <- inla.stack(stk.e, stk.p)
  # Model formula
  formula <- y ~ 0 + b0 + f(s, model = spde_t)
  # inla call
  res <- inla(formula,
              data = inla.stack.data(stk.full),
              control.predictor = list(
                compute = TRUE,
                A = inla.stack.A(stk.full)
              ), control.compute=list(config = TRUE)
  )
  # result
  samplefun <- function(a) {
    return(inla.rmarginal(100, res$marginals.fitted.values[[a+193]]))
  }
  #mean(samplefun(1))
  sample <- apply(matrix(c(1:nrow(coop)), ncol = 1), 1, samplefun) 
  for (j in 1:100){
    print(paste0(t, ":", j))
    new <- location
    new$temperature <- 0
    distcoop <- function(a, b) {
      return(dist(rbind(a, new[b, c("long", "lat")]))[1])
    }
    for (i in 1:nrow(new)) {
      if (! is.na(mean(sample[j, gridregionname == new$oldLOC[i]]))) { 
        new$temperature[i] <- mean(sample[j, gridregionname == new$oldLOC[i]]) 
      } else {
        distance <- apply(coop, 1, distcoop, b = i)
        new$temperature[i] <- sample[j, which(distance == min(distance))]  
      }
    }
    result_t[(((t-1)*328 + 1):(t*328)),j] <- new$temperature
  }
}



result_rh <- data.frame(matrix(0, ncol = 100, nrow = 52 * 328))# Relative humidity

for (t in 1:52) {
  # stack for estimation stk.e
  stk.e <- inla.stack(
    tag = "est",
    data = list(y = rh[,(t+2)]), 
    A = list(1, A_rh),
    effects = list(data.frame(b0 = rep(1, nrow(coo_rh))), s = indexs_rh)
  )
  # stack for prediction stk.p
  stk.p <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_rh),
    effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs_rh)
  )
  # stk.full has stk.e and stk.p
  stk.full <- inla.stack(stk.e, stk.p)
  # Model formula
  formula <- y ~ 0 + b0 + f(s, model = spde_rh)
  # inla call
  res <- inla(formula,
              data = inla.stack.data(stk.full),
              control.predictor = list(
                compute = TRUE,
                A = inla.stack.A(stk.full)
              ), control.compute=list(config = TRUE)
  )
  # result
  samplefun <- function(a) {
    return(inla.rmarginal(100, res$marginals.fitted.values[[a+113]]))
  }
  #mean(samplefun(1))
  sample <- apply(matrix(c(1:nrow(coop)), ncol = 1), 1, samplefun) 
  for (j in 1:100){
    print(paste0(t, ":", j))
    new <- location
    new$rh <- 0
    distcoop <- function(a, b) {
      return(dist(rbind(a, new[b, c("long", "lat")]))[1])
    }
    for (i in 1:nrow(new)) {
      if (! is.na(mean(sample[j, gridregionname == new$oldLOC[i]]))) { 
        new$rh[i] <- mean(sample[j, gridregionname == new$oldLOC[i]]) 
      } else {
        distance <- apply(coop, 1, distcoop, b = i)
        new$rh[i] <- sample[j, which(distance == min(distance))]  
      }
    }
    result_rh[(((t-1)*328 + 1):(t*328)),j] <- new$rh
  }
}



#save(A_p,coo_p,indexs_p,spde_p,coop,Ap_p,gridregionname,result_p,
#     A_t,coo_t,indexs_t,spde_t,Ap_t,result_t,
#     A_rh,coo_rh,indexs_rh,spde_rh,Ap_rh,result_rh,
#     file="./RData/SPDE.RData")






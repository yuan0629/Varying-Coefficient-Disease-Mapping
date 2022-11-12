# Varying-Coefficient-Disease-Mapping
All the code and data used for simulation and real data analysis for paper "Spatiotemporal Varying Coefficient Model for Respiratory Disease Mapping in Taiwan".

## Code 
realdata.R: real data analysis using varying coefficient model, fixed coefficient model and bakar's model respectively.(table D.1, table C.3, table C.4) \

SPDE.R: the INLA-SPDE method used to transform point data to regional data. \

simulation.R: \

simulation_misaligned.R: (table B.1) \

prediction: prediction performance of varying coefficient model, fixed coefficient model, Bakar's model and the disease mapping model without spatial random effects.(table B.2) \

MCMC.R: varying coefficient model using MCMC(table C.2) \

SPDE_kriging: compare the performance of SPDE and kriging in handling change of support problem.(appendix B.1) 

## Data
Folder "RData" contains all the data used for analysis. \
*data.RData: the original PM2.5, temperature and humidity data. \
*fulldata.RData: the full data used for real data analysis. Due to data privacy, "respiratory" and "total" are simulated data similar to the original ones, while the others are real original data. \
*map.RData: the adjacency matrix of 328 counties in Taiwan.
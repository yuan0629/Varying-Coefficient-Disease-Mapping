# Varying-Coefficient-Disease-Mapping
All the code and data used for simulation and real data analysis for paper "Spatiotemporal Varying Coefficient Model for Respiratory Disease Mapping in Taiwan".

## Code 
**realdata.R**: real data analysis using varying coefficient model, fixed coefficient model and bakar's model respectively.(table D.1, table C.3, table C.4) 

**SPDE.R**: the INLA-SPDE method used to transform point-level data PM2.5, temperature and relative humidity to areal-level data. 

**simulation.R**:  simulation under the scenarios of 1)fixed $\beta_1$ and varying $\beta_2$ 2)varing $\beta_1$ and $\beta_2$.(table 1, table2)

**simulation_misaligned.R**: simulation with real point-level PM2.5 concentration.(table B.1) 

**prediction**: prediction performance of varying coefficient model, fixed coefficient model, Bakar's model and the disease mapping model without spatial random effects.(table B.2) 

**MCMC.R**: varying coefficient model using MCMC(table C.2) 

**SPDE_kriging.R**: compare the performance of SPDE and kriging in handling change of support problem.(appendix B.1) 

## Data
Folder "RData" contains all the data used for analysis. \
**data.RData**: the original PM2.5, temperature and humidity data. 

**fulldata.RData**: the full data used for real data analysis. "respiratory" is visit counts for respiratory diseases and "total" is the visit counts for all diseases. Due to data privacy, these two are simulated data similar to the original ones, while the others are real original data. 

**map.RData**: the adjacency matrix of 328 counties in Taiwan.

## Reference
Spatiotemporal Varying Coefficient Model for Respiratory Disease Mapping in Taiwan, Feifei Wang, Congyuan Duan, Yang Li, Hui Huang, and Ben-Chang Shia. _Biostatistics(To appear)_
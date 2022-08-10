####Example run using Serotype 4 as an example
library(tidyverse)
source("./R/PneuVE.R")

ser4sim <- VEpred(nsim=100,seed=realization,obs.VEs=c(97, 65, 100),
                  nSubPlacebo=189,nSubCurrentPCV=123,nSubNewPCV=235,
                  serotype="pn4",GMC.Placebo=0.03,upperCI.Placebo=0.04,
                  GMC.CurrentPCV=1.61,upperCI.CurrentPCV=1.84,
                  GMC.NewPCV=1.55,upperCI.NewPCV=1.70)

#Extract summary of estimates protective threshold IgG concentration
Ser4Cp <- ser4sim[[2]]

#Extract summary of predicted effectiveness
Ser4VE <- ser4sim[[1]]

#Extract Simulation dataframe
dfSer4 <- ser4sim[[3]]

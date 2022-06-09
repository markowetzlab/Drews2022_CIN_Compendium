# Monte Carlo simulation of signature activities

The scripts of this directory run a MC simulation to test the variance in signature activities when adding noise to the features extracted from the samples. The code was optimised for running on the HPC infrastructure slurm. It is recommended to use a HPC environment as the 1000 runs can take a while on a normal laptop or computer.

Currently the scripts model 10% noise added and subtracted to each sample's feature distribution. These values can be changed with the parameters `FINALWIDTH` and `RANGECNAS` of the `addNoiseToFeatures` function. 

```
addNoiseToFeatures(lOriECNF, allFeatures = c("changepoint", "segsize", "bpchrarm", "osCN", "bp10MB"),
                   SDPROP = 20, FINALWIDTH = 0.1, RANGECNAS = 0.1, NUMSIMS = NUMSIMS) 
```

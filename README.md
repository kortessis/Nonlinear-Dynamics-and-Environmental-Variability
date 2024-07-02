# Increasing environmental fluctuations can dampen variability of endogenously cycling populations

## Paper and Code Summary
This repository has the code needed to reproduce the figures and dynamics presented in the manuscript "Increasing environmental fluctuations can dampen variability of endogenously cycling populations". The study investigates how environmental variability impacts population fluctuations in populations that cycle. Relationships between environmental variability and population size variability are investigated using the Ricker model. The Ricker model can exhibit a stable equilibrium, exogenous cycles, and endogenous cycles, where the endogenous cycles correspond to 2-point (and higher) cycles. The study finds that environmental variability elevates population size variation for parameter combinations exhibiting stable equilibria and exogenous cycles. However, the study finds environmental variability dampens population size fluctuations for the 2-pt cycle. The study shows that this reduction in population flucutations can be quantified using a variance partition. The variance partition is used to estimate the variance dampening effect in time series data on Lynx-Hare cycles and nonlinear dynamics of flour beetle populations, showing good support for variance dampening.



## Ricker Model Folder Contents
This folder contains code to analyze the Ricker model in accordance with analytical and comnputational methods described in the paper. The code in this folder is for use in Matlab.

Approximations2.m <- Compares approximation of the variance in the exogenous cycle Ricker model with a simulation. Used to generate figure S4 of the supplementary material.

History_Coeff.m <- Directly calculates the coefficient c in the approximation given by eqn. (2) of the main text. This coefficient summarizes the effect of autocorrelated environmental variation. It has inputs a_bar and cycle_length, being respectively, the average value of the environmental driver and the period of cycle of the environmental driver. 

Exogenous_approximation_only.m <- Same as Approximations2.m, except without the stochastic component. Makes figure S3 of the supplement.

EndogenousVsExogenousCycles.m <- Estimates the variance scaling relationship using Markov simulations of the model and illustrated in figure 2. 

viridis.m <- The "viridis" colormap used in figure creation.

## SupplementalCode Folder Contents
Files in the folder SupplementalCode are used to reproduce empirical analyses of the variance partition. The code in this folder is for use in the R programming language.

1. Folder list
  A.  FlourBeetle
    Short description: Contains .R files and data to analyze Desharnais (1979) flour beetle experimental data.
    
  B. LynxHare
    Short description: Contains .R files and data to analyze Lynx dataset.

2. File list
   A. Filename: /FlourBeetle/FB_MLE.R    
      Short description: Run this file to reproduce analysis of flour beetle dataset. Estiamtes model parameters and estimates variance components.

   B. Filename: /FlourBeetle/ModelFuncs.R
      Short description: Support functions for FB_MLE.R

   C. Filename: /FlourBeetle/FlourBeetle_LPA.csv
      Short description: Dataset from Desharnais (1979).

   D. Filename: /LynxHare/LV_MLE.R    
      Short description: Run this file to reproduce analysis of Lynx dataset. Estimates model parameters and estimates variance components.

   E. Filename: /LynxHare/ModelFuncs.R
      Short description: Support functions for LV_MLE.R

   F. Filename: /LynxHare/LynxHare.txt
      Short description: Dataset from http://people.whitman.edu/~hundledr/courses/M250F03/.
        

2. Relationship between files:        

To run the flour beetle analysis run /FlourBeetle/FB_MLE.R
To run the Lynx analysis run /LynxHare/LV_MLE.R



[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12627725.svg)](https://doi.org/10.5281/zenodo.12627725)


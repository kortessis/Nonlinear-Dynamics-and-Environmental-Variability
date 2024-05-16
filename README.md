# Nonlinear-Dynamics-and-Environmental-Variability

This repository has the code needed to reproduce the figures and dynamics presented in the manuscript "Increasing environmental fluctuations can dampen variability of endogenously cycling populations" posted on BioRxiv (https://doi.org/10.1101/2023.05.10.531506). The functions here are broken into two groups: the polynomial model and the ricker model. The polynomial model refers to a set of ordinary and stochastic differential equations (ode and sde) used to model repeated outbreaks. The ricker model is used to model stable equilibrium dynamics, exogenous cycles, and endogenous cycles. In all cases, the focus is how stochastic environmental variation influences population variation. 



Ricker Model

Approximations2.m <- Compares approximation of the variance in the exogenous cycle Ricker model with a simulation. Used to generate figure S4 of the supplementary material.

History_Coeff.m <- Directly calculates the coefficient c in the approximation given by eqn. (2) of the main text. This coefficient summarizes the effect of autocorrelated environmental variation. It has inputs a_bar and cycle_length, being respectively, the average value of the environmental driver and the period of cycle of the environmental driver. 

Exogenous_approximation_only.m <- Same as Approximations2.m, except without the stochastic component. Makes figure S3 of the supplement.

EndogenousVsExogenousCycles.m <- Estimates the variance scaling relationship using Markov simulations of the model and illustrated in figure 2. 



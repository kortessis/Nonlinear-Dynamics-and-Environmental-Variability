# Nonlinear-Dynamics-and-Environmental-Variability



Polynomial Model (figures 4-5)

AltStatePotential.m <- calculates the potential function

Polynomial_Growth_and_Potential.m <- a function to plot the per-capita growth rate and the potential function for the outbreak model. Requires: AltStatePotential.m

Stochastic_Polynomial_Sim_Updated.m <- Simulates the SDE model, calculates variance and plots the output. Requires: PolynomialODE.m

PolynomialODE.m <- a function that puts the ode into an ODE solver, specifically ode15s.


Ricker Model

Approximations2.m <- Compares approximation of the variance in the exogenous cycle Ricker model with a simulation. Used to generate figure S4 of the supplementary material.

History_Coeff.m <- Directly calculates the coefficient c in the approximation given by eqn. (5) of the main text. This coefficient summarizes the effect of autocorrelated environmental variation. It has inputs bar and cycle_length, being respectively, the average value of the environmental driver and the period of cycle of the environmental driver. 

Exogenous_approximation_only.m <- Same as Approximations2.m, except without the stochastic component. Makes figure S3 of the supplement.

EndogenousVsExogenousCycles.m <- Estimates the variance scaling relationship using Markov simulations of the model. It also uses these simulations to calculate the interaction effect. 




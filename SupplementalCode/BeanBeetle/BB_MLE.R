##analyze the bean beetle data from Cushing et al.
##Data are the control cultures from Table A.1
##Originally in THesis of Desharnais (1979)
##Analysis in Dennis et al 1995 Ecological Monographs
source("ModelFuncs.R")

library(mvtnorm) 
library(MASS)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

beetle.dat <- read.csv(file="BeanBeatle_LPA.csv")

beetle.dat$Total <- beetle.dat$L + beetle.dat$P + beetle.dat$A

beetle.long <- beetle.dat %>% pivot_longer(cols=c('L', 'P', 'A', 'Total'), names_to="Stage", values_to="Count") #format for plots

p1 <- ggplot(beetle.long) + 
  geom_line(data=beetle.long[beetle.long$Stage != 'Total',], aes(x=Week, y=Count, col=Stage, alpha=0.5)) +
  geom_line(data=beetle.long[beetle.long$Stage == 'Total',], aes(x=Week, y=Count), linewidth=1) +
  facet_wrap(~Replicate) + 
  scale_color_brewer(palette='Set2') +
  theme(legend.position="bottom")
  
  

#format data for analysis
state_obs <- as.matrix(beetle.dat[, 2:4])
state_lag  <- as.matrix(beetle.dat[, 5:7])

init_par <- c(b=log(11.68), ua=log(0.1108), ul=log(0.5129), cea=log(0.01097), cel=log(0.009264), cpa=log(0.01779), sig1=log(0.3), sig2=log(0.4), sig3=log(0.01), sig12=0.03, sig13=0.01, sig23=-0.01)

tmax <- dim(state_obs)[1]

#run optimization
LPA_opt <- optim(par=init_par, fn=LPA_nll, tmax=tmax, state_obs=state_obs, state_lag=state_lag, method="BFGS", control=list(maxit=1e4), hessian=TRUE)
parEst <- LPA_opt$par
parEst.cov <- solve(LPA_opt$hessian)

##get variance components
int_draws <- int_drawsA <- int_drawsB <- int_drawsC <- int_drawsD <- matrix(NA, 1e3, 5)
par.mat <- matrix(NA, dim(int_draws)[1], 12)

#deterministic simulation
for(i in 1:dim(int_draws)[1]) {

  par.draw      <- mvrnorm(n=1, mu=parEst, Sigma=parEst.cov)
  par.mat[i,]   <- par.draw
  Allstages.sim <- apply(LPA.sim(pars=par.draw, N=1e4, state_init=state_obs[1,]), 1, sum)
  
  fixed_point <- mean(Allstages.sim[-(1:100)]) #get the predicted fixed point Not sure if theis is correct.
  
  obs_stateA <- apply(state_obs[1:19,], 1, sum)
  obs_stateB <- apply(state_obs[20:38,], 1, sum)
  obs_stateC <- apply(state_obs[39:57,], 1, sum)
  obs_stateD <- apply(state_obs[58:76,], 1, sum)
  
  #set.seed(i)
  int_drawsA[i,]   <- int_estimate(x=Allstages.sim[2:19], X=obs_stateA[-1], mean_x=fixed_point, mean_X=mean(obs_stateA[-1]))
  
  #set.seed(i)
  int_drawsB[i,]   <- int_estimate(x=Allstages.sim[2:19], X=obs_stateB[-1], mean_x=fixed_point, mean_X=mean(obs_stateB[-1]))
  
  #set.seed(i)
  int_drawsC[i,]   <- int_estimate(x=Allstages.sim[2:19], X=obs_stateC[-1], mean_x=fixed_point, mean_X=mean(obs_stateC[-1]))
  
  #set.seed(i)
  int_drawsD[i,]   <- int_estimate(x=Allstages.sim[2:19], X=obs_stateD[-1], mean_x=fixed_point, mean_X=mean(obs_stateD[-1]))
  
  int_draws[i,] <- (int_drawsA[i,] + int_drawsB[i,] + int_drawsC[i,] + int_drawsD[i,])/4 #does this make sense?
}


dimnames(int_drawsA)[[2]] <- dimnames(int_drawsB)[[2]] <- dimnames(int_drawsC)[[2]] <- dimnames(int_drawsD)[[2]] <- dimnames(int_draws)[[2]] <- c("Interaction", "tot", "Deterministic", "Stochastic", "Nonlinear")
  
#weird distribution shape. what is drivign bimodality and the spike near 0?
hist(int_drawsA[,1], breaks=100)
hist(int_drawsB[,1], breaks=100)
hist(int_drawsC[,1], breaks=100)
hist(int_drawsD[,1], breaks=100)
hist(int_draws[,1], breaks=100)

print(quantile(int_drawsA[,1], probs=c(.025, 0.5, 0.975)))
print(quantile(int_drawsB[,1], probs=c(.025, 0.5, 0.975)))
print(quantile(int_drawsC[,1], probs=c(.025, 0.5, 0.975)))
print(quantile(int_drawsD[,1], probs=c(.025, 0.5, 0.975)))
print(quantile(int_draws[,1], probs=c(.025, 0.5, 0.975)))


int.draws.rescale <- int_draws[,-2]
int.draws.rescale[,1] <- int.draws.rescale[,1]/int_draws[,2]
int.draws.rescale[,2] <- int.draws.rescale[,2]/int_draws[,2]
int.draws.rescale[,3] <- int.draws.rescale[,3]/int_draws[,2]
int.draws.rescale[,4] <- int.draws.rescale[,4]/int_draws[,2]
#int.draws.rescale[,5] <- int.draws.rescale[,1]/int.draws[,-2]

apply(int.draws.rescale, 2, mean)
apply(int.draws.rescale, 2, quantile, c(0.025, 0.975))


var.df <- pivot_longer(data.frame(int.draws.rescale), everything(), values_to="Value", names_to="Component")

var.df$Component <- factor(var.df$Component , levels = c("Deterministic", "Stochastic", "Nonlinear", "Interaction"))
p2 <- ggplot(var.df)   + 
  geom_hline(yintercept=0, col="gray") +
  geom_violin(aes(x=Component, y=Value, fill=Component)) +
  lims(y=c(-5, 5)) + 
  scale_fill_brewer(palette="Set1") +
  labs(x="Variance component (scaled)") +
  theme(legend.position="none")

p3 <- plot_grid(p1, p2)

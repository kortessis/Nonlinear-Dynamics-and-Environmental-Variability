# TMB Lotka Volterra
source("ModelFuncs.R") #functions for variance decomposition

library(MASS)
library(TMB)
library(tidyverse)
library(deSolve)
library(cowplot)
theme_set(theme_cowplot())


LVlnorm <- LVlnorm <-"
  #include <TMB.hpp>
  #include <cppad/runge_45.hpp>      // for CppAD::Runge45
  
  template<class Type>
  class lv {
  private:
    Type alpha_;
    Type beta_;
    Type gamma_;
    Type delta_;
  public:
    //set f=x'(t)
    void Ode(const Type &t, const CppAD::vector<Type> &x, CppAD::vector<Type> &f)
    {
        f[0] = x[0]*alpha_ - x[0]*x[1]*beta_;
        f[1] = -gamma_*x[1] + delta_*x[0]*x[1];
        return;
    }
    void setpars(const Type& alpha, const Type& beta, const Type& gamma, const Type& delta)
    {
        alpha_ = alpha;
        beta_ = beta;
        gamma_ = gamma;
        delta_ = delta;
    }
    
  };

  /// Likelihood function
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    /// data and parameters
    DATA_VECTOR(X); // PREY
    DATA_VECTOR(Y); // PREDATORS
    DATA_VECTOR(times);
    DATA_INTEGER(m);
    
    PARAMETER_VECTOR(log_x0); 
    PARAMETER(log_alpha); Type alpha = exp(log_alpha);
    PARAMETER(log_beta); Type beta = exp(log_beta);
    PARAMETER(log_gamma); Type gamma = exp(log_gamma);
    PARAMETER(log_delta); Type delta = exp(log_delta);
    PARAMETER(log_sigmaPrey); Type sigmasqPrey = exp(log_sigmaPrey);
    PARAMETER(log_sigmaPred); Type sigmasqPred = exp(log_sigmaPred);

    int n = times.size();

    vector<Type> x0 = exp(log_x0);

    matrix<Type> x(n,m);
    
    x.row(0) = x0;
    lv<Type> lv_onestep;
    lv_onestep.setpars( alpha, beta, gamma, delta);
    CppAD::vector<Type> xi(m);
    Type ti; // current time step
    Type tf; // next time step
    
    for( int i=0; i<n-1; i++){
      xi = vector<Type>(x.row(i));
      ti = times(i);
      tf = times(i+1);
      x.row(i+1) = vector<Type>(CppAD::Runge45( lv_onestep, 1, ti, tf, xi ));
    }
    
    vector<Type> prey = x.col(0);
    vector<Type> pred = x.col(1);
  
    
    Type nll_prey = 0;
    Type nll_pred = 0;

    nll_prey = -dnorm(log(X(0)), log(prey(0)), sigmasqPrey, true);
    nll_pred = -dnorm(log(Y(0)), log(pred(0)), sigmasqPred, true);
    
    for(int i = 1; i<n; i++){
      
      nll_prey -= dnorm( log(X(i)), log(prey(i)), sigmasqPrey, true);
      nll_pred -= dnorm( log(Y(i)), log(pred(i)), sigmasqPred, true);
    }
    Type nll = nll_pred + nll_prey;
    
    
    REPORT(x);
    
    return nll;

  }
"


write(LVlnorm, file = "LVlnorm.cpp")
compile("LVlnorm.cpp") ### Compiles
dyn.load(dynlib("LVlnorm"))



#### Test run
lh.dat <- read.csv(file="LynxHare.txt", comment.char="#", sep=" ", header=F)
colnames(lh.dat) <- c("Year","Hare","Lynx")

#### Run model
dat      <- list(X=lh.dat$Hare, Y=lh.dat$Lynx, times=1:length(lh.dat$Year), m=2)

pars     <- list( log_alpha =log(0.55), log_beta = log(0.028), log_gamma = log(0.84), log_delta = log(0.026), log_sigmaPred = log(1.1), log_sigmaPrey = log(1.1), log_x0 = log(c(40,28)))

obj       <- MakeADFun(dat, pars, DLL="LVlnorm", hessian=TRUE)
opt       <- nlminb(obj$par, obj$fn, obj$gr)
exp(opt$par)
names(opt$par) <- c("x0", "x0", "alpha", "beta", "gamma", "delta", "sigmaPrey", "sigmaPred")

(2*opt$objective + 2*length(opt$par)) #AIC

det.pred <- LV.sim.tractory(pars=exp(opt$par[3:6]), y_init=matrix(exp(opt$par[1:2]), ncol=2), N=91,  Time=1:91)

summary( sdreport(obj), p.value = T ) 

lynx.pred <- obj$simulate()$x[,2]

pop.df <- data.frame(Year=lh.dat[,1], N=lh.dat[,3], predict=lynx.pred)

colors <- c("Data" = "black", "Model"="cornflowerblue")


pop.plot <- ggplot(data=pop.df) +
  geom_line(aes(x=Year, y=N, color="Data")) +
  geom_point(aes(x=Year, y=N, color="Data")) +
  geom_line(aes(x=Year, y=predict, color="Model")) +
  scale_color_manual(values = colors) +
  labs(y="Abundance", color="")

pop.plot <- ggplot(data=pop.df) +
  geom_line(aes(x=Year, y=N), linewidth=1) +
  labs(y="Pelts", color="")

plot(pop.plot)


X           <- dat$Y
N           <- length(X)
mean_X      <- mean(X)
tot_var     <- sum((dat$X - mean_X)^2)/(N)
x           <- obj$simulate()$x

parEst.mean <- opt$par
parEst.cov <- sdreport(obj)$cov.fixed

int_draws <- int_draws_prey <- matrix(NA, 1e4, 5)
Time <- seq(1, 91, by=1)

for(i in 1:dim(int_draws)[1]) {
  par.draw <- exp(mvrnorm(n=1, mu=parEst.mean, Sigma=parEst.cov))

  lv.sim <- LV.sim.tractory(pars=par.draw[3:6], y_init=matrix(par.draw[1:2], ncol=2), N=91,  Time=Time)

  prey.sim <- lv.sim[,2]
  pred.sim <- lv.sim[,3]

  mean_x <- mean(pred.sim) 
  int_draws[i,]   <- int_estimate(x=pred.sim, X=dat$Y, mean_x=mean(pred.sim), mean_X=mean(dat$Y))

  mean_x <- mean(prey.sim)
  int_draws_prey[i,]   <- int_estimate(x=prey.sim, X=dat$X, mean_x=mean(prey.sim), mean_X=mean(dat$X))
  
}

hist(int_draws[,1])
hist(int_draws_prey[,1])

dimnames(int_draws)[[2]] <- c('Interaction', 'Total', 'Deterministic', 'Stochastic', 'Nonlinear')

quantile(int_draws[,1], probs=c(0.025, 0.975))

int_draws.scaled <- int_draws[,-2]
int_draws.scaled[,1] <- int_draws.scaled[,1]/int_draws[,2]
int_draws.scaled[,2] <- int_draws.scaled[,2]/int_draws[,2]
int_draws.scaled[,3] <- int_draws.scaled[,3]/int_draws[,2]
int_draws.scaled[,4] <- int_draws.scaled[,4]/int_draws[,2]

quantile(int_draws.scaled[,1], probs=c(0.025, 0.975))

var.df <- pivot_longer(data.frame(int_draws.scaled), everything(), values_to="Value", names_to="Component")
var.df$Component <- factor(var.df$Component , levels = c("Deterministic", "Stochastic", "Nonlinear", "Interaction"))
p2 <- ggplot(var.df)   + 
  geom_hline(yintercept=0, col="gray") +
  geom_violin(aes(x=Component, y=Value, fill=Component)) +
  scale_fill_brewer(palette="Set1") +
  labs(x="Variance component (scaled)") +
  theme(legend.position="none")
p2

p3 <- plot_grid(pop.plot, p2, ncol=2)



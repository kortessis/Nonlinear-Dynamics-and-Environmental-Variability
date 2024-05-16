
##functions for analysis
#negative log-likelihood function
LPA_nll <- function(pars, state_obs, state_lag, tmax) {
  
  Sigma <- cov.construct(var.vec = exp(pars[7:9]), cov.vec = pars[10:12]) #diagonal is strictly positive, offdiagonal is not.
  
  temp <- try(solve(Sigma), silent=TRUE)
  
  if(class(temp)[1] == "try-error") {return(Inf)}
  
  nll <- numeric(tmax)
  state_pred <- matrix(NA, tmax,3)
  
  for(i in 1:tmax) {
    state_pred[i,] <- LPA.onestep(pars, state_lag[i,])
    nll[i] <- -sum(dmvnorm(x=log(state_obs[i,]), mean=log(state_pred[i,]), sigma=Sigma, log=TRUE))
  }
  
  return(sum(nll[-c(1,20,39,58)]))
}

#construct a covariance matrix
cov.construct <- function(var.vec, cov.vec) {
  Sigma <- matrix(NA, 3, 3)
  diag(Sigma) <- var.vec
  Sigma[1,2] <- Sigma[2,1] <- cov.vec[1]
  Sigma[1,3] <- Sigma[3,1] <- cov.vec[2]
  Sigma[2,3] <- Sigma[3,2] <- cov.vec[3]
  
  return(Sigma)
}
#from initial conditions predict N timesteps into the future
LPA.sim <- function(pars, N, state_init) {
  
  state_pred <- matrix(NA, N, 3)
  state_pred[1,] <- state_init
  
  for(i in 2:N) {
    state_pred[i,] <- LPA.onestep(pars, state_pred[i-1,])
  }
  
  return(state_pred)
}

#variance decomposition
int_estimate <- function(x, X, mean_x, mean_X) {
  
  N         <- length(X)
  mean_x      <- mean_x #mean(x)
  tot_var     <- sum((X - mean_X)^2)/(N)
  det_var     <- sum((x - mean_x)^2)/(N)
  stoch_var   <- sum((X - x)^2)/(N)
  del         <- mean_x - mean_X
  
  non_var     <- del^2 + 2*del*sum(x - mean_x)/(N) + 2*del*sum(X-x)/(N)
  int_var     <- 2*sum((X - x)*(x - mean_x))/N
  
  return(c(int=int_var, tot=tot_var, det=det_var, stoch=stoch_var, nonlin=non_var))
  
}

#predict one timestep given a past timestep
LPA.onestep <- function(pars, state_lag) {
  
  b   <- exp(pars[1])
  ua  <- exp(pars[2])  
  ul  <- exp(pars[3])
  cea <- exp(pars[4])
  cel <- exp(pars[5])
  cpa <- exp(pars[6])
  
  #from eqn 2.12 in cushing et al. 
  L <- b*state_lag[3]*exp(-cel*state_lag[1] - cea*state_lag[3])
  P <- (1-ul)*state_lag[1]
  A <- state_lag[2]*exp(-cpa*state_lag[3]) + (1-ua)*state_lag[3]
  
  return(c(L_pred=L, P_pred=P, A_pred=A))
}

#multivariate normal densnity fucntion
mvnormDensity <- function(x, mean, Sigma) {
  -3/2*log(2*pi) - 1/2*as.numeric(determinant(Sigma, logarithm=TRUE)$modulus) - 1/2*((x - mean)%*%solve(Sigma)%*%(x - mean))[1,1]
}

##Functions and code to run 

#function to obtain  variance components
#x is the deterministic prediction
#X is the stochastic realization
#mean_x is the deterministic mean
#mean_X is the stochastic mean
int_estimate <- function(x, X, mean_x, mean_X) {
  
  N         <- length(X)-1
  mean_x      <- mean_x 
  tot_var     <- sum((X - mean_X)^2)/(N)
  det_var     <- sum((x - mean_x)^2)/(N)
  stoch_var   <- sum((X - x)^2)/(N)
  del         <- mean_x - mean_X
  
  non_var     <- del^2 + 2*del*sum(x - mean_x)/(N) + 2*del*sum(X-x)/(N)
  int_var     <- 2*sum((X - x)*(x - mean_x))/N
  
  return(c(int=int_var, tot=tot_var, det=det_var, stoch=stoch_var, nonlin=non_var))
  
}

#Lotka-Volterra model
LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*(alpha - beta*y)
    dy = y*(-gamma + delta*x)
    return(list(c(dx, dy)))
  })
}

#simulate a trajectory from the Lotka-Volterra model
LV.sim.tractory <- function(pars, y_init, N, Time) {
  
  State <- c(x=y_init[1,1], y=y_init[1,2])
  
  lynx.prediction <- ode(func = LotVmod, y = State, parms = pars, times = Time)
  
  return(lynx.prediction) 
}
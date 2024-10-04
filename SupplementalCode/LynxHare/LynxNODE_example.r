
## author: Willem Bonnaffe (w.bonnaffe@gmail.com)
## 09-06-2022 - created v0_0
##m modifed by Jake Ferguson 10-03-2024

library(tidyverse)
###################
## Load the data ##
###################

## goal: initiate the NODE 

## load data
TS = read.csv(file="LynxHare.txt", comment.char="#", sep=" ", header=F)
TS = TS[,1:3]
colnames(TS) <- c("Year","Hare","Lynx")


## make out directory
pathToOut = "out/"
system(paste("mkdir",pathToOut))


#
###

##################
## PREPARE DATA ##
##################

## goal: prepare data for training of observation model

## data specs
N       = ncol(TS) - 1
n       = nrow(TS)

## prepare data
X_o = TS[,1]
Y_o = TS[,-1]

## predictive and response variable
t   = X_o 
Y   = Y_o
nt = seq(min(t),max(t),(t[2]-t[1])/1)

## standardise time steps
t_0 = min(t)
t_f = max(t)
dt  = diff(t[1:2])
t_  = 1:length(t)
nt_ = seq(min(t_),max(t_),(t_[2]-t_[1])/1)

## standardise data
Y_     = log(Y)
mean_y = apply(Y_,2,mean)
sd_y   = apply(Y_,2,sd)
Y_     = t((t(Y_)-mean_y)/sd_y)

#
###

## load NODE functions
#source("f_NODE_GM.r")
## argmax.logMarPost ##
## goal: compute parameter vector that maximises the log marginal density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - initial parameters 
# c     - scalar - parameter to control strength of regularisation
argmax.logMarPost = function(X,Y,f,df,Omega,c=1)
{
  error_     = function(x) -logMarPost(X,Y,f,x,c)
  graderror_ = function(x) -ddOmega.logMarPost(X,Y,f,df,x,c)
  Omega      = optim(par    = Omega,
                     fn     = error_,
                     gr     = graderror_,
                     method = "BFGS"# ,
                     # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
  )$par
  return(Omega)
}


## logMarPost ##
## goal: compute the log marginal posterior density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
# c     - scalar - parameter to control strength of regularisation
logMarPost = function(X,Y,f,Omega,c=1)
{
  res       = Y - f(X,Omega)
  logMarLik = - 0.5 * length(Y)     * log(0.5 * sum(res^2)   + 1)
  logMarPri = - 0.5 * length(Omega) * log(0.5 * sum(Omega^2) + 1)
  logMarPos = logMarLik + c*logMarPri
  return(logMarPos)
}



## ddOmega.logMarPost ##
## goal: compute derivate of log marginal posterior density wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
# c     - scalar - parameter to control strength of regularisation
ddOmega.logMarPost = function(X,Y,f,df,Omega,c=1)
{
  res                = Y - f(X,Omega)
  ddOmega.res        =   - df(X,Omega)
  ddOmega.logMarLik  = - 0.5 * length(Y)     * 1/(0.5 * sum(res^2)   + 1) * 0.5 * ddOmega.res%*%res
  ddOmega.logMarPri  = - 0.5 * length(Omega) * 1/(0.5 * sum(Omega^2) + 1) * Omega
  ddOmega.logMarPos  = ddOmega.logMarLik + c*ddOmega.logMarPri ## divide by number of neurones in the network
  return(ddOmega.logMarPos)
}



## parameters observation model
K_o      = 100
W_o      = rep(100,N)
N_o      = W_o * 3
sd1_o    = 0.1
sd2_o    = rep(0.01,N) 
rho      = 1 

#
###

#################################
## FUNCTIONS OBSERVATION MODEL ##
#################################

## goal: functions for the observation model 

## f_o ##
## goal: compute predicted values of response variable at time step t
# t     - float  - time step 
# Omega - vector - parameters 
f_o = function(t, Omega)
{	
  Omega = matrix(Omega,ncol=3)
  return(t(Omega[,1])%*%sin(pi*(t*Omega[,2] + Omega[,3])))
}

## ddt.f_o ##
## goal: compute time derivative of the predicted response t time step t
# t     - float  - time step 
# Omega - vector - parameters 
ddt.f_o = function(t,Omega)
{	
  Omega = matrix(Omega,ncol=3)
  return(t(Omega[,1])%*%(pi*Omega[,2]*cos(pi*(t*Omega[,2] + Omega[,3]))))
}

## ddOmega.f_o ##
## goal: compute derivative of the predicted response wtr to each network parameter
# t     - float  - time step
# Omega - vector - parameters 
ddOmega.f_o = function(t,Omega)
{	
  Omega      = matrix(Omega,ncol=3)
  dfdOmega_1 =                sin(pi*(t*Omega[,2] + Omega[,3]))
  dfdOmega_2 = Omega[,1]*pi*t*cos(pi*(t*Omega[,2] + Omega[,3]))
  dfdOmega_3 = Omega[,1]*pi*1*cos(pi*(t*Omega[,2] + Omega[,3]))
  return(c(dfdOmega_1,dfdOmega_2,dfdOmega_3))
}

## *.eval ##
## goal: compute functions across multiple time steps
# t     - vector - time steps in arbitrary units
# Omega - vector - parameters 
f_o.eval         = function(t,Omega) apply(t(t),2,function(x) f_o(x,Omega))
ddt.f_o.eval     = function(t,Omega) apply(t(t),2,function(x) ddt.f_o(x,Omega))
ddOmega.f_o.eval = function(t,Omega) apply(t(t),2,function(x) ddOmega.f_o(x,Omega))


## logMarLik ##
## goal: compute the log marginal likelihood 
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
logMarLik = function(X,Y,f,Omega)
{
  res      = Y - f(X,Omega)
  logMarLik   = - 0.5 * length(Y) * log(0.5 * sum(res^2)   + 1)
  return(logMarLik)
}



## int_estimate ##
## goal: compute the log marginal likelihood 
# x     - vector - deterministic model trajectory
# X     - vector - stochastic model trajectory
# mean_x     - numeric   - mean of the deterministic trajectory
# mean_X    - numeric   - mean of the stochastic trajectory
int_estimate <- function(x, X, mean_x, mean_X) {
  
  N         <- length(X)
  mean_x      <- mean_x 
  tot_var     <- sum((X - mean_X)^2)/(N)
  det_var     <- sum((x - mean_x)^2)/(N)
  stoch_var   <- sum((X - x)^2)/(N)
  del         <- mean_x - mean_X
  
  non_var     <- del^2 + 2*del*sum(x - mean_x)/(N) + 2*del*sum(X-x)/(N)
  int_var     <- 2*sum((X - x)*(x - mean_x))/N
  
  return(c(int=int_var, tot=tot_var, det=det_var, stoch=stoch_var, nonlin=non_var))
  
}

#
###


#
###

###########################
## FIT OBSERVATION MODEL ##
###########################

## goal: fit observation model 

## interpolate each column in the time series
Omega_o    = list()
Yhat_o     = list()
ddt.Yhat_o = list()
logPost_o  = list()

for (i in 2:N)
{
  
  counter <- 0
    ## iterator
    print(paste("fitting: ",i,"/",ncol(Y_o),sep=""))

    ## get ensemble
    Omega_o_i = NULL
    logPost_o_i = NULL
    for(k in 1:K_o)
    {   
        check = F
        while(check == F)
        {
          counter <- counter + 1
            ## fit
            Omega_0   = rnorm(N_o[i], 0, sd2_o[i]) #draw from the priors
            Omega_f   = argmax.logMarPost(t_, Y_[,i], f_o.eval, ddOmega.f_o.eval, Omega_0, 1/W_o[i])
            
            ## update
            logPost_0 = logMarPost(t_, Y_[,i], f_o.eval, Omega_0, 1/W_o[i])
            logPost_f = logMarPost(t_, Y_[,i], f_o.eval, Omega_f, 1/W_o[i])
            print(paste(sprintf("%03d",k),"/",K_o,"    ",
                    format(round(logPost_0,2),nsmall=2),"    ","-->","    ",
                    format(round(logPost_f,2),nsmall=2),sep=""))
            
            ## reject or keep sample
            check = (logPost_f >= logPost_0 + 40)
            if(check == T)
            {
                Omega_o_i = rbind(Omega_o_i,Omega_f)
                # logPost_o_i = c(logPost_o_i,logPost_f)
                logPost_o_i = c(logPost_o_i,logMarLik(t_, Y_[,i], f_o.eval, Omega_f))
            }
        }
    }

    ## sort by best models
    s = rev(order(logPost_o_i))[1:round(rho*K_o)]
    Omega_o_i = Omega_o_i[s,]

    ## store ensemble
    Omega_o[[i]] = Omega_o_i

    ## compute predictions
    Yhat_o[[i]]     = t(apply(Omega_o[[i]],1,function(x)     f_o.eval(nt_,x))) 
    ddt.Yhat_o[[i]] = t(apply(Omega_o[[i]],1,function(x) ddt.f_o.eval(nt_,x)))

    ## de-standardise
    Yhat_o[[i]]     = exp(mean_y[i] + sd_y[i] * Yhat_o[[i]])
    ddt.Yhat_o[[i]] = 1/dt * sd_y[i] * Yhat_o[[i]] * ddt.Yhat_o[[i]]
 }

## store results
names(Yhat_o)     = colnames(TS[,-1])
names(ddt.Yhat_o) = colnames(TS[,-1])
save(Yhat_o,    file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
save(ddt.Yhat_o,file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
save(Omega_o,   file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))

#
###



############################################
##Run the variance decomposition bootstrap##
############################################



var.est <- matrix(NA, dim(Yhat_o[[2]])[1], 5)
lynx.dat <- TS$Lynx

for(i in 1:dim(Yhat_o[[2]])[1]) {
  
  lynx.pred <- Yhat_o[[2]][i,]
  
  var.est[i,] <- int_estimate(X=lynx.dat, x=lynx.pred, mean_X=mean(lynx.dat), mean_x=mean(lynx.pred))
  
}

dimnames(var.est)[[2]] <- c("Interaction", "Total", "Deterministic", "Stochastic", "Nonlinear")

var.est.rescaled <- var.est[,-2] / var.est[,2]
var.df <- pivot_longer(data.frame(var.est.rescaled), everything(), values_to="Value", names_to="Component")

var.df$Component <- factor(var.df$Component , levels = c("Deterministic", "Stochastic", "Nonlinear", "Interaction"))


pop.plot <- ggplot(data=TS) +
  geom_line(aes(x=Year, y=Lynx)) +
  labs(x="Year", y="Pelts")

var.plot <- ggplot(var.df) + 
  geom_hline(yintercept=0, col="gray") +
  geom_violin(aes(x=Component, y=Value, fill=Component)) +
  lims(y=c(-2, 2)) + 
  scale_fill_brewer(palette="Set1") +
  labs(x="Variance component (scaled)") +
  theme(legend.position="none")


pop.plot
var.plot

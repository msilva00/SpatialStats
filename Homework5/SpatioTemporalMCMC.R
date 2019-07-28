#######################################################
#
# INPUTS
#
# Y        := ns x nt matrix of observations
# s        := nsx2 matrix of spatial coordinates
# X        := nsxpxnt matrix of covariates
#
# sp,Xp    := s and X for prediction locations
#
# mean_nu  := log(nu) ~ N(mean_nu,sd_nu)
# sd_nu
# range_nu := log(rangee) ~ N(mean_range,sd_range)
# sd_range
# a,b      := variance ~ InvGamma(a,b)
# sd_beta  := beta[j] ~ N(0,sd_beta)
#
# iters    := number of MCMC iteration
# burn     := length of burn-in
#
#######################################################


ST_Krige<-function(Y,s,X,
                   sp=NULL,Xp=NULL,
                   mean_nu=-1,sd_nu=1,
                   mean_range=0,sd_range=10,
                   a=.1,b=.01,
                   sd_beta=1000,
                   iters=25000,burn=5000){
  
  library(emulator)
  library(geoR)
  
  
  # Bookkeeping
  
  miss     <- is.na(Y)
  nm       <- sum(miss)
  ns       <- nrow(Y)
  nt       <- ncol(Y)
  p        <- dim(X)[2]
  d        <- rdist(s,s)
  diag(d)  <- 0
  theta.mn <- c(mean_range,mean_nu)
  theta.sd <- c(sd_range,sd_nu)
  
  # Initial values
  
  Y[miss] <- mean(Y,na.rm=TRUE)
  beta    <- rep(0,p)
  beta[1] <- mean(Y)
  mu      <- matrix(0,ns,nt)
  taue    <- 1
  taus    <- 1
  rho     <- 0.5
  theta   <- log(c(0.1,0.5))   
  Q       <- solve(corfx(d,theta))
  Xb      <- matrix(0,ns,nt)
  for(t in 1:nt){
    Xb[,t] <- X[,,t]%*%beta
  }
  
  # Keep track of stuff
  
  keep.beta <- matrix(0,iters,p)
  keepers   <- matrix(0,iters,5)
  Y1<-Y2<-P1<-P2<-0
  
  colnames(keepers)   <- c("sigma_e","sigma_s","rho","range","nu")
  colnames(keep.beta) <- colnames(X)
  
  # Precompute a few things
  
  Ptheta <- t(chol(1.8*diag(2)-.8))
  tXX <- 0
  for(t in 1:nt){
    tXX <- tXX + t(X[,,t])%*%X[,,t]
  }
  
  # GO!!!
  
  for(i in 1:iters){
    
    ##############################################:
    #####         MISSING Y (Gibbs)        #######:
    ##############################################:
    
    if(nm>0){
      Y[miss]<-rnorm(nm,Xb[miss]+mu[miss],1/sqrt(taue))
    }        
    
    ##############################################:
    #####       MEAN PARAMETERS (Gibbs)    #######:
    ##############################################:
    
    VVV  <- solve(taue*tXX + diag(p)/sd_beta^2)
    R    <- Y-mu
    MMM  <- 0
    for(t in 1:nt){
      MMM <- MMM + taue*t(X[,,t])%*%R[,t]
    }
    beta <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p)
    for(t in 1:nt){Xb[,t] <- X[,,t]%*%beta}
    
    ##############################################:
    #####       SPATIAL TERMS (Gibbs)      #######:
    ##############################################:
    
    R    <- Y-Xb
    rho2 <- rho^2
    V1   <- solve(taue*diag(ns) + taus*Q*(1+rho2/(1-rho2)))
    V2   <- solve(taue*diag(ns) + taus*Q*(1+rho2)/(1-rho2))
    V3   <- solve(taue*diag(ns) + taus*Q/(1-rho2))
    P1   <- t(chol(V1))
    P2   <- t(chol(V2))
    P3   <- t(chol(V3))
    
    for(t in 1:nt){
      if(t==1){
        VVV <- V1
        PPP <- P1
      }
      if(t>1 & t<nt){
        VVV <- V2
        PPP <- P2
      }
      if(t==nt){
        VVV <- V3
        PPP <- P3
      }
      
      MMM          <- taue*R[,t]
      if(t> 1){MMM <- MMM + rho*taus*Q%*%mu[,t-1]/(1-rho2)}
      if(t<nt){MMM <- MMM + rho*taus*Q%*%mu[,t+1]/(1-rho2)}
      
      mu[,t] <- VVV%*%MMM + PPP%*%rnorm(ns)
    }
    
    ##############################################:
    #####          VARIANCE (Gibbs)        #######:
    ##############################################:
    
    R      <- mu[,-1]-rho*mu[,-nt]
    R      <- sum(R*(Q%*%R))
    SS     <- quad.form(Q,mu[,1]) + R/(1-rho^2)
    taus   <- rgamma(1,ns*nt/2+a,SS/2+b)
    taue   <- rgamma(1,ns*nt/2+a,sum((Y-Xb-mu)^2)/2+b)
    
    ##############################################:
    #### CORRELATION PARAMETERS (Metropolis) #####:
    ##############################################:
    
    # Temporal
    
    R      <- mu[,-1]-rho*mu[,-nt]
    R      <- sum(R*(Q%*%R))
    curll  <- -0.5*ns*(nt-1)*log(1-rho*rho) - 
      0.5*taus*R/(1-rho*rho) 
    can    <- pnorm(rnorm(1,qnorm(rho),.02))
    R      <- mu[,-1]-can*mu[,-nt]
    R      <- sum(R*(Q%*%R))
    canll  <- -0.5*ns*(nt-1)*log(1-can*can) - 
      0.5*taus*R/(1-can*can) 
    R      <- canll-curll+
      dnorm(qnorm(can),log=TRUE)-
      dnorm(qnorm(rho),log=TRUE)
    if(log(runif(1))<R){
      rho <- can
    }
    
    # Spatial
    
    R <- (mu[,-1]-rho*mu[,-nt])/sqrt(1-rho^2)
    R <- sqrt(taus)*cbind(mu[,1],R)
    
    can   <- theta+Ptheta%*%rnorm(2,0,0.1)
    canQ  <- solve(corfx(d,can))
    curll <- 0.5*nt*determinant(Q)$modulus-0.5*sum(R*(Q%*%R))
    canll <- 0.5*nt*determinant(canQ)$modulus-0.5*sum(R*(canQ%*%R))
    MH    <- canll-curll+
      sum(dnorm(  can,theta.mn,theta.sd,log=TRUE))-
      sum(dnorm(theta,theta.mn,theta.sd,log=TRUE))
    if(log(runif(1))<MH){
      theta <- can
      Q     <- canQ
    }  
    
    ##############################################:
    #####        KEEP TRACK OF STUFF       #######:
    ##############################################:
    
    keep.beta[i,]  <- beta
    keepers[i,1]   <- 1/sqrt(taue)
    keepers[i,2]   <- 1/sqrt(taus)
    keepers[i,3]   <- rho
    keepers[i,4:5] <- exp(theta)
    
    ##############################################:
    #####           PREDICTIONS            #######:
    ##############################################:
    
    if(i>burn){
      Y1  <- Y1 + Y/(iters-burn)
      Y2  <- Y2 + Y*Y/(iters-burn)
    }
    
  }   
  
  output <- list(beta=keep.beta,keepers=keepers,Y.mn=Y1,Y.var=Y2-Y1^2)
  
  return(output)}


corfx <- function(d,theta){
  matern(d,exp(theta[1]),exp(theta[2]))
}


#### Simulate a fake dataset to try it out ####
library(fields)
library(geoR)

ns    <- 50
nt    <- 100
p     <- 2
sige  <- 1
sigs  <- 4
rho   <- 0.8
beta  <- c(10,1)
theta <- c(log(0.1),log(0.5))

S      <- cbind(runif(ns),runif(ns))
X      <- array(rnorm(ns*nt*p),c(ns,p,nt))
X[,1,] <- 1   
Cor    <- corfx(rdist(S),theta)
P      <- t(chol(Cor))

mu     <- sigs*P%*%rnorm(ns)
Xb     <- X[,,1]%*%beta
for(t in 2:nt){
  mu     <- cbind(mu,rho*mu[,t-1]+sigs*sqrt(1-rho*rho)*P%*%rnorm(ns))
  Xb     <- cbind(Xb,X[,,t]%*%beta)
}
Y    <- mu + Xb + rnorm(ns*nt,0,sige)

image.plot(1:nt,1:ns,t(Y),xlab="Time",ylab="Space",main="Simulated data")


#### Split the data into training and testing and fit the model ####
miss  <- matrix(runif(ns*nt),ns,nt)<.2
Yo    <- ifelse(miss,NA,Y)

burn  <- 5000
iters <- 25000
fit   <- ST_Krige(Yo,S,X,burn=burn,iters=iters)


#### Check predictions ####
Yp    <- Y[miss]
Y.mn  <- fit$Y.mn[miss]
Y.sd  <- sqrt(fit$Y.var[miss])

cov   <- mean(abs(Yp-Y.mn)<2*Y.sd)
plot(Yp,Y.mn,
     xlab="Test set value",ylab="Predicted value",
     main = paste0("Coverage = ",round(100*cov,1),"%"))
abline(0,1)  


#### Check convergence and fit ####
burnin = 1000:25000
plot(fit$keepers[burnin,1],type="l",main="Nugget SD")
abline(sige,0,lwd=2,col=4)

plot(fit$keepers[burnin,2],type="l",main="Partial sill (SD)")
abline(sigs,0,lwd=2,col=4)

plot(fit$keepers[burnin,3],type="l",main="AR coefficient")
abline(rho,0,lwd=2,col=4)

plot(fit$keepers[burnin,4],type="l",main="Matern range")
abline(exp(theta[1]),0,lwd=2,col=4)

plot(fit$keepers[burnin,5],type="l",main="Matern smoothness")
abline(exp(theta[2]),0,lwd=2,col=4)

pairs(fit$keepers[burn:iters,])

par(mfrow = c(2,1))
for(j in 1:p){
  plot(fit$beta[burnin,j],type="l",main=paste0("beta",j))
  abline(beta[j],0,lwd=2,col=4)
}

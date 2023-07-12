# ------------------------------------------------------------------------------------
# Title: R code for "An unbiased predictor for skewed response variable with  
# measurement error in covariate"
# Note: This code illustrates the coverage probabilities and expected lengths of 
# prediction intervals in Table 4.
# Author: Sepideh Mosaferi
# Date: March 2023
# ------------------------------------------------------------------------------------

setwd("/Users/sepidehmosaferi/Desktop/Code/Simulations")
rm(list = ls(all = TRUE))

source("Jackknife.R")
source("pars_estimation.R")
source("Bootstrap_CI.R")
require("sae"); require("dplyr"); require("infer"); require("gridExtra"); require("tidyverse") 

iter <- 0

sig2bhatmes <- c()
eblupmes <- c()
thetais <- c()
betahats <- c()
mseest <- c()
newparests <- c()
ebmeyls <- c()
m <- 50  #number of small areas  

# names of estimators: 
thetaiexps <- c()
predexp0s <- c()
predexp1s <- c()
predexpdirs <- c()
predexp2s <- c()
predexp3s <- c()
var_directs <- c()
R1is <- c()
mse_Js <- c()
CI_Boot_lows_90 <- c()
CI_Boot_ups_90 <- c()
CI_Boot_lows_95 <- c()
CI_Boot_ups_95 <- c()
CI_Boot_lows_99 <- c()
CI_Boot_ups_99 <- c()

# values of parameters:
mux <- 5
sdx <-  3
sdv <- sqrt(2)

set.seed(1200) #for the sake of reproducibility
psi <- rgamma(m,shape=4.5,scale=2)
beta0 <- 0
beta1 <- 3
sdu <- sqrt(sample(c(0,2),size=m,prob=c(1-50/100,50/100),replace=TRUE)) #percentage of measurement error
xi <- rnorm(m, mean = mux, sd = sdx)


repeat{
  
  iter <- iter + 1
  vi <- rnorm(m, mean  = 0, sd = sdv)
  ei <- rnorm(m, mean = 0, sd = sqrt(psi))
  thetai <- beta0 + beta1*xi + vi
  yi <- thetai + ei
  ui <- rnorm(m, mean = 0, sd = sdu)
  Xi <- xi + ui
  
  saedat <- data.frame(yi, Xi, psi)
  
  efh <- eblupFH(yi~Xi, vardir = psi, data = saedat, method = "ML")
  
  sig2bhatmes <- c(sig2bhatmes, efh$fit$refvar)
  eblupmes <- rbind(eblupmes, as.vector(efh$eblup))
  thetais <- rbind(thetais, thetai)
  betahats <- rbind(betahats, efh$fit$estcoef$beta)
  
  mseest <- rbind(mseest, mseFH(yi~Xi, vardir = psi, data = saedat, method = "ML")$mse)
  
  Cis <- sdu^2
  
  # using delta method:
  var_direct <- (exp(yi)^2)*psi
  var_directs <- rbind(var_directs,var_direct)
  
  ## Estimating unknown parameters based on iterative algorithm
  newpar.est <- FHME_estimation(yi, Xi, psi, Cis)$para
  
  newparests <- rbind(newparests, newpar.est) #(beta0,beta1,sigma^2v)
  
  gammameyl <- (newpar.est[3] + newpar.est[2]^2*Cis)/(newpar.est[3] + newpar.est[2]^2*Cis+psi)
  
  thetahatmeyl <- gammameyl*yi + (1-gammameyl)*as.vector(cbind(1, Xi)%*%newpar.est[c(1,2)])
  ebmeyls <- rbind(ebmeyls, thetahatmeyl)
  
  gammaeyl2 <- newpar.est[3]/(newpar.est[2]^2*Cis + psi + newpar.est[3])
  nuil <- (yi - newpar.est[1] - newpar.est[2]*Xi)
  xdot <- Xi + (newpar.est[2]*Cis/(newpar.est[2]^2*Cis + psi + newpar.est[3]))*nuil
  
  #### Prediction on the exponential scale:
  
  ####  True parameter
  thetaiexps <- rbind( thetaiexps, exp(thetai))
  
  ####  EBLUP ignoring measurement error
  gammameyl0 <- (newpar.est[3])/(newpar.est[3]+psi)
  thetahatmeyl0 <- gammameyl0*yi + (1-gammameyl0)*as.vector(cbind(1, xi)%*%newpar.est[c(1,2)])
  predexp0s <- rbind(predexp0s, exp(thetahatmeyl0 + 0.5*psi*gammameyl0))
  
  ####  EBLUP Mosaferi et al. (predictor A: biased predictor)
  predexp1s <- rbind(predexp1s, exp(thetahatmeyl + 0.5*psi*gammameyl))
  
  #### Direct estimator  
  predexpdirs <- rbind(predexpdirs, exp(yi))
  
  ####  New predictor estimating non-constant error variance (predictor B)
  
  sigmatilde2 <- Cis-(newpar.est[2]^2*Cis^2/(newpar.est[2]^2*Cis+psi+newpar.est[3]))
  predexp2a <- exp(newpar.est[1] + newpar.est[2]*xdot)
  predexp2b <- exp(-newpar.est[2]^2*sigmatilde2/2)
  predexp2c <- exp( gammaeyl2*nuil + 0.5*(gammaeyl2*(newpar.est[2]^2*Cis + psi)))
  
  predexp2 <- predexp2a*predexp2b*predexp2c
  predexp2s <-rbind(predexp2s, predexp2) 
  
  ####  New predictor assuming known parameters:
  
  gammaeyl22 <- sdv^2/(beta1^2*Cis + psi + sdv^2)
  nuil2 <- (yi - beta0 - beta1*Xi)
  xdot2 <- Xi - -beta1*Cis/(beta1^2*Cis + psi + sdv^2)*nuil2
  predexp3a <- exp(beta0 + beta1*xdot2)
  sigmatilde2_2 <- Cis- (beta1^2*Cis^2/(beta1^2*Cis + psi + sdv^2))
  predexp3b <- exp(-beta1^2*sigmatilde2_2/2)
  predexp3c <- exp( gammaeyl22*nuil2 + 0.5*(gammaeyl22*(beta1^2*Cis + psi) ))
  predexp3 <- predexp3a*predexp3b*predexp3c
  predexp3s <-rbind(predexp3s, predexp3) 
  
  ## MSE development
  
  ## First term of MSE
  di <- 2*psi*newpar.est[2]^2*Cis/(newpar.est[2]^2*Cis+newpar.est[3]+psi)
  M2i <- 1-2*exp(1.5*di-gammameyl*psi)+exp(di-gammameyl*psi)
  
  R1i <- (M2i^2)*exp(4*(yi-psi))*(1-exp(-4*newpar.est[3]-4*psi))
  R1is <- rbind(R1is,R1i)
  
  ## Jackknife estimator of MSE 
  
  R1i_del <- matrix(0,nrow=m,ncol=m)   #each row belongs to one area
  thetaB_del <- matrix(0,nrow=m,ncol=m) 
  
  Jack <- Jackknife(yi,Xi,psi,Cis,R1i_del,thetaB_del)
  
  # reassigning names
  thetaB_del <- Jack[,1:m]  #each row belongs to one area
  R1i_del <- Jack[,(m+1):(2*m)]
  
  R1i_Jack <- rep(0,m)
  for(u in 1:m){
    R1i_Jack[u] <- R1i[u]-((m-1)/m)*sum(R1i_del[u,]-R1i[u]) 
  }
  
  R2i_Jack <- rep(0,m)
  for(u in 1:m){
    R2i_Jack[u] <- ((m-1)/m)*sum(thetaB_del[u,]-predexp2[u])^2
  }
  
  mse_J <- R1i_Jack+R2i_Jack   #Final Jackknife estimator for MSE 
  mse_Js <- rbind(mse_Js,mse_J)
  
  ## Bootstrap confidence interval
  CI_Boot <- Bootstrap(yi,Xi,psi,Cis)
  
  CI_Boot_low_90 <- CI_Boot[1,1:m]
  CI_Boot_up_90 <- CI_Boot[2,1:m]
  
  CI_Boot_low_95 <- CI_Boot[1,(m+1):(2*m)]
  CI_Boot_up_95 <- CI_Boot[2,(m+1):(2*m)]
  
  CI_Boot_low_99 <- CI_Boot[1,(2*m+1):(3*m)]
  CI_Boot_up_99 <- CI_Boot[2,(2*m+1):(3*m)] 
  
  CI_Boot_lows_90 <- rbind(CI_Boot_lows_90,CI_Boot_low_90)
  CI_Boot_ups_90 <- rbind(CI_Boot_ups_90,CI_Boot_up_90)
  
  CI_Boot_lows_95 <- rbind(CI_Boot_lows_95,CI_Boot_low_95)
  CI_Boot_ups_95 <- rbind(CI_Boot_ups_95,CI_Boot_up_95)
  
  CI_Boot_lows_99 <- rbind(CI_Boot_lows_99,CI_Boot_low_99)
  CI_Boot_ups_99 <- rbind(CI_Boot_ups_99,CI_Boot_up_99)
  
  if(iter == 2000){break}
  
}

## prediction interval study
CI_direct_lows_90 <- predexpdirs-1.64*sqrt(var_directs)   
CI_direct_ups_90 <- predexpdirs+1.64*sqrt(var_directs) 
CI_direct_lows_95 <- predexpdirs-1.96*sqrt(var_directs)   
CI_direct_ups_95 <- predexpdirs+1.96*sqrt(var_directs) 
CI_direct_lows_99 <- predexpdirs-2.58*sqrt(var_directs)   
CI_direct_ups_99 <- predexpdirs+2.58*sqrt(var_directs) 

CI_R1i_lows_90 <- predexp2s-1.64*sqrt(R1is)   
CI_R1i_ups_90 <- predexp2s+1.64*sqrt(R1is) 
CI_R1i_lows_95 <- predexp2s-1.96*sqrt(R1is)   
CI_R1i_ups_95 <- predexp2s+1.96*sqrt(R1is) 
CI_R1i_lows_99 <- predexp2s-2.58*sqrt(R1is)   
CI_R1i_ups_99 <- predexp2s+2.58*sqrt(R1is) 

CI_Jack_lows_90 <- predexp2s-1.64*sqrt(mse_Js)   
CI_Jack_ups_90 <- predexp2s+1.64*sqrt(mse_Js) 
CI_Jack_lows_95 <- predexp2s-1.96*sqrt(mse_Js)   
CI_Jack_ups_95 <- predexp2s+1.96*sqrt(mse_Js) 
CI_Jack_lows_99 <- predexp2s-2.58*sqrt(mse_Js)   
CI_Jack_ups_99 <- predexp2s+2.58*sqrt(mse_Js) 

## empirical CI length
logCI_direct_length_90 <- apply(CI_direct_ups_90-CI_direct_lows_90,2,log); mean(logCI_direct_length_90)
logCI_direct_length_95 <- apply(CI_direct_ups_95-CI_direct_lows_95,2,log); mean(logCI_direct_length_95)
logCI_direct_length_99 <- apply(CI_direct_ups_99-CI_direct_lows_99,2,log); mean(logCI_direct_length_99) 

logCI_R1i_length_90 <- apply(CI_R1i_ups_90-CI_R1i_lows_90,2,log); mean(logCI_R1i_length_90) 
logCI_R1i_length_95 <- apply(CI_R1i_ups_95-CI_R1i_lows_95,2,log); mean(logCI_R1i_length_95)
logCI_R1i_length_99 <- apply(CI_R1i_ups_99-CI_R1i_lows_99,2,log); mean(logCI_R1i_length_99)

logCI_Jack_length_90 <- apply(CI_Jack_ups_90-CI_Jack_lows_90,2,log)
logCI_Jack_length_90_nona <- logCI_Jack_length_90[!is.na(logCI_Jack_length_90)]; mean(logCI_Jack_length_90_nona) 
logCI_Jack_length_95 <- apply(CI_Jack_ups_95-CI_Jack_lows_95,2,log)
logCI_Jack_length_95_nona <- logCI_Jack_length_95[!is.na(logCI_Jack_length_95)]; mean(logCI_Jack_length_95_nona) 
logCI_Jack_length_99 <- apply(CI_Jack_ups_99-CI_Jack_lows_99,2,log)
logCI_Jack_length_99_nona <- logCI_Jack_length_99[!is.na(logCI_Jack_length_99)]; mean(logCI_Jack_length_99_nona) 

logCI_Boot_length_90 <- apply(CI_Boot_ups_90-CI_Boot_lows_90,2,log); mean(logCI_Boot_length_90)
logCI_Boot_length_95 <- apply(CI_Boot_ups_95-CI_Boot_lows_95,2,log); mean(logCI_Boot_length_95)
logCI_Boot_length_99 <- apply(CI_Boot_ups_99-CI_Boot_lows_99,2,log); mean(logCI_Boot_length_99)

## empirical coverage

Coverage_90_direct <- rep(NA,m)
for(i in 1:m){
  Coverage_90_direct[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_direct_lows_90[u,i] && thetaiexps[u,i]<=CI_direct_ups_90[u,i])})) 
}
mean(Coverage_90_direct)

Coverage_95_direct <- rep(NA,m)
for(i in 1:m){
  Coverage_95_direct[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_direct_lows_95[u,i] && thetaiexps[u,i]<=CI_direct_ups_95[u,i])})) 
}
mean(Coverage_95_direct)

Coverage_99_direct <- rep(NA,m)
for(i in 1:m){
  Coverage_99_direct[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_direct_lows_99[u,i] && thetaiexps[u,i]<=CI_direct_ups_99[u,i])})) 
}
mean(Coverage_99_direct)

Coverage_90_R1i <- rep(NA,m)
for(i in 1:m){
  Coverage_90_R1i[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_R1i_lows_90[u,i] && thetaiexps[u,i]<=CI_R1i_ups_90[u,i])})) 
}
mean(Coverage_90_R1i)

Coverage_95_R1i <- rep(NA,m)
for(i in 1:m){
  Coverage_95_R1i[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_R1i_lows_95[u,i] && thetaiexps[u,i]<=CI_R1i_ups_95[u,i])})) 
}
mean(Coverage_95_R1i)

Coverage_99_R1i <- rep(NA,m)
for(i in 1:m){
  Coverage_99_R1i[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_R1i_lows_99[u,i] && thetaiexps[u,i]<=CI_R1i_ups_99[u,i])})) 
}
mean(Coverage_99_R1i)

Coverage_90_Jack <- rep(NA,m)
for(i in 1:m){
  Coverage_90_Jack[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_Jack_lows_90[u,i] && thetaiexps[u,i]<=CI_Jack_ups_90[u,i])})) 
}
mean(Coverage_90_Jack,na.rm=TRUE)

Coverage_95_Jack <- rep(NA,m)
for(i in 1:m){
  Coverage_95_Jack[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_Jack_lows_95[u,i] && thetaiexps[u,i]<=CI_Jack_ups_95[u,i])})) 
}
mean(Coverage_95_Jack,na.rm=TRUE)

Coverage_99_Jack <- rep(NA,m)
for(i in 1:m){
  Coverage_99_Jack[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_Jack_lows_99[u,i] && thetaiexps[u,i]<=CI_Jack_ups_99[u,i])})) 
}
mean(Coverage_99_Jack,na.rm=TRUE)

Coverage_90_Boot <- rep(NA,m)
for(i in 1:m){
  Coverage_90_Boot[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_Boot_lows_90[u,i] && thetaiexps[u,i]<=CI_Boot_ups_90[u,i])})) 
}
mean(Coverage_90_Boot)

Coverage_95_Boot <- rep(NA,m)
for(i in 1:m){
  Coverage_95_Boot[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_Boot_lows_95[u,i] && thetaiexps[u,i]<=CI_Boot_ups_95[u,i])})) 
}
mean(Coverage_95_Boot)

Coverage_99_Boot <- rep(NA,m)
for(i in 1:m){
  Coverage_99_Boot[i] <- mean(sapply(1:5,function(u){
    as.integer(thetaiexps[u,i]>=CI_Boot_lows_99[u,i] && thetaiexps[u,i]<=CI_Boot_ups_99[u,i])})) 
}
mean(Coverage_99_Boot)



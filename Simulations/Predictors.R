# ------------------------------------------------------------------------------------
# Title: R code for "An unbiased predictor for skewed response variable with  
# measurement error in covariate"
# Note: This code illustrates the performance of predictors in Tables 1--3
# and Figure 1.
# Author: Sepideh Mosaferi
# Date: March 2023
# ------------------------------------------------------------------------------------

setwd("/Users/sepidehmosaferi/Desktop/Code/Simulations")
rm(list = ls(all = TRUE))

source("pars_estimation.R")
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

# names of predictors: 
thetaiexps <- c()
predexp0s <- c()
predexp1s <- c()
predexpdirs <- c()
predexp2s <- c()
predexp3s <- c()
R1is <- c()
mse_Js <- c()

# values of parameters:
mux <- 5
sdx <-  3
sdv <- sqrt(2)

set.seed(1200) #for the sake of reproducibility
psi <- rgamma(m,shape=4.5,scale=2)
beta0 <- 0
beta1 <- 3
sdu <- sqrt(sample(c(0,2),size=m,prob=c(1-100/100,100/100),replace=TRUE)) #percentage of measurement error
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
  
  if(iter == 2000){break}
  
}

apply((eblupmes - thetais)^2, 2, mean)

mean(apply((eblupmes - thetais)^2, 2, mean))
mean(apply(mseest, 2, mean))
mean(apply((ebmeyls - thetais)^2, 2, mean))

apply(newparests, 2, mean)

apply(betahats, 2, mean)
mean(sig2bhatmes)

tapply(apply((eblupmes - thetais)^2, 2, mean), Cis, mean)
tapply(apply((ebmeyls  - thetais)^2, 2, mean), Cis, mean)

#####  Average MSE by value of SDU: 
round(log(tapply(apply((predexpdirs - thetaiexps)^2, 2, mean),sdu,mean)),3)
round(log(tapply(apply((predexp0s - thetaiexps)^2, 2, mean),sdu,mean)),3)
round(log(tapply(apply((predexp1s - thetaiexps)^2, 2, mean),sdu,mean)),3)
round(log(tapply(apply((predexp2s - thetaiexps)^2, 2, mean),sdu,mean)),3)

#####  Predictors per small area and their overall mean:
apply(thetaiexps,2,mean); mean(apply(thetaiexps,2,mean))
apply(predexp1s,2,mean); mean(apply(predexp1s,2,mean))
apply(predexp0s,2,mean); mean(apply(predexp0s,2,mean))
apply(predexpdirs,2,mean); mean(apply(predexpdirs,2,mean))
apply(predexp2s,2,mean); mean(apply(predexp2s,2,mean))

## Table 1 of Manuscript:
Result <- as.matrix(cbind(sdu,log(apply(thetaiexps,2,mean)),log(apply(predexpdirs,2,mean)),
                          log(apply(predexp0s,2,mean)),log(apply(predexp1s,2,mean)),
                          log(apply(predexp2s,2,mean))))
colnames(Result) <- c("sdu","True","Direct","No-ME","Pred A","Pred B")
apply(Result[,-1],2,mean) #Average of results

#####  empirical MSE per small area and their overall mean:
mean(apply((predexp1s - thetaiexps)^2, 2, mean))
mean(apply((predexp0s - thetaiexps)^2, 2, mean))
mean(apply((predexpdirs - thetaiexps)^2, 2, mean))
mean(apply((predexp2s - thetaiexps)^2, 2, mean))

## Table 2 of Manuscript:
logmean_mse_Js <- rep(0,m) #jackknife
for(i in 1:m){
  logmean_mse_Js[i] <- log(mean(mse_Js[,i][mse_Js[,i] >= 0]))
}

Result2 <- as.matrix(cbind(sdu,log(apply((predexpdirs-thetaiexps)^2,2,mean)),
                           log(apply((predexp0s - thetaiexps)^2, 2, mean)),
                           log(apply((predexp1s - thetaiexps)^2, 2, mean)),
                           log(apply((predexp2s - thetaiexps)^2, 2, mean)),
                           log(apply(R1is,2,mean)),logmean_mse_Js))
colnames(Result2) <- c("sdu","EMSE Direct","EMSE No-ME","EMSE Pred A","EMSE Pred B","R1is","Jakknife")
apply(Result2[,-1],2,mean) #Average of results



##### EMSE averaged by value of sdu
Result2 <- as.data.frame(Result2)
apply(Result2[Result2$sdu==unique(Result2$sdu)[1],],2,mean)
apply(Result2[Result2$sdu==unique(Result2$sdu)[2],],2,mean)

##### Ratio of average MSE predictor to average MSE direct 
mean(apply((predexp1s - thetaiexps)^2, 2, mean))/mean(apply((predexpdirs - thetaiexps)^2, 2, mean))
mean(apply((predexp0s - thetaiexps)^2, 2, mean))/mean(apply((predexpdirs - thetaiexps)^2, 2, mean))
mean(apply((predexp2s - thetaiexps)^2, 2, mean))/mean(apply((predexpdirs - thetaiexps)^2, 2, mean))

## Table 3 of Manuscript: 
##### Ratio of average MSE predictor to average MSE direct by values of sdu
tapply(apply((predexp1s - thetaiexps)^2, 2, mean),sdu,mean)/
  tapply(apply((predexpdirs - thetaiexps)^2, 2, mean),sdu,mean)

tapply(apply((predexp0s - thetaiexps)^2, 2, mean),sdu,mean)/
  tapply(apply((predexpdirs - thetaiexps)^2, 2, mean),sdu,mean)

tapply(apply((predexp2s - thetaiexps)^2, 2, mean),sdu,mean)/
  tapply(apply((predexpdirs - thetaiexps)^2, 2, mean),sdu,mean)


## Figure 1 of Manuscript: 
## Plots of relative bias and relative root MSE

Yi <- apply(thetaiexps,2,mean)
PredA.RB <- (apply(predexp1s,2,mean)-Yi)/Yi #predictor A
PredB.RB <- (apply(predexp2s,2,mean)-Yi)/Yi #predictor B
PredA.RRMSE <- sqrt(apply(predexp1s^2,2,mean)+Yi^2-2*apply(predexp1s,2,mean)*Yi)/Yi #predictor A
PredB.RRMSE <- sqrt(apply(predexp2s^2,2,mean)+Yi^2-2*apply(predexp2s,2,mean)*Yi)/Yi #predictor B

small_area <- c(1:m,1:m)
Predictor <- c(rep("A",m),rep("B",m))
RB <- c(PredA.RB,PredB.RB)
RRMSE <- c(PredA.RRMSE,PredB.RRMSE)

DATA_RB <- as_tibble(data.frame(cbind(Predictor,RB)))
DATA_RB <- DATA_RB %>% mutate(RB=as.numeric(RB))
DATA_RB <- DATA_RB %>% arrange(desc(RB)) 
DATA_RB$small_area <- small_area
DATA_RRMSE <- as_tibble(data.frame(cbind(Predictor,RRMSE)))
DATA_RRMSE <- DATA_RRMSE %>% mutate(RRMSE=as.numeric(RRMSE),samll_area=small_area)
DATA_RRMSE <- DATA_RRMSE %>% arrange(desc(RRMSE)) 
DATA_RRMSE$small_area <- small_area

Plot1 <- ggplot(DATA_RB, aes(x = small_area, y = RB, shape= Predictor, color = Predictor)) +
  geom_point(size=2.5)+scale_shape_manual(values = c(15, 16))+geom_line()+ 
  scale_color_manual(values = c(15, 16))+
  labs(x = "small area index (descending order based on RB)", y = "RB", 
       title="(a) Comparison of Relative Bias",color = "Predictor")+theme_light()

Plot2 <- ggplot(DATA_RRMSE, aes(x = small_area, y = RRMSE, shape=Predictor, color = Predictor)) +
  geom_point(size=2.5)+scale_shape_manual(values = c(15, 16))+geom_line()+ 
  scale_color_manual(values = c(15, 16))+
  labs(x = "small area index (descending order based on RRMSE)", y = "RRMSE", 
       title="(b) Comparison of Relative Root MSE",color = "Predictor")+theme_light()

grid.arrange(Plot1, Plot2, ncol=2)


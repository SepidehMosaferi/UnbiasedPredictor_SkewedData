#------------------------------------------------------------------------------------
# Title: R code for "An unbiased predictor for skewed response variable with  
# measurement error in covariate"
# Code for Application of Census of Governments (COG)
# Author: Sepideh Mosaferi
# Date: March 2023
#------------------------------------------------------------------------------------
setwd("/Users/sepidehmosaferi/Desktop/Code/Applications")
rm(list = ls(all = TRUE))

require("PracTools"); require("sae"); require("dplyr"); require("infer")
require("tidyverse"); require("gridExtra") 
source("Jackknife.R"); source("pars_estimation.R"); source("Bootstrap_CI.R")

# reading data
load(file="COG.RData")

set.seed(1200) #for the sake of reproducibility
m <- 49  #total number of small areas

## True values
MEANy <- tapply(FINALDATA$TOTEMP12rep,INDEX=FINALDATA$STATE,FUN=mean)

## Sample Selection for X 
Nh <- as.vector(table(FINALDATA$STATE))  #STATEs as STRATA
N <- sum(Nh)
nx <- 80000  #Overall sample size
ALLOC <- strAlloc(n.tot=nx,Nh=Nh,alloc="prop")  #Proportional Allocation
nhx <- floor(ALLOC$nh) #final selected sample size

SMP.IDs_x <- sapply(1:m,function(i){sample(FINALDATA$id[FINALDATA$STATE==unique(FINALDATA$STATE)[i]],
                                           size=nhx[i],replace=FALSE)})
SMP.IDs_x <- unlist(SMP.IDs_x)
SAMPLE_x <- FINALDATA[FINALDATA$id %in% SMP.IDs_x,] 

## mean and sample variance
SAMPLEMEANx <- tapply(SAMPLE_x$TOTEMP07rep,INDEX=SAMPLE_x$STATE,FUN=mean)
SAMPLEMEANx <- ifelse(SAMPLEMEANx==0,1,SAMPLEMEANx)
SAMPLEMEANx <- as.vector(SAMPLEMEANx)
SAMPLESDx <- tapply(SAMPLE_x$TOTEMP07rep,INDEX=SAMPLE_x$STATE,FUN=sd)
SAMPLEVARx <- as.vector(SAMPLESDx^2)

## Final Point Estimates
Xmatrix <- as.vector(SAMPLEMEANx)
LOGX <- as.vector(log(Xmatrix))
VARLOGX <- SAMPLEVARx/(Xmatrix^2)

## Sample Selection for Y
Nh <- as.vector(table(FINALDATA$STATE))
N <- sum(Nh)
ny <- 8000  #Overall sample size
ALLOC <- strAlloc(n.tot=ny,Nh=Nh,alloc="prop")  #Proportional Allocation
nhy <- floor(ALLOC$nh) #final selected sample size

SMP.IDs_y <- sapply(1:m,function(i){sample(FINALDATA$id[FINALDATA$STATE==unique(FINALDATA$STATE)[i]],
                                           size=nhy[i],replace=FALSE)})
SMP.IDs_y <- unlist(SMP.IDs_y)
SAMPLE_y <- FINALDATA[FINALDATA$id %in% SMP.IDs_y,] 

## Direct Estimator
SAMPLEMEANy <- tapply(SAMPLE_y$TOTEMP12rep,INDEX=SAMPLE_y$STATE,FUN=mean)
SAMPLEMEANy <- ifelse(SAMPLEMEANy==0,1,SAMPLEMEANy)
SAMPLEMEANy <- as.vector(SAMPLEMEANy)
SAMPLESDy <- tapply(SAMPLE_y$TOTEMP12rep,INDEX=SAMPLE_y$STATE,FUN=sd)
SAMPLEVARy <- as.vector(SAMPLESDy^2)

#Log transformation
LOGY <- as.vector(log(SAMPLEMEANy)) 
VARLOGY <- SAMPLEVARy/(SAMPLEMEANy^2)

## plots
DATA <- as_tibble(cbind(SAMPLEMEANy,SAMPLEMEANx,log(SAMPLEMEANy),log(SAMPLEMEANx),1:m))
DATA <- DATA %>% rename(log_SAMPLEMEANy="V3",log_SAMPLEMEANx="V4", area="V5")

Plot1 <- ggplot(DATA, aes(x=SAMPLEMEANy)) + 
  geom_histogram(color="black", fill="white")+
  labs(title="(a) Sample Mean from 2012",x="y",y="count")+theme_light()

Plot2 <- ggplot(DATA, aes(x=SAMPLEMEANx)) + 
  geom_histogram(color="black", fill="white")+
  labs(title="(b) Sample Mean from 2007",x="W",y="count")+theme_light()

Plot3 <- ggplot(DATA, aes(x = SAMPLEMEANx, y = SAMPLEMEANy)) +
  geom_point(size=2.5)+
  labs(x = "W", y = "y",title="(c) Scatterplot of y v.s. W")+theme_light()+
  geom_abline(intercept = 0, slope = 1, size = 0.25,color="grey",linetype = "dashed")

grid.arrange(Plot1, Plot2, Plot3, ncol=3)

Plot4 <- ggplot(DATA, aes(x=log_SAMPLEMEANy)) + 
  geom_histogram(color="black", fill="white",bins=20)+
  labs(title="(d) Sample Mean from 2012",x="log(y)",y="count")+theme_light()

Plot5 <- ggplot(DATA, aes(x=log_SAMPLEMEANx)) + 
  geom_histogram(color="black", fill="white",bins=20)+
  labs(title="(e) Sample Mean from 2007",x="log(W)",y="count")+theme_light()

Plot6 <- ggplot(DATA, aes(x = log_SAMPLEMEANx, y = log_SAMPLEMEANy)) +
  geom_point(size=2.5)+
  labs(x = "log(W)", y = "log(y)",title="(f) Scatterplot of log(y) v.s. log(W)")+theme_light()+
  geom_abline(intercept = 0, slope = 1, size = 0.25,color="grey",linetype = "dashed")

grid.arrange(Plot4, Plot5, Plot6, ncol=3)

## Fitting regression
FINALFIT <- lm(LOGY~LOGX)

## Variance Study
# for direct estimator:
direct <- SAMPLEMEANy 
Var_y <- SAMPLEVARy
CV_direct <- sqrt(Var_y)/direct

#### Modeling ####
D <- m
psi <- VARLOGY
Cis <- VARLOGX 
sdu <- sqrt(Cis)
Xi <- LOGX
yi <- LOGY 
beta0 <- as.vector(FINALFIT$coefficients[1]) 
beta1 <- as.vector(FINALFIT$coefficients[2]) 
sdv <- sqrt(mean(FINALFIT$residuals^2))

saedat <- data.frame(yi, Xi, psi) 
efh <- eblupFH(yi~Xi, vardir = psi, data = saedat, method = "ML") 

## Estimating unknown parameters based on iterative algorithm
newpar.est <- FHME_estimation(yi, Xi, psi, Cis)$para

# EB predictor A:
gammaeyl <- (newpar.est[3] + newpar.est[2]^2*Cis)/(newpar.est[3] + newpar.est[2]^2*Cis+psi)
thetahatmeyl <- gammaeyl*yi + (1-gammaeyl)*as.vector(cbind(1, Xi)%*%newpar.est[c(1,2)]) 
EB.A <- as.vector(exp(thetahatmeyl + 0.5*psi*gammaeyl)) 
var_EB.A <- as.vector(exp(psi*gammaeyl)*(exp(psi*gammaeyl)-1)*exp(2*thetahatmeyl)) 

# EB predictor B:
gammaeyl2 <- newpar.est[3]/(newpar.est[2]^2*Cis + psi + newpar.est[3])
nuil <- (yi - newpar.est[1] - newpar.est[2]*Xi)
xdot <- Xi+(newpar.est[2]*Cis/(newpar.est[2]^2*Cis + psi + newpar.est[3]))*nuil
sigmatilde2 <- Cis-(newpar.est[2]^2*Cis^2/(newpar.est[2]^2*Cis+psi+newpar.est[3]))

predexp_a <- exp(newpar.est[1] + newpar.est[2]*xdot)
predexp_b <- exp(-newpar.est[2]^2*sigmatilde2/2)
predexp_c <- exp( gammaeyl2*nuil + 0.5*(gammaeyl2*(newpar.est[2]^2*Cis + psi)))
EB.B <- as.vector(predexp_a*predexp_b*predexp_c) 

di <- 2*psi*newpar.est[2]^2*Cis/(newpar.est[2]^2*Cis+newpar.est[3]+psi)
M2i <- 1-2*exp(1.5*di-gammaeyl*psi)+exp(di-gammaeyl*psi)
R1i <- (M2i^2)*exp(4*(yi-psi))*(1-exp(-4*newpar.est[3]-4*psi))

## MSE of predictor B based on Jackknife

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
  R2i_Jack[u] <- ((m-1)/m)*sum(thetaB_del[u,]-EB.B[u])^2
}

mse_J <- R1i_Jack+R2i_Jack   #Final Jackknife MSE estimator 

## Bootstrap Prediction Interval for Predictor B
CI_Boot <- Bootstrap(yi,Xi,psi,Cis)

CI_Boot_90 <- CI_Boot[,1:m]
CI_Boot_90_low <- CI_Boot_90[1,] 
CI_Boot_90_up <- CI_Boot_90[2,] 

CI_Boot_95 <- CI_Boot[,(m+1):(2*m)]
CI_Boot_95_low <- CI_Boot_95[1,] 
CI_Boot_95_up <- CI_Boot_95[2,] 

CI_Boot_99 <- CI_Boot[,(2*m+1):(3*m)]
CI_Boot_99_low <- CI_Boot_99[1,] 
CI_Boot_99_up <- CI_Boot_99[2,] 

## Plots of Prediction Intervals

x <- 1:m
direct <- direct 
CI_direct_lower_90 <- direct-1.64*sqrt(Var_y) 
CI_direct_upper_90 <- direct+1.64*sqrt(Var_y)  
CI_direct_lower_95 <- direct-1.96*sqrt(Var_y) 
CI_direct_upper_95 <- direct+1.96*sqrt(Var_y)  
CI_direct_lower_99 <- direct-2.58*sqrt(Var_y) 
CI_direct_upper_99 <- direct+2.58*sqrt(Var_y)  

EB.B <- EB.B 
CI_B_lower_90 <- EB.B-1.64*sqrt(mse_J)  
CI_B_upper_90 <- EB.B+1.64*sqrt(mse_J)
CI_B_lower_95 <- EB.B-1.96*sqrt(mse_J)  
CI_B_upper_95 <- EB.B+1.96*sqrt(mse_J)
CI_B_lower_99 <- EB.B-2.58*sqrt(mse_J) 
CI_B_upper_99 <- EB.B+2.58*sqrt(mse_J)

method <- c(rep("Direct",m),rep("Bootstrap",m),rep("Jackknife",m))
length_direct_95 <- CI_direct_upper_95-CI_direct_lower_95 
length_Boot_95 <- as.numeric(CI_Boot_95_up-CI_Boot_95_low)
length_Jack_95 <- CI_B_upper_95-CI_B_lower_95 
CI_length <- c(length_direct_95,length_Boot_95,length_Jack_95)
DATA_length <- as_tibble(cbind(method,CI_length))
DATA_length <- DATA_length %>% mutate(CI_length=as.numeric(CI_length),
                                      log_CI_length=log(as.numeric(CI_length)))

ggplot(DATA_length, aes(x=factor(method), y=log_CI_length,color=method))+
  geom_boxplot()+
  labs(x="method",y="length",title="(b) Distribution of the Length of Prediction Intervals in Logarithmic Scale for n=8000")+theme_light()

Summary <- DATA_length %>% group_by(method) %>% summarize(min=min(log_CI_length,na.rm=TRUE),Q1=quantile(log_CI_length,0.25,na.rm=TRUE),Q2=quantile(log_CI_length,0.5,na.rm=TRUE),
                                                          mean=mean(log_CI_length,na.rm=TRUE),Q3=quantile(log_CI_length,0.75,na.rm=TRUE),max=max(log_CI_length,na.rm=TRUE))



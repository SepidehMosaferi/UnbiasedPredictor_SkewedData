#------------------------------------------------------------------------------------
# Title: R code for "An unbiased predictor for skewed response variable with  
# measurement error in covariate"
# Code for Application of SAIPE project
# Author: Sepideh Mosaferi
# Date: March 2023
#------------------------------------------------------------------------------------
getwd()
setwd("/Users/sepidehmosaferi/Desktop/Code/Applications")
rm(list = ls(all = TRUE))

require("PracTools"); require("sae"); require("dplyr"); require("infer") 
require("tidyverse"); require("gridExtra") 
source("Jackknife.R"); source("pars_estimation.R"); source("Bootstrap_CI.R")

# reading data
load(file="SAIPE.RData")

y2018 <- ACS2018$Estimate..Total..Population.for.whom.poverty.status.is.determined..AGE..Under.18.years..5.to.17.years
sd2018 <- ACS2018$Margin.of.Error..Total.MOE..Population.for.whom.poverty.status.is.determined..AGE..Under.18.years..5.to.17.years/1.645
var2018 <- sd2018^2 

x2017_5 <- ACS2017_5years$Estimate..Total..Population.for.whom.poverty.status.is.determined..AGE..Under.18.years..5.to.17.years 
sd2017_5 <- ACS2017_5years$Margin.of.Error..Total.MOE..Population.for.whom.poverty.status.is.determined..AGE..Under.18.years..5.to.17.years/1.645 
var2017_5 <- sd2017_5^2

# Predictors
logy <- log(y2018)
varlogy <- var2018/(y2018^2)
logx <- log(x2017_5)
varlogx <- var2017_5/(x2017_5^2)
m <- length(logx) 

# Fitting regression
FINALFIT <- lm(logy~logx)

#### Modeling ####
D <- m
psi <- varlogy 
Cis <- varlogx
sdu <- sqrt(Cis)
Xi <- logx
yi <- logy 
beta0 <- as.vector(FINALFIT$coefficients[1]) 
beta1 <- as.vector(FINALFIT$coefficients[2]) 
sdv <- sqrt(mean(FINALFIT$residuals^2))

saedat <- data.frame(yi, Xi, psi) 
efh <- eblupFH(yi~Xi, vardir = psi, data = saedat, method = "ML") 

## Estimating unknown parameters based on iterative algorithm
newpar.est <- FHME_estimation(yi, Xi, psi, Cis)$para

CV_direct <- sd2018/y2018 
x <- 1:D

#EB predictor A:
gammaeyl <- (newpar.est[3] + newpar.est[2]^2*Cis)/(newpar.est[3] + newpar.est[2]^2*Cis+psi)
thetahatmeyl <- gammaeyl*yi + (1-gammaeyl)*as.vector(cbind(1, Xi)%*%newpar.est[c(1,2)]) 
EB.A <- as.vector(exp(thetahatmeyl + 0.5*psi*gammaeyl)) #predictor 
var_EB.A <- exp(psi*gammaeyl)*(exp(psi*gammaeyl)-1)*exp(2*thetahatmeyl) #variance

# EB predictor B:
gammaeyl2 <- newpar.est[3]/(newpar.est[2]^2*Cis + psi + newpar.est[3])
nuil <- (yi - newpar.est[1] - newpar.est[2]*Xi) 
xdot <- Xi+(newpar.est[2]*Cis/(newpar.est[2]^2*Cis + psi + newpar.est[3]))*nuil
sigmatilde2 <- Cis-(newpar.est[2]^2*Cis^2/(newpar.est[2]^2*Cis+psi+newpar.est[3]))

predexp_a <- exp(newpar.est[1] + newpar.est[2]*xdot)
predexp_b <- exp(-newpar.est[2]^2*sigmatilde2/2)
predexp_c <- exp( gammaeyl2*nuil + 0.5*(gammaeyl2*(newpar.est[2]^2*Cis + psi)))
EB.B <- as.vector(predexp_a*predexp_b*predexp_c)  #predictor

di <- 2*psi*newpar.est[2]^2*Cis/(newpar.est[2]^2*Cis+newpar.est[3]+psi)
M2i <- 1-2*exp(1.5*di-gammaeyl*psi)+exp(di-gammaeyl*psi)
R1i <- (M2i^2)*exp(4*(yi-psi))*(1-exp(-4*newpar.est[3]-4*psi)) #variance

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


# Predictors FH Modeling
X_SNAP <- log(SNAP$`July 2018`)
saedat2 <- data.frame(yi, X_SNAP, psi) 
efh2 <- eblupFH(yi~X_SNAP, vardir = psi, data = saedat2, method = "ML") 
taui0 <- yi-efh2$fit$estcoef$beta[1]-efh2$fit$estcoef$beta[2]*X_SNAP
gammai0 <- efh2$fit$refvar/(efh2$fit$refvar+psi) 
Theta_SNAP <- exp(efh2$fit$estcoef$beta[1]+efh2$fit$estcoef$beta[2]*X_SNAP+
                    gammai0*taui0+(gammai0*psi)/2)
MSE_SNAP <- Theta_SNAP^2*(exp(gammai0*psi)-1) 
CV_SNAP <- sqrt(MSE_SNAP)/Theta_SNAP

## Plots of data sets

as_tibble(cbind(x2017_5,y2018,SNAP$`July 2018`,log(x2017_5),log(y2018),log(SNAP$`July 2018`)))
area <- c(rep(1:m,2))
y <- rep(y2018,2)
x <- c(x2017_5,SNAP$`July 2018`)
variable <- c(rep("ACS5yrs",m),rep("SNAP",m))
DATA_bftrans <- as_tibble(cbind(area,x,y,variable)) 
DATA_bftrans <- DATA_bftrans %>% mutate(x=as.numeric(x),y=as.numeric(y),area=as.numeric(area))
DATA_aftrans <- DATA_bftrans %>% mutate(logx=log(x),logy=log(y),
                                        variable=recode(variable,ACS5yrs="log(ACS5yrs)",SNAP="log(SNAP)"))

Plot1 <- ggplot(DATA_bftrans, aes(x = x, y = y, color = variable,shape=variable)) +
  geom_point(size=2.5)+scale_shape_manual(values = c(4, 5)) +
  scale_color_manual(values = c("black", "grey"))+
  geom_abline(intercept = 0, slope = 1, size = 0.25,color="grey",linetype = "dashed")+
  labs(x = "W", y = "y", 
       title="(a) Before Transformation")+theme_light()

Plot2 <- ggplot(DATA_aftrans, aes(x = logx, y = logy, color = variable,shape=variable)) +
  geom_point(size=2.5)+scale_shape_manual(values = c(4, 5)) +
  scale_color_manual(values = c("black", "grey"))+
  geom_abline(intercept = 0, slope = 1, size = 0.25,color="grey",linetype = "dashed")+
  labs(x = "log(W)", y = "log(y)", 
       title="(b) After Transformation")+theme_light()

grid.arrange(Plot1, Plot2, ncol=2)

cov <- c(rep("log(ACS5yrs)",m),rep("log(SNAP)",m))
Covariate <- as_tibble(cbind(rep(1:m,2),c(Xi,X_SNAP),cov)) 
Covariate2 <- Covariate %>% rename(area=V1,value=V2) %>% mutate(area=as.numeric(area),value=as.numeric(value))

ggplot(Covariate2, aes(x=factor(cov), y=value,color=cov))+
  geom_boxplot(col="black")+guides(col = FALSE)+theme_light()+labs(x = "covariates") 

## Plots of Prediction Intervals

direct <- y2018
CI_Direct_90_low <- direct-1.64*sqrt(var2018) 
CI_Direct_90_up <- direct+1.64*sqrt(var2018)  
CI_Direct_95_low <- direct-1.96*sqrt(var2018) 
CI_Direct_95_up <- direct+1.96*sqrt(var2018)  
CI_Direct_99_low <- direct-2.58*sqrt(var2018) 
CI_Direct_99_up <- direct+2.58*sqrt(var2018)  

EB.B <- EB.B 
CI_Jack_90_low <- EB.B-1.64*sqrt(mse_J)  
CI_Jack_90_up <- EB.B+1.64*sqrt(mse_J)
CI_Jack_95_low <- EB.B-1.96*sqrt(mse_J)  
CI_Jack_95_up <- EB.B+1.96*sqrt(mse_J)
CI_Jack_99_low <- EB.B-2.58*sqrt(mse_J) 
CI_Jack_99_up <- EB.B+2.58*sqrt(mse_J)

Theta_SNAP <- Theta_SNAP 
CI_SNAP_90_low <- Theta_SNAP-1.64*sqrt(MSE_SNAP)  
CI_SNAP_90_up <- Theta_SNAP+1.64*sqrt(MSE_SNAP)
CI_SNAP_95_low <- Theta_SNAP-1.96*sqrt(MSE_SNAP)  
CI_SNAP_95_up <- Theta_SNAP+1.96*sqrt(MSE_SNAP)
CI_SNAP_99_low <- Theta_SNAP-2.58*sqrt(MSE_SNAP) 
CI_SNAP_99_up <- Theta_SNAP+2.58*sqrt(MSE_SNAP)

method <- c(rep("Direct",m),rep("Bootstrap",m),rep("Jackknife",m),rep("FHeblup",m))
length_Direct_95 <- CI_Direct_95_up-CI_Direct_95_low 
length_Boot_95 <- as.numeric(CI_Boot_95_up-CI_Boot_95_low)
length_Jack_95 <- CI_Jack_95_up-CI_Jack_95_low 
length_SNAP_95 <- CI_SNAP_95_up-CI_SNAP_95_low 
CI_length <- c(length_Direct_95,length_Boot_95,length_Jack_95,length_SNAP_95)
DATA_length <- as_tibble(cbind(method,CI_length))
DATA_length <- DATA_length %>% mutate(CI_length=as.numeric(CI_length),
                                      log_CI_length=log(as.numeric(CI_length)))

ggplot(DATA_length, aes(x=factor(method), y=log_CI_length,color=method))+
  geom_boxplot()+guides(col = FALSE)+ 
  labs(x="method",y="length",title="Distribution of the Length of Prediction Intervals in Logarithmic Scale")+theme_light()


Summary <- DATA_length %>% group_by(method) %>% summarize(Q2=quantile(log_CI_length,0.5,na.rm=TRUE),mean=mean(log_CI_length,na.rm=TRUE),SD=sd(log_CI_length,na.rm=TRUE))


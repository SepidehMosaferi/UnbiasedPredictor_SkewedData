# ------------------------------------------------------------------------------
#  R code for implementing an iterative algorithm for parameter estimation 
# Author: Shonosuke Sugasawa
# ------------------------------------------------------------------------------

### INPUT 
# yi: response variable (direct estimate)
# Xi: single covariate (measured with error)  
# psi: sampling variance 
# Cis: MSE of Xi

FHME_estimation <- function(yi, Xi, psi, Cis, maxitr=1000, tolerance=10^(-5)){
  # initial value
  saedat <- data.frame(yi, Xi, psi)
  efh <- eblupFH(yi~Xi, vardir = psi, data = saedat, method = "ML")
  if(is.na(mean(efh$eblup))){
    lm.fit <- lm(yi~Xi)
    beta0 <- coef(lm.fit)[1]
    beta1 <- coef(lm.fit)[2]
    sigv <- mean(lm.fit$residuals^2) 
  }else{
    beta0 <- efh$fit$estcoef$beta[1]
    beta1 <- efh$fit$estcoef$beta[2] 
    sigv <- efh$fit$refvar
  }

  # iteration 
  for(itr in 1:maxitr){
    para <- c(beta0, beta1, sigv)
    taui <- yi-beta0-beta1*Xi
    sigma2bi <- beta1^2*Cis+sigv^2+psi
    # update beta0
    beta0 <- sum((yi-beta1*Xi)/sigma2bi) / sum(1/sigma2bi)
    # update beta1
    taui <- yi-beta0-beta1*Xi
    num <- sum(Xi*(yi-beta0)/sigma2bi)
    denom <- sum(Xi^2/sigma2bi) - sum(Cis*taui^2/sigma2bi^2)
    beta1 <- num/denom
    # update sigv
    taui <- yi-beta0-beta1*Xi
    sigma2bi <- beta1^2*Cis+sigv^2+psi
    eq <- function(sigv){
      taui <- yi-beta0-beta1*Xi
      sigma2bi <- beta1^2*Cis+sigv^2+psi
      U <- -0.5*sum(1/sigma2bi) + 0.5*sum(taui^2/sigma2bi^2)
      return(U)
    }
    if(eq(10^(-5))*eq(10^5)<0){
      sigv <- uniroot(f=eq, interval=c(10^(-5), 10^5))$root	
    }else{
      sigv <- sigv
    }
    para.new <- c(beta0, beta1, sigv) 
    # convergence check
    dd <- sum( abs(para.new - para) )
    if(dd<tolerance){ break }
  }
  
  # output
  result <- list(itr=itr, para=para.new)
  return(result)
}
  
  
  

# ------------------------------------------------------------------------------------
# R code for non-parametric bootstrap for the prediction interval
# Author: Sepideh Mosaferi
# ------------------------------------------------------------------------------------

Bootstrap <- function(yi,Xi,psi,Cis){
  
  saedat <- data.frame(yi, Xi, psi,Cis)
  
  nreps = 2000
  set.seed(2012)
  resamples <- saedat %>% rep_sample_n(size = nrow(saedat), replace = TRUE, reps = nreps)
  pred_B <- matrix(NA,nrow=nrow(saedat),ncol=nreps) #each row is one small area
  
  for(i in 1:nreps){
    
    saedat_i <- data.frame(resamples$yi[resamples$replicate==i], resamples$Xi[resamples$replicate==i], 
                           resamples$psi[resamples$replicate==i], resamples$Cis[resamples$replicate==i])
    colnames(saedat_i) <- c("yi","Xi","psi","Cis")
    
    newpar.est <- FHME_estimation(saedat_i$yi, saedat_i$Xi, saedat_i$psi, saedat_i$Cis)$para
    
    sigmatilde2 <- saedat_i$Cis-(newpar.est[2]^2*saedat_i$Cis^2/(newpar.est[2]^2*saedat_i$Cis+saedat_i$psi+newpar.est[3]))
    gammaeyl2 <- newpar.est[3]/(newpar.est[2]^2*saedat_i$Cis + saedat_i$psi + newpar.est[3])
    nuil <- (saedat_i$yi - newpar.est[1] - newpar.est[2]*saedat_i$Xi)
    xdot <- saedat_i$Xi + (newpar.est[2]*saedat_i$Cis/(newpar.est[2]^2*saedat_i$Cis + saedat_i$psi + newpar.est[3]))*nuil
    
    predexp2a <- exp(newpar.est[1] + newpar.est[2]*xdot)
    predexp2b <- exp(-newpar.est[2]^2*sigmatilde2/2)
    predexp2c <- exp(gammaeyl2*nuil + 0.5*(gammaeyl2*(newpar.est[2]^2*saedat_i$Cis + saedat_i$psi)))
    
    # output
    pred_B[,i] <- predexp2a*predexp2b*predexp2c
    
  }
  CI_90 <- sapply(1:nrow(saedat),function(i){as.vector(quantile(pred_B[i,],probs=c(0.05,0.95)))})
  CI_95 <- sapply(1:nrow(saedat),function(i){as.vector(quantile(pred_B[i,],probs=c(0.025,0.975)))})
  CI_99 <- sapply(1:nrow(saedat),function(i){as.vector(quantile(pred_B[i,],probs=c(0.005,0.995)))})
  return(data.frame(CI_90,CI_95,CI_99)) 
}


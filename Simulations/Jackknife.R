# ------------------------------------------------------------------------------
## R code for Jackknife estimator for the MSE
# Author: Sepideh Mosaferi
# ------------------------------------------------------------------------------------

Jackknife <- function(yi,Xi,psi,Cis,R1i_del,thetaB_del){
  
  pars_del <- array(0,c(3,1,m))
  for(u in 1:m){  
    yi_del <- yi[-u]; Xi_del <- Xi[-u]; psi_del <- psi[-u]; Cis_del <- Cis[-u]   
    efh_del <- eblupFH(yi_del~Xi_del, vardir = psi_del, method = "ML")
    
    ## Estimating unknown parameters based on iterative algorithm
    pars_del[,,u] <- FHME_estimation(yi_del, Xi_del, psi_del, Cis_del)$para
    }
  
  gammaeyl_del <- matrix(0,m,m)  #each row belongs to one area
  di_del <- matrix(0,m,m) 
  for(u in 1:m){
    gammaeyl_del[u,] <- sapply(1:m,function(i){
      (pars_del[,,i][3] + pars_del[,,i][2]^2*Cis[u])/(pars_del[,,i][3] + pars_del[,,i][2]^2*Cis[u]+psi[u])})
    di_del[u,] <- sapply(1:m,function(i){2*psi[u]*pars_del[,,i][2]^2*Cis[u]/(pars_del[,,i][2]^2*Cis[u]+pars_del[,,i][3]+psi[u])})
  }
  
  M2i_del <- matrix(0,m,m) #each row belongs to one area
  for(u in 1:m){
    M2i_del[u,] <- 1-2*exp(1.5*di_del[u,]-gammaeyl_del[u,]*psi[u])+exp(di_del[u,]-gammaeyl_del[u,]*psi[u])
    }
  
  for(u in 1:m){
    R1i_del[u,] <- sapply(1:m,function(i){(M2i_del[u,i]^2)*exp(4*(yi[u]-psi[u]))*(1-exp(-4*pars_del[,,i][3]-4*psi[u]))})
    }
  
  gammaeyl2_del <- matrix(0,m,m)  #each row belongs to one area
  nuil_del <- matrix(0,m,m)
  xdot_del <- matrix(0,m,m)
  sigmatilde2_del <- matrix(0,m,m)
  for(u in 1:m){
    gammaeyl2_del[u,] <- sapply(1:m,function(i){pars_del[,,i][3]/(pars_del[,,i][2]^2*Cis[u] + psi[u] + pars_del[,,i][3])})
    nuil_del[u,] <- sapply(1:m,function(i){(yi[u] - pars_del[,,i][1] - pars_del[,,i][2]*Xi[u])})
    xdot_del[u,] <- sapply(1:m,function(i){Xi[u]+(pars_del[,,i][2]*Cis[u]/(pars_del[,,i][2]^2*Cis[u] + psi[u] + pars_del[,,i][3]))*nuil_del[u,i]})
    sigmatilde2_del[u,] <- sapply(1:m,function(i){Cis[u]-(pars_del[,,i][2]^2*Cis[u]^2/(pars_del[,,i][2]^2*Cis[u]+psi[u]+pars_del[,,i][3]))})
    }
  
  predexp_a_del <- matrix(0,m,m)  #each row belongs to one area
  predexp_b_del <- matrix(0,m,m)
  predexp_c_del <- matrix(0,m,m)
  for(u in 1:m){
    predexp_a_del[u,] <- sapply(1:m,function(i){exp(pars_del[,,i][1] + pars_del[,,i][2]*xdot_del[u,i])})
    predexp_b_del[u,] <- sapply(1:m,function(i){exp(-pars_del[,,i][2]^2*sigmatilde2_del[u,i]/2)}) 
    predexp_c_del[u,] <- sapply(1:m,function(i){
      exp(gammaeyl2_del[u,i]*nuil_del[u,i]+0.5*(gammaeyl2_del[u,i]*(pars_del[,,i][2]^2*Cis[u] + psi[u])))}) 
  }

  for(u in 1:m){
    thetaB_del[u,] <- as.vector(predexp_a_del[u,]*predexp_b_del[u,]*predexp_c_del[u,])  
  }
  
  return(data.frame(thetaB_del,R1i_del))
}


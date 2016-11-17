gammakernel <- function(y, x, z, betas.ini,gammas.now,gammas.old,gpri,Gpri,model,m,ni) {
  
  if (model=="NB1"){
    Etas_1 = x%*%betas.ini
    Mu = exp(Etas_1)
    Etas_2N = z%*%gammas.old
    Alpha_N = exp(Etas_2N)
    
    Ym2_N = Etas_2N+y/Mu-1
    SIGMA2_N=diag(as.vector((Mu+Alpha_N)/(Alpha_N*Mu)))
    
  } else if (model=="NB2") {
    
    Etas_1 = x%*%betas.ini
    Mu = exp(Etas_1)
    Etas_2N = z%*%gammas.old
    sigma2 = exp(Etas_2N)
    
    Ym2_N = Etas_2N+y/Mu-1
    SIGMA2_N=diag(as.vector(sigma2/Mu^2))
    
  } else {
    
    Etas_1 = x%*%betas.ini
    Etas_2N = z%*%gammas.old
    
    Mu=inv.logit(Etas_1)
    phi=exp(Etas_2N)
    
    if (is.null(m)){
      m=length(y)
    }
    
    if (is.null(ni)){
      
      nro <- length(y)
      
      ni <- rep(m,nro) 
      
    }
    
    Ym2_N = Etas_2N+(y/(Mu*ni))-1
    
    Q2=as.vector(((phi+ni)/(phi+1))*((1-Mu)/(ni+Mu)))
    SIGMA2_N=diag(Q2)
    
  }
  
  Gpos <-qr.solve(qr.solve(Gpri,tol = 1e-100)+ t(z)%*%solve(SIGMA2_N,tol = 1e-100)%*%z,tol = 1e-100)
  Gpos <- as.matrix(forceSymmetric(as.matrix(Gpos)))
  gpos <-Gpos%*%(qr.solve(Gpri,tol = 1e-100)%*%gpri + t(z)%*%solve(SIGMA2_N,tol = 1e-100)%*%Ym2_N)
  dmvnorm(t(gammas.now),gpos,Gpos) #These functions provide the density function for the multivariate normal
}
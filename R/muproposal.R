muproposal<-function(y, x, z, betas.ini,gammas.ini,bpri,Bpri,model,m,ni) {
  
  if (model=="NB1"){
    
    Etas_1 = x%*%betas.ini
    Etas_2 = z%*%gammas.ini
    Mu = exp(Etas_1)
    Alpha = exp(Etas_2)
    
    Ym_1 = Etas_1 + y/Mu-1 #Var de trabajo
    
    Q1=as.vector((Mu+Alpha)/Alpha*1/Mu)
    SIGMA1=diag(Q1)
    
  } else if (model=="NB2") {
    
    Etas_1 = x%*%betas.ini
    Etas_2 = z%*%gammas.ini
    Mu = exp(Etas_1)
    sigma2 = exp(Etas_2)
    
    Ym_1 = Etas_1 + y/Mu-1 #Var de trabajo
    
    Q1=as.vector(sigma2/Mu^2)
    SIGMA1=diag(Q1)       
    
  } else {
    
    Etas_1 = x%*%betas.ini
    Etas_2 = z%*%gammas.ini
    
    Mu=inv.logit(Etas_1)
    phi=exp(Etas_2)
    
    if (is.null(m)){
      m=length(y)
    }
    
    if (is.null(ni)){
      
      nro <- length(y)
      
      ni <- rep(m,nro) 
      
    }
    
    Ym_1 = Etas_1 + ((y/ni)-Mu)/(Mu*(1-Mu)) #Var de trabajo
    
    Q1=as.vector((ni*Mu*(1-Mu))^(-1)*(phi+ni)/(phi+1))
    SIGMA1=diag(Q1)          
    
  }
  
  Bpos <- qr.solve(qr.solve(Bpri,tol = 1e-100)+ t(x)%*%qr.solve(SIGMA1,tol = 1e-100)%*%x,tol = 1e-100)
  Bpos <- as.matrix(forceSymmetric(as.matrix(Bpos)))
  bpos <- Bpos%*%(qr.solve(Bpri,tol = 1e-100)%*%bpri + t(x)%*%qr.solve(SIGMA1,tol = 1e-100)%*%Ym_1)
  betas.pro <- rmvnorm(1,bpos,Bpos)
  betas.pro
}

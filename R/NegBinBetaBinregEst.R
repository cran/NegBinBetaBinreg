NegBinBetaBinregEst <-
  function (y,x,z,nsim,bpri,Bpri,gpri,Gpri,burn,jump,bini,gini,model,m,ni,graph1,graph2){
    
    
    ### Cepeda - Metropolis - Hastings
    
    Y=as.matrix(y)
    
    if(is.null(x)|is.null(z)|is.null(y)){
      stop("There is no data")
    }
    
    if(burn> 1 | burn < 0){
      stop("Burn must be a proportion between 0 and 1")
    }
    
    if(nsim <= 0){
      stop("the number of simulations must be greater than 0")
    }
    
    if(jump < 0|jump > nsim){
      stop("Jumper must be a positive number lesser than nsim")
    }
    
    
    ind1<-rep(0,nsim)
    ind2<-rep(0,nsim)
    
    if (is.null(bini)){
      betas.ind <- matrix(bpri,nrow=ncol(x))
    }else{
      betas.ind <- matrix(bini,nrow=ncol(x))
    }
    
    
    if (is.null(gini)){
      gammas.ind <- matrix(gpri,nrow=ncol(z))
    }else{
      gammas.ind <- matrix(gini,nrow=ncol(z))
    }
    
    beta.mcmc<-matrix(NA,nrow=nsim,ncol=ncol(x))
    gamma.mcmc<-matrix(NA,nrow=nsim,ncol=ncol(z))
    
    for(i in 1:nsim) {
      
      #Betas
      
      betas.sim <- matrix(muproposal(y,x,z,betas.ind,gammas.ind,bpri,Bpri,model,m,ni),nrow=ncol(x))
      gammas.sim <- matrix(gammaproposal(y,x,z,betas.ind,gammas.ind,gpri,Gpri,model,m,ni),nrow=ncol(z))
      
      q1.mu <- mukernel(y,x,z,betas.sim,betas.ind,gammas.sim,bpri,Bpri,model,m,ni)
      q2.mu <- mukernel(y,x,z,betas.ind,betas.sim,gammas.ind,bpri,Bpri,model,m,ni)
      p1.mu <- dpostb(y,x,z,betas.sim,gammas.ind,bpri,Bpri,model,m)
      p2.mu <- dpostb(y,x,z,betas.ind,gammas.ind,bpri,Bpri,model,m)
      
      alfa1<-((p1.mu/p2.mu)*(q1.mu/q2.mu))
      
      if(is.na(alfa1)==T |alfa1==Inf){
        alfa1=10
      }
      
      Mu.val<-min(1,alfa1)
      u<-runif(1)
      if (u <=Mu.val) {
        betas.ind <- betas.sim
        ind1[i] = 1
      }
      
      beta.mcmc[i, ] <- betas.ind
      beta.mcmc <- as.ts(beta.mcmc)
      
      q1.gamma <- gammakernel(y,x,z,betas.sim,gammas.sim,gammas.ind,gpri,Gpri,model,m,ni)
      q2.gamma <- gammakernel(y,x,z,betas.ind,gammas.ind,gammas.sim,gpri,Gpri,model,m,ni)
      p1.gamma <- dpostg(y,x,z,betas.ind,gammas.sim,gpri,Gpri,model,m)
      p2.gamma <- dpostg(y,x,z,betas.ind,gammas.ind,gpri,Gpri,model,m)  
      
      
      alfa2<-((p1.gamma/p2.gamma)*(q1.gamma/q2.gamma))
      
      if(is.na(alfa2)==T |alfa2==Inf){
        alfa2=10
      }
      
      Gamma.val<-min(1,alfa2)
      u<-runif(1)
      if (u <=Gamma.val) {
        gammas.ind <- gammas.sim
        ind2[i] = 1
      }
      gamma.mcmc[i,]<-gammas.ind
      gamma.mcmc <- as.ts(gamma.mcmc)
      
      if (i%%1000 == 0)
        cat("Burn-in iteration : ", i, "\n")
    }
    
    
    tburn <- nsim*burn
    extr <- seq(0,(nsim-tburn),jump)
    
    betas.burn <-as.matrix(beta.mcmc[(tburn+1):nrow(beta.mcmc),])
    gammas.burn <-as.matrix(gamma.mcmc[(tburn+1):nrow(gamma.mcmc),])
    
    beta.mcmc.auto <- as.matrix(betas.burn[extr,])
    beta.mcmc.auto <- as.ts(beta.mcmc.auto)
    gamma.mcmc.auto <- as.matrix(gammas.burn[extr,])
    gamma.mcmc.auto <- as.ts(gamma.mcmc.auto)
    
    
    if (graph1==TRUE) {
      
      for(i in 1:ncol(x)){
        dev.new()
        ts.plot(beta.mcmc[,i], main=paste("Complete chain for beta",i), xlab="number of iterations", ylab=paste("parameter beta",i))
      }
      
      for(i in 1:ncol(z)){
        dev.new()
        ts.plot(gamma.mcmc[,i], main=paste("Complete chain for gamma",i), xlab="number of iterations", ylab=paste("parameter gamma",i))
        
      }
      
    } else{
    }
    
    
    if (graph2==TRUE) {
      
      for(i in 1:ncol(x)){
        dev.new()
        ts.plot(beta.mcmc.auto[,i], main=paste("Burn chain for beta",i), xlab="number of iterations", ylab=paste("parameter beta",i))
        
      }
      
      for(i in 1:ncol(z)){
        dev.new()
        ts.plot(gamma.mcmc.auto[,i], main=paste("Burn chain for gamma",i), xlab="number of iterations", ylab=paste("parameter gamma",i))
        
      }
      
    } else{
    }
    
    #Beta y Gamma estimations
    Bestimado <- colMeans(beta.mcmc.auto)
    Gammaest <- colMeans(gamma.mcmc.auto)
    
    #estandar errors of beta and gamma
    DesvBeta <- matrix(apply(beta.mcmc.auto,2,sd))
    DesvGamma <- matrix(apply(gamma.mcmc.auto,2,sd))
    
    #estimate values of the dependent variable
    
    if (model=="NB1"){
      yestimado = exp(x%*%Bestimado)
    } else if (model=="NB2") {
      yestimado = exp(x%*%Bestimado)
    } else {
      yestimado = exp(x%*%Bestimado)/(1+exp(x%*%Bestimado))
    }
    
    #Approximate varaince
    varianza <- exp(z%*%Gammaest)
    
    
    #Residuals 
    residuales = as.matrix(y) - yestimado
    
    #Standardized Weighted Residual 1
    swr1 <- residuales/sqrt(varianza)
    
    #Credibility intervals for beta
    B1 <- matrix(0, ncol(x),1)
    B2 <- matrix(0, ncol(x),1)
    
      
    for(i in 1:ncol(x)){
      B1[i,]<-quantile(beta.mcmc.auto[,i],0.025)
      B2[i,]<-quantile(beta.mcmc.auto[,i],0.975)
      B <- cbind(B1,B2)
    }
    
    # Credibility intervals for gamma
    
    G1 <- matrix(0, ncol(z),1)
    G2 <- matrix(0, ncol(z),1)
    
    for(i in 1:ncol(z)){
      G1[i,]<-quantile(gamma.mcmc.auto[,i],0.025)
      G2[i,]<-quantile(gamma.mcmc.auto[,i],0.975)
      G <- cbind(G1,G2)
    }
    
    aceptbeta <- sum(ind1)/nsim
    aceptgamma <- sum(ind2)/nsim
    
    list(Bestimado=Bestimado,Gammaest=Gammaest,X=x,Z=z,DesvBeta=DesvBeta, DesvGamma=DesvGamma, B=B, G=G, yestimado=yestimado, residuales=residuales, estresiduals=swr1, beta.mcmc=beta.mcmc, gamma.mcmc=gamma.mcmc, beta.mcmc.auto=beta.mcmc.auto, gamma.mcmc.auto=gamma.mcmc.auto, Y = y, aceptbeta=aceptbeta, aceptgamma=aceptgamma, model=model, m=m)
  }      
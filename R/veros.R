veros <- function(y,x,z,betas,gammas,model,m){
  
  if (model=="NB1"){
    
    Etas_1 = x%*%betas
    Etas_2 = z%*%gammas
    mu = exp(Etas_1)
    alpha=exp(Etas_2)
    
    a = sum(alpha*(log(alpha)-log(mu)))
    b = sum(lgamma(y+alpha))
    c = sum(lfactorial(y)+lgamma(alpha)+(y+alpha)*log(1+alpha/mu))
    pot = a+b-c
    vero = exp(pot)
    
  } else if (model=="NB2") {
    
    Etas_1 = x%*%betas
    Etas_2 = z%*%gammas
    mu = exp(Etas_1)
    sigma2=exp(Etas_2)
    
    a = (mu/sigma2)^(mu^2/(sigma2-mu))
    b = gamma(y+(mu^2/(sigma2-mu)))
    c = factorial(y)*gamma(mu^2/(sigma2-mu))*(sigma2/(sigma2-mu))^y
    vero = prod(a*b/c)
    
  } else {
    
    Etas_1 = x%*%betas
    Etas_2 = z%*%gammas
    
    mu=inv.logit(Etas_1)
    phi=exp(Etas_2)
    
    if (is.null(m)){
      m=length(y)
    }
    
    combinatorias <- matrix(0,length(y),1)
    
    for (i in 1:nrow(combinatorias)){
      
      combinatorias[i,] <- choose(m,y[i])
      
    }
    
    b <- beta(y+(mu*phi),(m-y)+(phi*(1-mu)))/beta(mu*phi,phi*(1-mu))
    
    vero = prod(combinatorias*b)
    
  }
  
  vero
}

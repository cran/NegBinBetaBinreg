criteria <-
function(objeto){
  if(is.null(objeto)){
    stop("you have to include a NegBinBetaBinreg object")
  }
  
  y <- objeto$Y
  x <- objeto$X
  z <- objeto$Z
  betas <- objeto$Bestimado
  gammas <- objeto$Gammaest
  m <- objeto$m
  model <- objeto$model
  n <- nrow(x)
  
  residuales <- model$residuals
  variance <- model$variance
  phi <- model$precision
  yestimado <-  model$fitted.values
  
  
  L <- veros(y,x,z,betas,gammas,model,m)
  
  if (model=="NB1"){
    AIC <- 2*(length(betas)+length(gammas))-2*L
    BIC <- log(n)*(length(betas)+length(gammas))-2*L
  } else {
    AIC <- 2*(length(betas)+length(gammas))-2*log(L)
    BIC <- log(n)*(length(betas)+length(gammas))-2*log(L)
    
  }
  

  criteria <- list()
  
  criteria$AIC <- AIC
  criteria$BIC <- BIC
  
  criteria
}

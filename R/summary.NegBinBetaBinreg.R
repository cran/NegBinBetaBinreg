summary.NegBinBetaBinreg <-
function(object, ...){
  
  se <- object$desv
  
  TAB <- cbind( Coefficient = coef(object),
                Desv. = se,
                L.CredIntv = object$interv[,1],
                U.CredIntv = object$interv[,2]
  )
  
  colnames(TAB) <- c("Estimate", "Est.Error", "L.CredIntv",  "U.CredIntv")
  
  
  criteria <- criteria(object$X)
  
  res <- list(call=object$call, coefficients=TAB, AIC=criteria$AIC, BIC=criteria$BIC)  
  
  class(res) <- "summary.NegBinBetaBinreg"
  res  
}

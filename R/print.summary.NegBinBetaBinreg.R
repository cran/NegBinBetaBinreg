print.summary.NegBinBetaBinreg <-
function(x, ...){
  
  cat (" \n            ################################################################
            ###     Negative Binomial or Beta Binomial Regression        ###
            ################################################################ \n")
  
  cat("\n Call: \n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, digits=4)
  
  cat("\n AIC: \n")
  print(x$AIC)
  
  cat("\n BIC: \n")
  print(x$BIC)
  
}

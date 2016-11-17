print.NegBinBetaBinreg <-
function(x, ...)
{
  
  cat (" \n            ################################################################
            ###    Negative Binomial or Beta Binomial Regression         ###
            ################################################################ \n")
  
  cat("\n Call: \n")
  
  print(x$call)
  cat("\n Coefficients: \n")
  print(x$coefficients)
  
}

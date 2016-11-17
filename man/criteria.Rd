\name{criteria}
\alias{criteria}
\title{criteria for comparison the Bayesian Negative Binomial regression models with mean and shape (or variance) regression structures, and  Beta Binomial regression with mean and dispersion regression structures.}
\description{Performs the comparison criterias for the Bayesian Negative Binomial regression models with mean and shape (or variance) regression structures, and  Beta Binomial regression with mean and dispersion regression structures.}
\usage{
criteria(objeto)
}
\arguments{

 \item{objeto}{object of class NegBinBetaBinreg}
}
\details{
This function calculate the information criteria for a Bayesian Negative Binomial regression with mean and shape modeling and mean and variance modeling and Beta Binomial regression with mean and dispersion modeling. 
}
\value{
 \item{AIC}{the AiC criteria}
 \item{BIC}{the BIC criteria}
}
\author{
 Edilberto Cepeda-Cuervo \email{ecepedac@unal.edu.co},
 Maria Victoria Cifuentes-Amado \email{mvcifuentesa@unal.edu.co},
 Margarita Marin \email{mmarinj@unal.edu.co}       
}

\keyword{Bayesian}
\keyword{Metropolis Hastings}
\keyword{Meancovariance modelling}
\keyword{Negative Binomial}
\keyword{Beta Binomial}

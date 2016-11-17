\name{NegBinBetaBinregEst}
\alias{NegBinBetaBinregEst}
\title{Negative Binomial and Beta Binomial regression}
\description{Function to estimate a Negative Binomial regression models with mean and shape (or variance) regression structures, and  Beta Binomial regression with mean and dispersion regression structures.}
\usage{
NegBinBetaBinregEst(y,x,z,nsim,bpri,Bpri,
gpri,Gpri,burn,jump,bini,gini,model,m,ni,graph1,graph2)
}
\arguments{
  \item{y}{object of class matrix or vector, with the dependent variable.}
  \item{x}{object of class matrix, with the variables for modelling the mean.}
  \item{z}{object of class matrix, with the variables for modelling the shape, variance or dispersion.}
  \item{nsim}{a number that indicate the number of iterations.}
  \item{bpri}{a vector with the prior values of beta.}
  \item{Bpri}{a matrix with the prior values of the variance of beta.}
  \item{gpri}{a vector with the prior values of gamma.}
  \item{Gpri}{a matrix with the prior values of the variance of gamma.}
  \item{burn}{a proportion that indicate the number of iterations to be burn at the beginning of the chain.}
  \item{jump}{a number that indicate the distance between samples of the autocorrelated the chain, to be excluded from the final chain.}
  \item{bini}{a vector with the initial values of beta.}
  \item{gini}{a vector with the initial values of gamma.}
  \item{model}{it indicates the model that will be used. By default, is the Beta Binomial model (BB), but it could also be the Negative Binomial with mean and shape (NB1) or the Negative Binomial with mean and variance (NB2).}
  \item{m}{Is positive integer that In the Beta Binomial model indicates the number of trials. By default, is the number of data}
  \item{ni}{Is a vector of positive integer that In the Beta Binomial model indicates the number of trials to each individual. By default, is a vector of m}
  \item{graph1}{if it is TRUE present the graph of the chains without jump and burn.}
  \item{graph2}{if it is TRUE present the graph of the chains with jump and burn.}
}
\details{
The Bayesian Negative Binomial regression allow the joint modelling of mean and shape or variance of a negative binomial distributed variable, as is proposed in Cepeda (2001), with exponential link for the mean and the shape or variance. 
The Bayesian Beta Binomial regression allow the joint modelling of mean and precision of a beta binomial distributed variable, as is proposed in Cepeda (2001), with logit link for the mean and exponential link for the precision.
}
\value{
object of class \code{bayesbetareg} with the following:
  \item{Bestimado}{object of class matrix with the estimated coefficients of beta}
  \item{Gammaest}{object of class matrix with the estimated coefficients of gamma}
  \item{X}{object of class matrix, with the variables for modelling the mean}
  \item{Z}{object of class matrix, with the variables for modelling the shape, variance or dispersion.}
  \item{DesvBeta}{object of class matrix with the estimated desviations of beta}
  \item{DesvGamma}{object of class matrix with the estimated desviations of gamma}
  \item{B}{object of class matrix with the B values of the confidence intervals for beta}
  \item{G}{object of class matrix with the G values of the confidence intervals for gamma}
  \item{yestimado}{object of class matrix with the fitted values of y}
  \item{residuales}{object of class matrix with the residuals of the regression}
  \item{residuales}{object of class matrix with the standardized residuals of the regression}
  \item{beta.mcmc}{object of class matrix with the complete chains for beta}
  \item{gamma.mcmc}{object of class matrix with the complete chains for gamma}
  \item{beta.mcmc.auto}{object of class matrix with the chains for beta after the burned process}
  \item{gamma.mcmc.auto}{object of class matrix with the chains for gamma after the burned process}
\item{aceptbeta}{object of class matrix with the acceptance rate for the betas}
\item{aceptgamma}{object of class matrix with the acceptance rate for the gammas}
}
\references{
1. Cepeda C. E. (2001). Modelagem da variabilidade em modelos lineares generalizados. Unpublished Ph.D. tesis. Instituto de Matematicas. 
Universidade Federal do Rio do Janeiro.
//http://www.docentes.unal.edu.co/ecepedac/docs/MODELAGEM20DA20VARIABILIDADE.pdf.
http://www.bdigital.unal.edu.co/9394/.
2.Cepeda, E. C. and Gamerman D. (2005). Bayesian Methodology for modeling parameters in the two-parameter exponential family. Estadistica 57, 93 105.
//
3.Cepeda, E. and Garrido, L. (2011). Bayesian beta regression models: joint mean and precision modeling. Universidad Nacional
//
4.Cepeda, E. and Migon, H. and Garrido, L. and Achcar, J. (2012) Generalized Linear models with random effects in the two parameter exponential family. 
Journal of Statistical Computation and Simulation. 1, 1 13.
//
5.Cepeda-Cuervo, E. and  Cifuentes-Amado, V. (2016) Double generalized beta-binomial and negative binomial regression. To appear.
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


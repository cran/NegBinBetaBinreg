\name{veros}
\alias{veros}
\title{Likelihood}
\description{calculate the likelihood value for the Negative Binomial regression models with mean and shape (or variance) regression structures, and  Beta Binomial regression with mean and dispersion regression structures.}
\usage{
veros(y,x,z,betas,gammas,model,m)
}
\arguments{
  \item{y}{object of class matrix, with the dependent variable}
  \item{x}{object of class matrix, with the variables for modelling the mean}
  \item{z}{object of class matrix, with the variables for modelling the variance}
  \item{betas}{a vector with the previous proposal beta parameters}
  \item{gammas}{a vector with the previous proposal gamma parameters}
  \item{model}{it indicates the model that will be used. By default, is the Beta Binomial model (BB), but it could also be the Negative Binomial with mean and shape (NB1) or the Negative Binomial with mean and variance (NB2).}
  \item{m}{It is positive integer that In the Beta Binomial model indicates the number of trials. By default, is the number of data}
}
\details{
calculate the likelihood value for the Negative Binomial regression with mean and shape modeling and mean and variance modeling and Beta Binomial regression with mean and dispersion modeling. 
}
\value{
\item{value}{a integer with the likelihood}
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


\name{muproposal}
\alias{muproposal}
\title{A proposal for beta parameter}
\description{Propose a value for the beta parameter}
\usage{
muproposal(y, x, z, betas.ini,gammas.ini,bpri,Bpri,model,m,ni)
}
\arguments{
  \item{y}{object of class matrix or vector, with the dependent variable.}
  \item{x}{object of class matrix, with the variables for modelling the mean.}
  \item{z}{object of class matrix, with the variables for modelling the shape, variance or dispersion.}
  \item{betas.ini}{a vector with the beta that define the old p.d.f}
  \item{gammas.ini}{a vector with the gamma that define the old p.d.f}
  \item{bpri}{a vector with the prior values of beta.}
  \item{Bpri}{a matrix with the prior values of the variance of beta.}
  \item{model}{it indicates the model that will be used. By default, is the Beta Binomial model (BB), but it could also be the Negative Binomial with mean and shape (NB1) or the Negative Binomial with mean and variance (NB2).}
  \item{m}{It is positive integer that In the Beta Binomial model indicates the number of trials. By default, is the number of data}
  \item{ni}{It is a vector of positive integer that In the Beta Binomial model indicates the number of trials to each individual. By default, is a vector of m}
}
\details{
Generate a proposal for the beta parameter according to the model proposed by Cepeda(2001) and Cepeda and Gamerman(2005). 
}
\value{
 \item{value}{a matrix with the proposal for beta}
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

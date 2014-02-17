\name{Derv1}
\alias{Derv1}
\title{Calculating the first derivative of the pencopula likelihood function w.r.t. parameter b}
\description{
 Calculating the first derivative of the pencopula likelihood function w.r.t. parameter b. 
}
\usage{
Derv1(penden.env)
}

\arguments{
	\item{penden.env}{Containing all information, environment of pencopula().}
}
\details{
The calculation of the first derivative of the pencopula likelihood function w.r.t. b equals
\deqn{s(b,\lambda)= {\partial l(b,\lambda)}/{\partial b}= \sum_{i=1}^n
  \Phi(u_i)/c(u_i,b) - P(\lambda)b}{%
  \eqn{s(b,lambda)=}}
  with \deqn{P(\lambda)}{%
    \eqn{P(lambda)}}
  is the penalty matrix, saved in the environment.
}

\value{
  \item{Derv1.pen}{first order derivation of the penalized likelihood
    function w.r.t. parameter b.}
  Derv1.pen is saved in the environment.
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2013), Scandinavian Journal of Statistics 40(4), 685-705.}
\keyword{math}


\name{penalty.matrix}
\alias{penalty.matrix}

\title{Calculating the penalty matrix P(lambda)}

\description{
Calculating the penalty matrix P depends on the number of
covariates 'p', the order of differences to be penalized 'pen.order',
the number of observations 'n' and the penalty parameters 'lambda''.
}
\usage{
penalty.matrix(penden.env, temp = FALSE)
}

\arguments{
  \item{penden.env}{Containing all information, environment of
    pencopula().}
  \item{temp}{If TRUE, the iteration for a new 'b' is not finished and a temporary penalty matrix is calculated, default = FALSE.}
}

\details{The penalty matrix is calculated as

  \eqn{P(\lambda)=\sum_{j=1}^{p} \lambda_j P_j}

  with

  \eqn{P_j=\left(\bigotimes_{l=1}^{j-1}{I}\right) \otimes \{({A^{-1})^T} P {A}^{-1}\} \otimes \left(\bigotimes_{l=j+1}^p {I}\right)}

  The needed values are calculated or saved in the environment 'penden.env'.

  There is no penalty used for the Bernstein polynomial basis, P=0.
}

\value{
  \item{DDD.sum}{Penalty matrix P}
  Matrix is saved in the environment.
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{nonparametric}

\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2013), Scandinavian Journal of Statistics 40(4), 685-705.}

\name{pencopula-package}
\alias{pencopula-package}
\docType{package}
\title{
The package 'pencopula' offers routines for estimating multivariate penalized copula densities and copula distribution.
}
\description{
The package 'pencopula' offers routines for estimating multivariate penalized copula densities and copula distribution. For details see the description of the function pencopula().
}
\details{
\tabular{ll}{
Package: \tab pencopula\cr
Type: \tab Package\cr
Version: \tab 0.3.4\cr
Date: \tab 2012-08-15\cr
License: GPL (>= 2)
LazyLoad: \tab yes\cr
}
The packages contributes the function 'pencopula()' for estimating copula densities and copula distributions using penalized splines techniques. 
}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2011), to appear.}
\keyword{nonparametric}
\examples{
Y <- cbind(runif(200),runif(200))
plot(pencopula(Y,d=4,D=4,lambda=rep(10,2),base="B-spline"))
}
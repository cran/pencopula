\name{Derv2}
\alias{Derv2}
\title{Calculating the second order derivative with and without penalty.}
\description{
Calculating the second order derivative of the likelihood function of
the pencopula approach w.r.t. the parameter b. Thereby, for later use,
the program calculates the second order derivative with and without the
penalty. Moreover, Derv2 seperates the calculation for temporary weights
b in iteration and final weights b.
}
\usage{
Derv2(penden.env, temp = FALSE)
}

\details{
  We approximate the second order derivative in this approach with the negative fisher information. 
}

\arguments{
  \item{penden.env}{Containing all information, environment of pendensity()}
  \item{temp}{smoothing parameter lambda}
}

\value{
  \item{Derv2.pen}{second order derivative w.r.t. beta with penalty}
  \item{Derv2.cal}{second order derivative w.r.t. beta without
    penalty. Needed for calculating of e.g. AIC.}
  Derv2.cal and Derv2.pen are saved in the environment.
}
\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2011), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}


\keyword{math}


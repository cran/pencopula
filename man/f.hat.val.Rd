\name{f.hat.val}
\alias{f.hat.val}
\title{Calculating the actual fitted values 'f.hat.val' of the estimated density function}

\description{
Calculating the actual fitted values of the response, depending on the actual parameter set b
}
\usage{
f.hat.val(penden.env,cal=FALSE,temp=FALSE)
}
\arguments{
  	\item{penden.env}{Containing all information, environment of pencopula()}
	\item{cal}{if TRUE, the final weights of one iteration are used for
    the calculation of the fitted values.}
	\item{temp}{if TRUE, the iteration for optimal weights is still in
    progress and the temporary weights are used for calculation of the
    fitted values.}
}
\details{
Calculating the actual fitted values of the response, depending on the actual parameter set b. Multiplying the actual set of parameters \eqn{b} with the base 'base.den' delivers the fitted values.
}
\value{
  \item{f.hat.val}{Fitted values for the current coefficents}.
  f.hat.val is saved in the environment.
}
\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2011), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\keyword{nonparametric}

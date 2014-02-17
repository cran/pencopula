\name{my.IC}
\alias{my.IC}

\title{Calculating the AIC- and BIC-value}

\description{
Calculating the AIC- and BIC- value of the copula density estimation. Therefore, we add the unpenalized log likelihood of the estimation and the degree of freedom, which are 	
}
\usage{
my.IC(penden.env)
}

\arguments{
  \item{penden.env}{Containing all information, environment of pencopula()}
}

\details{
AIC is calculated as
\eqn{AIC(\lambda)= - 2*l({\bf u},\hat{\bf{b}}) + 2*df(\lambda)}

BIC is calculated as
\eqn{BIC(\lambda)= 2*l({\bf u},\hat{\bf{b}}) + 2*df(\lambda)*log(n)}
}

\value{
  \item{AIC}{sum of twice the negative non-penalized log likelihood and mytrace}
  \item{trace}{calculated mytrace as the sum of the diagonal matrix
    df, which results as the product of the inverse of the penalized
    second order derivative of the log likelihood with the non-penalized
    second order derivative of the log likelihood}
  \item{BIC}{sum of twice the non-penalized log likelihood and log(n)}
  All values are saved in the environment.
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\keyword{math}

\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2013), Scandinavian Journal of Statistics 40(4), 685-705.}

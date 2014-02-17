\name{pen.log.like}
\alias{pen.log.like}
\title{Calculating the log likelihood}

\description{
Calculating the considered log likelihood.
}

\usage{
pen.log.like(penden.env, cal=FALSE, temp=FALSE)
}

\arguments{
  \item{penden.env}{Containing all information, environment of pencopula()}
  \item{cal}{if TRUE, the final weights of one iteration are used for
    the calculation of the penalized log likelihood.}
  \item{temp}{if TRUE, the iteration for optimal weights is still in
    progress and the temporary weights are used for calculation.}
}
\details{
The calculation depends on the estimated weights b, the penalized
hierarchical B-splines Phi and the penalty paramters lambda.\cr
\deqn{l(b,\lambda)=\sum_{i=1}^{n} \left[ \log \{\sum_{i=1}^n
  \boldsymbol\Phi(u_i)\} b\right]- \frac 12 b^T \boldsymbol{P}(\lambda) b}{%
  \eqn{l(beta,lambda)=sum(log(\Phi(u_i)b))-0.5*b^T \tilde{P}(\lambda) b}}
with
\deqn{\boldsymbol{P}(\lambda)=\sum_{j=1}{p}\lambda_j\boldsymbol{P}_j}

The needed values are saved in the environment.
}

\value{
  \item{pen.log.like}{Penalized log likelihood of the copula density.}
  \item{log.like}{Log-Likelihood of the copula density.}
The values are saved in the environment.
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{nonparametric}

\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2013), Scandinavian Journal of Statistics 40(4), 685-705.}

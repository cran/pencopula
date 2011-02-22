\name{start.valgrid}
\alias{start.valgrid}
\alias{grid.points}
\title{Calculating the start values 'b' for the first iteration of the
  quadratic program.}

\description{Calculating the start values 'b' for the first iteration of the
  quadratic program. Moreover, the grid of values for the side condition
  c>=0 of the quadratic program are calculated. If "adapt.grid", the
  number of grid points is reduced to speed up the quadratic program.
}

\usage{
start.valgrid(penden.env)
}

\arguments{
  \item{penden.env}{Containing all information, environment of pencopula()}
}

\details{The grid of values for the side conditions of the quadratic
  program c>=0 is constructed as the tensor product of all knots. If $p$
  and $d$ increase, the number of conditions and computational time of
  the quadratic programm increase enormously, e.g. a full tensor product
  $u$ for $p=4$ and $d=4$ contains 83521 entries. If the data
  $u$ is not high correlated, i.e. the data is not from a extreme
  value copula like a Clayton copula, one can reduce the full tensor
  product. In 'pencopula' one can choose the option
  'adapt.grid' which effects the following and may reduce the
  calculating time without any loss of accuracy. One can omit points in
  $u$ in sections of $[0,1]^p$ which are in the neighbourhood
  of many observations in $u$, because the data itself induces a
  positive density in these areas by construction. Therefore, we
  calculate the minimal $p$-dimensional euclidean distance $e_i$
  of each $u_i, i=1,\dots,(2^d+1)^p$ to the data $u$ and omit the points
  corresponding to the first quartile of minimal euclidean distance
  $e_i$ in $u$, we call this new set of points
  $u_{min}$. This amout of points is used in the first
  iteration step to estimate weights $b$ corresponding to a copula
  density.}


\value{
  \item{X.knots.g}{If adapt.grid=TRUE, set of reduced grid values, in which the side condition of
    the quadratic program c(u,b)>=0 will be postulated.}
  \item{X.knots.g.all}{Set of all grid values, in which the side condition of
    the quadratic program c(u,b)>=0 will be postulated.}
  The values are saved in the environment.
}

\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2011), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{nonparametric}

\name{my.bspline}
\alias{my.bspline}
\alias{int.my.bspline}
\title{my.bspline}

\description{
'my.bspline' Integrates the normal B-Spline basis to a B-spline density
basis. The dimension of the basis depends on the input of number of
knots 'k' and of the order of the B-spline basis 'q'. 'int.my.bspline'
is a function for transformation of open B-spline basis at the
boundary to become a B-spline basis density.
}
\usage{
my.bspline(h, q, knots, y, K, plot.bsp, typ)
int.my.bspline(help.env)
}

\arguments{
  \item{h}{if equidistant knots are used (default in pencopula()), h is the distance between two neighbouring knots}
  \item{q}{selected order of the B-spline basis}
  \item{knots}{selected values for the knots}
  \item{y}{values of the response variable}
  \item{K}{the number of knots for the construction of the base}
  \item{plot.bsp}{Indicator variable TRUE/FALSE if the integrated
    B-spline basis should be plotted}
  \item{typ}{typ==1 without open B-splines at the boundary
    typ==2 with open B-splines at the boundary}
  \item{help.env}{Internal environment of my.bspline().}
}

\details{
  Firstly, the function constructs the B-spline basis to the given number
  of knots 'K' and the given locations of the knots. Due
  to the recursive construction of the B-Spline, for all orders greater
  than 2, the dimension of the B-spline basis of given K grows up with
  help.degree=q-2.
  There are two typs of B-spline basis possible. First, a B-spline basis
  without open B-splines at the boundary (typ==1) and a regulat B-spline
  basis with open B-splines at the boundary (typ==2).
  Both typs are integrated to become B-spline density basis. To
  integrate a basis of typ 1 we use the well-known factor
  'q/(knots.val[i+q]-knots.val[i])'. For typ==2 we determine
  functions analytically for the integration.
  Moreover, one can draw the integrated basis and, if one calls this function with the argument 'plot.bsp=TRUE'.
}
\value{
  \item{base.den}{The integrated B-Spline base of order q}
  \item{stand.num}{The coefficients for standardization of the ordinary B-Spline basis}
  \item{knots.val}{This return is a list. It consider of the used knots 'knots.val\$val', the help knots 'knots.val\$help' and the additional knots 'knots.val\$all', used for the construction of the base and the calculation of the distribution function of each B-Spline.}
  \item{K}{The transformed value of K, due to used order 'q' and the input of 'K'}
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\note{
This functions uses the fda-package to build the B-Spline density basis.
}
\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2013), Scandinavian Journal of Statistics 40(4), 685-705.}

\keyword{math}


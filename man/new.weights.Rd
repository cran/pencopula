\name{new.weights}
\alias{new.weights}
\title{Calculating new weights b.}

\description{
Calculating new weights b using quadratic programing.
}

\usage{
new.weights(penden.env)
}

\arguments{
  \item{penden.env}{Containing all information, environment of pencopula()}
}

\details{
  The new weights are calculated solving a quadratic program. Therefore,
  the derivates of first and second order are needed, 'Derv1.pen' and
  'Derv2.pen'. Moreover, we have to fulfill the side conditions
  c(u,b)>=0, sum(b)=1 and that the marginal densities are uniform. All side
  conditions are saved as 'AA.help' in the environment. There exists two
  different algorithms, due to the fact if one works with or without an
  adaptive grid for the side condition c(u,b)>=0.

  If the quadratic program does not find a new feasible solution, the whole
  program terminates. For solving the quadratic program, we use the
  function 'solve.QP' from the  R-package 'quadprog'.

  Using a Bernstein polynomial basis, the side conditions are sum(b)=1
  ,b>=0 and that the marginal densities are uniform.
}

\note{For high dimensional data (p>=4) and high hierachy level 'd' of the
  marginal hierarchical B-spline, the computational time for solving a
  quadratic program increases.
  }

\value{
  \item{ck.val.temp}{Calculated new values for the weights 'b'. The add
    on 'temp' means, that there is a check in the next step if the
    weights 'b' have been converted or not. If converted, the new values
    'ck.val.temp' are unnoted. If not converted, 'ck.val.temp' become
    the ordinary 'ck.val' for the next iteration. This check is done in my.loop.}
  'ck.val.temp' is saved in the environment.
}

\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2011), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{nonparametric}
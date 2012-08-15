\name{pencopula}
\alias{pencopula}

\title{Calculating penalized copula density with penalized hierarchical B-splines}
\description{
  Main program for estimation of copula densities with penalized
  hierarchical B-splines. The estimation can be done for multivariate
  datasets. The response is called 'y', the covariates 'x'. We estimate
  densities using penalized hierarchical B-splines. This is done by
  choosing the polynomial degree of the univariate B-spline density
  basis, that is built upon 2^d+1 equidistant knots. 'd' refers to the
  hierachy level of the marginal hierarchical  B-spline and 'D' to the
  maximum hierarchical level of the hierarchical B-spline basis. We penalize the m-order differences of the coefficients 'b' to estimate new weights 'b'.
}
\usage{
pencopula(data, d=3, D=d, q=1, base = "B-spline", max.iter = 20,
plot.bsp = FALSE, lambda = NULL, pen.order=2, adapt.grid = FALSE, add =
TRUE, alpha = 0, symmetric= TRUE, data.frame=parent.frame())
}

\arguments{
  \item{data}{'data' contains the data. 'data' has to be a matrix or a
    data.frame. The number of columns of 'data' is p.}
  \item{d}{refers to the hierachy level of the marginal hierarchical
    B-spline, default is d=3.}
  \item{D}{referes to the maximum hierachy level, default is D=3. If
    D<d, it follows D<-d.}
  \item{q}{degree of the marginal hierarchical B-spline.}
  \item{base}{By default, the used marginal basis is a
    'B-spline'. Second possible option is 'Bernstein', using a Bernstein
    polynomial basis.}
  \item{max.iter}{maximum number of iteration, the default is max.iter=20.}
  \item{plot.bsp}{TRUE or FALSE, if TRUE the used B-Spline base should
    be plotted.}
  \item{lambda}{p-dimensional vector of penalty parameters, the values
    can be different. Default is rep(10000,p).}
  \item{pen.order}{The order of differences for the penalization,
    default is pen.order=2.}
  \item{adapt.grid}{Default = FALSE, if TRUE the used grid for the side
    condition c(u,b)>=0 is reduced.}
  \item{add}{Default = TRUE. Due to numerical rounding errors, some
    results in the quadratic programming are misleading. So we add a
    small epsilon in the side conditions c(u,b) >=0 to the actual weights 'b' in each iteration.}
  \item{alpha}{Default = 0, that results in equidistant knots. If alpha
    !=0 the knots are moved.}
  \item{symmetric}{Default = TRUE. Tuning parameter for the location of
    the knots if alpha != 0.}
  \item{data.frame}{reference to the data. Default reference is the parent.frame().}
}
\details{
  The estimation of the copula density is done in several steps.
  1. Preparations: Calculating the number of marginal knots 'ddb' and the
  number of spline coefficients 'DD' , depending on 'd' and
  'D'.'Index.Basis.D' is created, that is an index set,
  which marginal basis will belong to the hierarchical B-spline density
  basis 'tilde.Psi.d.D'.
  2. Creating marginal spline coefficients 'T.marg', depending on
  'Index.basis.D'.
  3. Building the hierarchical B-spline density basis 'tilde.PSI.d.D'
  4. Calculating the knots for the estimation, by default the knots are
  equidistant. Transformations of the knots are possible, see
  'knots.start'. Moreover a B-spline density basis 'tilde.Psi.d.knots'
  is created, with the knots in the hierachical order
  'knots.t'. 'tilde.Psi.d.knots' is needed for the restriction, that
  each marginal density has to be uniform.
  5. Depending on the knots, start values for the weights 'b' are
  calculated. This is done in 'startval.grid'. Moreover 'startval.grid'
  determines an grid for the side condition c(u,b)>=0 of the quadratic
  program.
  6. Each marginal basis is estimated to be uniform. Therefore a matrix
  'A.Restrict' is calculated, such that \deqn{A.Restrict \%*\% b = 1}{%
    \eqn{A.Restrict \%*\% b = 1}}
  for all p-covariates.
  7. The penalty matrix P is created, see 'penalty.matrix'.
  8. The first calculation of coefficients 'b' is done.
  9. 'my.loop' iterates the calculation of the optimal weights 'b' until
  some convergence, see 'my.loop'. The maximal number of iterations are
  limited, default is max.iter=20.
  10. If 'my.loop' terminates, the information criteria are calculated,
  see 'my.IC'.
}

\value{
Returning an object of class pencopula. The class pencopula consists of
the environment 'penden.env', which includes all calculated values of
the estimation approach.
For a fast overview of the main results, one can use the function 'print.pencopula()'.
}

\note{If the estimation with the current setup does not work, e.g. the
  quadratic program can not be solved, the routine returns log.like=0, pen.log.like=0,
  AIC=0 and BIC=0. Please restart the estimation, e.g. with a different
  penalty parameter.}

\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2011), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\examples{
library(MASS)
data(Lufthansa)
data(DeutscheBank)

data.dbank <- data.lufth <- c()

dim.data <- dim(DeutscheBank)[1]

for(i in 2:dim.data) data.dbank[i-1] <- log(DeutscheBank[i,2]/DeutscheBank[i-1,2])
for(i in 2:dim.data) data.lufth[i-1] <- log(Lufthansa[i,2]/Lufthansa[i-1,2])

dbank1 <- fitdistr(data.dbank,"t")
lufth1 <- fitdistr(data.lufth,"t")

mypt <- function(x, m, s, df) pt((x-m)/s, df)

Y <- cbind(mypt(data.dbank,dbank1$estimate[1],s=dbank1$estimate[2],df=dbank1$estimate[3]),
mypt(data.lufth,lufth1$estimate[1],s=lufth1$estimate[2],df=lufth1$estimate[3]))

cop <- pencopula(Y,d=4,D=4,lambda=rep(10,2))

}

\keyword{nonparametric}
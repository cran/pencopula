\name{plot.pencopula}
\alias{plot.pencopula}

\title{Plot the estimated copula density or copula distribution.}

\description{The function plots the estimated copula density or the
  copula distrubtion, using the R-package 'lattice'.}

\usage{
plot.pencopula(x, val = NULL, marg = TRUE, plot = TRUE, int = FALSE,
main.txt = NULL, sub.txt = NULL, contour = FALSE, cond = NULL, cuts =
20, cex = 1, cex.axes = 1, xlab = NULL, ylab = NULL, zlab=NULL,...)
}

\arguments{
  \item{x}{object of class 'pencopula'.}
  \item{val}{Default val = NULL, one can calculate the estimated density
    in for p-dimensional vector, e.g. val=c(0.5,1) for the two
    dimensional case.}
  \item{marg}{Default = TRUE, plotting the marginal densities.}
  \item{plot}{Default = TRUE, if 'FALSE' no plot is shown, e.g. for
    calculations with val != NULL.}
  \item{int}{Default = FALSE, if TRUE, the integral, i.e. the distribution of
    the copula density is plotted.}
  \item{main.txt}{Default = NULL shows 'd', 'D', the values of lambda,
    the penalty order and the degree of the B-splines.}
  \item{sub.txt}{Default = NULL shows the log-likelihood, the penalized
    log-likelihood and the AIC-value of the estimation.}
  \item{contour}{If TRUE, a contour plot is shown. Default = FALSE.}
  \item{cond}{Default = NULL, if the dimension of data 'p' is higher
    than 2, one can plot a two-dimensional conditional plot. The user
    specifies p-2 values for the plot, indicating with '-1'. So for a
    three-dimensional plot, cond=c(0,-1,-1) shows the
    density/distribution  ith fixed first covariate and the second
    and third covariates vary.}
  \item{cuts}{Number of cuts for the contour plots, if contour=TRUE.}
  \item{cex}{Default = 1, determing the size of the main of the plot.}
  \item{cex.axes}{Default = 1, determing the size of the labels at the axes.}
  \item{xlab}{Default = NULL and no text is printed at the xlab}
  \item{ylab}{Default = NULL and no text is printed at the ylab}
  \item{zlab}{Default = NULL and 'density' is printed at the zlab for
    int=FALSE and 'distribution' for int=TRUE.}
  \item{...}{further arguments}
}

\details{
  For the two dimensional plots, a equidistant grid of 21 values between
  0 and 1 is constructed. The plot consists of the density or
  distribution values in this grid points. For plots of high dimensional
  data (p>2), one has to fix p-2 covariates (see 'cond').
}

\value{If 'val' is not NULL, the function returns a matrix with the
  calculated density or distribution values for the set 'val'.
}

\references{Flexible Copula Density Estimation with Penalized
  Hierarchical B-Splines, Kauermann G., Schellhase C. and Ruppert, D. (2011), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{plot}


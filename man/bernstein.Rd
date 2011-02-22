\name{bernstein}
\alias{bernstein}
\title{Calculating a bernstein polynomial.}
\description{
 Calculating a bernstein polynomial. 
}
\usage{
bernstein(v,x,n)
}

\arguments{
	\item{v}{number of the current basis element, e.g. v is in the interval [0,2^K]}
        \item{x}{data points for constructing of the bernstein polynomial basis.}
        \item{n}{number of polynomials in the bernstein polynomial basis, e.g n = 2^K}    
}

\value{
  The bernstein polynomial basis is created using the function 'apply'
  for some data 'x'.

  K <- 3
  index.b <- matrix(0:2^K)
  ## Bernstein polynomial 
  B <- apply(index.b,1,bernstein,x,n=2^K)
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{math}

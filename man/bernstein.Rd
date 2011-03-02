\name{bernstein}
\alias{bernstein}
\alias{int.bernstein}
\title{Calculating a bernstein polynomial.}
\description{
 Calculating a bernstein polynomial, transformed to be a
 density. 'int.bernstein' calculates the integral of the bernstein
 density polynomial.
}
\usage{
bernstein(v,x,n)
int.bernstein(x,n)
}

\arguments{
	\item{v}{number of the current basis element, e.g. v is in the interval [0,2^K]}
        \item{x}{data points for constructing of the bernstein polynomial basis.}
        \item{n}{number of polynomials in the bernstein polynomial basis, e.g n = 2^K}    
}

\value{
  The bernstein polynomial density basis is created using the function 'apply'
  for some data 'x'.

  K <- 3
  index.b <- matrix(0:2^K)
  ## Bernstein polynomial 
  B <- apply(index.b,1,bernstein,x,n=2^K)

  The integral of the Basis B of degree n is:
  B <- int.bernstein(x,n=2^K)
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{math}

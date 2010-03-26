\name{Reciprocal Trigonometric Functions}
\alias{trig}
\alias{csc}
\alias{sec}
\alias{cot}
\title{Reciprocal Trigonometric Functions}
\description{
  Calculates the reciprocal trigonometric functions: cosecant (csc),
  secant (sec), and cotangent (cot).
}
\usage{
csc(x)
sec(x)
cot(x)
}
\arguments{
  \item{x}{the value of the angle in radians.}
}
\value{
  Reciprocal trigonometric functions:
  \tabular{l}{
    \code{csc(x) = 1 / sin(x)} \cr
    \code{sec(x) = 1 / cos(x)} \cr
    \code{cot(x) = 1 / tan(x)}
  }
}
\seealso{
  \code{\link{Rx.calc}}, \code{\link{Rrup.calc}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\keyword{internal}

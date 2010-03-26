\name{CB08 Period-Dependent Model Coefficients}
\alias{coefs.t.cb}
\alias{c0.cb}
\alias{c1.cb}
\alias{c2.cb}
\alias{c3.cb}
\alias{c4.cb}
\alias{c5.cb}
\alias{c6.cb}
\alias{c7.cb}
\alias{c8.cb}
\alias{c9.cb}
\alias{c10.cb}
\alias{c11.cb}
\alias{c12.cb}
\alias{k1.cb}
\alias{k2.cb}
\alias{k3.cb}
\alias{Sigma.LnY.cb}
\alias{Tau.LnY.cb}
\alias{Sigma.C.cb}
\alias{rho.cb}
\title{CB08 Period-Dependent Model Coefficients}
\description{
  Functions for computing the period-dependent coefficients of the
  Campbell & Bozorgnia NGA model.
}
\usage{
c0.cb(T)
c1.cb(T)
c2.cb(T)
c3.cb(T)
c4.cb(T)
c5.cb(T)
c6.cb(T)
c7.cb(T)
c8.cb(T)
c9.cb(T)
c10.cb(T)
c11.cb(T)
c12.cb(T)
k1.cb(T)
k2.cb(T)
k3.cb(T)
Sigma.LnY.cb(T)
Tau.LnY.cb(T)
Sigma.C.cb(T)
rho.cb(T)
}
\arguments{
  \item{T}{Spectral period (sec).  Use 0 for PGA, -1 for PGV, and -2 for
  PGD.}
}
\details{
  \describe{
    \item{Coefficients for median ground motion terms (Table 2 in
      Campbell and Bozorgnia (2008)):}{\code{c0.cb, c1.cb, c2.cb, c3.cb,
	c4.cb, c5.cb, c6.cb, c7.cb, c8.cb, c9.cb, c10.cb, c11.cb,
	c12.cb, k1.cb, k2.cb,} and \code{k3.cb}.}
    \item{Coefficients for standard deviation terms (Table
      3):}{\code{Sigma.LnY.cb, Tau.LnY.cb, Sigma.C.cb,} and
      \code{rho.cb}.}
  }
}
\value{
  Coefficients for the Campbell & Bozorgnia NGA model.
}
\seealso{
  \code{\link{coefs}}, \code{\link{periods.cb}}, \code{\link{subs.cb}},
  \code{\link{Sa.cb}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Campbell, K. W., and Y. Bozorgnia (2008). NGA Ground Motion Model for
  the Geometric Mean Horizontal Component of PGA, PGV, PGD, and 5\%
  Damped Linear Elastic Response Spectra for Periods Ranging from 0.01
  to 10 s. \emph{Earthquake Spectra} \bold{24,} 139--171.
  Boore DM, Atkinson, GM (2008).
}
\keyword{internal}

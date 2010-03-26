\name{AS08 Period-Dependent Model Coefficients}
\alias{coefs.t.as}
\alias{Vlin.as}
\alias{b.as}
\alias{a1.as}
\alias{a2.as}
\alias{a3.as}
\alias{a4.as}
\alias{a5.as}
\alias{a8.as}
\alias{a10.as}
\alias{a12.as}
\alias{a13.as}
\alias{a14.as}
\alias{a15.as}
\alias{a16.as}
\alias{a18.as}
\alias{s1.est.as}
\alias{s2.est.as}
\alias{s1.meas.as}
\alias{s2.meas.as}
\alias{s3.as}
\alias{s4.as}
\alias{rho.as}
\title{AS08 Period-Dependent Model Coefficients}
\description{
  Functions for computing the period-dependent coefficients of the
  Abrahamson & Silva NGA model.
}
\usage{
Vlin.as(T)
b.as(T)
a1.as(T)
a2.as(T)
a3.as(T)
a4.as(T)
a5.as(T)
a8.as(T)
a10.as(T)
a12.as(T)
a13.as(T)
a14.as(T)
a15.as(T)
a16.as(T)
a18.as(T)
s1.est.as(T)
s2.est.as(T)
s1.meas.as(T)
s2.meas.as(T)
s3.as(T)
s4.as(T)
rho.as(T)
}
\arguments{
  \item{T}{Spectral period (sec).  Use 0 for PGA and -1 for PGV.}
}
\details{
  \describe{
    \item{Coefficients for median ground motion terms
      (Tables 5a and 5b of Abrahamson and Silva (2008)):}{\code{Vlin.as,
	b.as, a1.as, a2.as, a3.as, a4.as, a5.as, a8.as, a10.as, a12.as,
	a13.as, a14.as, a15.as, a16.as,} and \code{a18.as}.}
    \item{Coefficients for standard deviation terms
      (Table 6):}{\code{s1.est.as, s2.est.as} (estimated / inferred
      Vs30)\code{; s1.meas.as, s2.meas.as} (measured Vs30)\code{;
	s3.as; s4.as;} and \code{rho.as}.}
  }
}
\value{
  Coefficients for the Abrahamson & Silva NGA model.
}
\seealso{
  \code{\link{coefs}}, \code{\link{periods.as}}, \code{\link{subs.as}},
  \code{\link{Sa.as}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Abrahamson, N., and W. Silva (2008). Summary of the Abrahamson & Silva
  NGA Ground-Motion Relations. \emph{Earthquake Spectra} \bold{24,}
  67--97.
}
\keyword{internal}

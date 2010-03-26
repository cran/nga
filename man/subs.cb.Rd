\name{CB08 Sub-Functions}
\alias{subs.cb}
\alias{Fmag.cb}
\alias{Fdist.cb}
\alias{Fflt.cb}
\alias{Fhng.cb}
\alias{Fsite.cb}
\alias{Fsed.cb}
\alias{Alpha.cb}
\alias{Sigma.cb}
\alias{Tau.cb}
\alias{SigmaTot.cb}
\alias{SigmaArb.cb}
\alias{A1100.cb}
\alias{SaMedian.cb}
\title{CB08 Sub-Functions}
\description{
  The sub-functions needed to estimate ground motions using the
  Campbell & Bozorgnia NGA model.  These functions are called by the main
  function used for calculations, \code{\link{Sa.cb}}.
}
\usage{
Fmag.cb(M, T)
Fdist.cb(M, Rrup, T)
Fflt.cb(Ztor, Frv, Fnm, T)
Fhng.cb(M, Rrup, Rjb, Ztor, dip, T)
Fsite.cb(Vs30, A1100, T)
Fsed.cb(Z2.5, T)
Alpha.cb(A1100, Vs30, T)
Sigma.cb(A1100, Vs30, T)
Tau.cb(T)
SigmaTot.cb(A1100, Vs30, T)
SigmaArb.cb(A1100, Vs30, T)
A1100.cb(M, Rrup, Rjb, Ztor, Frv, Fnm, dip, Z2.5)
SaMedian.cb(M, Rjb, Rrup, Ztor, Frv, Fnm, dip, Vs30, Z2.5, T)
}
\arguments{
  \item{M}{Moment magnitude of earthquake.}
  \item{Rjb}{Horizontal distance to the surface projection of the
    rupture plane; Joyner-Boore distance (km).}
  \item{Rrup}{Closest distance to the rupture plane; rupture distance
    (km).}
  \item{Ztor}{Depth to top of rupture (km).}
  \item{Frv}{Reverse style-of-faulting flag; equal to 1 for reverse
    faulting (30 <= rake <= 150 deg), and 0 otherwise.}
  \item{Fnm}{Normal style-of-faulting flag; equal to 1 for normal
    faulting (-150 <= rake <= -30 deg), and 0 otherwise.}
  \item{dip}{Fault dip angle (deg).}
  \item{Vs30}{Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s).}
  \item{Z2.5}{Depth to Vs = 2.5 km/s (km).}
  \item{A1100}{Median PGA when Vs30 = 1100 m/s (g).}
  \item{T}{Spectral period (sec). Use 0 for PGA, -1 for PGV, and -2 for PGD.}
}
\details{
  \tabular{ll}{
    \code{Fmag.cb} \tab Magnitude term (Eqn 2 of Campbell and Bozorgnia (2008)).\cr
    \code{Fdist.cb} \tab Distance term (Eqn 3).\cr
    \code{Fflt.cb} \tab Fault mechanism term (Eqns 4 and 5).\cr
    \code{Fhng.cb} \tab Hanging wall term (Eqns 6 to 10).\cr
    \code{Fsite.cb} \tab Shallow site response term (Eqn 11).\cr
    \code{Fsed.cb} \tab Basin response term (Eqn 12).\cr
    \code{Alpha.cb} \tab Partial derivative of site response function (\code{Fsite.cb})
    with respect to A1100 (Eqn 17).\cr
    \code{Sigma.cb} \tab Intra-event standard deviation (Eqn 15).\cr
    \code{Tau.cb} \tab Inter-event standard deviation (Eqn 14).\cr
    \code{SigmaTot.cb} \tab Total standard deviation of the geometric
    mean horizontal component of LnY (Eqn 16).\cr
    \code{SigmaArb.cb} \tab Total standard deviation of the arbitrary
    horizontal component of LnY (Eqn 18).\cr
    \code{A1100.cb} \tab Calculation of median PGA when Vs30 = 1100 m/s\cr
    \code{SaMedian.cb} \tab Median estimate of ground motion (Eqn 1).
  }
}
\value{
  Returns values of the sub-functions for the Campbell & Bozorgnia
  NGA model.  The outputs of the functions are in log space, with
  the exception of (1) \code{SaMedian.cb}, which outputs values in g
  (for Sa and PGA), cm/s (for PGV), or cm (for PGD); and (2)
 \code{A1100.cb}, which outputs values in g.
}
\seealso{
  \code{\link{coefs}}, \code{\link{coefs.t.cb}},
  \code{\link{periods.cb}}, \code{\link{Sa.cb}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Campbell, K. W., and Y. Bozorgnia (2008). NGA Ground Motion Model for
  the Geometric Mean Horizontal Component of PGA, PGV, PGD, and 5\%
  Damped Linear Elastic Response Spectra for Periods Ranging from 0.01
  to 10 s. \emph{Earthquake Spectra} \bold{24,} 139--171.
}
\keyword{internal}

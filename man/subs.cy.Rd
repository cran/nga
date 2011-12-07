\name{CY08 Sub-Functions}
\alias{subs.cy}
\alias{LnYref.cy}
\alias{NL.cy}
\alias{Sigma.cy}
\alias{Tau.cy}
\alias{SigmaTot.cy}
\alias{SaMedian.cy}
\title{CY08 Sub-Functions}
\description{
  The sub-functions needed to estimate ground motions using the
  Chiou & Youngs NGA model.  These functions are called by the main
  function used for calculations, \code{\link{Sa.cy}}.
}
\usage{
LnYref.cy(M, Rrup, Rjb, Rx, Ztor, dip, Frv, Fnm, Fhw, AS, T)
NL.cy(Yref, Vs30, T)
Sigma.cy(M, Yref, Vs30, VsFlag, AS, T)
Tau.cy(M, T)
SigmaTot.cy(M, Yref, Vs30, VsFlag, AS, T)
SaMedian.cy(M, Rjb, Rrup, Rx, Ztor, Frv, Fnm, AS, Fhw, dip, Vs30, Z1.0, T)
}
\arguments{
  \item{M}{Moment magnitude of earthquake.}
  \item{Rjb}{Horizontal distance to the surface projection of the
    rupture plane; Joyner-Boore distance (km).}
  \item{Rrup}{Closest distance to the rupture plane; rupture distance
    (km).}
  \item{Rx}{Horizontal distance to the surface projection of the
    top edge of the rupture plane, measured perpendicular to the
    strike; site coordinate (km).}
  \item{Ztor}{Depth to top of rupture (km).}
  \item{Frv}{Reverse style-of-faulting flag; equal to 1 for reverse
    faulting (\code{30 <= rake <= 150 deg}), and 0 otherwise.}
  \item{Fnm}{Normal style-of-faulting flag; equal to 1 for normal
    faulting (\code{-120 <= rake <= -60 deg}), and 0 otherwise.}
  \item{AS}{Aftershock flag; equal to 1 for aftershocks and 0 for
    mainshocks.}
  \item{Fhw}{Hanging wall flag; equal to 1 for sites on the hanging
    wall side of the fault (\code{Rx >= 0}) and 0 otherwise.}
  \item{dip}{Fault dip angle (deg).}
  \item{Vs30}{Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s).}
  \item{VsFlag}{Flag variable indicating how \code{Vs30} is obtained;
    equal to 1 if \code{Vs30} is measured, and 0 if \code{Vs30} is
    estimated or inferred.}
  \item{Z1.0}{Depth to Vs = 1.0 km/s (m).}
  \item{Yref}{Median ground motion for the reference rock site
    condition (g).}
  \item{T}{Spectral period (sec). Use 0 for PGA and -1 for PGV.}
}
\details{
  \describe{
    \item{\code{LnYref.cy}}{Median ground motion on reference rock
    site (Eqn 13a of Chiou and Youngs (2008)).}
    \item{\code{NL.cy}}{Nonlinear site response term (Eqns 10 and 21).}
    \item{\code{Sigma.cy}}{Intra-event standard deviation (Eqn 20).}
    \item{\code{Tau.cy}}{Inter-event standard deviation (Eqn 19).}
    \item{\code{SigmaTot.cy}}{Total standard deviation (Eqn 21).}
    \item{\code{SaMedian.cy}}{Median estimate of ground motion, with
    site response (Eqn 13b).}
  }
}
\value{
  Returns values of the sub-functions for the Chiou & Youngs
  NGA model.  The outputs of the functions are in log space, with
  the exception of \code{SaMedian.cy}, which outputs values in g
  (for Sa and PGA) or cm/s (for PGV).
}
\seealso{
  \code{\link{coefs}}, \code{\link{coefs.t.cy}},
  \code{\link{periods.cy}}, \code{\link{Sa.cy}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Chiou, B. S.-J., and R. R. Youngs (2008). An NGA Model for the Average
  Horizontal Component of Peak Ground Motion and Response Spectra.
  \emph{Earthquake Spectra} \bold{24,} 173--215.
}
\keyword{internal}

\name{BA08 Sub-Functions}
\alias{subs.ba}
\alias{Fd.ba}
\alias{Fm.ba}
\alias{Fs.ba}
\alias{SaMedian.ba}
\title{BA08 Sub-Functions}
\description{
  The sub-functions needed to estimate ground motions using the
  Boore & Atkinson NGA model.  These functions are called by the main
  function used for calculations, \code{\link{Sa.ba}}.
}
\usage{
Fd.ba(M, Rjb, T)
Fm.ba(M, U, SS, NS, RS, T)
Fs.ba(M, Rjb, Vs30, U, SS, NS, RS, T)
SaMedian.ba(M, Rjb, Vs30, U, SS, NS, RS, T)
}
\arguments{
  \item{M}{Moment magnitude of earthquake.}
  \item{Rjb}{Horizontal distance to the surface projection of the
    rupture plane; Joyner-Boore distance (km).}
  \item{Vs30}{Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s).}
  \item{U}{Unspecified style-of-faulting flag; equal to 1 if the user
    wishes to perform a generic ground motion calculation when the style
    of faulting is unspecified, and 0 otherwise.}
  \item{RS}{Reverse style-of-faulting flag; equal to 1 for reverse
    faulting (\code{30 <= rake <= 150 deg}), and 0 otherwise.}
  \item{NS}{Normal style-of-faulting flag; equal to 1 for normal
    faulting (\code{-150 <= rake <= -30 deg}), and 0 otherwise.}
  \item{SS}{Strike-slip style-of-faulting flag; equal to 1 for
    strike-slip faulting (when the rake is not in either of the ranges
    specified for \code{RS} or \code{NS}), and 0 otherwise.}
  \item{T}{Spectral period (sec). Use 0 for PGA and -1 for PGV.}
}
\details{
  \tabular{ll}{
    \code{Fd.ba} \tab Distance function (Eqns 3 and 4 of Boore and Atkinson
    (2008)).\cr
    \code{Fm.ba} \tab Magnitude function (Eqn 5).\cr
    \code{Fs.ba} \tab Site amplification function (Eqns 6 to 13).\cr
    \code{SaMedian.ba} \tab Median estimate of ground motion (Eqn 1).
  }
}
\value{
  Returns values of the sub-functions for the Boore & Atkinson
  NGA model.  The outputs of the functions are in log space, with
  the exception of \code{SaMedian.ba}, which outputs values in g
  (for Sa and PGA) or cm/s (for PGV).
}
\seealso{
  \code{\link{coefs}}, \code{\link{coefs.t.ba}},
  \code{\link{periods.ba}}, \code{\link{Sa.ba}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Boore, D. M., and G. M. Atkinson (2008). Ground-Motion Prediction
  Equations for the Average Horizontal Component of PGA, PGV, and
  5\%-Damped PSA at Spectral Periods between 0.01 s and 10.0
  s. \emph{Earthquake Spectra} \bold{24,} 99--138.
}
\keyword{internal}

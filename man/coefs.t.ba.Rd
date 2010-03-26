\name{BA08 Period-Dependent Model Coefficients}
\alias{coefs.t.ba}
\alias{blin.ba}
\alias{b1.ba}
\alias{b2.ba}
\alias{c1.ba}
\alias{c2.ba}
\alias{c3.ba}
\alias{h.ba}
\alias{e1.ba}
\alias{e2.ba}
\alias{e3.ba}
\alias{e4.ba}
\alias{e5.ba}
\alias{e6.ba}
\alias{e7.ba}
\alias{Mh.ba}
\alias{Sigma.ba}
\alias{TauU.ba}
\alias{TauM.ba}
\alias{SigmaTotU.ba}
\alias{SigmaTotM.ba}
\title{BA08 Period-Dependent Model Coefficients}
\description{
  Functions for computing the period-dependent coefficients of the
  Boore & Atkinson NGA model. 
}
\usage{
blin.ba(T)
b1.ba(T)
b2.ba(T)
c1.ba(T)
c2.ba(T)
c3.ba(T)
h.ba(T)
e1.ba(T)
e2.ba(T)
e3.ba(T)
e4.ba(T)
e5.ba(T)
e6.ba(T)
e7.ba(T)
Mh.ba(T)
Sigma.ba(T)
TauU.ba(T)
TauM.ba(T)
SigmaTotU.ba(T)
SigmaTotM.ba(T)
}
\arguments{
  \item{T}{Spectral period (sec).  Use 0 for PGA and -1 for PGV.}
}
\details{
  \describe{
    \item{Site-amplification coefficients (Table 3 in
      Boore and Atkinson (2008)):}{\code{blin.ba, b1.ba,} and
      \code{b2.ba}.}
    \item{Distance-scaling coefficients (Table 6):}{\code{c1.ba, c2.ba,
	c3.ba,} and \code{h.ba}.} 
    \item{Magnitude-scaling coefficients (Table 7):}{\code{e1.ba, e2.ba,
	e3.ba, e4.ba, e5.ba, e6.ba, e7.ba,} and \code{Mh.ba}.}
    \item{Coefficients for standard deviation term (Table
      8):}{\code{Sigma.ba; TauU.ba, SigTotU.ba} (unspecified fault
      mechanism) \code{; TauM.ba, SigTotM.ba} (specified fault
      mechanism)}
  }
}
\value{
  Coefficients for the Boore & Atkinson NGA model.
}
\seealso{
  \code{\link{coefs}}, \code{\link{periods.ba}}, \code{\link{subs.ba}},
  \code{\link{Sa.ba}}, \code{\link{Sa.nga}}.
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

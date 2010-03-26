\name{CY08 Period-Dependent Model Coefficients}
\alias{coefs.t.cy}
\alias{c1.cy}
\alias{c1a.cy}
\alias{c1b.cy}
\alias{cn.cy}
\alias{cM.cy}
\alias{c5.cy}
\alias{c6.cy}
\alias{c7.cy}
\alias{c7a.cy}
\alias{c9.cy}
\alias{c9a.cy}
\alias{c10.cy}
\alias{cGamma1.cy}
\alias{cGamma2.cy}
\alias{phi1.cy}
\alias{phi2.cy}
\alias{phi3.cy}
\alias{phi4.cy}
\alias{phi5.cy}
\alias{phi6.cy}
\alias{phi7.cy}
\alias{phi8.cy}
\alias{tau1.cy}
\alias{tau2.cy}
\alias{sigma1.cy}
\alias{sigma2.cy}
\alias{sigma3.cy}
\alias{sigma4.cy}
\title{CY08 Period-Dependent Model Coefficients}
\description{
  Functions for computing the period-dependent coefficients of the
  Chiou & Youngs NGA model. 
}
\usage{
c1.cy(T)
c1a.cy(T)
c1b.cy(T)
cn.cy(T)
cM.cy(T)
c5.cy(T)
c6.cy(T)
c7.cy(T)
c7a.cy(T)
c9.cy(T)
c9a.cy(T)
c10.cy(T)
cGamma1.cy(T)
cGamma2.cy(T)
phi1.cy(T)
phi2.cy(T)
phi3.cy(T)
phi4.cy(T)
phi5.cy(T)
phi6.cy(T)
phi7.cy(T)
phi8.cy(T)
tau1.cy(T)
tau2.cy(T)
sigma1.cy(T)
sigma2.cy(T)
sigma3.cy(T)
sigma4.cy(T)
}
\arguments{
  \item{T}{Spectral period (sec).  Use 0 for PGA and -1 for PGV.}
}
\details{
  \describe{
    \item{Coefficients of model for Ln(Yref) (Table 2 in Chiou and
      Youngs (2008)):}{\code{c1.cy, c1a.cy, c1b.cy, cn.cy, cM.cy, c5.cy,
	c6.cy, c7.cy, c7a.cy, c9.cy, c9a.cy, c10.cy, cGamma1.cy,} and
      \code{cGamma2.cy}.}
    \item{Coefficients of site response model for Ln(Y)
      (Table 3):}{\code{phi1.cy, phi2.cy, phi3.cy, phi4.cy, phi5.cy,
	phi6.cy, phi7.cy,} and \code{phi8.cy}.}
    \item{Coefficients of standard deviation model
      (Table 4):}{\code{tau1.cy, tau2.cy, sigma1.cy, sigma2.cy,
	sigma3.cy,} and \code{sigma4.cy}.}
  }
}
\value{
  Coefficients for the Chiou & Youngs NGA model.
}
\seealso{
  \code{\link{coefs}}, \code{\link{periods.cy}}, \code{\link{subs.cy}},
  \code{\link{Sa.cy}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Chiou, B. S.-J., and R. R. Youngs (2008). An NGA Model for the Average
  Horizontal Component of Peak Ground Motion and Response Spectra.
  \emph{Earthquake Spectra} \bold{24,} 173--215.
}
\keyword{internal}

\name{Period-Independent Model Coefficients}
\alias{coefs}
\alias{coefs.as}
\alias{coefs.ba}
\alias{coefs.cb}
\alias{coefs.cy}
\title{Period-Independent Model Coefficients}
\description{
  Defines the period-independent coefficients for the different NGA
  models. 
}
\usage{
coefs.as()
coefs.ba()
coefs.cb()
coefs.cy()
}
\details{
  Here is a list of the period-independent coefficients for each of the
  NGA models:

  \describe{
   \item{\code{coefs.as:}}{\code{c, c1, c2, c4,} and \code{n} (Table 4 in
           Abrahamson and Silva (2008)).}
  \item{\code{coefs.ba:}}{\code{a1, pga.low, a2, V1, V2,} and \code{Vref}
          (Site-amplification coefficients in Table 4 of Boore and Atkinson
  	  (2008));}
  \item{}{\code{Mref} and \code{Rref} (Distance-scaling coefficients in
            Table 6).}
  \item{\code{coefs.cb:}}{\code{c} and \code{n} (bottom of Table 2 in
          Campbell and Bozorgnia (2008)).}
  \item{\code{coefs.cy:}}{\code{c2, c3, c4, c4a, cRB, cHM,} and \code{cGamma3} (Table
          1 in Chiou and Youngs (2008)).}
   }
}
\value{
  Nothing is returned; variables are assigned in the parent
  environment. 
}
\seealso{
  \code{\link{coefs.t.as}}, \code{\link{coefs.t.ba}},
  \code{\link{coefs.t.cb}}, \code{\link{coefs.t.cy}},
  \code{\link{periods}}
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Abrahamson, N., and W. Silva (2008). Summary of the Abrahamson & Silva
  NGA Ground-Motion Relations. \emph{Earthquake Spectra} \bold{24,}
  67--97.
  
  Boore, D. M., and G. M. Atkinson (2008). Ground-Motion Prediction
  Equations for the Average Horizontal Component of PGA, PGV, and
  5\%-Damped PSA at Spectral Periods between 0.01 s and 10.0
  s. \emph{Earthquake Spectra} \bold{24,} 99--138.
  
  Campbell, K. W., and Y. Bozorgnia (2008). NGA Ground Motion Model for
  the Geometric Mean Horizontal Component of PGA, PGV, PGD, and 5\%
  Damped Linear Elastic Response Spectra for Periods Ranging from 0.01
  to 10 s. \emph{Earthquake Spectra} \bold{24,} 139--171.

  Chiou, B. S.-J., and R. R. Youngs (2008). An NGA Model for the Average
  Horizontal Component of Peak Ground Motion and Response Spectra.
  \emph{Earthquake Spectra} \bold{24,} 173--215.
}
\keyword{internal}

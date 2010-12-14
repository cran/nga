\name{Distance Calculations}
\alias{Rx.calc}
\alias{Rrup.calc}
\title{Calculation of Source-to-Site Distance Measures}
\description{
  Calculates the values of the rupture distance (\code{Rrup}) and
  site coordinate (\code{Rx}) from the other distance parameters and the
  geometric source characteristics of the fault rupture.  The equations
  for \code{Rx} and \code{Rrup} are derived and explained in Kaklamanos
  et al. (2011).
}
\usage{
Rx.calc(Rjb, Ztor, W, dip, azimuth, Rrup = NA)
Rrup.calc(Rx, Ztor, W, dip, azimuth, Rjb = NA)
}
\arguments{
  \item{Rjb}{Horizontal distance to the surface projection of the
    rupture plane; Joyner-Boore distance (km).}
  \item{Rrup}{Closest distance to the rupture plane; rupture distance
    (km).}
  \item{Rx}{Horizontal distance to the surface projection of the
    top edge of the rupture plane, measured perpendicular to the
    strike; site coordinate (km).}
  \item{Ztor}{Depth to top of rupture (km).}
  \item{W}{Down-dip rupture width (km).}
  \item{dip}{Fault dip angle (deg).}
  \item{azimuth}{source-to-site azimuth (deg); see Kaklamanos et
    al. (2011) for a description.}
}
\details{
  The distance functions for \code{Rx} and \code{Rrup} require that the
  Joyner-Boore distance (\code{Rjb}) be known.  The source-to-site
  azimuth is also a necessary argument; if the exact azimuth is unknown,
  assume a generic value of 50 degrees for sites on the hanging wall
  side of the fault and -50 degrees for sites on the footwall side of
  the fault.  An analysis of the database used to derive the NGA
  relations suggests that these values are reasonable.  The geometric
  source parameters \code{Ztor}, \code{W}, and \code{dip} are also
  required; for methods of estimating these source parameters when they
  are not known beforehand, see \code{\link{Ztor.calc},}
  \code{\link{W.calc},} and \code{\link{dip.calc},} respectively.
  
  A general strategy for calculating distances is to first calculate
  \code{Rx}, and then calculate \code{Rrup} using \code{Rx}.  In order
  to calculate \code{Rx} using the function \code{Rx.calc}, the argument
  \code{Rrup} is only necessary when the site is located directly over
  the ruptured area (\code{Rjb = 0}).  If \code{Rrup} is unknown in this
  case, then the function assumes that the site is located in the middle
  of the surface projection of the ruptured area.  In the function
  \code{Rrup.calc}, the argument \code{Rjb} is only necessary in the
  rare case that the site is located directly on the surface projection
  of fault strike (\code{azimuth = 0, 180,} or \code{-180}).
}
\value{
  \code{Rx.calc} outputs \code{Rx} (the \dQuote{site coordinate}),
  the horizontal distance to the surface projection of the top edge of
  the rupture plane, measured perpendicular to the strike (km).
  
  \code{Rrup.calc} outputs \code{Rrup} (the \dQuote{rupture distance}),
  the closest distance to the rupture plane (km).
}
\seealso{
  \code{\link{Ztor.calc}}, \code{\link{W.calc}}, \code{\link{dip.calc}},
  \code{\link{trig}}, \code{\link{Sa}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Kaklamanos, J., L. G. Baise, and D. M. Boore (2011).  Estimating
  Unknown Input Parameters when Implementing the NGA Ground-Motion
  Prediction Equations in Engineering Practice.
  \emph{Earthquake Spectra} (in press).
}
\examples{
#########################################################################
# Example 1:  Calculate the distance measures for a synthetic example,
#             with Rjb = 5

# Assumed source and location parameters
M <- 6
Rjb <- 5
azimuth <- 15
rake <- 90   # Reverse fault

# Estimate Ztor, W, and dip, using the respective functions

W <- W.calc(M = M, rake = rake)
W

dip <- dip.calc(rake = rake)
dip

Zhyp <- Zhyp.calc(M = M, rake = rake)
Zhyp
# Zhyp is needed in order to estimate Ztor

Ztor <- Ztor.calc(W = W, dip = dip, Zhyp = Zhyp)
Ztor


# Estimate Rx and Rrup
Rx <- Rx.calc(Rjb = Rjb, Ztor = Ztor, W = W, dip = dip,
              azimuth = azimuth, Rrup = NA)
Rx

Rrup <- Rrup.calc(Rx = Rx, Ztor = Ztor, W = W, dip = dip,
                  azimuth = azimuth, Rjb = Rjb)
Rrup



#########################################################################
# Example 2:  Calculate and plot the distance measures for a synthetic
#             example, for values of Rjb ranging from 0 to 100

# Redefine Rjb as a vector
Rjb <- seq(from = 0, to = 20, by = 0.5)

# Calculate Rx; vectorize the calculation using the intrinsic
# R "sapply" function
Rx <- sapply(Rjb, Rx.calc, Ztor = Ztor, W = W, dip = dip,
             azimuth = azimuth, Rrup = NA)

# Calculate Rrup, again using the "sapply" function
Rrup <- sapply(Rx, Rrup.calc, Ztor = Ztor, W = W, dip = dip,
               azimuth = azimuth, Rjb = Rjb)


# Plot the results against Rjb
# Make basic plot
plot(Rjb, Rjb, type = "l", xaxs = "i", yaxs = "i", xlab = "Rjb (km)",
     ylab = "Rjb, Rrup, and Rx (km)", main = "Comparison of Distance Measures",
     col = "black", lwd = 2)

# Add line for Rrup
lines(Rjb, Rrup, col = "red", lwd = 2)

# Add line for Rx
lines(Rjb, Rx, col = "blue", lwd = 2)

# Add legend
legend(x = "topleft", inset = 0.02, lwd = 2, bty = "n",
       col = c("black", "red", "blue"),
       legend = c("Rjb", "Rrup", "Rx"))

}


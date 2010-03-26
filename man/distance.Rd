\name{Distance Calculations}
\alias{Rx.calc}
\alias{Rrup.calc}
\title{Calculation of Source-to-Site Distance Measures}
\description{
  Calculates the values of the rupture distance (Rrup) and
  site coordinate (Rx) from the other distance parameters and the
  geometric source characteristics of the fault rupture.  The equations
  for Rx and Rrup are derived and explained in Kaklamanos and Baise (2010).
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
  \item{azimuth}{source-to-site azimuth (deg); see Figure 2 in
    Kaklamanos and Baise (2010).}
}
\details{
  The distance functions for Rx and Rrup require that the Joyner-Boore
  distance (Rjb) be known.  The source-to-site azimuth is also a
  necessary argument; if the exact azimuth is unknown, assume a
  generic value of 50 degrees for sites on the hanging wall side of
  the fault and -50 degrees for sites on the footwall side of the fault.
  The geometric source parameters Ztor, W, and dip are also required;
  for methods of estimating these source parameters
  when they are not known beforehand, see \code{\link{Ztor.calc},}
  \code{\link{W.calc},} and \code{\link{dip.calc},} respectively.
  
  A general strategy for calculating distances is to first calculate Rx,
  and then calculate Rrup using Rx.   In order to calculate Rx using
  the function \code{Rx.calc}, the argument Rrup is only necessary when
  the site is located directly over the ruptured area (Rjb = 0).  If
  Rrup is unknown in this case, then the function assumes that the site
  is located in the middle of the surface projection of the ruptured
  area.  In the function \code{Rrup.calc}, the argument Rjb is only
  necessary in the rare case that the site is located directly on the
  surface projection of fault strike (azimuth = 0, 180, or -180).
}
\value{
  \code{Rx.calc} outputs Rx (the \dQuote{site coordinate}),
  the horizontal distance to the surface projection of the top edge of
  the rupture plane, measured perpendicular to the strike (km).
  
  \code{Rrup.calc} outputs Rrup (the \dQuote{rupture distance}), the closest
  distance to the rupture plane (km).
}
\seealso{
  \code{\link{Ztor.calc}}, \code{\link{W.calc}}, \code{\link{dip.calc}},
  \code{\link{trig}}, \code{\link{Sa}}, \code{\link{Sa.nga}}.
}
\author{James Kaklamanos <james.kaklamanos@tufts.edu> and
  Eric M. Thompson <eric.thompson@tufts.edu>}
\references{
  Kaklamanos, J., and L. G. Baise (2010).  Relationships between
  Distance Measures in Recent Ground Motion Prediction Equations.
  \emph{Seismological Research Letters} (in review).
}
\examples{
# Example 1:  Calculate the distance measures for the synthetic example
#             presented in Kaklamanos and Baise (2010), Rjb = 5

# Assumed source and location parameters
M <- 6
Rjb <- 5
# Rjb <- seq(from = 0, to = 200, by = 1)
azimuth <- 15
rake <- 90   # Reverse fault

# Estimate Ztor, W, and dip, using the respective functions

W <- W.calc(M = M, rake = rake)
W

dip <- dip.calc(rake = rake)
dip

Ztor <- Ztor.calc(W = W, dip = dip, M = M, rake = rake)
Ztor


# Estimate Rx and Rrup
Rx <- Rx.calc(Rjb = Rjb, Ztor = Ztor, W = W, dip = dip,
              azimuth = azimuth, Rrup = NA)
Rx

Rrup <- Rrup.calc(Rx = Rx, Ztor = Ztor, W = W, dip = dip,
                  azimuth = azimuth, Rjb = Rjb)
Rrup




# Example 2:  Calculate the distance measures for the synthetic example
#             presented in Kaklamanos and Baise (2010), for values of
#             Rjb ranging from 0 to 100

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


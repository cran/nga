\name{Ground Motion Predictions for all NGA Models}
\alias{Sa.nga}
\alias{Sa.ngaR}
\alias{nga}
\title{Ground Motion Predictions for all NGA Models}
\description{
  Comprehensive functions that estimate ground motion parameters using the
  AS08, BA08, CB08, and CY08 models from the Next Generation Attenuation
  of Ground Motions (NGA) project in 2008.  The function \code{Sa.nga}
  is designed to mimic the output from Boore and Campbell's FORTRAN
  output files; the function \code{Sa.ngaR} has reduced output.
}
\usage{
Sa.nga(M, Rjb, Rrup = NA, Rx = NA, rake = NA, Frv = NA,
       Fnm = NA, Fhw = NA, dip = NA, W = NA, Ztor = NA, Vs30,
       Z1.0 = NA, Z1.5 = NA, Z2.5 = NA, VsFlag = 0, Fas = 0,
       Zhyp = NA, azimuth = NA, U = 0, arb = 0, epsilon, T)
Sa.ngaR(M, Rjb, Rrup = NA, Rx = NA, rake = NA, Frv = NA,
        Fnm = NA, Fhw = NA, dip = NA, W = NA, Ztor = NA, Vs30,
        Z1.0 = NA, Z1.5 = NA, Z2.5 = NA,  VsFlag, Fas = 0,
        Zhyp = NA, azimuth = NA, U = 0, arb = 0, epsilon, T)
}
\arguments{
  \item{M}{Moment magnitude of earthquake.}
  \item{Rjb}{Joyner-Boore distance (km):  the horizontal distance to the
    surface projection of the rupture plane.}
  \item{Rrup}{Rupture distance (km):  the closest distance to the
    rupture plane; if left empty (or set negative), Rrup is calculated
    from Rx, the source-to-site azimuth, and the geometric rupture
    parameters (Ztor, W, and dip) using \code{\link{Rrup.calc}}.}
  \item{Rx}{Site coordinate (km):  The horizontal distance to the
    surface projection of the top edge of the rupture plane, measured
    perpendicular to the strike.  If left empty, Rx is calculated from
    Rjb, the source-to-site azimuth, and the geometric rupture
    parameters (Ztor, W, and dip) using \code{\link{Rx.calc}}.  When
    only Rjb and the azimuth are assumed, Rjb is used to calculate Rx,
    which is then used to calculate Rrup.}
  \item{rake}{Rake angle of fault movement (deg).  Either the rake angle
    or the style-of-faulting flag variables (Frv and Fnm) must be specified.}
  \item{Frv}{Reverse style-of-faulting flag (1 for reverse faulting, 0
    otherwise).  Either (a) the rake angle, or (b) both Frv and Fnm,
    must be specified.  Reverse faulting is characterized by rake angles
    in the range 30 <= rake <= 150 deg.}
  \item{Fnm}{Normal style-of-faulting flag (1 for normal faulting, 0
    otherwise).  Either (a) the rake angle, or (b) both Frv and Fnm,
    must be specified.  Normal faulting is characterized by rake angles
    in the range -120 <= rake <= -60 deg for the AS08 and CY08 models
    and -150 <= rake <= -30 deg for the BA08 and CB08 models.}
  \item{Fhw}{Hanging wall flag; equal to 1 for sites on the hanging
    wall side of the fault (Rx >= 0; azimuth >= 0), and 0 otherwise.} 
  \item{dip}{Dip angle of the rupture plane (deg).  If left empty (or set
    negative), the dip is estimated using \code{\link{dip.calc}}.}
  \item{W}{Down-dip width of rupture plane (km).  If left empty (or set
    negative), W is estimated using \code{\link{W.calc}}.}
  \item{Ztor}{Depth to top of rupture (km).  If left empty (or set
    negative), Ztor is estimated using \code{\link{Ztor.calc}}.}
  \item{Vs30}{Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s).}
  \item{Z1.0}{Depth to Vs = 1.0 km/s (m).  If left empty (or set
    negative), Z1.0 is estimated using \code{\link{Z1.calc.as}} for the
    AS08 model and \code{\link{Z1.calc.cy}} for the CY08 model.}
  \item{Z1.5}{Depth to Vs = 1.5 km/s (m).  Z1.5 is not utilized in
    ground motion calculations, but if available, it may be used to
    estimate Z2.5 for the CB08 model.}
  \item{Z2.5}{Depth to Vs = 2.5 km/s (km).  If left empty (or set negative),
    Z2.5 is estimated from Z1.5 or Z1.0 if available, using the
    recommendations in Campbell and Bozorgnia (2007).  If neither Z1.5 nor
    Z1.0 is available, then Vs30 is used to estimate Z1.0 using
    \code{\link{Z1.calc.as}}, which is in turn used to estimate Z2.5.}
  \item{VsFlag}{Flag variable indicating how Vs30 is obtained; equal
    to 1 if Vs30 is measured, and 0 if Vs30 is estimated or inferred.
    Necessary for \code{Sa.ngaR}, but not necessary for \code{Sa.nga}
    because ground motions for both cases are output.}
  \item{Fas}{Aftershock flag; equal to 1 for aftershocks and 0 for
    mainshocks (the default)}
  \item{Zhyp}{Hypocentral depth of the earthquake (km).  Zhyp is not
    utilized in ground motion calculations, but it may be used to estimate
    Ztor.  See \code{\link{Ztor.calc}} for details.}
  \item{azimuth}{Source-to-site azimuth (deg); see Figure 2 in
    Kaklamanos and Baise (2010) for details.  Used by
    \code{\link{Rx.calc}} and \code{\link{Rrup.calc}} for distance
    calculations.}
  \item{U}{Unspecified style-of-faulting flag for the BA08 model; equal
    to 1 if the user wishes to perform a generic ground motion
    calculation when the style of faulting is unspecified, and 0
    otherwise (the default).  Not considered for \code{Sa.nga}
    because ground motions for both cases are output.}
  \item{arb}{Flag variable indicating the method of determining aleatory
    uncertainty for the CB08 model; equal to 1 if the standard deviation
    should be calculated for the arbitrary horizontal component of ground
    motion, and 0 if the standard deviation should be calculated for the
    geometric mean horizontal ground motion (the default).  Not
    considered for \code{Sa.nga} because ground motions for both cases
    are output.}
  \item{epsilon}{number of standard deviations to be considered in the
    calculations.}
  \item{T}{Spectral period (sec).  Use 0 for PGA, -1 for PGV, and -2 for
    PGD (CB08 only).}
}
\value{
  The function \code{Sa.nga} outputs a list composed of the following
  elements:

  \bold{INPUT VARIABLES:}
  \tabular{ll}{
    \code{T} \tab Spectral period, sec [input] \cr
    \code{epsilon} \tab number of standard deviations considered in the
    calculations [input] \cr
    \code{M} \tab Moment magnitude [input] \cr
    \code{Rjb} \tab Joyner-Boore distance (km) [input] \cr
    \code{Rrup.in} \tab Rupture distance (km) [input] \cr
    \code{Rrup.out} \tab Rupture distance (km) [calculated if \code{Rrup.in} is
    not specified] \cr
    \code{Rx.in} \tab Site coordinate (km) [input] \cr
    \code{Rx.out} \tab Site coordinate (km) [calculated if \code{Rx.in} is not
    specified] \cr
    \code{rake.in} \tab Rake angle of fault movement (deg) [input] \cr
    \code{rake.out} \tab Rake angle of fault movement (deg) [calculated
    if \code{rake.in} is not specified] \cr
    \code{Frv} \tab Reverse style-of-faulting flag [input] \cr
    \code{Fnm1} \tab Normal style-of-faulting flag for AB08 and CY08 \cr
    \code{Fnm2} \tab Normal style-of-faulting flag for BA08 and CB08 \cr
    \code{Fhw} \tab Hanging wall flag \cr
    \code{dip.in} \tab Fault dip angle (deg) [input] \cr
    \code{dip.out} \tab Fault dip angle (deg) [calculated if \code{dip.in} is
    not specified] \cr
    \code{W.in} \tab Down-dip rupture width (km) [input] \cr
    \code{W.out} \tab Down-dip rupture width (km) [calculated if
    \code{W.in} is not specified] \cr
    \code{Ztor.in} \tab Depth to top of rupture (km) [input] \cr
    \code{Ztor.out} \tab Depth to top of rupture (km) [calculated if
    \code{Ztor.in} is not specified] \cr
    \code{Vs30} \tab Time-averaged shear wave velocity over 30 m
    subsurface depth  (m/sec) [input] \cr
    \code{Z1.0in} \tab Depth to Vs of 1.0 km/sec  (m) [input] \cr
    \code{Z1.0as} \tab Depth to Vs of 1.0 km/sec  (m) [calculated for
    use in AS08 model] \cr
    \code{Z1.0cy} \tab Depth to Vs of 1.0 km/sec  (m) [calculated for
    use in CY08 model] \cr
    \code{Z1.5in} \tab Depth to Vs of 1.5 km/sec  (m) [input] \cr
    \code{Z2.5in} \tab Depth to Vs of 2.5 km/sec  (km) [input] \cr
    \code{Z2.5out} \tab Depth to Vs of 2.5 km/sec  (km) [calculated from
    Z1.0 for use in CB08 model] \cr
    \code{VsFlag} \tab Flag variable indicating how Vs30 is obtained
    [input] \cr
    \code{Fas} \tab Aftershock flag [input] \cr
    \code{Zhyp} \tab hypocentral depth (km) [input] \cr
    \code{azimuth.in} \tab source-to-site azimuth (deg) [input] \cr
    \code{azimuth.out} \tab source-to-site azimuth (deg) [calculated if
    \code{azimuth.in} is not specified] \cr
    \code{U} \tab Unspecified style-of-faulting flag for the BA08 model
    [input] \cr
    \code{arb} \tab Flag variable indicating the method of determining aleatory
    uncertainty for the CB08 model [input]
  }
  \bold{OUTPUT VARIABLES:}
  \tabular{ll}{
    \bold{AS08 Model:} \tab  \cr
    \code{Y50.as} \tab Median ground motion estimate using AS08 (epsilon
    = 0) \cr
    \code{YplusEpsilon.meas.as} \tab Upper ground motion estimate
    using AS08, for measured Vs30 (VsFlag = 1) \cr
    \code{YplusEpsilon.est.as} \tab Upper ground motion estimate
    using AS08, for estimated Vs30 (VsFlag = 0) \cr
    \code{YminusEpsilon.meas.as} \tab Lower ground motion estimate
    using AS08, for measured Vs30 (VsFlag = 1) \cr
    \code{YminusEpsilon.est.as} \tab Lower ground motion estimate
    using AS08, for estimated Vs30 (VsFlag = 0) \cr
    \code{sdMeas.as} \tab total standard deviation using AS08, for
    measured Vs30 (VsFlag = 1) \cr
    \code{sdEst.as} \tab total standard deviation using AS08, for
    estimated Vs30 (VsFlag = 0) \cr
    \bold{BA08 Model:} \tab  \cr
    \code{Y50M.ba} \tab Median ground motion estimate using BA08, when
    fault type is specified (U = 0) \cr
    \code{Y50U.ba} \tab Median ground motion estimate using BA08, when
    fault type is unspecified (U = 1) \cr
    \code{YplusEpsilon.M.ba} \tab Upper ground motion estimate
    using BA08, when fault type is specified (U = 0) \cr
    \code{YplusEpsilon.U.ba} \tab Upper ground motion estimate
    using BA08, when fault type is unspecified (U = 1) \cr
    \code{YminusEpsilon.M.ba} \tab Lower ground motion estimate
    using BA08, when fault type is specified (U = 0) \cr
    \code{YminusEpsilon.U.ba} \tab Lower ground motion estimate
    using BA08, when fault type is unspecified (U = 1) \cr
    \code{sdM.ba} \tab total standard deviation using BA08, when
    fault type is specified (U = 0) \cr
    \code{sdU.ba} \tab total standard deviation using BA08, when
    fault type is unspecified (U = 1) \cr
    \bold{CB08 Model:} \tab  \cr
    \code{Y50.cb} \tab Median ground motion estimate using CB08 (epsilon
    = 0) \cr
    \code{YplusEpsilon.GM.cb} \tab Upper CB08 estimate for the
    geometric mean horizontal component (arb = 0) \cr
    \code{YplusEpsilon.arb.cb} \tab Upper CB08 estimate for the
    arbitrary horizontal component (arb = 1) \cr
    \code{YminusEpsilon.GM.cb} \tab Lower CB08 estimate for the
    geometric mean horizontal component (arb = 0) \cr
    \code{YminusEpsilon.arb.cb} \tab Lower CB08 estimate for the
    arbitrary horizontal component (arb = 1) \cr
    \code{sdGM.cb} \tab CB08 total standard deviation
    for the geometric mean horizontal component (arb = 0) \cr
    \code{sdArb.cb} \tab CB08 total standard deviation for the
    arbitrary horizontal component (arb = 1) \cr
    \bold{CY08 Model:} \tab  \cr
    \code{Y50.cy} \tab Median ground motion estimate using CY08 (epsilon
    = 0) \cr
    \code{YplusEpsilon.meas.cy} \tab Upper ground motion estimate
    using CY08, for measured Vs30 (VsFlag = 1) \cr
    \code{YplusEpsilon.est.cy} \tab Upper ground motion estimate
    using CY08, for estimated Vs30 (VsFlag = 0) \cr
    \code{YminusEpsilon.meas.cy} \tab Lower ground motion estimate
    using CY08, for measured Vs30 (VsFlag = 1) \cr
    \code{YminusEpsilon.est.cy} \tab Lower ground motion estimate
    using CY08, for estimated Vs30 (VsFlag = 0) \cr
    \code{sdMeas.cy} \tab total standard deviation using CY08, for
    measured Vs30 (VsFlag = 1) \cr
    \code{sdEst.cy} \tab total standard deviation using CY08, for
    estimated Vs30 (VsFlag = 0) \cr
     \tab
  }
  
  The function \code{Sa.ngaR} takes advantage of the input parameters
  \code{VsFlag} (AS08 and CY08), \code{U} (BA08), and \code{arb} (CB08),
  and returns a list of reduced output.  The \dQuote{Input Variables}
  section of the list is the same, except \code{VsFlag}, \code{U}, and
  \code{arb} are also included in the list.  The simplified \dQuote{Output
    Variables} section of \code{Sa.ngaR} is:
  \tabular{ll}{
    \bold{AS08 Model:} \tab  \cr
    \code{Y50.as} \tab Median ground motion estimate using AS08 (epsilon
    = 0) \cr
    \code{YplusEpsilon.as} \tab Upper bound estimate of ground
    motion using AS08 \cr
    \code{YminusEpsilon.as} \tab Lower bound estimate of ground
    motion using AS08 \cr
    \bold{BA08 Model:} \tab  \cr   
    \code{Y50.ba} \tab  Median ground motion estimate using BA08
    (epsilon = 0) \cr
    \code{YplusEpsilon.ba} \tab Upper bound estimate of ground motion
    using BA08 \cr
    \code{YminusEpsilon.ba} \tab Lower bound estimate of ground motion
    using BA08 \cr
    \bold{CB08 Model:} \tab  \cr
    \code{Y50.cb} \tab Median ground motion estimate using CB08 (epsilon
    = 0) \cr
    \code{YplusEpsilon.cb} \tab Upper bound estimate of ground motion
    using CB08 \cr
    \code{YminusEpsilon.cb} \tab Lower bound estimate of ground
    motion using CB08 \cr
    \bold{CY08 Model:} \tab  \cr    
    \code{Y50.cy} \tab Median ground motion estimate using CY08 (epsilon
    = 0) \cr
    \code{YplusEpsilon.cy} \tab Upper bound estimate of ground
    motion using CY08 \cr
    \code{YminusEpsilon.cy} \tab Lower bound estimate of ground
    motion using CY08 \cr
     \tab
   }
   
   In the \dQuote{Output Section} of the lists returned by these
   functions, \dQuote{Y} refers to the ground motion parameter of
   interest, which can be:
   
   Sa = Spectral acceleration (g);
   
   PGA = Peak ground acceleration (g), calculated by evaluating Sa at T
   = 0;
   
   PGV = Peak ground velocity (cm/sec), calculated by evaluating Sa at T
   = -1; or
   
   PGD = Peak ground displacement (cm), calculated by evaluating Sa at T = -2 (CB08 only).

   
   \dQuote{sd} refers to the standard deviation of the ground motion
   estimate, which is presented in log space.
}
\details{
Note that T (spectral period) can be a vector, while all other arguments are 
scalars.
}
\seealso{
  See \code{\link{Sa.as}}, \code{\link{Sa.ba}}, \code{\link{Sa.cb}},
  \code{\link{Sa.cy}} for separate functions that compute ground motion
  parameters using the individual NGA models. See
  \code{\link{mainshocks}} for an example of inputting and outputting
  earthquake data and predictions.
    
  For details on the sub procedures used for the individual NGA models,
  see \code{\link{subs.as}}, \code{\link{subs.ba}},
  \code{\link{subs.cb}}, and \code{\link{subs.cy}}.  See
  \code{\link{coefs}} for details on the period-independent model
  coefficients, and \code{\link{coefs.t.as}}, \code{\link{coefs.t.ba}},
  \code{\link{coefs.t.cb}}, and \code{\link{coefs.t.cy}} for details on
  the period-dependent model coefficients.
   
  For procedures on estimating input variables when they are not known,
  see \code{\link{Rx.calc}}, \code{\link{Rrup.calc}},
  \code{\link{dip.calc}}, \code{\link{W.calc}}, \code{\link{Ztor.calc}},
  and \code{\link{Z1.calc}}.

  For details on the spectral periods and ground motion parameters
  defined for each of the models, see \code{\link{modelPeriods}} or
  \code{\link{periods}}.  The functions \code{\link{getPeriod}} and
  \code{\link{interpolate}} provide methods for estimating spectral
  accelerations at intermediate periods between those with defined model
  coefficients. 
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

  Campbell, K. W., and Y. Bozorgnia (2007). Campbell-Bozorgnia NGA
  Ground Motion Relations for the Geometric Mean Horizontal Component of
  Peak and Spectral Ground Motion Parameters, \emph{PEER Report
    No. 2007/02}, Pacific Earthquake Engineering Research Center,
  University of California, Berkeley.
   
  Campbell, K. W., and Y. Bozorgnia (2008). NGA Ground Motion Model for
  the Geometric Mean Horizontal Component of PGA, PGV, PGD, and 5\%
  Damped Linear Elastic Response Spectra for Periods Ranging from 0.01
  to 10 s. \emph{Earthquake Spectra} \bold{24,} 139--171.

  Chiou, B. S.-J., and R. R. Youngs (2008). An NGA Model for the Average
  Horizontal Component of Peak Ground Motion and Response Spectra.
  \emph{Earthquake Spectra} \bold{24,} 173--215. 

  Kaklamanos, J., and L. G. Baise (2010).  Relationships between
  Distance Measures in Recent Ground Motion Prediction Equations.
  \emph{Seismological Research Letters} (in review).
}
\examples{

# Assumed earthquake parameters for these examples:
M <- 7
Rjb <- 50
Rrup <- 51
Vs30 <- 300
Ztor <- 2
W <- 20
dip <- 80
VsFlag <- 0
Fhw <- 0
rake <- 180
Fas <- 0
T.list <- c(0, 0.1, 0.5, 1)


# Example 1:  Comparison of output from Sa.nga and Sa.ngaR

Sa.nga(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, Fhw = Fhw, dip = dip,
       W = W, Ztor = Ztor, Vs30 = Vs30, VsFlag = VsFlag, Fas = Fas,
       epsilon = 1, T = T.list)

Sa.ngaR(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, Fhw = Fhw, dip = dip,
        W = W, Ztor = Ztor, Vs30 = Vs30, VsFlag = VsFlag, Fas = Fas,
        epsilon = 1, T = T.list)



# Example 2:  Illustration of the versatility of input

# First calculate ground motions using the known input variables
# Some of the variables (such as Z1.0 and Rx) are unknown, and will
# be calculated by the program
Sa.nga(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, Fhw = Fhw, dip = dip,
       W = W, Ztor = Ztor, Vs30 = Vs30, VsFlag = VsFlag, Fas = Fas,
       epsilon = 1, T = T.list)

# Repeat the ground motion calculation the bare minimum necessary requirements
Sa.nga(M = M, Rjb = Rjb, rake = rake, Fhw = Fhw, Vs30 = Vs30,
       epsilon = 1, T = T.list)

# Note that the style-of-faulting flag variables may be used in place
# of the rake, and that the azimuth (if known) may be used instead of Fhw
Sa.nga(M = M, Rjb = Rjb, Frv = 0, Fnm = 0, azimuth = -30, Vs30 = Vs30,
       epsilon = 1, T = T.list)



# Example 3: Generate a plot of the predicted response spectrum (and
#            uncertainty) for a hypothetical earthquake using the BA08
#            model

# Redefine T to be a vector
# We only desire T >= 0.01 for plotting
T.list <- modelPeriods(model = "BA08", positive = TRUE) 

# Ground motion calculations
ResultsMatrix <- Sa.ngaR(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake,
                         Fhw = Fhw, dip = dip, W = W, Ztor = Ztor,
                         Vs30 = Vs30, VsFlag = VsFlag, Fas = Fas,
                         epsilon = 1, T = T.list)

# To see the names of all the columns in the data frame, use the "names"
# function on a column of the matrix:
names(ResultsMatrix)

# To access individual columns of the data frame, use the "$" operator:
SaMedian <- ResultsMatrix$Y50.ba
SaPlusEpsilon <- ResultsMatrix$YplusEpsilon.ba
SaMinusEpsilon <- ResultsMatrix$YminusEpsilon.ba

# Plot
plot(T.list, SaMedian, type = "p", log = "xy", col = "blue", pch = 19,
     xlim = c(0.01, 10), ylim = c(0.001, 1), xaxs = "i", yaxs = "i",
     xlab = "Spectral Period, T [sec]", ylab = "Spectral Acceleration, Sa [g]",
     main = "BA08 Ground Motion Predictions:  Median +/- 1 SD")
points(T.list, SaMedian, pch = 19, col = "blue")
points(T.list, SaPlusEpsilon, pch = 19, col = "red")
points(T.list, SaMinusEpsilon, pch = 19, col = "red")
lines(T.list, SaMedian, lwd = 3, col = "blue")
lines(T.list, SaPlusEpsilon, lwd = 1, col = "red")
lines(T.list, SaMinusEpsilon, lwd = 1, col = "red")



# Example 4:  Generate a plot of the median response spectra for the
#             same hypothetical earthquake, comparing the different
#             NGA models

# Redefine T to be a vector
# We only desire T >= 0.01 for plotting
T.list <- modelPeriods(model = "BA08", positive = TRUE) 

# Ground motion calculations
ResultsMatrix <- Sa.ngaR(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake,
                         Fhw = Fhw, dip = dip, W = W, Ztor = Ztor,
                         Vs30 = Vs30, VsFlag = VsFlag, Fas = Fas,
                         epsilon = 1, T = T.list)

# Access individual columns of the data frame using the "$" operator:
SaAS08 <- ResultsMatrix$Y50.as
SaBA08 <- ResultsMatrix$Y50.ba
SaCB08 <- ResultsMatrix$Y50.cb
SaCY08 <- ResultsMatrix$Y50.cy

# Plot
plot(T.list, SaAS08, type = "l", log = "xy", col = "blue", pch = 19, lwd = 2,
     xlim = c(0.01, 10), ylim = c(0.001, 1), xaxs = "i", yaxs = "i",
     xlab = "Spectral Period, T [sec]", ylab = "Spectral Acceleration, Sa [g]",
     main = "Comparison of NGA Ground Motion Predictions")
lines(T.list, SaBA08, lwd = 2, col = "red")
lines(T.list, SaCB08, lwd = 2, col = "darkgreen")
lines(T.list, SaCY08, lwd = 2, col = "black")
legend(x = "bottomleft", inset = 0.02, lwd = 2, bty = "n",
       col = c("blue", "red", "darkgreen", "black"),
       legend = c("AS08", "BA08", "CB08", "CY08"))
}

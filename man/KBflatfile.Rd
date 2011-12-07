%---------------------------------------------------
\name{Example Data Analysis Using the nga Package:  KB Flatfile Data}
\docType{data}
\alias{KBflatfile}
\title{Example Earthquake Records from Recent California Earthquakes}
\description{This data set contains 1060 ground motion records from
  seven recent earthquakes recorded in California: the (1) 2003 \bold{M}
  6.5 San Simeon, (2) 2004 \bold{M} 6.0 Parkfield, (3) 2005 \bold{M} 5.2
  Anza, (4) 2007 \bold{M} 5.4 Alum Rock, (5) 2008 \bold{M} 5.4 Chino
  Hills, (6) 2010 \bold{M} 7.2 Baja, and (7) 2010 \bold{M} 5.7 Ocotillo
  earthquakes.  None of these earthquakes were present in the database
  used to develop the NGA models (the NGA flatfile), and thus these
  records were used in a blind comparison test of the models in
  Kaklamanos and Baise (2011).  The headers of this data frame are
  designed to be similar to those in the NGA flatfile; this data frame
  is termed the \dQuote{KB flatfile} (\dQuote{KB} stands for
  \dQuote{Kaklamanos and Baise}).  For further details on this dataset,
  please refer to Kaklamanos and Baise (2011) and the electronic
  supplement available at
  \url{http://www.seismosoc.org/publications/BSSA_html/bssa_101-1/2010038-esupp/index.html}.}
\usage{
data(KBflatfile)
}
\format{A dataframe containing 1060 rows and 45 columns.  For further
  details about these columns, see the documentation for the electronic
  supplement of Kaklamanos and Baise (2011).  The ground motion
  parameters at the bottom of the list are comprised of the geometric
  mean of the as-recorded horizontal components, and are presented in
  units of \emph{g}.
  \describe{
    \item{1. \code{RecNum}}{Record sequence number in the KB flatfile}
    \item{2. \code{EQID}}{Earthquake identification number in the KB
    flatfile}
    \item{3. \code{EQName}}{Earthquake name}
    \item{4. \code{Month}}{Month of the earthquake}
    \item{5. \code{Day}}{Day of the earthquake}
    \item{6. \code{Year}}{Year of the earthquake}
    \item{7. \code{StationName}}{Name of the strong-motion station}
    \item{8. \code{StaID}}{Identification number of the strong-motion
    station}
    \item{9. \code{StaNetwork}}{Network code of the strong-motion
    station}
    \item{10. \code{StaSeqNum}}{Sequence number of the strong-motion
    station in the KB flatfile}
    \item{11. \code{StaLat}}{Latitude of the strong-motion station
    (deg)}
    \item{12. \code{StaLong}}{Longitude of the strong-motion station
    (deg)}
    \item{13. \code{M}}{Moment magnitude of earthquake}
    \item{14. \code{Strike}}{Strike of the rupture plane (deg)}
    \item{15. \code{Dip}}{Dip angle of the rupture plane (deg)}
    \item{16. \code{Rake}}{Rake angle of fault movement (deg)}
    \item{17. \code{EQmechanism}}{Earthquake mechanism defined by
    rake angle}
    \item{18. \code{HypocenterLat}}{Hypocenter latitude (deg)}
    \item{19. \code{HypocenterLong}}{Hypocenter longitude (deg)}
    \item{20. \code{Zhyp}}{Depth of hypocenter (km)}
    \item{21. \code{FiniteFaultModelFlag}}{Flag variable indicating
    if a finite fault model was used (1 = Yes, 0 = No)}
    \item{22. \code{Source_of_SourceParameters}}{Reference for
    source parameters (finite fault model / moment tensor soln.)}
    \item{23. \code{Ztor}}{Depth to top of rupture (km)}
    \item{24. \code{L}}{Length of rupture plane (km)}
    \item{25. \code{W}}{Down-dip width of rupture plane (km)}
    \item{26. \code{Repi}}{Epicentral distance (km)}
    \item{27. \code{Rhyp}}{Hypocentral distance (km)}
    \item{28. \code{Rjb}}{Joyner-Boore distance (km)}
    \item{29. \code{Rrup}}{Rupture distance (km)}
    \item{30. \code{Rseis}}{Seismogenic distance (km)}
    \item{31. \code{Rx}}{Site coordinate (km)}
    \item{32. \code{Azimuth}}{Source-to-site azimuth (deg)}
    \item{33. \code{Geology}}{Surficial geologic unit}
    \item{34. \code{Vs30}}{Time-averaged shear wave velocity over a
    subsurface depth of 30 meters (m/s)}
    \item{35. \code{VsFlag}}{Vs flag variable: 1 for measured Vs, 0
    for inferred Vs}
    \item{36. \code{VsSource}}{Source of Vs / geology data}
    \item{37. \code{VsReference}}{Reference for Vs / geology data}
    \item{38. \code{GroundMotionDataSource}}{Source of ground motion
    data}
    \item{39. \code{PGA}}{Observed peak ground acceleration}
    \item{40. \code{T0.1S}}{Observed spectral acceleration (Sa) at T
    = 0.1 sec}
    \item{41. \code{T0.2S}}{Observed Sa at T = 0.2 sec}
    \item{42. \code{T0.3S}}{Observed Sa at T = 0.3 sec}
    \item{43. \code{T0.5S}}{Observed Sa at T = 0.5 sec}
    \item{44. \code{T1.0S}}{Observed Sa at T = 1.0 sec}
    \item{45. \code{T2.0S}}{Observed Sa at T = 2.0 sec}
  }
}
\source{Electronic supplement of Kaklamanos and Baise (2011), available
  at \url{http://www.seismosoc.org/publications/BSSA_html/bssa_101-1/2010038-esupp/index.html}}
\seealso{
  \code{\link{Sa}}, \code{\link{Sa.nga}}
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

  Kaklamanos, J., and L. G. Baise (2011).  Model Validations and
  Comparisons of the Next Generation Attenuation of Ground Motions
  (NGA-West) Project.
  \emph{Bulletin of the Seismological Society of America}, \bold{101,}
  160--175. 

  Kaklamanos, J., L. G. Baise, and D. M. Boore (2011).  Estimating
  Unknown Input Parameters when Implementing the NGA Ground-Motion
  Prediction Equations in Engineering Practice.
  \emph{Earthquake Spectra} \bold{27,} 1219--1235. 
}
\keyword{datasets}
\examples{

# Load dataset (this command MUST be typed prior to using the dataset)
data(KBflatfile)

# See the column names of the dataset
names(KBflatfile)


######################################################################
# Example 1:  Generate a plot of observed versus predicted response
#             spectrum for a ground motion record in the database

# Use Rec No. 824, the first ground motion record for the Baja
# earthquake of 2010 listed in the dataset.


# Read data from the 824th row
# Only read columns that are necessary for ground motion calculations

# Input variables:
n <- 824
M <- KBflatfile$M[n]
dip <- KBflatfile$Dip[n]
rake <- KBflatfile$Rake[n]
Ztor <- KBflatfile$Ztor[n]
W <- KBflatfile$W[n]
Rjb <- KBflatfile$Rjb[n]
Rrup <- KBflatfile$Rrup[n]
Rx <- KBflatfile$Rx[n]
azimuth <- KBflatfile$Azimuth[n]
Vs30 <- KBflatfile$Vs30[n]
# VsFlag is not read, because it is only necessary for standard
# deviation calculations (i.e., epsilon != 0)

# Observed response spectral values:
PGA <- KBflatfile$PGA[n]
Sa0.1 <- KBflatfile$T0.1S[n]
Sa0.2 <- KBflatfile$T0.2S[n]
Sa0.3 <- KBflatfile$T0.3S[n]
Sa0.5 <- KBflatfile$T0.5S[n]
Sa1.0 <- KBflatfile$T1.0S[n]
Sa2.0 <- KBflatfile$T2.0S[n]

# Vectorize the observed spectral acceleration and corresponding periods
# NOTE: Observed PGA is assumed to have a spectral period of T = 0.01 sec
T.obs <- c(0.01, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0)
Sa.obs <- c(PGA, Sa0.1, Sa0.2, Sa0.3, Sa0.5, Sa1.0, Sa2.0)

# Define the periods at which ground motion calculations will be performed
# NOTE: the same could be achieved by using the function call
# modelPeriods(model = "AS08", positive = TRUE).
T.list <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20,
            0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0,
            5.0, 7.5, 10.0)

# Compute ground motion predictions
ResultsMatrix <- Sa.nga(M = M, dip = dip, rake = rake, Ztor = Ztor,
                        W = W, Rjb = Rjb, Rrup = Rrup, Rx = Rx,
                        azimuth = azimuth, Vs30 = Vs30,
                        epsilon = 1, T = T.list)

# Access individual columns of the data frame using the "$" operator:
SaAS08 <- ResultsMatrix$Y50.as
SaBA08 <- ResultsMatrix$Y50.ba
SaCB08 <- ResultsMatrix$Y50.cb
SaCY08 <- ResultsMatrix$Y50.cy

# Plot the results
plot(T.obs, Sa.obs, type = "p", log = "xy", col = "black", pch = 19, lwd = 4,
     xlim = c(0.01, 10), ylim = c(0.001, 1), xaxs = "i", yaxs = "i",
     xlab = "Spectral Period, T [sec]", ylab = "Spectral Acceleration, Sa [g]",
     main = paste("Comparison of NGA Ground Motion Predictions:", "\n",
                  "Record No. 824; Baja Earthquake of 2010"))
lines(T.list, SaAS08, lwd = 2, col = "blue")
lines(T.list, SaBA08, lwd = 2, col = "red")
lines(T.list, SaCB08, lwd = 2, col = "darkgreen")
lines(T.list, SaCY08, lwd = 2, col = "purple")
legend(x = "bottomleft", inset = 0.02, lwd = c(-1,2,2,2,2),
       lty = c(-1,1,1,1,1), pch = c(19,-1,-1,-1,-1), bty = "n",
       col = c("black", "blue", "red", "darkgreen", "purple"),
       legend = c("Observed", "AS08", "BA08", "CB08", "CY08"))



############################################################################
# Example 2:  Generate a plot of peak ground acceleration versus distance
#             for the Chino Hills earthquake of 2008

# The relevant ground motion records are present in rows 447 to 823 of
# the KB flatfile.  Note that because a finite fault model was not
# developed for this earthquake, some of the source and distance
# parameters are unknown and must be estimated by the program.

# Read data
start <- 447
end <- 823
n <- seq(from = start, to = end, by = 1)
M <- KBflatfile$M[n]
rake <- KBflatfile$Rake[n]
dip <- KBflatfile$Dip[n]
Zhyp <- KBflatfile$Zhyp[n]
Repi <- KBflatfile$Repi[n]
Vs30 <- KBflatfile$Vs30[n]
PGA <- KBflatfile$PGA[n]


# Generate NGA ground motion predictions versus distance

# Extract source parameters from the vectors.
# These are constants for each of the 337 ground motion records in the
# subset, so it does not matter which row we extract.
M.value <- M[1]
rake.value <- rake[1]
dip.value <- dip[1]
Zhyp.value <- Zhyp[1]

# Assume an average Vs30 for the purpose of drawing the graphs
Vs30.value <- mean(Vs30)

# Assume site is on footwall (since the earthquake is low-magnitude,
# the hanging wall effects are not likely to be significant).
Fhw <- 0

# First, illustrate the calculation for one point:
ResultsMatrix1 <- Sa.nga(M = M.value, dip = dip.value, rake = rake.value,
                         Rjb = 0, Fhw = 0, Vs30 = Vs30.value, epsilon = 0, T = 0)

# Generate a vector of Rjb values from 0 to 200 to be used for
# plotting and for generating ground motion predictions
Rjb.plot <- seq(from = 0, to = 200, by = 4)

# Perform ground motion calculations for all points.
# Define ResultsMatrix2; use the column names of ResultsMatrix1
ResultsMatrix2 <- matrix(nrow = length(Rjb.plot), ncol = length(ResultsMatrix1))
ResultsMatrix2 <- as.data.frame(ResultsMatrix2)
names(ResultsMatrix2) <- names(ResultsMatrix1)
# It is necessary to place the calculation in a loop since we are varying Rjb.
for(i in 1:length(Rjb.plot)){
  ResultsMatrix2[i,] <- Sa.nga(M = M.value, dip = dip.value, rake = rake.value,
                               Rjb = Rjb.plot[i], Fhw = 0, Vs30 = Vs30.value,
                               epsilon = 0, T = 0)
}

# Access individual columns of the data frame using the "$" operator:
pgaAS08 <- ResultsMatrix2$Y50.as
pgaBA08 <- ResultsMatrix2$Y50.ba
pgaCB08 <- ResultsMatrix2$Y50.cb
pgaCY08 <- ResultsMatrix2$Y50.cy


# Plot the results.
# For the purpose of generating the plot, Repi is used in place of
# Rjb. For small-magnitude events, the area of fault rupture is
# small, and the assumption Repi = Rjb is not unreasonable.
plot(Repi, PGA, type = "p", log = "y", pch = 1,
     xlab = "Joyner-Boore Distance, Rjb [km]",
     ylab = "Peak Ground Acceleration, PGA [g]",
     main = paste("Comparison of NGA Ground Motion Predictions:", "\n",
     "PGA versus Rjb for the Chino Hills Earthquake of 2008"))
lines(Rjb.plot, pgaAS08, lwd = 2, col = "blue")
lines(Rjb.plot, pgaBA08, lwd = 2, col = "red")
lines(Rjb.plot, pgaCB08, lwd = 2, col = "darkgreen")
lines(Rjb.plot, pgaCY08, lwd = 2, col = "purple")
legend(x = "bottomleft", inset = 0.02, pch = c(1,-1,-1,-1,-1),
       lwd = c(-1,2,2,2,2), lty = c(-1,1,1,1,1), bty = "n",
       col = c("black", "blue", "red", "darkgreen", "purple"),
       legend = c("Observed", "AS08", "BA08", "CB08", "CY08"))



########################################################################
# Example 3:  Tabulate predicted versus observed peak ground acceleration
#             for the ground motion records of the San Simeon earthquake
#             of 2003

# The relevant ground motion records are present in rows 1 to 30 of
# the KB flatfile:

# Read data
start <- 1
end <- 30
n <- seq(from = start, to = end, by = 1)
M <- KBflatfile$M[n]
dip <- KBflatfile$Dip[n]
rake <- KBflatfile$Rake[n]
Ztor <- KBflatfile$Ztor[n]
W <- KBflatfile$W[n]
Rjb <- KBflatfile$Rjb[n]
Rrup <- KBflatfile$Rrup[n]
Rx <- KBflatfile$Rx[n]
azimuth <- KBflatfile$Azimuth[n]
Vs30 <- KBflatfile$Vs30[n]
PGA.obs <- KBflatfile$PGA[n]

# Create matrices to store the calculated values
pgaAS08 <- matrix(nrow = length(n), ncol = 1)
pgaBA08 <- matrix(nrow = length(n), ncol = 1)
pgaCB08 <- matrix(nrow = length(n), ncol = 1)
pgaCY08 <- matrix(nrow = length(n), ncol = 1)

# Perform ground motion predictions
for(i in 1:length(n)){
  ResultsMatrix <- Sa.nga(M = M[i], dip = dip[i], rake = rake[i],
                          Ztor = Ztor[i], W = W[i], Rjb = Rjb[i],
                          Rrup = Rrup[i], Rx = Rx[i],
                          azimuth = azimuth[i], Vs30 = Vs30[i],
                          epsilon = 0, T = 0)
  pgaAS08[i] <- ResultsMatrix$Y50.as
  pgaBA08[i] <- ResultsMatrix$Y50.ba
  pgaCB08[i] <- ResultsMatrix$Y50.cb
  pgaCY08[i] <- ResultsMatrix$Y50.cy
}

# Combine the results into a data frame
Ex3 <- cbind(PGA.obs, pgaAS08, pgaBA08, pgaCB08, pgaCY08)
colnames(Ex3) <- c("pgaObs", "pgaAS08", "pgaBA08", "pgaCB08", "pgaCY08")

# Display results
Ex3

# You could now use a function such as "write.csv" or "write.table" to export Ex3




########################################################################
# Example 4:  Generate matrices of median predicted response spectra
#             for the San Simeon earthquake of 2003

# The relevant ground motion records are present in rows 1 to 30 of
# the KB flatfile (same as example 3)

# Read data
start <- 1
end <- 30
n <- seq(from = start, to = end, by = 1)
M <- KBflatfile$M[n]
dip <- KBflatfile$Dip[n]
rake <- KBflatfile$Rake[n]
Ztor <- KBflatfile$Ztor[n]
W <- KBflatfile$W[n]
Rjb <- KBflatfile$Rjb[n]
Rrup <- KBflatfile$Rrup[n]
Rx <- KBflatfile$Rx[n]
azimuth <- KBflatfile$Azimuth[n]
Vs30 <- KBflatfile$Vs30[n]
VsFlag <- KBflatfile$VsFlag[n]

# Create matrix of observed response spectra

# Read observed data as vectors
PGA <- KBflatfile$PGA[n]
Sa0.1 <- KBflatfile$T0.1S[n]
Sa0.2 <- KBflatfile$T0.2S[n]
Sa0.3 <- KBflatfile$T0.3S[n]
Sa0.5 <- KBflatfile$T0.5S[n]
Sa1.0 <- KBflatfile$T1.0S[n]
Sa2.0 <- KBflatfile$T2.0S[n]

# Combine the individual vectors into a matrix using the
# "cbind" function
Obs <- cbind(PGA, Sa0.1, Sa0.2, Sa0.3, Sa0.5, Sa1.0, Sa2.0)

# Periods for analysis
T.list <- c(0, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0)


# Create matrices to store the calculated values
PredAS08 <- matrix(nrow = length(n), ncol = length(T.list))
PredBA08 <- matrix(nrow = length(n), ncol = length(T.list))
PredCB08 <- matrix(nrow = length(n), ncol = length(T.list))
PredCY08 <- matrix(nrow = length(n), ncol = length(T.list))
colnames(PredAS08) <- colnames(Obs)
colnames(PredBA08) <- colnames(Obs)
colnames(PredCB08) <- colnames(Obs)
colnames(PredCY08) <- colnames(Obs)


# Perform ground motion predictions (this example illustrates the
# use of the individual functions Sa.as, Sa.ba, Sa.cb, and Sa.cy,
# which are faster and generate less output than Sa.nga)

# Ground motion calculations
for(i in 1:length(n)){
  PredAS08[i,] <- Sa.as(M = M[i], dip = dip[i], rake = rake[i],
                        Ztor = Ztor[i], W = W[i], Rjb = Rjb[i],
                        Rrup = Rrup[i], Rx = Rx[i],
                        azimuth = azimuth[i], Vs30 = Vs30[i],
                        VsFlag = VsFlag[i], Fas = 0, epsilon = 0,
                        T = T.list)
  PredBA08[i,] <- Sa.ba(M = M[i], rake = rake[i], Rjb = Rjb[i],
                        Vs30 = Vs30[i], epsilon = 0, T = T.list)
  PredCB08[i,] <- Sa.cb(M = M[i], dip = dip[i], rake = rake[i],
                        Ztor = Ztor[i], Rjb = Rjb[i], Rrup = Rrup[i],
                        Vs30 = Vs30[i], epsilon = 0, T = T.list)
  PredCY08[i,] <- Sa.cy(M = M[i], dip = dip[i], rake = rake[i],
                        Ztor = Ztor[i], W = W[i], Rjb = Rjb[i],
                        Rrup = Rrup[i], Rx = Rx[i],
                        azimuth = azimuth[i], Vs30 = Vs30[i],
                        VsFlag = VsFlag[i], AS = 0, epsilon = 0,
                        T = T.list)
}

# Display results
Obs
PredAS08
PredBA08
PredCB08
PredCY08

# Now each of the matrices may be used in later calculations, or
# written to a text or csv file.

}

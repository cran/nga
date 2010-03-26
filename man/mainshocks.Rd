%---------------------------------------------------
\name{Example Data Analysis Using the nga Package:  Mainshocks Data}
\docType{data}
\alias{mainshocks}
\title{Example Earthquake Records from the NGA Flatfile}
\description{This data set contains 1231 ground motion records from
  mainshocks in the NGA flatfile.  These records, which are a subset of
  the 3551 ground motion records in the NGA flatfile, meet the criteria of
  applicability of the NGA models, have finite fault models (i.e.,
  populated values for the key earthquake parameters), are
  representative of free-field conditions, and do not have identified
  problems.  For further details on the selection process of this
  dataset, see Kaklamanos and Baise (2010a).}
\usage{
data(mainshocks)
}
\format{A dataframe containing 1231 rows and 46 columns.  For further
  details about these columns, see the documentation for the NGA
  flatfile.
  \tabular{lll}{
    1. \tab RecNum \tab Record sequence number in the NGA flatfile  \cr
    2. \tab EQID \tab Earthquake identification number in the NGA flatfile  \cr
    3. \tab EQName \tab Earthquake name \cr
    4. \tab Year \tab Year of the earthquake \cr
    5. \tab MoDay \tab Month and day of the earthquake \cr
    6. \tab StationName \tab Name of the strong-motion station \cr
    7. \tab StaSeqNum \tab Sequence number of the strong-motion station \cr
    8. \tab StaID \tab Identification number of the strong-motion
    station (provided by data provider) \cr
    9. \tab M \tab Moment magnitude of earthquake \cr
    10. \tab Dip \tab Dip angle of the rupture plane (deg) \cr
    11. \tab Rake \tab Rake angle of fault movement (deg) \cr
    12. \tab Ztor \tab Depth to top of rupture (km) \cr
    13. \tab W \tab Down-dip width of rupture plane (km) \cr
    14. \tab Rjb \tab Joyner-Boore distance (km):  the horizontal
    distance to the surface projection of the rupture plane \cr
    15. \tab Rrup \tab Rupture distance (km):  the closest distance to the
    rupture plane \cr
    16. \tab Azimuth \tab Source-to-site azimuth (deg) \cr
    17. \tab Vs30 \tab Time-averaged shear wave velocity over a subsurface
    depth of 30 meters (m/s) \cr
    18. \tab VsFlag \tab Vs Flag variable: 1 for measured Vs, 0 for inferred Vs \cr
    19. \tab Z1.0 \tab Depth to Vs = 1.0 km/s (m) \cr
    20. \tab Z1.5 \tab Depth to Vs = 1.5 km/s (m) \cr
    21. \tab Z2.5 \tab Depth to Vs = 2.5 km/s (m) \cr
    22. \tab PGA \tab Observed peak ground acceleration (g) \cr
    23. \tab PGV \tab Observed peak ground velocity (cm/s) \cr
    24. \tab PGD \tab Observed peak ground displacement (cm) \cr
    25. \tab T0.010S \tab Observed spectral acceleration (Sa) at T =
    0.010 sec \cr
    26. \tab T0.020S \tab Observed Sa at T = 0.020 sec \cr
    27. \tab T0.030S \tab Observed Sa at T = 0.030 sec \cr
    28. \tab T0.040S \tab Observed Sa at T = 0.040 sec \cr
    29. \tab T0.050S \tab Observed Sa at T = 0.050 sec \cr
    30. \tab T0.075S \tab Observed Sa at T = 0.750 sec \cr
    31. \tab T0.100S \tab Observed Sa at T = 0.100 sec \cr
    32. \tab T0.150S \tab Observed Sa at T = 0.150 sec \cr
    33. \tab T0.200S \tab Observed Sa at T = 0.200 sec \cr
    34. \tab T0.250S \tab Observed Sa at T = 0.250 sec \cr
    35. \tab T0.300S \tab Observed Sa at T = 0.300 sec \cr
    36. \tab T0.400S \tab Observed Sa at T = 0.400 sec \cr
    37. \tab T0.500S \tab Observed Sa at T = 0.500 sec \cr
    38. \tab T0.750S \tab Observed Sa at T = 0.750 sec \cr
    39. \tab T1.000S \tab Observed Sa at T = 1.000 sec \cr
    40. \tab T1.500S \tab Observed Sa at T = 1.500 sec \cr
    41. \tab T2.000S \tab Observed Sa at T = 2.000 sec \cr
    42. \tab T3.000S \tab Observed Sa at T = 3.000 sec \cr
    43. \tab T4.000S \tab Observed Sa at T = 4.000 sec \cr
    44. \tab T5.000S \tab Observed Sa at T = 5.000 sec \cr
    45. \tab T7.500S \tab Observed Sa at T = 7.500 sec \cr
    46. \tab T10.000S \tab Observed Sa at T = 10.000 sec \cr
     \tab \tab
  }
}
\source{Pacific Earthquake Engineering Research Center. Next Generation
  Attenuation of Ground Motions (NGA) Project,
  \url{http://peer.berkeley.edu/products/nga_project.html}}
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

  Kaklamanos, J., and L. G. Baise (2010a).  Model Validations and
  Comparisons of the Next Generation Attenuation of Ground Motions
  (NGA-West) Project.
  \emph{Bulletin of the Seismological Society of America} (in review).

  Kaklamanos, J., and L. G. Baise (2010b).  Relationships between
  Distance Measures in Recent Ground Motion Prediction Equations.
  \emph{Seismological Research Letters} (in review).  
}
\keyword{datasets}
\examples{

# Load dataset
data(mainshocks)

# See the column names of the dataset
names(mainshocks)


# Example 1:  Generate a plot of observed versus predicted response
#             spectrum for a ground motion record in the database

# Use Rec No. 731, the first ground motion record for the Loma Prieta
# EQ of 1989 listed in the dataset.  This record corresponds to
# the 364th row in the mainshocks dataset.

# Read data from the 364th row
# Input variables:
n <- 364
M <- mainshocks$M[n]
dip <- mainshocks$Dip[n]
rake <- mainshocks$Rake[n]
Ztor <- mainshocks$Ztor[n]
W <- mainshocks$W[n]
Rjb <- mainshocks$Rjb[n]
Rrup <- mainshocks$Rrup[n]
azimuth <- mainshocks$Azimuth[n]
Vs30 <- mainshocks$Vs30[n]
VsFlag <- mainshocks$VsFlag[n]
Z1.0 <- mainshocks$Z1.0[n]
Z1.5 <- mainshocks$Z1.5[n]
Z2.5 <- mainshocks$Z2.5[n]
# Observed response spectral values:
Sa0.01 <- mainshocks$T0.010S[n]
Sa0.02 <- mainshocks$T0.020S[n]
Sa0.03 <- mainshocks$T0.030S[n]
Sa0.04 <- mainshocks$T0.040S[n]
Sa0.05 <- mainshocks$T0.050S[n]
Sa0.075 <- mainshocks$T0.075S[n]
Sa0.10 <- mainshocks$T0.100S[n]
Sa0.15 <- mainshocks$T0.150S[n]
Sa0.20 <- mainshocks$T0.200S[n]
Sa0.25 <- mainshocks$T0.250S[n]
Sa0.30 <- mainshocks$T0.300S[n]
Sa0.40 <- mainshocks$T0.400S[n]
Sa0.50 <- mainshocks$T0.500S[n]
Sa0.75 <- mainshocks$T0.750S[n]
Sa1.00 <- mainshocks$T1.000S[n]
Sa1.50 <- mainshocks$T1.500S[n]
Sa2.00 <- mainshocks$T2.000S[n]
Sa3.00 <- mainshocks$T3.000S[n]
Sa4.00 <- mainshocks$T4.000S[n]
Sa5.00 <- mainshocks$T5.000S[n]
Sa7.50 <- mainshocks$T7.500S[n]
Sa10.0 <- mainshocks$T10.000S[n]

# Vectorize the spectral period and observed spectral acceleration
# NOTE: the same could be achieved by using the function
# modelPeriods(model = "AS08", positive = TRUE).
T.list <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20,
            0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0,
            5.0, 7.5, 10.0)
Sa.obs <- c(Sa0.01, Sa0.02, Sa0.03, Sa0.04, Sa0.05, Sa0.075, Sa0.10,
            Sa0.15, Sa0.20, Sa0.25, Sa0.30, Sa0.40, Sa0.50, Sa0.75,
            Sa1.00, Sa1.50, Sa2.00, Sa3.00, Sa4.00, Sa5.00, Sa7.50, Sa10.0)

# Compute ground motion predictions
ResultsMatrix <- Sa.ngaR(M = M, dip = dip, rake = rake, Ztor = Ztor,
                         W = W, Rjb = Rjb, Rrup = Rrup, azimuth = azimuth,
                         Vs30 = Vs30, Z1.0 = Z1.0, Z1.5 = Z1.5, Z2.5 = Z2.5,
                         VsFlag = VsFlag, epsilon = 1, T = T.list)

# Access individual columns of the data frame using the "$" operator:
SaAS08 <- ResultsMatrix$Y50.as
SaBA08 <- ResultsMatrix$Y50.ba
SaCB08 <- ResultsMatrix$Y50.cb
SaCY08 <- ResultsMatrix$Y50.cy

# Plot
plot(T.list, Sa.obs, type = "l", log = "xy", col = "black", pch = 19, lwd = 4,
     xlim = c(0.01, 10), ylim = c(0.001, 1), xaxs = "i", yaxs = "i",
     xlab = "Spectral Period, T [sec]", ylab = "Spectral Acceleration, Sa [g]",
     main = paste("Comparison of NGA Ground Motion Predictions:", "\n",
                  "Flatfile Record No. 731; Loma Prieta Earthquake of 1989"))
lines(T.list, SaAS08, lwd = 2, col = "blue")
lines(T.list, SaBA08, lwd = 2, col = "red")
lines(T.list, SaCB08, lwd = 2, col = "darkgreen")
lines(T.list, SaCY08, lwd = 2, col = "purple")
legend(x = "bottomleft", inset = 0.02, lwd = c(4,2,2,2,2), bty = "n",
       col = c("black", "blue", "red", "darkgreen", "purple"),
       legend = c("Observed", "AS08", "BA08", "CB08", "CY08"))




# Example 2:  Generate a plot of peak ground acceleration versus distance
#             for the Loma Prieta earthquake of 1989

# The relevant ground motion records are present in rows 364 to 430 of
# the mainshocks database:

# Read data
start <- 364
end <- 430
n <- seq(from = start, to = end, by = 1)
M <- mainshocks$M[n]
dip <- mainshocks$Dip[n]
rake <- mainshocks$Rake[n]
Ztor <- mainshocks$Ztor[n]
W <- mainshocks$W[n]
Rjb <- mainshocks$Rjb[n]
Rrup <- mainshocks$Rrup[n]
azimuth <- mainshocks$Azimuth[n]
Vs30 <- mainshocks$Vs30[n]
VsFlag <- mainshocks$VsFlag[n]
Z1.0 <- mainshocks$Z1.0[n]
Z1.5 <- mainshocks$Z1.5[n]
Z2.5 <- mainshocks$Z2.5[n]
PGA <- mainshocks$PGA[n]


# Generate NGA ground motion predictions versus distance

# Extract source parameters from the vectors
#   These are constants for each of the 67 ground motion records in the
#   subset, so it does not matter which value we extract
M.value <- M[1]
dip.value <- dip[1]
rake.value <- rake[1]
Ztor.value <- Ztor[1]
W.value <- W[1]

# Assume an average Vs30 for the purpose of drawing the graphs
Vs30.value <- mean(Vs30)

# Generate predictions for value of Rjb from 0 to 120
Rjb.plot <- seq(from = 0, to = 120, by = 1)
Fhw <- 0   # Assume site is on footwall

# First illustrate the calculation for one point
ResultsMatrix1 <- Sa.ngaR(M = M.value, dip = dip.value, rake = rake.value,
                          Ztor = Ztor.value, W = W.value, Rjb = 0, Fhw = 0,
			  Vs30 = Vs30.value, VsFlag = 0, epsilon = 0, T = 0)		  

# Perform ground motion calculations for all points
#   It is necessary to place this in a loop since we are varying Rjb			  
ResultsMatrix2 <- matrix(nrow = length(Rjb), ncol = dim(ResultsMatrix1)[2])
ResultsMatrix2 <- as.data.frame(ResultsMatrix2)
names(ResultsMatrix2) <- names(ResultsMatrix1)
for(i in 1:length(Rjb.plot)){
  ResultsMatrix2[i,] <- Sa.ngaR(M = M.value, dip = dip.value, rake = rake.value,
                                Ztor = Ztor.value, W = W.value, Rjb = Rjb.plot[i],
                                Fhw = 0, Vs30 = Vs30.value, VsFlag = 0,
                                epsilon = 0, T = 0)
}

# Access individual columns of the data frame using the "$" operator:
pgaAS08 <- ResultsMatrix2$Y50.as
pgaBA08 <- ResultsMatrix2$Y50.ba
pgaCB08 <- ResultsMatrix2$Y50.cb
pgaCY08 <- ResultsMatrix2$Y50.cy

# Plot
plot(Rjb, PGA, type = "p", log = "y", pch = 1,
     xlim = c(0, 120), ylim = c(0.03, 1), xaxs = "i", yaxs = "i",
     xlab = "Joyner-Boore Distance, Rjb [km]",
     ylab = "Peak Ground Acceleration, PGA [g]",
     main = paste("Comparison of NGA Ground Motion Predictions", "\n",
     "PGA versus Rjb for the Loma Prieta Earthquake of 1989"))
lines(Rjb.plot, pgaAS08, lwd = 2, col = "blue")
lines(Rjb.plot, pgaBA08, lwd = 2, col = "red")
lines(Rjb.plot, pgaCB08, lwd = 2, col = "darkgreen")
lines(Rjb.plot, pgaCY08, lwd = 2, col = "purple")
legend(x = "bottomleft", inset = 0.02, pch = c(1,-1,-1,-1,-1),
       lwd = c(-1,2,2,2,2), lty = c(-1,1,1,1,1), bty = "n",
       col = c("black", "blue", "red", "darkgreen", "purple"),
       legend = c("Observed", "AS08", "BA08", "CB08", "CY08"))



# Example 3:  Tabulate predicted versus observed peak ground velocity
#             for the ground motion records in the Loma Prieta Earthquake

# Read data (same as example 2)
start <- 364
end <- 430
n <- seq(from = start, to = end, by = 1)
M <- mainshocks$M[n]
dip <- mainshocks$Dip[n]
rake <- mainshocks$Rake[n]
Ztor <- mainshocks$Ztor[n]
W <- mainshocks$W[n]
Rjb <- mainshocks$Rjb[n]
Rrup <- mainshocks$Rrup[n]
azimuth <- mainshocks$Azimuth[n]
Vs30 <- mainshocks$Vs30[n]
VsFlag <- mainshocks$VsFlag[n]
Z1.0 <- mainshocks$Z1.0[n]
Z1.5 <- mainshocks$Z1.5[n]
Z2.5 <- mainshocks$Z2.5[n]
PGV.obs <- mainshocks$PGV[n]

# Create matrices to store data
pgvAS08 <- matrix(nrow = length(n), ncol = 1)
pgvBA08 <- matrix(nrow = length(n), ncol = 1)
pgvCB08 <- matrix(nrow = length(n), ncol = 1)
pgvCY08 <- matrix(nrow = length(n), ncol = 1)

# Perform ground motion predictions
for(i in 1:length(n)){
  ResultsMatrix <- Sa.ngaR(M = M[i], dip = dip[i], rake = rake[i],
                           Ztor = Ztor[i], W = W[i], Rjb = Rjb[i],
                           Rrup = Rrup[i], azimuth = azimuth[i],
                           Vs30 = Vs30[i], Z1.0 = Z1.0[i], Z1.5 = Z1.5[i],
                           Z2.5 = Z2.5[i], VsFlag = VsFlag[i],
                           epsilon = 0, T = -1)
  pgvAS08[i] <- ResultsMatrix$Y50.as
  pgvBA08[i] <- ResultsMatrix$Y50.ba
  pgvCB08[i] <- ResultsMatrix$Y50.cb
  pgvCY08[i] <- ResultsMatrix$Y50.cy
}

# Display results
PGV.obs
pgvAS08
pgvBA08
pgvCB08
pgvCY08

}


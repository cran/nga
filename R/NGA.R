# GROUND MOTION COMPUTATIONS FOR ALL NGA MODELS
# Abrahamson and Silva (2008)
# Boore and Atkinson (2008)
# Campbell and Bozorgnia (2008)
# Chiou and Youngs (2008)


# FUNCTION I:  Sa.nga

# OUTLINE OF CODE
# 1. Check Input Parameters
# 2. Obtain Estimates of Unspecified Input Parameters
# 3. Calculation of Ground Motions
# 4. Return List of Inputs and Outputs


# INPUT PARAMETERS:
#   M = Moment magnitude
#   Rjb = Joyner-Boore distance (km)
#   Rrup = Rupture distance (km)
#   Rx = Site coordinate (km)
#   rake = Rake angle of fault movement (deg)
#   Frv = Reverse style-of-faulting flag (1 for reverse faulting,
#         0 otherwise)
#   Fnm = Normal style-of-faulting flag (1 for normal faulting,
#         0 otherwise)
#   Fhw = Hanging wall flag; 1 for site on the hanging wall side of the fault,
#         and 0 otherwise
#   dip = Fault dip (deg)
#   W = Down-dip rupture width (km)
#   Ztor = Depth to top of rupture (km)
#   Vs30 = Time-averaged shear wave velocity over 30 m subsurface depth  (m/sec)
#   Z1.0 = Depth to Vs = 1.0 km/sec  (m)
#   Z1.5 = Depth to Vs = 1.5 km/sec  (m)
#   Z2.5 = Depth to Vs = 2.5 km/sec  (km)
#   VsFlag = Flag variable indicating how Vs30 is obtained
#            (1 if measured, 0 if estimated/inferred)
#   Fas = Aftershock flag (1 for aftershocks, 0 for mainshocks)
#   Zhyp = hypocentral depth (km)
#   azimuth = source-to-site azimuth (deg); see Figure 2 in Kaklamanos and Baise (2010)
#   U = Unspecified style-of-faulting flag for the BA08 model;
#       equal to 1 if the user wishes to perform a generic ground motion
#       calculation when the style of faulting is unspecified, and 0 otherwise.
#   arb =  Flag variable indicating the method of determining aleatory uncertainty
#          for the CB08 model; equal to "1" if the standard deviation should be
#          calculated for the arbitrary horizontal component; 0 if the standard
#          deviation should be calculated for the geometric mean of LnY
#   epsilon = number of standard deviations to be considered in the calculations
#   T = Spectral period, sec (0 for PGA; -1 for PGV; -2 for PGD)

# OUTPUTS:  A list composed of the following elements:
#
#  INPUT VARIABLES:
#   T = Spectral period, sec [input]
#   epsilon = number of standard deviations considered in the calculations [input]
#   M = Moment magnitude [input]
#   Rjb = Joyner-Boore distance (km) [input]
#   Rrup.in = Rupture distance (km) [input]
#   Rrup.out = Rupture distance (km) [calculated if Rrup.in is not specified]
#   Rx.in = Site coordinate (km) [input]
#   Rx.out = Site coordinate (km) [calculated if Rx.in is not specified]
#   rake.in = Rake angle of fault movement (deg) [input]
#   rake.out = Rake angle of fault movement (deg) [calculated if rake.in is not specified]
#   Frv = Reverse style-of-faulting flag [input]
#   Fnm1 = Normal style-of-faulting flag for AB08 and CY08
#   Fnm2 = Normal style-of-faulting flag for BA08 and CB08
#   Fhw = Hanging wall flag
#   dip.in = Fault dip angle (deg) [input]
#   dip.out = Fault dip angle (deg) [calculated if dip.in is not specified]
#   W.in = Down-dip rupture width (km) [input]
#   W.out = Down-dip rupture width (km) [calculated if W.in is not specified]
#   Ztor.in = Depth to top of rupture (km) [input]
#   Ztor.out = Depth to top of rupture (km) [calculated if Ztor.in is not specified]
#   Vs30 = Time-averaged shear wave velocity over 30 m subsurface depth  (m/sec) [input]
#   Z1.0in = Depth to Vs of 1.0 km/sec  (m) [input]
#   Z1.0as = Depth to Vs of 1.0 km/sec  (m) [calculated for use in AS08 model]
#   Z1.0cy = Depth to Vs of 1.0 km/sec  (m) [calculated for use in CY08 model]
#   Z1.5in = Depth to Vs of 1.5 km/sec  (m) [input]
#   Z2.5in = Depth to Vs of 2.5 km/sec  (km) [input]
#   Z2.5out = Depth to Vs of 2.5 km/sec  (km) [calculated from Z1.0 for use in CB08 model]
#   VsFlag = Flag variable indicating how Vs30 is obtained [input]
#   Fas = Aftershock flag [input]
#   Zhyp = hypocentral depth (km) [input]
#   azimuth.in = source-to-site azimuth (deg) [input]
#   azimuth.out = source-to-site azimuth (deg) [calculated if azimuth.in is not specified]
#   U = Unspecified style-of-faulting flag for the BA08 model [input]
#   arb = "1" if the standard deviation should be calculated for the arbitrary
#         horizontal component; "0" if the standard deviation should be calculated
#         for the geometric mean of LnY [input]
#
#  OUTPUT VARIABLES:
#   NOTE:  "Y" refers to the ground motion parameter of interest, which can be
#             Sa = Spectral acceleration (g);
#             PGA = Peak ground acceleration (g), calculated by evaluating Sa at T = 0;
#             PGV = Peak ground velocity (cm/sec), calculated by evaluating Sa at T = -1; or
#             PGD = Peak ground displacement (cm), calculated by evaluating Sa at T = -2
#                   (CB08 only).
#          "sd" refers to the standard deviation of the ground motion estimate, which
#          is presented in log space.
#   AS08 Model:
#    Y50.as = Median ground motion estimate using AS08 (epsilon = 0)
#    YplusEpsilon.meas.as = Upper bound estimate of ground motion using AS08,
#                           for measured Vs30 (VsFlag = 1)
#    YplusEpsilon.est.as = Upper bound estimate of ground motion using AS08,
#                          for estimated Vs30 (VsFlag = 0)
#    YminusEpsilon.meas.as = Lower bound estimate of ground motion using AS08,
#                            for measured Vs30 (VsFlag = 1)
#    YminusEpsilon.est.as = Lower bound estimate of ground motion using AS08,
#                           for estimated Vs30 (VsFlag = 0)
#    sdMeas.as = total standard deviation using AS08, for measured Vs30 (VsFlag = 1)
#    sdEst.as = total standard deviation using AS08, for estimated Vs30 (VsFlag = 0)
#   BA08 Model:
#    Y50M.ba = Median ground motion estimate using BA08,
#              when the fault type is specified (U = 0)
#    Y50U.ba = Median ground motion estimate using BA08,
#              when the fault type is unspecified (U = 1)
#    YplusEpsilon.M.ba = Upper bound estimate of ground motion using BA08,
#                        when the fault type is specified (U = 0)
#    YplusEpsilon.U.ba = Upper bound estimate of ground motion using BA08,
#                        when the fault type is unspecified (U = 1)
#    YminusEpsilon.M.ba = Lower bound estimate of ground motion using BA08,
#                         when the fault type is specified (U = 0)
#    YminusEpsilon.U.ba = Lower bound estimate of ground motion using BA08,
#                         when the fault type is unspecified (U = 1)
#    sdM.ba = total standard deviation using BA08,
#             when the fault type is specified (U = 0)
#    sdU.ba = total standard deviation using BA08,
#             when the fault type is unspecified (U = 1)
#   CB08 Model:
#    Y50.cb = Median ground motion estimate using CB08 (epsilon = 0)
#    YplusEpsilon.GM.cb = Upper bound estimate of ground motion using CB08, for the 
#                         geometric mean horizontal component of ground motion (arb = 0)
#    YplusEpsilon.arb.cb = Upper bound estimate of ground motion using CB08, for the
#                          arbitrary horizontal component of ground motion (arb = 1)
#    YminusEpsilon.GM.cb = Lower bound estimate of ground motion using CB08, for the 
#                          geometric mean horizontal component of ground motion (arb = 0)
#    YminusEpsilon.arb.cb = Lower bound estimate of ground motion using CB08, for the 
#                           arbitrary horizontal component of ground motion (arb = 1)  
#    sdGM.cb =  total standard deviation using CB08, for the
#               geometric mean horizontal component of ground motion (arb = 0)
#    sdArb.cb = total standard deviation using CB08, for the
#               arbitrary horizontal component of ground motion (arb = 1)
#   CY08 Model:
#    Y50.cy = Median ground motion estimate using CY08 (epsilon = 0)
#    YplusEpsilon.meas.cy = Upper bound estimate of ground motion using CY08,
#                           for measured Vs30 (VsFlag = 1)
#    YplusEpsilon.est.cy = Upper bound estimate of ground motion using CY08,
#                          for estimated Vs30 (VsFlag = 0)
#    YminusEpsilon.meas.cy = Lower bound estimate of ground motion using CY08,
#                            for measured Vs30 (VsFlag = 1)
#    YminusEpsilon.est.cy = Lower bound estimate of ground motion using CY08,
#                           for estimated Vs30 (VsFlag = 0)
#    sdMeas.cy = total standard deviation using CY08, for measured Vs30 (VsFlag = 1)
#    sdEst.cy = total standard deviation using CY08, for estimated Vs30 (VsFlag = 0)


Sa.nga <- function(M, Rjb, Rrup = NA, Rx = NA, rake = NA, Frv = NA, Fnm = NA, Fhw = NA,
                   dip = NA, W = NA, Ztor = NA, Vs30, Z1.0 = NA, Z1.5 = NA, Z2.5 = NA,
                   VsFlag = 0, Fas = 0, Zhyp = NA, azimuth = NA, U = 0, arb = 0,
                   epsilon, T){

  # If T is a vector, perform calculation for each of the elements
  if(length(T) > 1){

    # Obtain original data frame
    OrigDF <- sapply(T, Sa.nga, M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake,
                  Frv = Frv, Fnm = Fnm, Fhw = Fhw, dip = dip, W = W, Ztor = Ztor,
                  Vs30 = Vs30, Z1.0 = Z1.0, Z1.5 = Z1.5, Z2.5 = Z2.5, VsFlag = VsFlag,
                  Fas = Fas, Zhyp = Zhyp, azimuth = azimuth, U = U, arb = arb,
                  epsilon = epsilon)
    OrigDims <- dim(OrigDF)

    # Fill new data frame (basically the transpose of the original data frame)
    NewDF <- matrix(nrow = OrigDims[2], ncol = OrigDims[1])
    ColNames <- rownames(OrigDF)
    for(i in 1:OrigDims[1])
      NewDF[ , i] <- as.numeric(OrigDF[i, ])
    NewDF <- as.data.frame(NewDF)
    names(NewDF) <- ColNames
    return(NewDF)


  # Perform calculation for single value of T:
  } else {

  
    # 1. CHECK INPUT PARAMETERS
  
    Rrup.in <- Rrup
    Rx.in <- Rx
    rake.in <- rake
    dip.in <- dip
    W.in <- W
    Ztor.in <- Ztor
    Z1.0in <- Z1.0
    Z2.5in <- Z2.5
    azimuth.in <- azimuth
    epsilon.in <- epsilon

    # Check mandatory input parameters
    if(is.na(M) == TRUE | M < 0)
      stop("M must be a positive number")
    if(is.na(Rjb) == TRUE | Rjb < 0)
      stop("Rjb must be a non-negative number")
    if(is.na(Vs30) == TRUE | Vs30 < 0)
      stop("Vs30 must be a positive number")
    if(is.na(VsFlag == TRUE) | !(VsFlag == 1 | VsFlag == 0))
      stop("VsFlag must be either 1 (for measured Vs) or 0 (for estimated / inferred Vs)")
    if(is.na(epsilon) == TRUE)
      stop("epsilon must be numeric")
    if(is.na(Fas) == TRUE | !(Fas == 1 | Fas == 0))
      stop("Fas must be either 1 (for aftershocks) or 0 (for mainshocks, the default)")
    if(Fas == 1)
      warning("The BA08 and CB08 models are not defined for aftershocks.")
  
    # Check style of faulting parameters
    if(is.na(rake.in) == TRUE & (is.na(Frv) == TRUE | is.na(Fnm) == TRUE))
      stop("either (1) the rake angle, or (2) both Frv and Fnm must be specified")
    if(is.na(rake.in) == FALSE & (is.na(Frv) == FALSE | is.na(Fnm) == FALSE) &
       (is.na(Frv) == TRUE | is.na(Fnm) == TRUE))
      stop("either (1) the rake angle, or (2) both Frv and Fnm must be specified")  
    if(is.na(rake.in) == FALSE & abs(rake.in) > 180)
      stop("rake angle must be between -180 and 180, inclusive")
    if(is.na(Frv) == FALSE & is.na(Fnm) == FALSE){
      if(Frv == 1 & Fnm == 1)
        stop("either Frv or Fnm may be equal to 1, but both flags cannot be equal to 1")
      if(!(Frv == 1 & Fnm == 0) & !(Frv == 0 & Fnm == 1) & !(Frv == 0 & Fnm == 0))
        stop("Frv must be 1 for reverse faulting and 0 otherwise.",
             "\n", "Fnm must be 1 for normal faulting and 0 otherwise.")
    }

    # Check hanging wall parameters
    if(is.na(azimuth.in) == TRUE & is.na(Fhw) == TRUE & is.na(Rx.in) == TRUE)
      stop("At least one of Rx, azimuth, and Fhw must be specified")
    if(is.na(azimuth.in) == FALSE & abs(azimuth.in) > 180)
      stop("the source-to-site azimuth must be between -180 and 180, inclusive")
    if(is.na(Fhw) == FALSE & !(Fhw == 1 | Fhw == 0))
      stop("Fhw must be either 1 (for sites on the hanging wall side of the fault)", "\n",
           "or 0 (for sites on the footwall side of the fault)")
    # Ensure consistency between Rx, azimuth, and Fhw
    if(is.na(azimuth.in) == FALSE & is.na(Fhw) == FALSE){
      if(azimuth.in < 0 & Fhw == 1)
        stop("Inconsistency between azimuth and Fhw. Fhw must be 0 when azimuth < 0.")
      if(azimuth.in > 0 & Fhw == 0)
        stop("Inconsistency between azimuth and Fhw. Fhw must be 1 when azimuth > 0.")
    }
    if(is.na(Rx.in) == FALSE & is.na(Fhw) == FALSE){
      if(Rx.in < 0 & Fhw == 1)
        stop("Inconsistency between Rx and Fhw. Fhw must be 0 when Rx < 0.")
      if(Rx.in > 0 & Fhw == 0)
        stop("Inconsistency between Rx and Fhw. Fhw must be 1 when Rx > 0.")
    }
    if(is.na(Rx.in) == FALSE & is.na(azimuth.in) == FALSE){
      if(!(Rx.in <= 0 & azimuth.in <= 0) & !(Rx.in >= 0 & azimuth.in >= 0))
        stop("Rx and azimuth must have the same sign.")
    }



    # 2. OBTAIN ESTIMATES OF UNSPECIFIED INPUT PARAMETERS
  
    # Assign generic rake angle if Fhw and Fnm are specified
    if(is.na(rake.in) == TRUE){
      if(Frv == 1 & Fnm == 0){  # Reverse faulting
        rake <- 90
      } else{
        if(Frv == 0 & Fnm == 1){  # Normal faulting
          rake <- -90
        } else{
          if(Frv == 0 & Fnm == 0){  # Strike-slip faulting
            rake <- 180
          }
        }
      }
      # Set Fnm1 and Fnm2 to be equal to Fnm
      Fnm1 <- Fnm
      Fnm2 <- Fnm
    }
    
    # Convert rake to fault type if rake is provided
    if(is.na(Frv) == TRUE & is.na(Fnm) == TRUE){
      # Frv
      if(rake.in >= 30 & rake.in <= 150)
        Frv <- 1
      else
        Frv <- 0
      # Fnm1 (AS08 and CY08)
      if(rake.in >= -120 & rake.in <= -60)
        Fnm1 <- 1
      else
        Fnm1 <- 0
      # Fnm2 (BA08 and CB08)
      if(rake.in >= -150 & rake.in <= -30)
        Fnm2 <- 1
      else
        Fnm2 <- 0 
    }
    
    # Dip angle
    if(is.na(dip.in) == TRUE | dip.in < 0)
      dip <- dip.calc(rake)

    # Down-dip rupture width, W
    if(is.na(W.in) == TRUE | W.in < 0)
      W <- W.calc(M, rake)
    
    # Depth to top of rupture, Ztor
    if(is.na(Ztor.in) == TRUE | Ztor.in < 0)
      Ztor <- Ztor.calc(Zhyp, W, dip, M, rake)

    # Determine hanging wall flag, Fhw
    # Calculate from azimuth if provided
    if(is.na(azimuth.in) == FALSE){
      if(azimuth.in > 0 & azimuth.in < 180 & dip != 90)
        Fhw <- 1
      else
        Fhw <- 0
     # Calculate from Rx if azimuth is not provided
    } else{
      if(is.na(Rx.in) == FALSE){
        if(Rx.in >= 0)
          Fhw <- 1
        else
          Fhw <- 0
      }
    }

    # Azimuth angle
    if(is.na(azimuth.in) == TRUE){
      # If hanging wall site, assume azimuth = 50 deg
      if(Fhw == 1)
        azimuth <- 50
      # If footwall site, assume azimuth = -50 deg  
      else if(Fhw == 0)
        azimuth <- -50
    }
    
    # Site coordinate, Rx
    if(is.na(Rx.in) == TRUE)
      Rx <- Rx.calc(Rjb = Rjb, Ztor = Ztor, W = W, dip = dip,
                    azimuth = azimuth, Rrup = Rrup.in)
  
    # Rupture distance, Rrup
    if(is.na(Rrup.in) == TRUE | Rrup.in < 0)
      Rrup <- Rrup.calc(Rx, Ztor, W, dip, azimuth, Rjb)

    # Depth parameter, Z1.0
    if(is.na(Z1.0in) == TRUE | Z1.0in < 0){
      Z1.0as <- Z1.calc.as(Vs30)
      Z1.0cy <- Z1.calc.cy(Vs30)
    } else{
      Z1.0as <- Z1.0
      Z1.0cy <- Z1.0
    }
     
    # Depth parameter, Z2.5
    if(is.na(Z2.5in) == TRUE | Z2.5in < 0){
      # Calculate from Z1.5 if provided
      # (Eqn 6.4 in Campbell and Bozorgnia (2007); final report to PEER)
      if(is.na(Z1.5) == FALSE & Z1.5 > 0){
        Z2.5 <- 0.636 + 0.001549*Z1.5
        # Calculate from Z1.0 if provided
        # (Eqn 6.3 in Campbell and Bozorgnia (2007); final report to PEER)
      } else{
        if(is.na(Z1.0) == FALSE & Z1.0 > 0){
          Z2.5 <- 0.519 + 0.003595*Z1.0
          # If neither Z1.0 nor Z1.5 is provided, estimate from
          # Vs30 using the AS08 relation for Z1.0 = f(Vs30)
        } else{
          Z1.0cb <- Z1.calc.as(Vs30)
          Z2.5 <- 0.519 + 0.003595*Z1.0cb
        }
      }
    }


  
    # 3. CALCULATION OF GROUND MOTIONS

    # Median Ground Motions (50th percentile)
  
    Y50.as <-  Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake, Frv = Frv,
                     Fnm = Fnm1, Fhw = Fhw, dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30,
                     Z1.0 = Z1.0as, VsFlag = VsFlag, Fas = Fas, Zhyp = Zhyp,
                     azimuth = azimuth, epsilon = 0, T = T)
  
    # Specified fault mechanism
    Y50M.ba <- Sa.ba(M = M, Rjb = Rjb, Vs30 = Vs30, rake = rake, U = U,
                     SS = if(Frv == 0 & Fnm2 == 0) 1 else 0, NS = Fnm2,
                     RS = Frv, epsilon = 0, T = T)

    # Unspecified fault mechanism
    Y50U.ba <- Sa.ba(M = M, Rjb = Rjb, Vs30 = Vs30, rake = NA, U = 1,
                     SS = 0, NS = 0, RS = 0, epsilon = 0, T = T)
    
    Y50.cb <- Sa.cb(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, Frv = Frv, Fnm = Fnm2,
                    dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30, Z1.0 = Z1.0cb,
                    Z1.5 = Z1.5, Z2.5 = Z2.5, Zhyp = Zhyp, Fhw = Fhw, azimuth = azimuth,
                    arb = arb, epsilon = 0, T = T)
    
    Y50.cy <- Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake, Frv = Frv,
                    Fnm = Fnm1, Fhw = Fhw, dip = dip, W = W, Ztor = Ztor,
                    Vs30 = Vs30, Z1.0 = Z1.0cy, VsFlag = VsFlag, AS = Fas, Zhyp = Zhyp,
                    azimuth = azimuth, epsilon = 0, T = T)
    

    # Median + epsilon*SD

    # Special case of epsilon = 0
    # This is necessary to avoid division by 0;
    # correction occurs later in code
    if(epsilon.in == 0)
      epsilon <- 1

    # AS08
    # Measured Vs30
    YplusEpsilon.meas.as <-  Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake,
                                   Frv = Frv, Fnm = Fnm1, Fhw = Fhw, dip = dip, W = W,
                                   Ztor = Ztor, Vs30 = Vs30, Z1.0 = Z1.0as, VsFlag = 1,
                                   Fas = Fas, Zhyp = Zhyp, azimuth = azimuth,
                                   epsilon = epsilon, T = T)
    # Estimated Vs30
    YplusEpsilon.est.as <-  Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake,
                                  Frv = Frv, Fnm = Fnm1, Fhw = Fhw, dip = dip, W = W,
                                  Ztor = Ztor, Vs30 = Vs30, Z1.0 = Z1.0as, VsFlag = 0,
                                  Fas = Fas, Zhyp = Zhyp, azimuth = azimuth,
                                  epsilon = epsilon, T = T) 

    # BA08
    # Specified fault type
    YplusEpsilon.M.ba <- Sa.ba(M = M, Rjb = Rjb, Vs30 = Vs30, rake = rake, U = 0,
                               SS = if(Frv == 0 & Fnm2 == 0) 1 else 0, NS = Fnm2,
                               RS = Frv, epsilon = epsilon, T = T)
    # Unspecified fault type
    YplusEpsilon.U.ba <- Sa.ba(M = M, Rjb = Rjb, Vs30 = Vs30, rake = NA, U = 1,
                               SS = 0, NS = 0, RS = 0, epsilon = epsilon, T = T)

    # CB08
    # Geometric mean horizontal component
    YplusEpsilon.GM.cb <- Sa.cb(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, Frv = Frv,
                                Fnm = Fnm2, dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30,
                                Z1.0 = Z1.0cb, Z1.5 = Z1.5, Z2.5 = Z2.5, Zhyp = Zhyp,
                                Fhw = Fhw, azimuth = azimuth, arb = 0,
                                epsilon = epsilon, T = T) 
    # Arbitrary horizontal component
    YplusEpsilon.arb.cb <- Sa.cb(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, Frv = Frv,
                                 Fnm = Fnm2, dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30,
                                 Z1.0 = Z1.0cb, Z1.5 = Z1.5, Z2.5 = Z2.5, Zhyp = Zhyp,
                                 Fhw = Fhw, azimuth = azimuth, arb = 1,
                                 epsilon = epsilon, T = T) 

    # CY08
    # Measured Vs30 
    YplusEpsilon.meas.cy <- Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake,
                                  Frv = Frv, Fnm = Fnm1, Fhw = Fhw, dip = dip, W = W,
                                  Ztor = Ztor, Vs30 = Vs30, Z1.0 = Z1.0cy, VsFlag = 1,
                                  AS = Fas, Zhyp = Zhyp, azimuth = azimuth,
                                  epsilon = epsilon, T = T) 
    # Estimated Vs30  
    YplusEpsilon.est.cy <- Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake,
                                 Frv = Frv, Fnm = Fnm1, Fhw = Fhw, dip = dip, W = W,
                                 Ztor = Ztor, Vs30 = Vs30, Z1.0 = Z1.0cy, VsFlag = 0,
                                 AS = Fas, Zhyp = Zhyp, azimuth = azimuth,
                                 epsilon = epsilon, T = T)


    # Median - epsilon*SD

    # AS08
    # Measured Vs30
    YminusEpsilon.meas.as <-  Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake,
                                    Frv = Frv, Fnm = Fnm1, Fhw = Fhw, dip = dip, W = W,
                                    Ztor = Ztor, Vs30 = Vs30, Z1.0 = Z1.0as, VsFlag = 1,
                                    Fas = Fas, Zhyp = Zhyp, azimuth = azimuth,
                                    epsilon = -epsilon, T = T)
    # Estimated Vs30
    YminusEpsilon.est.as <-  Sa.as(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake,
                                   Frv = Frv, Fnm = Fnm1, Fhw = Fhw, dip = dip, W = W,
                                   Ztor = Ztor, Vs30 = Vs30, Z1.0 = Z1.0as, VsFlag = 0,
                                   Fas = Fas, Zhyp = Zhyp, azimuth = azimuth,
                                   epsilon = -epsilon, T = T) 

    # BA08
    # Specified fault type
    YminusEpsilon.M.ba <- Sa.ba(M = M, Rjb = Rjb, Vs30 = Vs30, rake = rake, U = 0,
                                SS = if(Frv == 0 & Fnm2 == 0) 1 else 0, NS = Fnm2,
                                RS = Frv, epsilon = -epsilon, T = T)
    # Unspecified fault type
    YminusEpsilon.U.ba <- Sa.ba(M = M, Rjb = Rjb, Vs30 = Vs30, rake = NA, U = 1,
                                SS = 0, NS = 0, RS = 0, epsilon = -epsilon, T = T)

    # CB08
    # Geometric mean horizontal component
    YminusEpsilon.GM.cb <- Sa.cb(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, Frv = Frv,
                                 Fnm = Fnm2, dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30,
                                 Z1.0 = Z1.0cb, Z1.5 = Z1.5, Z2.5 = Z2.5, Zhyp = Zhyp,
                                 Fhw = Fhw, azimuth = azimuth, arb = 0,
                                 epsilon = -epsilon, T = T) 
    # Arbitrary horizontal component
    YminusEpsilon.arb.cb <- Sa.cb(M = M, Rjb = Rjb, Rrup = Rrup, rake = rake, Frv = Frv,
                                  Fnm = Fnm2, dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30,
                                  Z1.0 = Z1.0cb, Z1.5 = Z1.5, Z2.5 = Z2.5, Zhyp = Zhyp,
                                  Fhw = Fhw, azimuth = azimuth, arb = 1,
                                  epsilon = -epsilon, T = T) 

    # CY08
    # Measured Vs30 
    YminusEpsilon.meas.cy <- Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake,
                                   Frv = Frv, Fnm = Fnm1, Fhw = Fhw, dip = dip, W = W,
                                   Ztor = Ztor, Vs30 = Vs30, Z1.0 = Z1.0cy, VsFlag = 1,
                                   AS = Fas, Zhyp = Zhyp, azimuth = azimuth,
                                   epsilon = -epsilon, T = T) 
    # Estimated Vs30  
    YminusEpsilon.est.cy <- Sa.cy(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake,
                                  Frv = Frv, Fnm = Fnm1, Fhw = Fhw, dip = dip, W = W,
                                  Ztor = Ztor, Vs30 = Vs30, Z1.0 = Z1.0cy, VsFlag = 0,
                                  AS = Fas, Zhyp = Zhyp, azimuth = azimuth,
                                  epsilon = -epsilon, T = T)


    # Standard Deviations (log space)
    sdMeas.as <- (log(YplusEpsilon.meas.as) - log(Y50.as)) / epsilon
    sdEst.as <- (log(YplusEpsilon.est.as) - log(Y50.as)) / epsilon
    sdM.ba <- (log(YplusEpsilon.M.ba) - log(Y50M.ba)) / epsilon
    sdU.ba <- (log(YplusEpsilon.U.ba) - log(Y50U.ba)) / epsilon
    sdGM.cb <- (log(YplusEpsilon.GM.cb) - log(Y50.cb)) / epsilon
    sdArb.cb <- (log(YplusEpsilon.arb.cb) - log(Y50.cb)) / epsilon
    sdMeas.cy <- (log(YplusEpsilon.meas.cy) - log(Y50.cy)) / epsilon
    sdEst.cy <- (log(YplusEpsilon.est.cy) - log(Y50.cy)) / epsilon

    # Correction for zero epsilon
    if(epsilon.in == 0){
      YplusEpsilon.meas.as = Y50.as
      YplusEpsilon.est.as = Y50.as
      YminusEpsilon.meas.as = Y50.as
      YminusEpsilon.est.as = Y50.as
      YplusEpsilon.M.ba = Y50M.ba
      YplusEpsilon.U.ba = Y50U.ba
      YminusEpsilon.M.ba = Y50M.ba
      YminusEpsilon.U.ba = Y50U.ba
      YplusEpsilon.GM.cb = Y50.cb
      YplusEpsilon.arb.cb = Y50.cb
      YminusEpsilon.GM.cb = Y50.cb
      YminusEpsilon.arb.cb = Y50.cb
      YplusEpsilon.meas.cy = Y50.cy
      YplusEpsilon.est.cy = Y50.cy 
      YminusEpsilon.meas.cy = Y50.cy
      YminusEpsilon.est.cy = Y50.cy
    }



    # 4. RETURN LIST OF INPUTS AND OUTPUTS
    return(list(T = T,
                epsilon = epsilon.in,
                M = M,
                Rjb = Rjb,
                Rrup.in = Rrup.in,
                Rrup.out = Rrup,
                Rx.in = Rx.in,
                Rx.out = Rx,
                rake.in = rake.in,
                rake.out = rake,
                Frv = Frv,
                Fnm1 = Fnm1,
                Fnm2 = Fnm2,
                Fhw = Fhw,
                dip.in = dip.in,
                dip.out = dip,
                W.in = W.in,
                W.out = W,
                Ztor.in = Ztor.in,
                Ztor.out = Ztor,
                Vs30 = Vs30,
                Z1.0in = Z1.0in,
                Z1.0as = Z1.0as,
                Z1.0cy = Z1.0cy,
                Z1.5in = Z1.5,
                Z2.5in = Z2.5in,
                Z2.5out = Z2.5,
                Fas = Fas,
                Zhyp = Zhyp,
                azimuth.in = azimuth.in,
                azimuth.out = azimuth,
                Y50.as = Y50.as,
                YplusEpsilon.meas.as = YplusEpsilon.meas.as,
                YplusEpsilon.est.as = YplusEpsilon.est.as,
                YminusEpsilon.meas.as = YminusEpsilon.meas.as,
                YminusEpsilon.est.as = YminusEpsilon.est.as,
                sdMeas.as = sdMeas.as,
                sdEst.as = sdEst.as,
                Y50M.ba = Y50M.ba,
                Y50U.ba = Y50U.ba,
                YplusEpsilon.M.ba = YplusEpsilon.M.ba,
                YplusEpsilon.U.ba = YplusEpsilon.U.ba,
                YminusEpsilon.M.ba = YminusEpsilon.M.ba,
                YminusEpsilon.U.ba = YminusEpsilon.U.ba,
                sdM.ba = sdM.ba,
                sdU.ba = sdU.ba,
                Y50.cb = Y50.cb,
                YplusEpsilon.GM.cb = YplusEpsilon.GM.cb,
                YplusEpsilon.arb.cb = YplusEpsilon.arb.cb,            
                YminusEpsilon.GM.cb = YminusEpsilon.GM.cb,
                YminusEpsilon.arb.cb = YminusEpsilon.arb.cb,
                sdGM.cb =  sdGM.cb,
                sdArb.cb = sdArb.cb,
                Y50.cy = Y50.cy,
                YplusEpsilon.meas.cy = YplusEpsilon.meas.cy,
                YplusEpsilon.est.cy = YplusEpsilon.est.cy,         
                YminusEpsilon.meas.cy = YminusEpsilon.meas.cy,
                YminusEpsilon.est.cy = YminusEpsilon.est.cy,
                sdMeas.cy = sdMeas.cy,
                sdEst.cy = sdEst.cy))
  }
}




# FUNCTION II:  Sa.ngaR (reduced output)

# OUTLINE OF CODE
# 1. Obtain extended data frame of ground motion predictions
# 2. Obtain elements of reduced data frame
# 3. Return reduced list of values

# This function advantage of the input parameters VsFlag (AS08 and CY08),
# U (BA08), and arb (CB08), and returns a list of reduced output.  The
# "Input Variables" section of the list is the same, except VsFlag, U,
# and arb are also included in the list.  The simplified "Output Variables"
# section of "Sa.ngaR" is:
# AS08 Model:
#   Y50.as = Median ground motion estimate using AS08 (epsilon = 0)
#   YplusEpsilon.as = Upper bound estimate of ground motion using AS08
#   YminusEpsilon.as = Lower bound estimate of ground motion using AS08
# BA08 Model:
#   Y50.ba =  Median ground motion estimate using BA08 (epsilon = 0)
#   YplusEpsilon.ba = Upper bound estimate of ground motion using BA08
#   YminusEpsilon.ba = Lower bound estimate of ground motion using BA08
# CB08 Model:
#   Y50.cb = Median ground motion estimate using CB08 (epsilon = 0)
#   YplusEpsilon.cb = Upper bound estimate of ground motion using CB08
#   YminusEpsilon.cb = Lower bound estimate of ground motion using CB08
# CY08 Model:
#   Y50.cy = Median ground motion estimate using CY08 (epsilon = 0)
#   YplusEpsilon.cy = Upper bound estimate of ground motion using CY08
#   YminusEpsilon.cy = Lower bound estimate of ground motion using CY08


Sa.ngaR <- function(M, Rjb, Rrup = NA, Rx = NA, rake = NA, Frv = NA, Fnm = NA, Fhw = NA,
                    dip = NA, W = NA, Ztor = NA, Vs30, Z1.0 = NA, Z1.5 = NA, Z2.5 = NA,
                    VsFlag, Fas = 0, Zhyp = NA, azimuth = NA, U = 0, arb = 0,
                    epsilon, T){

  # 1.  Obtain extended data frame of ground motion predictions
  #     "E" stands for extended, as in "ngaE"
  #     "R" stands for reduced, as in "ngaR"
  
  ngaE <- Sa.nga(M = M, Rjb = Rjb, Rrup = Rrup, Rx = Rx, rake = rake, Frv = Frv, Fnm = Fnm,
                 Fhw = Fhw, dip = dip, W = W, Ztor = Ztor, Vs30 = Vs30, Z1.0 = Z1.0,
                 Z1.5 = Z1.5, Z2.5 = Z2.5, VsFlag = VsFlag, Fas = Fas, Zhyp = Zhyp,
                 azimuth = azimuth, U = U, arb = arb, epsilon = epsilon, T = T)
  
  # Check other flag variables
  if(is.na(VsFlag) == TRUE | !(VsFlag == 1 | VsFlag == 0))
    stop("VsFlag must be either 1 (for measured Vs) or 0 (for estimated Vs)")
  if(is.na(U) == TRUE | !(U == 1 | U == 0))
    stop("U must be either 0 or 1.")
  if(is.na(arb) == TRUE | !(arb == 1 | arb == 0))
    stop("arb must be either 0 or 1.")


  # 2.  Obtain elements of reduced data frame

  # Select values corresponding to VsFlag for AS08 and CY08
  # Estimated Vs30 (VsFlag = 0)
  if(VsFlag == 0){
    YplusEpsilon.as <- ngaE$YplusEpsilon.est.as
    YminusEpsilon.as <- ngaE$YminusEpsilon.est.as
    YplusEpsilon.cy <- ngaE$YplusEpsilon.est.cy
    YminusEpsilon.cy <- ngaE$YminusEpsilon.est.cy
  # Measured Vs30 (VsFlag = 1)
  } else{
    if(VsFlag == 1){
      YplusEpsilon.as <- ngaE$YplusEpsilon.meas.as
      YminusEpsilon.as <- ngaE$YminusEpsilon.meas.as
      YplusEpsilon.cy <- ngaE$YplusEpsilon.meas.cy
      YminusEpsilon.cy <- ngaE$YminusEpsilon.meas.cy
    }
  }
  
  # Select values corresponding to U for BA08
  # Specified fault mechanism (U = 0)
  if(U == 0){
    Y50.ba <- ngaE$Y50M.ba
    YplusEpsilon.ba <- ngaE$YplusEpsilon.M.ba
    YminusEpsilon.ba <- ngaE$YminusEpsilon.U.ba
  # Unspecified fault mechanism (U = 1)
  } else{
    if(U == 1){
      Y50.ba <- ngaE$Y50U.ba
      YplusEpsilon.ba <- ngaE$YplusEpsilon.U.ba
      YminusEpsilon.ba <- ngaE$YminusEpsilon.U.ba
    }
  }

  # Select values corresponding to arb for CB08
  # Geometric mean horizontal component (arb = 0)
  if(arb == 0){
    YplusEpsilon.cb <- ngaE$YplusEpsilon.GM.cb
    YminusEpsilon.cb <- ngaE$YminusEpsilon.GM.cb
  # Arbitrary horizontal component (arb = 1)
  } else{
    if(arb == 1){
      YplusEpsilon.cb <- ngaE$YplusEpsilon.arb.cb
      YminusEpsilon.cb <- ngaE$YminusEpsilon.arb.cb
    }
  }
  

  # 3. Develop data frame to return

  # Extract relevant columns from ngaE
  T <- ngaE$T
  epsilon <- ngaE$epsilon
  M <- ngaE$M
  Rjb <- ngaE$Rjb
  Rrup.in <- ngaE$Rrup.in
  Rrup.out <- ngaE$Rrup.out
  Rx.in <- ngaE$Rx.in
  Rx.out <- ngaE$Rx.out
  rake.in <- ngaE$rake.in
  rake.out <- ngaE$rake.out
  Frv <- ngaE$Frv
  Fnm1 <- ngaE$Fnm1
  Fnm2 <- ngaE$Fnm2
  Fhw <- ngaE$Fhw
  dip.in <- ngaE$dip.in
  dip.out <- ngaE$dip.out
  W.in <- ngaE$W.in
  W.out <- ngaE$W.out
  Ztor.in <- ngaE$Ztor.in
  Ztor.out <- ngaE$Ztor.out
  Vs30 <- ngaE$Vs30
  Z1.0in <- ngaE$Z1.0in
  Z1.0as <- ngaE$Z1.0as
  Z1.0cy <- ngaE$Z1.0cy
  Z1.5in <- ngaE$Z1.5
  Z2.5in <- ngaE$Z2.5in
  Z2.5out <- ngaE$Z2.5out
  Fas <- ngaE$Fas
  Zhyp <- ngaE$Zhyp
  azimuth.in <- ngaE$azimuth.in
  azimuth.out <- ngaE$azimuth.out
  Y50.as <- ngaE$Y50.as
  Y50.cb <- ngaE$Y50.cb
  Y50.cy <- ngaE$Y50.cy

  # Generate column vectors for VsFlag, U, and arb
  # Store values of VsFlag, U, and arb in temporary variables
  VsFlag.value <- VsFlag
  U.value <- U
  arb.value <- arb
  # Name vectors
  VsFlag <- matrix(nrow = length(T), ncol = 1)
  U <- matrix(nrow = length(T), ncol = 1)
  arb <- matrix(nrow = length(T), ncol = 1)
  for(i in 1:length(T)){
    VsFlag[i] <- VsFlag.value
    U[i] <- U.value
    arb[i] <- arb.value
  }
  
  # Define data frame with relevant vectors
  ngaR <- data.frame(T, epsilon, M, Rjb, Rrup.in, Rrup.out, Rx.in, Rx.out,
                     rake.in, rake.out, Frv, Fnm1, Fnm2, Fhw, dip.in, dip.out,
                     W.in, W.out, Ztor.in, Ztor.out, Vs30, Z1.0in, Z1.0as,
                     Z1.0cy, Z1.5in, Z2.5in, Z2.5out, Fas, Zhyp, azimuth.in,
                     azimuth.out, VsFlag, U, arb, Y50.as, YplusEpsilon.as,
                     YminusEpsilon.as, Y50.ba, YplusEpsilon.ba, YminusEpsilon.ba,
                     Y50.cb, YplusEpsilon.cb, YminusEpsilon.cb, Y50.cy,
                     YplusEpsilon.cy, YminusEpsilon.cy)
  return(ngaR)
}

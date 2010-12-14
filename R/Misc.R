# IMPLEMENTATION OF THE NEXT GENERATION ATTENUATION (NGA)
# GROUND-MOTION PREDICTION EQUATIONS IN R


# Misc.R:  MISCELLANEOUS FUNCTIONS NEEDED FOR MULTIPLE NGA MODELS


# James Kaklamanos
# Tufts University
# Department of Civil and Environmental Engineering
# James.Kaklamanos@tufts.edu
# http://geohazards.cee.tufts.edu/people/jkakla01
# December 14, 2010

# For further details, see the following papers:
# o  Kaklamanos, J., D. M. Boore, E. M. Thompson, and K. W. Campbell (2010).
#    Implementation of the Next Generation Attenuation (NGA) ground-motion
#    prediction equations in Fortran and R, U.S. Geological Survey Open-File
#    Report 2010-1296.
# o  Kaklamanos, J., L. G. Baise, and D. M. Boore (2011). Estimating unknown input
#    parameters when implementing the NGA ground-motion prediction equations in
#    engineering practice, Earthquake Spectra, in press.




# OUTLINE OF CODE:
# This code contains functions that are needed for multiple NGA models, including:
# 1. Interpolation of spectral acceleration
# 2. Calculation of distance measures (Rx and Rrup)
# 3. Calculation of depth parameter, Z1.0
# 4. Calculation of hypocentral depth, Zhyp
# 5. Calculation of depth to top of rupture, Ztor
# 6. Calculation of dip angle
# 7. Calculation of down-dip width, W
# 8. Calculation of depth parameter, Z2.5




# 1. INTERPOLATION OF SPECTRAL ACCELERATION

# This function returns the list of model periods for the specified model,
# and whether or not PGA, PGV, and (in the case of CB08, PGD) should be included

modelPeriods <- function(model, positive = FALSE){
  if(toupper(model) == "AS" | toupper(model) == "AS08")
    periods.as(positive)
  else{
    if(toupper(model) == "BA" | toupper(model) == "BA08")
      periods.ba(positive)
    else{
      if(toupper(model) == "CB" | toupper(model) == "CB08")
        periods.cb(positive)
      else{
        if(toupper(model) == "CY" | toupper(model) == "CY08")
          periods.cy(positive)
        else{
          stop("please enter a valid NGA model")
        }
      }
    }
  }
}


# This function receives a period T and returns a three-element list
# "interp" states whether or not interpolation is necessary
# "lower" gives the period less than T that will be used for interpolation
#         (if interpolation is necessary)
# "upper" gives the period greater than T that will be used for interpolation
#         (if interpolation is necessary)

getPeriod <- function(T, model){

  # Check to be sure the model is valid
  if(toupper(model) != "AS" & toupper(model) != "AS08" &
     toupper(model) != "BA" & toupper(model) != "BA08" &
     toupper(model) != "CB" & toupper(model) != "CB08" &
     toupper(model) != "CY" & toupper(model) != "CY08")
    stop("please enter a valid NGA model name")
    
  # See if the period matches one of those with defined coefficients
  if(T %in% modelPeriods(model, 0)){
    interp <- FALSE  # interpolation not necessary
  } else {
    if(T < 0.01 | T > 10){  # period not within allowable range
      stop("spectral period not within allowable range.")
    } else{
      interp <- TRUE  # period OK; interpolation is necessary
    }
  }
  
  # If interpolation is necessary, determine the periods with defined
  # coefficients directly above and below the period of interest
 
  if(interp == TRUE){
    index.lower <- findInterval(T, modelPeriods(model, 1))
    lower <- modelPeriods(model, 1)[index.lower]
    upper <- modelPeriods(model, 1)[index.lower + 1]
  } else {
    lower <- NA
    upper <- NA
  }
  
  return(list(interp = interp, lower = lower, upper = upper))
}


# Calculate Sa using interpolation
interpolate <- function(x, x1, x2, y1, y2){
  approx(c(x1, x2), c(y1, y2), x)$y
}




# 2. CALCULATION OF DISTANCE MEASURES (Rx and Rrup)
# For details on these derived equations (and on the estimation of other input
# parameters), please see the following paper:
# Kaklamanos, J., L. G. Baise, and D. M. Boore (2011).  Estimating
# Unknown Input Parameters when Implementing the NGA Ground Motion Prediction
# Equations in Engineering Practice. Earthquake Spectra (in press).
# The equation numbers and figures in this code refer to the equation numbers
# and figures in Kaklamanos et al (2011).

# Necessary trigonometric functions for Rx calculations
csc <- function(x) {1/sin(x)}
sec <- function(x) {1/cos(x)}
cot <- function(x) {1/tan(x)}


# Function for Site Coordinate, Rx
Rx.calc <- function(Rjb, Ztor, W, dip, azimuth, Rrup = NA){

  # Check input
  if(is.na(Rjb) == TRUE | Rjb < 0)
    stop("Rjb must be a non-negative number")
  if(is.na(Ztor) == TRUE | Ztor < 0)
    stop("Ztor must be a non-negative number")
  if(is.na(W) == TRUE | W <= 0)
    stop("W must be a positive number")
  if(is.na(dip) == TRUE | dip <= 0 | dip > 90)
    stop("dip angle must be positive and no greater than 90")
  if(is.na(azimuth) == TRUE | abs(azimuth) > 180)
    stop("azimuth must be between -180 and 180, inclusive")
  
  
  # Define angles in terms of radians
  d <- dip*pi/180
  a <- azimuth*pi/180


  # Non-vertical faults
  if(dip != 90){

    if(azimuth == 90 | Rjb == 0){

      # Rjb > 0
      # Case 6 in Figure 3 (Eqn 9)
      if(Rjb > 0)
        Rx <- Rjb + W*cos(d)

      # Rjb = 0
      else if(Rjb == 0){
        
        # If Rrup is known, then calculate Rx from Rrup
        if(is.na(Rrup) == FALSE){
          # Case 5A in Figure 3 (Eqn 10)
          if(Rrup < Ztor*sec(d))
            Rx <- sqrt(Rrup^2 - Ztor^2)
          # Case 5B in Figure 3 (Eqn 11)
          else{
            if(Rrup >= Ztor*sec(d))
              Rx <- Rrup*csc(d) - Ztor*cot(d)
          }
        }

        # If Rrup is unknown, assume that the site is located at the
        # center of the surface projection of the ruptured area (Eqn 5)
        else
          Rx <- W*cos(d)/2
      }

    } else{
      if(azimuth >= 0 & azimuth != 90) {
        
        # Cases 2 and 8 (Eqn 7)
        if(Rjb*abs(tan(a)) <= W*cos(d)){
          Rx <- Rjb*abs(tan(a))

        # Cases 3 and 9 (Eqn 8)
        } else{
          if(Rjb*abs(tan(a)) > W*cos(d)){
            Rx <- Rjb*tan(a)*cos(a - asin(W*cos(d)*cos(a)/Rjb)) }
        }
      
      # Cases 1, 4, and 7 (Eqn 12)
      } else{
        if(azimuth >= -180 & azimuth < 0)
          Rx <- Rjb*sin(a)
      }
    }

  # Special case for vertical faults (Eqn 13)
  } else{
    if(dip == 90)
      Rx <- Rjb*sin(a)
  }
    
  return(Rx)
}



# Function for Rupture Distance, Rrup
Rrup.calc <- function(Rx, Ztor, W, dip, azimuth, Rjb = NA) {

  # Check input
  if(is.na(Rx) == TRUE)
    stop("Rx must be numeric")
  if(is.na(Ztor) == TRUE | Ztor < 0)
    stop("Ztor must be a non-negative number")
  if(is.na(W) == TRUE | W <= 0)
    stop("W must be a positive number")
  if(is.na(dip) == TRUE | dip <= 0 | dip > 90)
    stop("dip angle must be positive and no greater than 90")
  if(is.na(azimuth) == TRUE | abs(azimuth) > 180)
    stop("azimuth must be between -180 and 180, inclusive")
  
  # Define angles in terms of radians
  d <- dip*pi/180
  a <- azimuth*pi/180

  # Shortcut calculation for vertical faults, if Rjb is entered (Eqn 21)
  if(dip == 90 & is.na(Rjb) == FALSE){
    return(sqrt(Rjb^2 + Ztor^2))
  }


  # Calculate Rrup'
  if(dip != 90){

    # Zone A in Figure 4 (Eqn 15)
    if(Rx < Ztor*tan(d)){
      Rrup.prime <- sqrt(Rx^2 + Ztor^2)
      
    # Zone B in Figure 4 (Eqn 16)
    } else{
      if(Rx >= Ztor*tan(d) & Rx <= Ztor*tan(d) + W*sec(d)){
        Rrup.prime <- Rx*sin(d) + Ztor*cos(d)
        
      # Zone C in Figure 4 (Eqn 17)
      } else{
        if(Rx > Ztor*tan(d) + W*sec(d)){
          Rrup.prime <- sqrt((Rx - W*cos(d))^2 + (Ztor + W*sin(d))^2)
        }
      }
    }
  } else if(dip == 90)  # Vertical faults, if Rjb is not entered
                        # (equivalent to Eqn 21)    
    {Rrup.prime <- sqrt(Rx^2 + Ztor^2)}


  # Calculate Ry

  # Eqn 18
  if(azimuth == 90 | azimuth == -90){
    Ry <- 0
  
  # Eqn 19
  } else{
    if(azimuth == 0 | azimuth == 180 | azimuth == -180){
      if(is.na(Rjb) == TRUE)
        stop("Rjb must be specified")
      else
        Ry <- Rjb
      
    # Eqn 20
    } else{
      Ry <- abs(Rx*cot(a))
    }
  }

  # Calculate Rrup (Eqn 14)
  Rrup <- sqrt(Rrup.prime^2 + Ry^2)
  
  return(Rrup)
}




# 3. CALCULATION OF DEPTH PARAMETER, Z1.0

# AS08 [Eqn 17 in Abrahamson and Silva (2008)]
Z1.calc.as <- function(Vs30) {
  
  # Check input
  if(is.na(Vs30) == TRUE | Vs30 < 0)
    stop("Vs30 must be a positive number")
  
  # Calculate Ln(Z1.0)
  if(Vs30 < 180) {
    LnZ1.0 <- 6.745
  } else { 
    if(Vs30 >= 180 & Vs30 <= 500) {
      LnZ1.0 <- 6.745 - 1.35*log(Vs30/180)
    } else { 
      if(Vs30 > 500) {
        LnZ1.0 <- 5.394 - 4.48*log(Vs30/500)
      }
    }
  }

  # Calculate Z1.0
  return(exp(LnZ1.0))
}


# CY08 [Eqn 1 in Chiou and Youngs (2008)]
Z1.calc.cy <- function(Vs30) {

  # Check input
  if(is.na(Vs30) == TRUE | Vs30 < 0)
    stop("Vs30 must be a positive number")
  
  # Calculate Z1.0
  LnZ1.0 <- 28.5 - (3.82/8)*log(Vs30^8 + 378.7^8)
  return(exp(LnZ1.0))
}




# 4. CALCULATION OF HYPOCENTRAL DEPTH, Zhyp

Zhyp.calc <- function(M, rake) {

  # Check input
  if(is.na(M) == TRUE | M < 0)
    stop("M must be a positive number")
  if(is.na(rake) == TRUE | abs(rake) > 180)
    stop("rake angle must be between -180 and 180, inclusive")

  # Estimate Zhyp from M and rake, using correlations in Scherbaum et al (2004)
  # Strike-slip fault (SS)
  if(abs(rake) < 30 | abs(rake) > 150){
    Zhyp <- 5.63 + 0.68*M
  # Shallow-dipping fault (SHD)  
  } else{
    Zhyp <- 11.24 - 0.20*M
  }
  return(Zhyp)
}


 

# 5. CALCULATION OF DEPTH TO TOP OF RUPTURE, Ztor

Ztor.calc <- function(W, dip, Zhyp) {

  # Check input
  if(is.na(W) == TRUE | W < 0)
    stop("W must be a non-negative number")
  if(is.na(dip) == TRUE | dip <= 0 | dip > 90)
      stop("dip must be be greater than 0 and less than or equal to 90 deg")
  if(is.na(Zhyp) == TRUE | Zhyp < 0)
    stop("Zhyp must be a non-negative number; use NA if unknown.")

  # Calculate Ztor from Zhyp, by assuming that the
  # hypocenter is located 60% down the fault width
  Ztor <- max(Zhyp - 0.6*W*sin(dip*pi/180), 0)
  return(Ztor)
}





# 6. DETERMINATION OF DIP ANGLE

dip.calc <- function(rake) {

  # Check input
  if(is.na(rake) == TRUE | abs(rake) > 180)
    stop("rake angle must be between -180 and 180, inclusive")

  # Strike-slip fault
  if(abs(rake) < 30 | abs(rake) > 150){
    dip <- 90

  # Normal fault
  } else{
    if(rake <= -30 & rake >= -150){
      dip <- 50

    # Reverse fault
    } else{
      if(rake >= 30 & rake <= 150){
        dip <- 40
      }
    }
  }
  return(dip)
}




# 7. CALCULATION OF DOWN-DIP RUPTURE WIDTH, W

W.calc <- function(M, rake){

  # Check input
  if(is.na(rake) == TRUE | abs(rake) > 180)
    stop("rake angle must be between -180 and 180, inclusive")
  if(is.na(M) == TRUE | M < 0)
    stop("M must be a positive number.")
  
  # Strike-slip fault
  if(abs(rake) < 30 | abs(rake) > 150){
    W <- 10^(-0.76+0.27*M)

  # Normal fault
  } else{
    if(rake <= -30 & rake >= -150){
      W <- 10^(-1.14+0.35*M)

    # Reverse fault
    } else{
      if(rake >= 30 & rake <= 150){
        W <- 10^(-1.61+0.41*M)
      }
    }
  }
  return(W)
}





# 8. CALCULATION OF DEPTH PARAMETER, Z2.5

Z2.5.calc <- function(Vs30 = NA, Z1.0 = NA, Z1.5 = NA){

  # Check input
  if(is.na(Vs30) == FALSE & Vs30 < 0)
    stop("Vs30 must be a positive number")
  if(is.na(Z1.0) == FALSE & Z1.0 < 0)
    stop("Z1.0 must be a positive number")
  if(is.na(Z1.5) == FALSE & Z1.5 < 0)
    stop("Z1.5 must be a positive number")
  if(is.na(Vs30) == TRUE & is.na(Z1.0) == TRUE & is.na(Z1.5) == TRUE)
    stop("Either Vs30, Z1.0, or Z1.5 must be specified.")

  # First choice: calculate from Z1.5 if available
  if(is.na(Z1.5) == FALSE){
    Z2.5 <- 636 + 1.549*Z1.5

  # Second choice: calculate from Z1.0
  } else{
    if(is.na(Z1.0) == FALSE){
      Z2.5 <- 519 + 3.595*Z1.0

  # Third choice: calculate from Vs30
    } else{
      if(is.na(Vs30) == FALSE){
        Z1.0 <- Z1.calc.as(Vs30 = Vs30)
        Z2.5 <- 519 + 3.595*Z1.0
      }
    }
  }

  return(Z2.5)
}


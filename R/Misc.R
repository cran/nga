# MISCELLANEOUS FUNCTIONS NEEDED FOR MULTIPLE NGA MODELS

# OUTLINE OF CODE:
# This code contains functions that are needed for multiple NGA models, including:
# 1. Interpolation of spectral acceleration
# 2. Calculation of distance measures (Rx and Rrup)
# 3. Calculation of depth parameter, Z1.0
# 4. Calculation of depth to top of rupture, Ztor
# 5. Calculation of dip angle
# 6. Calculation of down-dip width, W




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
    if(T <= 0.01 | T > 10){  # period not within allowable range
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
# For details on these derived equations, please see the following paper:
# Kaklamanos, J., and L. G. Baise (2010).  Relationships between Distance
# Measures in Recent Ground Motion Prediction Equations, Seismological
# Research Letters (in review).
# The equation numbers and figures in this code refer to the equation numbers
# and figures in Kaklamanos and Baise (2010).

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
  
  if(azimuth == 90){

      # Rjb > 0
      # Case 6 in Figure 3 (Eqn 5)
      if(Rjb > 0)
        Rx <- Rjb + W*cos(d)

      # Rjb = 0
      else if(Rjb == 0){
             
        # If Rrup is known...
        if(!(is.na(Rrup) == TRUE | Rrup < 0)){
          # Case 5A in Figure 3 (Eqn 6)
          if(Rrup < Ztor*sec(d))
             Rx <- sqrt(Rrup^2 - Ztor^2)
          # Case 5B (Eqn 7)
          else{
            if(Rrup >= Ztor*sec(d))
              Rx <- Rrup*csc(d) - Ztor*cot(d)
          }
        }

        # If Rrup is unknown, assume that the site is located at
        # the center of the surface projection of the ruptured area
        else
          Rx <- W*cos(d)/2
      }

    } else{
      if(azimuth >= 0 & azimuth <= 180 & azimuth != 90) {
        
        # Cases 2 and 8 (Eqn 3)
        if(Rjb*abs(tan(a)) <= W*cos(d)){
          Rx <- Rjb*abs(tan(a))

        # Cases 3 and 9 (Eqn 4)
        }else{
          if(Rjb*abs(tan(a)) > W*cos(d)){
            Rx <- Rjb*tan(a)*cos(a - asin(W*cos(d)*cos(a)/Rjb)) }
        }
      
      # Cases 1, 4, and 7 (Eqn 8)
      } else{
        if(azimuth >= -180 & azimuth < 0){
          Rx <- Rjb*sin(a) }
      }
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
    
  # Calculate Rrup'
  if(dip != 90)
    {
      # Zone A in Figure 5 (Eqn 10)
      if(Rx < Ztor*tan(d)){
        Rrup.prime <- sqrt(Rx^2 + Ztor^2)

      # Zone B (Eqn 11)
      } else{
        if(Rx >= Ztor*tan(d) & Rx <= Ztor*tan(d) + W*sec(d)){
          Rrup.prime <- Rx*sin(d) + Ztor*cos(d)

        # Zone C (Eqn 12)
        } else{
          if(Rx > Ztor*tan(d) + W*sec(d)){
            Rrup.prime <- sqrt((Rx - W*cos(d))^2 + (Ztor + W*sin(d))^2)
          }
        }
      }
    }
  else if(dip == 90) # Eqn 16
    {Rrup.prime <- sqrt(Rx^2 + Ztor^2)}


  # Calculate Ry

  # Eqn 13
  if(azimuth == 90 | azimuth == -90){
    Ry <- 0
  
  # Eqn 14
  } else{
    if(azimuth == 0 | azimuth == 180 | azimuth == -180){
      if(is.na(Rjb) == TRUE | Rjb < 0)
        stop("Rjb must be a positive number")
      else
        Ry <- Rjb
      
    # Eqn 15  
    } else{
      Ry <- abs(Rx*cot(a))
    }
  }

  # Calculate Rrup (Eqn 9)
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
  exp(LnZ1.0)
}


# CY08 [Eqn 1 in Chiou and Youngs (2008)]
Z1.calc.cy <- function(Vs30) {

  # Check input
  if(is.na(Vs30) == TRUE | Vs30 < 0)
    stop("Vs30 must be a positive number")
  
  # Calculate Z1.0
  LnZ1.0 <- 28.5 - (3.82/8)*log(Vs30^8 + 378.7^8)
  exp(LnZ1.0)
}





# 4. CALCULATION OF DEPTH TO TOP OF RUPTURE, Ztor

Ztor.calc <- function(Zhyp = NA, W, dip, M = NA, rake = NA) {

  # Check input
  if(is.na(W) == TRUE | is.na(dip) == TRUE |
     (is.na(Zhyp) == TRUE & (is.na(M) == TRUE | is.na(rake) == TRUE)))
    stop("W and dip must be specified, as well as either (1) Zhyp, or (2) M and rake,",
         "\n", "which are used to estimate Zhyp if Zhyp is unspecified.")

  # Estimate Zhyp from M and rake if Zhyp is unspecified, using
  # correlations in Scherbaum et al (2004)
  if(is.na(Zhyp) == TRUE){

    # Strike-slip fault (SS)
    if(abs(rake) < 30 | abs(rake) > 150){
      Zhyp <- 5.63 + 0.68*M

    # Shallow-dipping fault (SHD)  
    } else{
      Zhyp <- 11.24 - 0.2*M
    }
  }

  # Calculate Ztor from Zhyp, by assuming that the
  # hypocenter is located 60% down the fault width
  Ztor <- max(Zhyp - 0.6*W*sin(dip*pi/180), 0)
  

  # ALTERNATE METHOD (currently not implemented)
  # If Zhyp has not been provided, estimate Ztor from M using the basic
  # equation derived from the correlation suggested by Campbell et al. (2009)
  #  } else {
  #    if(M < 8){
  #      Ztor <- 0.5*M^2 - 8.5*M + 36
  #    } else{
  #      if(M >= 8){
  #        Ztor <- 0
  #      }
  #    }
  #  }
    
  return(Ztor)
}




  


# 5. DETERMINATION OF DIP ANGLE

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
      dip <- 55

    # Reverse fault
    } else{
      if(rake >= 30 & rake <= 150){
        dip <- 40
      }
    }
  }
  return(dip)
}




# 6. CALCULATION OF DOWN-DIP RUPTURE WIDTH, W

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


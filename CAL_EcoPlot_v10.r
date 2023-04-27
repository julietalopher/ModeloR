# Soil Water model function for automatic calibration#
SWMc <- function (         # Model parameter definition
  rE = -0.19,           # extinction factor
  alpha = 4.41,          # interception threshold parameter, default set to 4.5 mm (Simunek et al, 2008)
  Smax = 586,             # maximum measured soil moisture content upper horizon
  Ic = 57.6,              # max Infiltration capacity, if exceeding surface runoff Qs generation
  ks1 = 11.7,						# saturated hydraulic conductivity upper soil horizon
  ks2 = 18.8,					# saturated hydraulic conductivity lower horizon
  ks3 = 25,             # saturated hydraulic conductivity deepest horizon
  GWmax = 918,					# maximum measured soil moisture content lower horizon		
  Lmax = 2200,           # maximum measured soil moisture content deeper horizon
  g1 = 2.8,            # nonlinear scaling parameter of upper horizon, if g1 = 1, linear case	
  g2 = 4.6,       # nonlinear scaling parameter of lower horizon, if g2 = 1, linear case
  g3 = 6,             # nonlinear scaling parameter of deeper horizon, if g2 = 1, linear case
  PF_Scale = 0.40     # preferential flowpath parameter
)            				

{  ##### BEGIN OF FUNCTION BODY #####
  
  nrun <- nrow (inp)     # model time step equal to input
  
  ## Dynamic meteorological and interception variables
  
  PN <- rep(NA,nrun)	            # net precipitation minus interception threshold
  I <- rep(NA,nrun)               # interception = 1 - PN
  D <- rep(NA,nrun)               # interception threshold 
  IMax <- rep(NA,nrun)            # time-varying maximum canopy storage capacity - AJN
  SCF <- rep(NA,nrun)             # surface cover fraction 
  Ei	<- rep(NA,nrun)							# interception evaporation
  Tr <- rep(NA,nrun)              # transpiration rate
  Es	<- rep(NA,nrun)							# soil evaporation
  Ep	<- rep(NA,nrun)							# potential evaporation
  Tp	<- rep(NA,nrun)							# potential transpiration
  Th <- rep(NA, nrun)              # throughfall
  Tr_Upper <- rep(NA, nrun) 
  Tr_Lower <- rep(NA, nrun)        # Transpiration taken from the lower store
  Tr_Deep <- rep(NA, nrun) 
  
  ## Iteratively filled storage and flux variables
  STO <- rep(NA,nrun)					# upper soil horizon storage in mm
  GW <- rep(NA,nrun)					# lower soil horizon storage in mm
  Sdeep <- rep(NA,nrun)       # deep soil horizon storage in mm
  GWflow <- rep(NA,nrun)			# recharge from lower soil horizon in mm
  Recharge <- rep(NA,nrun)		# recharge loss from Sdeep
  Perc <- rep(NA,nrun)        # vertical flux from STO into GW in mm
  Qs  <- rep(NA,nrun)         # Overland flow if infiltration capacity parameter Ic is exceeded in mm
  API  <- rep(NA,nrun)         # API vector
  Pref_Flow<- rep(NA, nrun)   # Preferential flow occuring in larger Net Precipitation events

  #~~~~~~~~~~~~~~ calculate iterative dependent time series ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for (i in 1:nrun) {					# the loop iterates through the length of the input time series
    
    #~~~~~~~~~~~~~~ calculate dynamic dependent time series ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ## storage initialization
    
    if (i==1) {							# conditional: filling of storage variables with initial volumes...
      API[i] <- 0.995*2.6      # set the first api value to the last rainfall 
      I[i] <- (alpha*LAI[i]*0.5) # interception store initialised half full for current LAI??
      STO[i] <- 500					# measured start value				
      GW[i] <- 700            # measured start value
      Sdeep[i] <- 3000
    }
    else {							# ... or i-1 storage fillings
      API[i] <- API[i-1]
      I[i] <- I[i-1]
      STO[i] <- STO[i-1]
      GW[i] <- GW[i-1]
      Sdeep[i] <- Sdeep[i-1]
    }
    
    API[i] <- 0.995*API[i] + P[i]        # daily decay constant of 0.95 by Hill et al. (2014)
    
    if ( is.na(P[i]) || is.na(PET[i]) || is.na(LAI[i]))   # no input data NAs allowed! - AJN updated; only need P and PET here now, but LAI added too...
    {
      P[i] <- 0
      PET[i] <- 0
      LAI[i] <- 0
    }
    
    # taken outside else statement so NA not assigned to Ep and TP - PN also defined later now... 
    SCF[i] <- 1 - exp(rE * LAI[i])     # surface cover fraction, rExtinct = -0.463 Rutter (1972)
    
    Tp[i] <- (PET[i] * SCF[i])            # fraction potential canopy ET (Ei plus Tr when I = 0) AJN
    Ep[i] <- (PET[i] * (1-SCF[i]))        # fraction potential soil evaporation only AJN
    
    #~~ Interception, simplified Simunek et al (2008) model ##~~#~~#~~
    # Available canopy storage as maximum storage for the current LAI (alpha * LAI[i]) minus the depth already filled
    IMax[i] <- (alpha * LAI[i])
    
    IAvail <- IMax[i] - I[i] 
    
    # Check if change in LAI reduces available storage below current depth in storage (i.e. IAvail < 0)
    if(IAvail<0){
      
      # Excess water due to change in LAI added to throughfall (MIXING!)
      Th[i] <- I[i] - IMax[i]
      
      # Then no available interception storage, interception store is set to max for this timestep, no interception
      IAvail <- 0
      I[i] <- IMax[i]
      D[i] <- 0
      
    } else if(IAvail == 0){
      
      # If there's exactly 0 available storage - can happen if LAI constant and PET = 0...
      Th[i] <- 0        # No throughfall from max storage capacity change
      D[i] <- 0         # No interception
      
      I[i] <- IMax[i]   # Interception store is at maximum (already was the case!)
      
    } else{
      
      # Calculate the interception based on available storage
      D[i] <- (IAvail) * (1 - (1 / (1 + ((SCF[i] * P[i])/(IAvail)))))    
      
      # Update interception store
      I[i] <- I[i] + D[i]
      
      # No throughfall
      Th[i] <- 0
      
    }
    
    # PN will be P plus excess storage Th if LAI reduces available storage below depth already stored (MIXING IMPLICATIONS)
    # PN will equal P minus interception D if there is available storage, and Th will be 0.
    
    PN[i] <- (P[i] - D[i]) + Th[i]            # net Precip input to storages 
    
    ## water balance calculations
    if(Ep[i]<= I[i]){
      Ei[i]<- Ep[i]           # If sufficient water in store to meet full potential canopy ET (Tr), full volume taken from canopy
      I[i]<- I[i] - Ei[i]     # Interception store reduced by evaporated volume
      Ep[i]<- 0}              # Potential canopy ET now 0
    else{
      Ei[i]<- I[i]            # If potential canopy ET greater than canopy storage all of interception store evaporated
      I[i]<- 0                # Interception store therefore empty
      Tp[i]<- Tp[i] - Ei[i]   # Potential canopy  ET updated to reflect used volume - remaining volume available for transpiration
    }
    
    # Hortonion OLF
    if (PN[i] > Ic) {          #surface runoff Qs if Ic is exceeded
      Qs[i] <- PN[i] - Ic
      PN[i] <- Ic}
    else{
      Qs[i] <- 0
      PN[i] <- PN[i]}
    
    if(PN[i] > 6){                        # In larger net precipitation events (>P90) preferential flow allowed
      Pref_Flow[i]<- PN[i] * PF_Scale     # Fraction of net precipitation constituting preferential flow determined by calibrated parameter - added below to Recharge[i]
      PN[i] <- PN[i] * (1 - PF_Scale)}    # Net precipitation to upper soil box reduced accordingly
    else{
      Pref_Flow[i]<- 0
      PN[i]<- PN[i] }
    
    STO[i] <- STO[i] + PN[i]				# upper soil horizon storage plus net precipitation
    
    if(is.na(STO[i]))
    {STO[i] <- 0}
    
    if (STO[i] < Smax) {
      abs(Tr_Upper[i] <- (STO[i] / Smax) * Tp[i]) # Transpiration from upper store calculated
      STO[i] <- STO[i] - Tr_Upper[i]              # Upper store reduced accordingly
      Tp[i]<- Tp[i] - Tr_Upper[i]}                # Potential transpiration reduced accordingly
    else {                                        # If store equalling or greater than Smax then all potenital transpiration used.
      Tr_Upper[i] <- Tp[i]                        # Cannot use (STO[i] / Smax) as this will create number >1 and so transpiration would be greater than potential
      STO[i]<- STO[i] - Tr_Upper[i]
      Tp[i]<- 0}   
    
    if (is.na(Tr_Upper[i]))                       # Code to ensure no NA or negative transpiration values
    {Tr_Upper[i] <- 0}
    if (Tr_Upper[i] < 0) 
    { Tr_Upper[i] <- Tr_Upper[i] * -1} 
    else {Tr_Upper[i] <- Tr_Upper[i]}
    
    if(Ep[i] == 0){                     # If no potential evaporation voulme, soil evaporation cannot occur
      Es[i] <- 0}
    else{
      if(STO[i] > Ep[i]){                 # If greater soil storage than potential evporation:
        Es[i]<- (STO[i] / Smax) * Ep[i]   # Soil evaporation is calculated using full potential evaporation volume as multiplier
        STO[i] <- STO[i] - Es[i]
        Ep[i]<- Ep[i] - Es[i]}	
      else{                               # If less soil storage than potential evaporation:
        Es[i]<- (STO[i] / Smax) * STO[i]  # Soil evaporation is calculated using a potential evaporation volume reduced to be identical to the volume of water in storage (therefore == STO[i])
        STO[i] <- STO[i] - Es[i]
        Ep[i]<- Ep[i] - Es[i]}
      if (API[i] < 6.71)   # only Es if median API is exceeded, else Es becomes Tr
      {Tr_Upper[i] <- Tr_Upper[i] + Es[i]
      Es[i] <- 0}
      else {Tr_Upper[i] <- Tr_Upper[i]
      Es[i] <- Es[i]}
    }
    
    if (PN[i] > 0)
    { Perc[i] <- (ks1 * (STO[i] / Smax) ^ g1)				# nonlinear calculation of percolation to lower soil storage 
    STO[i] <- STO[i] - Perc[i]	}
    else {Perc[i] <- 0}
    
    if(is.na(GW[i]))
    {GW[i] <- 0}
    if (is.na(Perc[i]))
    {Perc[i] <- 0}
    
    GW[i] <- GW[i] + Perc[i] + Pref_Flow[i]       # fill lower soil horizon
    
    if(Tp[i] <= 0){                                    # If all potential transpiration taken from upper store no volume can be taken from lower store (<= to be on safe side)
      Tr_Lower[i]<- 0}
    else{                                              
      if (GW[i] < GWmax) {
        abs(Tr_Lower[i] <- (GW[i] / GWmax) * Tp[i])      # If potenital transpiration volume available the lower store transpiration volume calculated as per upper store
        GW[i] <- GW[i] - Tr_Lower[i]
        Tp[i]<- Tp[i] - Tr_Lower[i]}                     
      else {                                           
        Tr_Lower[i] <- Tp[i]                           # If GW[i] exceeds GWmax then full potential volume removed (using GW[i]/ GWmax would create value > potenital)
        GW[i]<- GW[i] - Tr_Lower[i]
        Tp[i]<- 0}   
      
      if (is.na(Tr_Lower[i]))                          # Code to ensure no NA or negative transpiration values
      {Tr_Lower[i] <- 0}
      if (Tr_Lower[i] < 0) 
      { Tr_Lower[i] <- Tr_Lower[i] * -1} 
      else {Tr_Lower[i] <- Tr_Lower[i]}
    }

    if ( Perc[i] > 0 ) {
      GWflow[i] <- ks2 * (GW[i] / GWmax) ^ g2					# nonlinear recharge calculation based on soil storage volume
      GW[i] <- GW[i] - GWflow[i]}
    
    if ( is.na(GWflow[i]))  						
    {
      GWflow[i] <- 0
    }
    else {GWflow[i] <-  GWflow[i]} 			# simulated percolation
    
    if(is.na(Sdeep[i]))
    {Sdeep[i] <- 0}
    
    Sdeep[i] <- Sdeep[i] + GWflow[i]        # fill lower soil horizon
    
    if(Tp[i] <= 0){                                    # If all potential transpiration taken from upper store no volume can be taken from lower store (<= to be on safe side)
      Tr_Deep[i]<- 0}
    else{                                              
      if (Sdeep[i] < Lmax) {
        abs(Tr_Deep[i] <- (Sdeep[i] / Lmax) * Tp[i])      # If potenital transpiration volume available the lower store transpiration volume calculated as per upper store
        Sdeep[i] <- Sdeep[i] - Tr_Deep[i]
        Tp[i]<- Tp[i] - Tr_Deep[i]}                     
      else {                                           
        Tr_Deep[i] <- Tp[i]                           # If GW[i] exceeds GWmax then full potential volume removed (using GW[i]/ GWmax would create value > potenital)
        Sdeep[i]<- Sdeep[i] - Tr_Deep[i]
        Tp[i]<- 0}   
      
      if (is.na(Tr_Deep[i]))                          # Code to ensure no NA or negative transpiration values
      {Tr_Deep[i] <- 0}
      if (Tr_Deep[i] < 0) 
      { Tr_Deep[i] <- Tr_Deep[i] * -1} 
      else {Tr_Deep[i] <- Tr_Deep[i]}
    }
    Tr[i]<- Tr_Upper[i] + Tr_Lower[i] + Tr_Deep[i]                # Total transpiration volume 
    
    if ( GWflow[i] > 0 ) {
      Recharge[i] <- ks3 * (Sdeep[i] / Lmax) ^ g3					# nonlinear recharge calculation based on soil storage volume
      Sdeep[i] <- Sdeep[i] - Recharge[i]}
    
    if ( is.na(Recharge[i]))  						
    {
      Recharge[i] <- 0
    }
    else {Recharge[i] <- Recharge[i]} 			# simulated recharge
  }  
  # objective functions for calibration
  Qtemp1  <- na.omit(data.frame(SW_30[597:2920], STO[597:2920]))								        # create time frame without missing discharge time steps
  Q_Rsqu <- cor(Qtemp1$SW_30,Qtemp1$STO)			                                          # calculate the Q correlation coefficient r
  Q_var  <-  (sd(Qtemp1$STO)/mean(Qtemp1$STO)) / (sd(Qtemp1$SW_30)/mean(Qtemp1$SW_30))    # variability
  Q_b    <-  mean(Qtemp1$STO)/mean(Qtemp1$SW_30)                                        # bias
  KGE1  <- (1- (sqrt((Q_Rsqu-1)^2 + (Q_var-1)^2 + (Q_b -1)^2)))                    # modified Kling-Gupta efficiency (Gupta et al, 2009; Kling et al., 2012) used to minimize
  MAE1 <- mean(sum(abs(Qtemp1$STO - Qtemp1$SW_30)))
  
  Qtemp2  <- na.omit(data.frame(SW_120[597:2920], GW[597:2920]))								        # create time frame without missing discharge time steps
  Q2_Rsqu <- cor(Qtemp2$SW_120,Qtemp2$GW)			                                          # calculate the Q correlation coefficient r
  Q2_var  <-  (sd(Qtemp2$GW)/mean(Qtemp2$GW)) / (sd(Qtemp2$SW_120)/mean(Qtemp2$SW_120))    # variability
  Q2_b    <-  mean(Qtemp2$GW)/mean(Qtemp2$SW_120)                                        # bias
  KGE2  <- (1- (sqrt((Q2_Rsqu-1)^2 + (Q2_var-1)^2 + (Q2_b -1)^2)))                    # modified Kling-Gupta efficiency (Gupta et al, 2009; Kling et al., 2012) used to minimize
  MAE2 <- mean(sum(abs(Qtemp2$GW - Qtemp2$SW_120)))
  
  Qtemp3  <- na.omit(data.frame(SW_570[597:2920], Sdeep[597:2920]))								        # create time frame without missing discharge time steps
  Q3_Rsqu <- cor(Qtemp3$SW_570,Qtemp3$Sdeep)			                                          # calculate the Q correlation coefficient r
  Q3_var  <-  (sd(Qtemp3$Sdeep)/mean(Qtemp3$Sdeep)) / (sd(Qtemp3$SW_570)/mean(Qtemp3$SW_570))    # variability
  Q3_b    <-  mean(Qtemp3$Sdeep)/mean(Qtemp3$SW_570)                                        # bias
  KGE3  <- (1- (sqrt((Q3_Rsqu-1)^2 + (Q3_var-1)^2 + (Q3_b -1)^2)))                    # modified Kling-Gupta efficiency (Gupta et al, 2009; Kling et al., 2012) used to minimize
  MAE3 <- mean(sum(abs(Qtemp3$Sdeep - Qtemp3$SW_570)))
  
  Qtemp4  <- na.omit(data.frame(Tobs[597:2920], Tr[597:2920]))								        # create time frame without missing discharge time steps
  Q4_Rsqu <- cor(Qtemp4$Tobs,Qtemp4$Tr)			                                          # calculate the Q correlation coefficient r
  Q4_var  <-  (sd(Qtemp4$Tr)/mean(Qtemp4$Tr)) / (sd(Qtemp4$Tobs)/mean(Qtemp4$Tobs))    # variability
  Q4_b    <-  mean(Qtemp4$Tr)/mean(Qtemp4$Tobs)                                        # bias
  KGE4  <- (1- (sqrt((Q4_Rsqu-1)^2 + (Q4_var-1)^2 + (Q4_b -1)^2)))                    # modified Kling-Gupta efficiency (Gupta et al, 2009; Kling et al., 2012) used to minimize
  MAE4 <- mean(sum(abs(Qtemp4$Tr - Qtemp4$Tobs)))
  
  Qtemp5  <- na.omit(data.frame(Es_obs[597:2920], Es[597:2920]))								        # create time frame without missing discharge time steps
#  Q5_Rsqu <- cor(Qtemp5$Es_obs,Qtemp5$Es)			                                          # calculate the Q correlation coefficient r
  Q5_var  <-  (sd(Qtemp5$Es)/mean(Qtemp5$Es)) / (sd(Qtemp5$Es_obs)/mean(Qtemp5$Es_obs))    # variability
  Q5_b    <-  mean(Qtemp5$Es)/mean(Qtemp5$Es_obs)                                        # bias
#  KGE5  <- (1- (sqrt((Q5_Rsqu-1)^2 + (Q5_var-1)^2 + (Q5_b -1)^2)))                    # modified Kling-Gupta efficiency (Gupta et al, 2009; Kling et al., 2012) used to minimize
  MAE5 <- mean(sum(abs(Qtemp5$Es - Qtemp5$Es_obs)))
  
  Qtemp6  <- na.omit(data.frame(Robs[597:2920], Recharge[597:2920]))								        # create time frame without missing discharge time steps
  Q6_Rsqu <- cor(Qtemp6$Robs,Qtemp6$Recharge)			                                          # calculate the Q correlation coefficient r
  Q6_var  <-  (sd(Qtemp6$Recharge)/mean(Qtemp6$Recharge)) / (sd(Qtemp6$Robs)/mean(Qtemp6$Robs))    # variability
  Q6_b    <-  mean(Qtemp6$Recharge)/mean(Qtemp6$Robs)                                        # bias
  KGE6  <- (1- (sqrt((Q6_Rsqu-1)^2 + (Q6_var-1)^2 + (Q6_b -1)^2)))                    # modified Kling-Gupta efficiency (Gupta et al, 2009; Kling et al., 2012) used to minimize
  MAE6 <- mean(sum(abs(Qtemp6$Recharge - Qtemp6$Robs)))
  
  Qtemp7  <- na.omit(data.frame(PNobs[597:2920], PN[597:2920]))								        # create time frame without missing discharge time steps
#  Q7_Rsqu <- cor(Qtemp7$PNobs,Qtemp7$PN)			                                          # calculate the Q correlation coefficient r
  Q7_var  <-  (sd(Qtemp7$PN)/mean(Qtemp7$PN)) / (sd(Qtemp7$PNobs)/mean(Qtemp7$PNobs))    # variability
  Q7_b    <-  mean(Qtemp7$PN)/mean(Qtemp7$PNobs)                                        # bias
#  KGE7  <- (1- (sqrt((Q7_Rsqu-1)^2 + (Q7_var-1)^2 + (Q7_b -1)^2)))                    # modified Kling-Gupta efficiency (Gupta et al, 2009; Kling et al., 2012) used to minimize
  MAE7 <- mean(sum(abs(Qtemp7$PN - Qtemp7$PNobs)))
   
  return (cbind(KGE1, Q_Rsqu, MAE1,KGE2, Q2_Rsqu, MAE2, Q3_Rsqu, KGE3, MAE3, 
                Q4_Rsqu, KGE4, MAE4, MAE5,Q6_Rsqu, KGE6, MAE6, MAE7))                  # return evaluation criteria for multi-objective calibration
}     ##### END OF FUNCTION BODY #####

#    Embedding RCT Health Economic Analysis using the Sheffield Type 2 Diabetes Treatment Model - version 3
#    Copyright (C) 2023   Pollard, Pidd, Breeze, Brennan, Thomas

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#    Contact person: Dan Pollard, Email: d.j.pollard@sheffield.ac.uk, 
#    Address: Regent Court, 30 Regent Court, Sheffield, United Kingdom, S1 4DA


####'This function runs the model for a given set of patients and parameters
####'@param population_, is the population matrix
####'@param parameters_, is the full parameters matrix
####'@param endtime_, is the number of years to run the simulation for
####'@param GlobalVars_, is the matrix giving the global variables
####'@param random_numbs_, is an array of common random numbers giving a random number
####'@param LifeTables_, is a dataframe containing life table information in formate
####'that can easily be matched to the population matrix
####'@param SOUR_, is the current second order uncertainty run
####'draw for each patient in each year for every event in which a random number 
####'is required 
####'@return results, is the results matrix
####'@return psaresults, is a summary set of results to produce in the PSA

run_simulation <- function(population_, parameters_, endtime_, treatment_, GlobalVars_, random_numbs_,LifeTables_, SOUR_){
  
  ##reduce the parameters matrix down to the correct row
  if(GlobalVars_["run_psa","Value"]==F){
    parameters_ <- parameters_[1,]
  }else{
    parameters_ <- parameters_[SOUR_+1,]
  }
  
  #Use the correct trajectories
  if(GlobalVars_["Trajectory","Value"]=="Constant"){
  
  ##Create the underlying trajectory matrices - for constant values
  HBA1c_underlying  <- constant_A1c(population_,endtime_)
  BMI_underlying    <- constant_BMI(population_,endtime_)
  SBP_underlying    <- constant_SBP(population_,endtime_)
  HDL_underlying    <- constant_HDL(population_,endtime_)
  LDL_underlying    <- constant_LDL(population_,endtime_)
  HEARTR_underlying <- constant_HEARTRATE(population_,endtime_)
  WBC_underlying    <- constant_WBC(population_,endtime_)
  HAEM_underlying   <- constant_HAEM(population_,endtime_)  
  eGFR_underlying   <- constant_eGFR(population_, endtime_)
    
  }else{
  ##Create the underlying UKPDS trajectory matrices
  HBA1c_underlying  <- UKPDS_90_contrisk_A1c(population_,parameters_,endtime_)
  BMI_underlying    <- UKPDS_90_contrisk_BMI(population_,parameters_,endtime_)
  SBP_underlying    <- UKPDS_90_contrisk_SBP(population_,parameters_,endtime_)
  HDL_underlying    <- UKPDS_90_contrisk_HDL(population_,parameters_,endtime_)
  LDL_underlying    <- UKPDS_90_contrisk_LDL(population_,parameters_,endtime_)
  HEARTR_underlying <- UKPDS_90_HEARTR(population_,parameters_,endtime_)
  WBC_underlying    <- UKPDS_90_WBC(population_,parameters_,endtime_)
  HAEM_underlying   <- UKPDS_90_HAEM(population_,parameters_,endtime_)
  
  }
  
  #Place to add Intervention effects
  attend_se <- matrix(data = 0, nrow = length(population_[,"ID"]), ncol = 3 ) # not embedding, haven't yet edited costs, so just set attend_se to 0
  attend_se[,1] <- 1:length(population_[,"ID"]) #change the first column to a row ID
  HBA1c_INTV        <- initialise_intervention_dt_HbA1c(length(population_[,"ID"]),treatment_,parameters_,endtime_,GlobalVars_,attend_se)
  BMI_INTV          <- initialise_intervention_dt_BMI(length(population_[,"ID"]),treatment_,parameters_,endtime_,GlobalVars_,attend_se)
  SBP_INTV          <- initialise_intervention_dt_SBP(length(population_[,"ID"]),treatment_,parameters_,endtime_,GlobalVars_,attend_se)
  HDL_INTV          <- initialise_intervention_dt_HDL(length(population_[,"ID"]),treatment_,parameters_,endtime_,GlobalVars_,attend_se)
  LDL_INTV          <- initialise_intervention_dt_LDL(length(population_[,"ID"]),treatment_,parameters_,endtime_,GlobalVars_,attend_se)
  
  #start year at 0
  year <- 0
  
  #initialise the results matrix
  results <- GenerateResultsMatrix(GlobalVars_, endtime_)
  #initialise a PSA results matrix, if it is a PSA
  if(GlobalVars_["run_psa","Value"]==T){
  psaresults <- matrix(data=NA, nrow = as.numeric(GlobalVars_["psa_count", "Value"]), ncol = 5)  
  colnames(psaresults) <- c("Life Years per Patient", "Undiscounted QALYs per Patient", 
                            "QALYs per Patient", "Undiscounted Costs per Patient",
                            "Discounted Costs per Patient")
  }
  
  #run the model up to the specified end time or so long as at least one person is alive
  while (year < endtime_ & 
         (sum(is.na(population_[,"F_ALLCAUSE"]))) >= 1){
  #Get a logical vector indicating if people are alive, use the fact that FOTH
  #in the population matrix is NA if alive
  alive <- is.na(population_[,"F_ALLCAUSE"])
  #create matrix of ture = dead too for unit tests
  dead <- is.na(population_[,"F_ALLCAUSE"])==F
  
  #Estimate Diabetes Related complications and all cause deaths for this year
  population_ <- update_events_UKPDS82(population_,parameters_, treatment_, year, alive, random_numbs_, LifeTables_)
  #Estimate PVD,ATFIB,MMALB
  population_ <- update_events_UKPDS90(population_,parameters_, year, alive,random_numbs_, GlobalVars_)
  #Estimate Depression
  population_ <- update_events_SPHR_depression(population_,parameters_,year,alive,random_numbs_)
  #Estimate Osteoarthritis
  population_ <- update_events_SPHR_osteoarthritis(population_,parameters_,year,alive,random_numbs_)
  #Estimate Cancer incidence
  population_ <- update_events_SPHR_cancer(population_, parameters_,year, alive,random_numbs_)
  
  #Stop the model if dead people have events
  if(sum(population_[,"MI_E"][dead])  !=0|
     sum(population_[,"MI2_E"][dead]) != 0|
     sum(population_[,"AMP_E"][dead]) != 0| 
     sum(population_[,"AMP2_E"][dead]) != 0|
     sum(population_[,"BLIND_E"][dead]) != 0|
     sum(population_[,"ULCER_E"][dead]) != 0|
     sum(population_[,"RENAL_E"][dead]) != 0|
     sum(population_[,"CHF_E"][dead]) != 0|
     sum(population_[,"IHD_E"][dead]) != 0|
     sum(population_[,"STRO_E"][dead]) != 0|
     sum(population_[,"STRO2_E"][dead]) != 0|
     sum(population_[,"ATFIB_E"][dead]) != 0|
     sum(population_[,"PVD_E"][dead]) != 0|
     sum(population_[,"MMALB_E"][dead]) != 0|
     sum(population_[,"CANB_E"][dead]) != 0|
     sum(population_[,"CANC_E"][dead]) != 0|
     sum(population_[,"DEP_E"][dead]) != 0|
     sum(population_[,"OST_E"][dead]) != 0){
    #stop the model if the number of people with dead people are recorded as
    #having an event
    #push everything to the global enviroment for bug checking
    for (variable in ls()) {
      assign(variable, get(variable), envir = .GlobalEnv)
    }
    
    stop("dead people are getting events")
  }
  
  ##QALYs
  population_ <- calculate_QALYs(population_, parameters_,  year, alive, GlobalVars_)
  ##Costs
  population_ <- calculate_costs(population_, parameters_, year, alive, GlobalVars_,treatment_,attend_se)
  
  #Record results
  results <- GenerateDetailedresults(results,population_, year, alive, GlobalVars_)
  
  
  #update histories
  population_ <- update_history(population_,
                                HBA1c_underlying,
                                BMI_underlying,
                                SBP_underlying,    
                                HDL_underlying,
                                LDL_underlying, 
                                HEARTR_underlying,
                                WBC_underlying,
                                HAEM_underlying,
                                HBA1c_INTV,
                                BMI_INTV,
                                SBP_INTV,
                                HDL_INTV,
                                LDL_INTV,
                                year,
                                GlobalVars_,
                                eGFR_underlying)
  
  population_ <- update_patchars(population_, parameters_, alive)
  
  #Unit test
  #Stop the model if there are events still in the population matrix
 if(sum(population_[,"MI_E"])  !=0|
     sum(population_[,"MI2_E"]) != 0|
     sum(population_[,"AMP_E"]) != 0| 
     sum(population_[,"AMP2_E"]) != 0|
     sum(population_[,"BLIND_E"]) != 0|
     sum(population_[,"ULCER_E"]) != 0|
     sum(population_[,"RENAL_E"]) != 0|
     sum(population_[,"CHF_E"]) != 0|
     sum(population_[,"IHD_E"]) != 0|
     sum(population_[,"STRO_E"]) != 0|
     sum(population_[,"STRO2_E"]) != 0){
   
     #push everything to the global enviroment for bug checking
     for (variable in ls()) {
       assign(variable, get(variable), envir = .GlobalEnv)
     }
   #stop the model if the events are not reset to 0
    stop("events are not set to zero before reset")
  }
  
  year <- year + 1
  
  }

  
  #For now return the Detailed results table or the population matrix if the run is deterministic
  if(GlobalVars_["Results_output", "Value"] == "Summary"&
     GlobalVars_["run_psa", "Value"]==T){
    psaresults <- matrix(data=NA,nrow=1,ncol=24)
    #Life Years
    psaresults[,1] <- sum(results["Undiscounted life years accrued",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,2] <- sum(results["Discounted life years accrued",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,3] <- sum(results["Undiscounted QALYs",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,4] <- sum(results["Discounted QALYs",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,5] <- sum(results["Undiscounted Costs",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,6] <- sum(results["Discounted Costs",],na.rm=TRUE)/length(population_[,"ID"])
    #10 year histories
    psaresults[,7] <- results["1st MI Hist",11]
    psaresults[,8] <- results["2nd MI Hist",11]
    psaresults[,9] <- results["1st Stroke Hist",11]
    psaresults[,10] <- results["2nd Stroke Hist",11]
    psaresults[,11] <- results["CHF Hist",11]
    psaresults[,12] <- results["IHD Hist",11]
    psaresults[,13] <- results["Blindness Hist",11]
    psaresults[,14] <- results["Ulcer Hist",11]
    psaresults[,15] <- results["1st Amputation Hist",11]
    psaresults[,16] <- results["2nd Amputation Hist",11]
    psaresults[,17] <- results["Renal Failure Hist",11]
    psaresults[,18] <- results["PVD Hist",11]
    psaresults[,19] <- results["MMALB Hist",11]
    psaresults[,20] <- results["ATFIB Hist",11]
    psaresults[,21] <- results["Breast Cancer Hist",11]
    psaresults[,22] <- results["Colorectal Cancer Hist",11]
    psaresults[,23] <- results["Depression Hist",11]
    psaresults[,24] <- results["Osteoarthritis Hist",11]
    #delete the original results matrix
    rm(results)
    return(psaresults)
  }else if(GlobalVars_["Results_output","Value"]=="Patient Level"&
     GlobalVars_["run_psa","Value"]==F){#option to produce the patient 
    #charateristics matrix for checking results stability for number of patients
    return(population_)
  }else{#default is the 
  return(results)
  }
}

####'This function runs the model for a given set of patients and parameters
####'It is designed to run the simulation in MtHood2025 Challenge 1, where every patient
####'follows the same specified time path
####'@param population_, is the population matrix
####'@param parameters_, is the full parameters matrix
####'@param endtime_, is the number of years to run the simulation for
####'@param GlobalVars_, is the matrix giving the global variables
####'@param random_numbs_, is an array of common random numbers giving a random number
####'@param LifeTables_, is a dataframe containing life table information in formate
####'that can easily be matched to the population matrix
####'@param SOUR_, is the current second order uncertainty run
####'draw for each patient in each year for every event in which a random number 
####'is required 
####'@return results, is the results matrix
####'@return psaresults, is a summary set of results to produce in the PSA

run_simulation_MtHood2025_C1_v2 <- function(population_, parameters_, endtime_, treatment_, GlobalVars_, random_numbs_,LifeTables_, SOUR_,MtHood2025C1_ ){
  
  ##reduce the parameters matrix down to the correct row
  if(GlobalVars_["run_psa","Value"]==F){
    parameters_ <- parameters_[1,]
  }else{
    parameters_ <- parameters_[SOUR_+1,]
  }
  
  #Use the trajectories
  A1c         <- MtHood2025C1_$A1c.csv
  BMI         <- MtHood2025C1_$BMI.csv
  SBP         <- MtHood2025C1_$SBP.csv
  SMO         <- MtHood2025C1_$SMO.csv
  HDL         <- MtHood2025C1_$HDL.csv
  LDL         <- MtHood2025C1_$LDL.csv
  PVD         <- MtHood2025C1_$PVD.csv
  ATF         <- MtHood2025C1_$AF.csv
  MMALB       <- MtHood2025C1_$MMALB.csv
  eGFR        <- MtHood2025C1_$eGFR.csv
  HEAM        <- MtHood2025C1_$HEAM.csv
  WBC         <- MtHood2025C1_$WBC.csv
  HR          <- MtHood2025C1_$HR.csv
  PVD         <- MtHood2025C1_$PVD.csv
  
  #If control we only want group 1 people, if intervention, we want group 2 and to change their ID to start at 1 again
  if(mean(population_[,"MALE"])==1){
    #If MALE only get the rows associated with being Male
    #remove the 1st row, which is text
    A1c       <- A1c[1:6,-1]
    BMI       <- BMI[1:6,-1]
    SBP       <- SBP[1:6,-1]
    SMO       <- SMO[1:6,-1]
    HDL       <- HDL[1:6,-1]
    LDL       <- LDL[1:6,-1]
    ATF       <- ATF[1:6,-1]
    MMALB     <- MMALB[1:6,-1]
    eGFR      <- eGFR[1:6,-1]
    HEAM      <- HEAM[1:6,-1]
    WBC       <- WBC[1:6,-1]
    HR        <- HR[1:6,-1]
    PVD       <- PVD[1:6,-1]
  }else if(mean(population_[,"MALE"])==0){
    #If FEMALE only get the rows associated with being Female
    A1c       <- A1c[7:12,-1]
    BMI       <- BMI[7:12,-1]
    SBP       <- SBP[7:12,-1]
    SMO       <- SMO[7:12,-1]
    HDL       <- HDL[7:12,-1]
    LDL       <- LDL[7:12,-1]
    ATF       <- ATF[7:12,-1]
    MMALB     <- MMALB[7:12,-1]
    eGFR      <- eGFR[7:12,-1]
    HEAM      <- HEAM[7:12,-1]
    WBC       <- WBC[7:12,-1]
    HR        <- HR[7:12,-1]
    PVD       <- PVD[7:12,-1]
  }else{
    stop("the population is neither entirely Male or Female for Mt Hood 2025 Challenge 1")
  } 
  
  
  ##select the right row of each matrix depending on the treatment arm
  rowlookup <- ifelse(treatment_=="Baseline",1,#Baseline use the 1st row
                      ifelse(treatment_=="Mt_HOOD_RS_A1c",2, #A1 effect only, use the 2nd row
                             ifelse(treatment_=="Mt_HOOD_RS_SBP",3, #SBP effect only use the 3rd row
                                    ifelse(treatment_=="Mt_HOOD_RS_LDL",4, #LDL effect only use the 4th row
                                           ifelse(treatment_=="Mt_HOOD_RS_BMI",5, #BMI effect only use the 5th row
                                                  6))))) #otherwise, use all effects
  
  #Get the right row for the intervention effects
  A1c       <- as.matrix(A1c[rowlookup,])
  BMI       <- as.matrix(BMI[rowlookup,])
  SBP       <- as.matrix(SBP[rowlookup,])
  SMO       <- as.matrix(SMO[rowlookup,])
  HDL       <- as.matrix(HDL[rowlookup,])
  LDL       <- as.matrix(LDL[rowlookup,])
  ATF       <- as.matrix(ATF[rowlookup,])
  MMALB     <- as.matrix(MMALB[rowlookup,])
  eGFR      <- as.matrix(eGFR[rowlookup,])
  HEAM      <- as.matrix(HEAM[rowlookup,])
  WBC       <- as.matrix(WBC[rowlookup,])
  HR        <- as.matrix(HR[rowlookup,])
  PVD       <- as.matrix(PVD[rowlookup,])
  
  #Expand these into matrices for each person in the simulation
  #Don't want the 1st data item in the vector, as it is description of the intervention, not data
  A1c <- matrix(A1c[1,], nrow = length(population_[,"ID"]), ncol = (length(A1c)), byrow = TRUE)
  BMI <- matrix(BMI[1,], nrow = length(population_[,"ID"]), ncol = (length(BMI)), byrow = TRUE)
  SBP <- matrix(SBP[1,], nrow = length(population_[,"ID"]), ncol = (length(SBP)), byrow = TRUE)
  SMO <- matrix(SMO[1,], nrow = length(population_[,"ID"]), ncol = (length(SMO)), byrow = TRUE)
  HDL <- matrix(HDL[1,], nrow = length(population_[,"ID"]), ncol = (length(HDL)), byrow = TRUE)
  LDL <- matrix(LDL[1,], nrow = length(population_[,"ID"]), ncol = (length(LDL)), byrow = TRUE)
  ATF <- matrix(ATF[1,], nrow = length(population_[,"ID"]), ncol = (length(ATF)), byrow = TRUE)
  MMALB <- matrix(MMALB[1,], nrow = length(population_[,"ID"]), ncol = (length(MMALB)), byrow = TRUE)
  eGFR <- matrix(eGFR[1,], nrow = length(population_[,"ID"]), ncol = (length(eGFR)), byrow = TRUE)
  HEAM <- matrix(HEAM[1,], nrow = length(population_[,"ID"]), ncol = (length(HEAM)), byrow = TRUE)
  WBC <- matrix(WBC[1,], nrow = length(population_[,"ID"]), ncol = (length(WBC)), byrow = TRUE)
  HR <- matrix(HR[1,], nrow = length(population_[,"ID"]), ncol = (length(HR)), byrow = TRUE)
  PVD <- matrix(PVD[1,], nrow = length(population_[,"ID"]), ncol = (length(PVD)-1), byrow = TRUE)
  
  #Dummy attend_SE matrix, where no one attends SE
  attend_se <-matrix(data = 0, nrow = length(population_[,"ID"]), ncol = 3)
  attend_se[,1] <- 1:length(population_[,"ID"])
  
  
  #start year at 0
  year <- 0
  
  #initialise the results matrix
  results <- GenerateResultsMatrix(GlobalVars_, endtime_)
  #initialise a PSA results matrix, if it is a PSA
  if(GlobalVars_["run_psa","Value"]==T){
    psaresults <- matrix(data=NA, nrow = as.numeric(GlobalVars_["psa_count", "Value"]), ncol = 5)  
    colnames(psaresults) <- c("Life Years per Patient", "Undiscounted QALYs per Patient", 
                              "QALYs per Patient", "Undiscounted Costs per Patient",
                              "Discounted Costs per Patient")
  }
  
  #run the model up to the specified end time or so long as at least one person is alive
  while (year < endtime_ & 
         (sum(is.na(population_[,"F_ALLCAUSE"]))) >= 1){
    #Get a logical vector indicating if people are alive, use the fact that FOTH
    #in the population matrix is NA if alive
    alive <- is.na(population_[,"F_ALLCAUSE"])
    #create matrix of ture = dead too for unit tests
    dead <- is.na(population_[,"F_ALLCAUSE"])==F
    
    #Estimate Diabetes Related complications and all cause deaths for this year
    population_ <- update_events_UKPDS82(population_,parameters_, treatment_, year, alive, random_numbs_, LifeTables_)
    #Estimate Depression
    population_ <- update_events_SPHR_depression(population_,parameters_,year,alive,random_numbs_)
    #Estimate Osteoarthritis
    population_ <- update_events_SPHR_osteoarthritis(population_,parameters_,year,alive,random_numbs_)
    #Estimate Cancer incidence
    population_ <- update_events_SPHR_cancer(population_, parameters_,year, alive,random_numbs_)
    
    #Stop the model if dead people have events
    if(sum(population_[,"MI_E"][dead])  !=0|
       sum(population_[,"MI2_E"][dead]) != 0|
       sum(population_[,"AMP_E"][dead]) != 0| 
       sum(population_[,"AMP2_E"][dead]) != 0|
       sum(population_[,"BLIND_E"][dead]) != 0|
       sum(population_[,"ULCER_E"][dead]) != 0|
       sum(population_[,"RENAL_E"][dead]) != 0|
       sum(population_[,"CHF_E"][dead]) != 0|
       sum(population_[,"IHD_E"][dead]) != 0|
       sum(population_[,"STRO_E"][dead]) != 0|
       sum(population_[,"STRO2_E"][dead]) != 0|
       sum(population_[,"ATFIB_E"][dead]) != 0|
       sum(population_[,"PVD_E"][dead]) != 0|
       sum(population_[,"MMALB_E"][dead]) != 0|
       sum(population_[,"CANB_E"][dead]) != 0|
       sum(population_[,"CANC_E"][dead]) != 0|
       sum(population_[,"DEP_E"][dead]) != 0|
       sum(population_[,"OST_E"][dead]) != 0){
      #stop the model if the number of people with dead people are recorded as
      #having an event
      #push everything to the global enviroment for bug checking
      for (variable in ls()) {
        assign(variable, get(variable), envir = .GlobalEnv)
      }
      
      stop("dead people are getting events")
    }
    
    ##QALYs
    population_ <- calculate_QALYs(population_, parameters_,  year, alive, GlobalVars_)
    ##Costs
    population_ <- calculate_costs(population_, parameters_, year, alive, GlobalVars_,treatment_,attend_se)
    
    #Record results
    results <- GenerateDetailedresults(results,population_, year, alive, GlobalVars_)
    
    
    #update histories
    #update histories
    population_ <- update_history_MtHood2025C1(population_,
                                               alive,
                                               A1c,
                                               BMI,
                                               SBP,
                                               SMO,       
                                               HDL,
                                               LDL,       
                                               ATF,
                                               MMALB,
                                               HEAM,
                                               WBC,
                                               HR,
                                               PVD,
                                               eGFR,
                                               treatment_,
                                               year,
                                               GlobalVars_)
    
    population_ <- update_patchars(population_, parameters_, alive)
    
    #Unit test
    #Stop the model if there are events still in the population matrix
    if(sum(population_[,"MI_E"])  !=0|
       sum(population_[,"MI2_E"]) != 0|
       sum(population_[,"AMP_E"]) != 0| 
       sum(population_[,"AMP2_E"]) != 0|
       sum(population_[,"BLIND_E"]) != 0|
       sum(population_[,"ULCER_E"]) != 0|
       sum(population_[,"RENAL_E"]) != 0|
       sum(population_[,"CHF_E"]) != 0|
       sum(population_[,"IHD_E"]) != 0|
       sum(population_[,"STRO_E"]) != 0|
       sum(population_[,"STRO2_E"]) != 0){
      
      #push everything to the global enviroment for bug checking
      for (variable in ls()) {
        assign(variable, get(variable), envir = .GlobalEnv)
      }
      #stop the model if the events are not reset to 0
      stop("events are not set to zero before reset")
    }
    
    year <- year + 1
    
  }
  
  
  #For now return the Detailed results table or the population matrix if the run is deterministic
  if(GlobalVars_["Results_output", "Value"] == "Summary"&
     GlobalVars_["run_psa", "Value"]==T){
    psaresults <- matrix(data=NA,nrow=1,ncol=24)
    #Life Years
    psaresults[,1] <- sum(results["Undiscounted life years accrued",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,2] <- sum(results["Discounted life years accrued",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,3] <- sum(results["Undiscounted QALYs",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,4] <- sum(results["Discounted QALYs",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,5] <- sum(results["Undiscounted Costs",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,6] <- sum(results["Discounted Costs",],na.rm=TRUE)/length(population_[,"ID"])
    #10 year histories
    psaresults[,7] <- results["1st MI Hist",11]
    psaresults[,8] <- results["2nd MI Hist",11]
    psaresults[,9] <- results["1st Stroke Hist",11]
    psaresults[,10] <- results["2nd Stroke Hist",11]
    psaresults[,11] <- results["CHF Hist",11]
    psaresults[,12] <- results["IHD Hist",11]
    psaresults[,13] <- results["Blindness Hist",11]
    psaresults[,14] <- results["Ulcer Hist",11]
    psaresults[,15] <- results["1st Amputation Hist",11]
    psaresults[,16] <- results["2nd Amputation Hist",11]
    psaresults[,17] <- results["Renal Failure Hist",11]
    psaresults[,18] <- results["PVD Hist",11]
    psaresults[,19] <- results["MMALB Hist",11]
    psaresults[,20] <- results["ATFIB Hist",11]
    psaresults[,21] <- results["Breast Cancer Hist",11]
    psaresults[,22] <- results["Colorectal Cancer Hist",11]
    psaresults[,23] <- results["Depression Hist",11]
    psaresults[,24] <- results["Osteoarthritis Hist",11]
    #delete the original results matrix
    rm(results)
    return(psaresults)
  }else if(GlobalVars_["Results_output","Value"]=="Patient Level"&
           GlobalVars_["run_psa","Value"]==F){#option to produce the patient 
    #charateristics matrix for checking results stability for number of patients
    return(population_)
  }else{#default is the 
    return(results)
  }
}



####'This function runs the model for a given set of patients and parameters
####'it is designed for the Mt Hood Diabetes Challenge 2025, Challenge 2 as this 
####'function uses predefined life paths for each simulated inidividual's trajectories
####'this function differs from the function for challenge 1, as each of the 1,000 individuals for this
####'challenge has a different life course
####'@param population_, is the population matrix
####'@param parameters_, is the full parameters matrix
####'@param endtime_, is the number of years to run the simulation for
####'@param GlobalVars_, is the matrix giving the global variables
####'@param LifeTables_, is a dataframe containing life table information in formate
####'that can easily be matched to the population matrix
####'@param SOUR_, is the current second order uncertainty run
####'draw for each patient in each year for every event in which a random number 
####'is required 
####'@param MtHood2025C2_ is a list of baseline characteristics and risk factor trajectories
####'for the Mt Hood 2025 challenge 2
####'@return results, is the results matrix
####'@return psaresults, is a summary set of results to produce in the PSA

run_simulation_MtHood2025_C2 <- function(population_, parameters_, endtime_, treatment_, GlobalVars_,LifeTables_, SOUR_,MtHood2025C2_,bootit_){
  
  ##reduce the parameters matrix down to the correct row
  if(GlobalVars_["run_psa","Value"]==F){
    parameters_ <- parameters_[1,]
  }else{
    parameters_ <- parameters_[SOUR_+1,]
  }
  
  attend_se <-matrix(data = 0, nrow = length(population_[,"ID"]), ncol = 3)
  attend_se[,1] <- 1:length(population_[,"ID"])
  
  ##Get the right data frames for risk factor / event progression by treatment
  #HbA1c, BMI, SBP
  if(treatment_=="BC_Control"|treatment_=="BC_INTV"){ # If running intervention or control for the basecase, load the basecase A1c trajectory
    A1c <- MtHood2025C2_$A1cBC.csv
    BMI <- MtHood2025C2_$BMIBC.csv
    SBP <- MtHood2025C2_$SBPBC.csv
  }else if (treatment_=="S1_Control"|treatment_=="S1_INTV"){ # If running intervention or control for the Scenario 1, load the Scenario 1 A1c trajectory
    A1c <- MtHood2025C2_$A1cS1.csv
    BMI <- MtHood2025C2_$BMIS1.csv
    SBP <- MtHood2025C2_$SBPS1.csv
    
  }else if (treatment_=="S2_Control"|treatment_=="S2_INTV"){ # If running intervention or control for the Scenario 2, load the Scenario 2 A1c trajectory
    A1c <- MtHood2025C2_$A1cS2.csv
    BMI <- MtHood2025C2_$BMIS2.csv
    SBP <- MtHood2025C2_$SBPS2.csv
    
  }else {
    A1c <- MtHood2025C2_$A1cNoEffect.csv #Otherwise load the no effect A1c trajectories
    BMI <- MtHood2025C2_$BMINoEffect.csv
    SBP <- MtHood2025C2_$SBPNoEffect.csv
    
    
  }
  #Extract the other risk factors that do not depend on treatment arm for progression
  SMO         <- MtHood2025C2_$Smoking.csv
  HDL         <- MtHood2025C2_$HDL.csv
  LDL         <- MtHood2025C2_$LDL.csv
  PVD         <- MtHood2025C2_$PVD.csv
  ATF         <- MtHood2025C2_$AF.csv
  MMALB       <- MtHood2025C2_$MMALB.csv
  eGFR        <- MtHood2025C2_$eGFR.csv
  HEAM        <- MtHood2025C2_$HEAM.csv
  WBC         <- MtHood2025C2_$WBC.csv
  HR          <- MtHood2025C2_$HR.csv
  PVD         <- MtHood2025C2_$PVD.csv
  
  #If control we only want group 1 people, if intervention, we want group 2 and to change their ID to start at 1 again
  if(treatment_=="BC_Control"|treatment_=="S1_Control"|treatment_=="S2_Control"|treatment_=="NoEffect_Control"){
    A1c       <- subset(A1c,Group==1)
    BMI       <- subset(BMI,Group==1)
    SBP       <- subset(SBP,Group==1)
    SMO       <- subset(SMO,Group==1)
    HDL       <- subset(HDL,Group==1)
    LDL       <- subset(LDL,Group==1)
    ATF       <- subset(ATF,Group==1)
    MMALB     <- subset(MMALB,Group==1)
    HEAM      <- subset(HEAM,Group==1)
    WBC       <- subset(WBC,Group==1)
    HR        <- subset(HR,Group==1)
    PVD       <- subset(PVD,Group==1)
    eGFR      <- subset(eGFR, Group==1)
  }else{
    A1c       <- subset(A1c,Group==2)
    BMI       <- subset(BMI,Group==2)
    SBP       <- subset(SBP,Group==2)
    SMO       <- subset(SMO,Group==2)
    HDL       <- subset(HDL,Group==2)
    LDL       <- subset(LDL,Group==2)
    ATF       <- subset(ATF,Group==2)
    MMALB     <- subset(MMALB,Group==2)
    HEAM      <- subset(HEAM,Group==2)
    WBC       <- subset(WBC,Group==2)
    HR        <- subset(HR,Group==2)
    PVD       <- subset(PVD,Group==2)
    eGFR      <- subset(eGFR, Group==2)
    
    
    #Change ID's back to 1 to 1000
    A1c$ID    <- 1:1000
    BMI$ID    <- 1:1000
    SBP$ID    <- 1:1000
    SMO$ID    <- 1:1000
    HDL$ID    <- 1:1000
    LDL$ID    <- 1:1000
    ATF$ID    <- 1:1000
    MMALB$ID  <- 1:1000
    HEAM$ID   <- 1:1000
    WBC$ID    <- 1:1000
    HR$ID     <- 1:1000
    PVD$ID    <- 1:1000
    eGFR$ID   <- 1:1000
  }
  
  
  #start year at 0
  year <- 0
  
  #As this function need random numbers on the fly generate them here
  random_numbs <- generate_random(length(population_[,"ID"]))
  
  
  #initialise the results matrix
  results <- GenerateResultsMatrix(GlobalVars_, endtime_)
  #initialise a PSA results matrix, if it is a PSA
  if(GlobalVars_["run_psa","Value"]==T){
    psaresults <- matrix(data=NA, nrow = as.numeric(GlobalVars_["psa_count", "Value"]), ncol = 5)  
    colnames(psaresults) <- c("Life Years per Patient", "Undiscounted QALYs per Patient", 
                              "QALYs per Patient", "Undiscounted Costs per Patient",
                              "Discounted Costs per Patient")
  }
  
  #run the model up to the specified end time or so long as at least one person is alive
  while (year < endtime_ & 
         (sum(is.na(population_[,"F_ALLCAUSE"]))) >= 1){
    #Get a logical vector indicating if people are alive, use the fact that FOTH
    #in the population matrix is NA if alive
    alive <- is.na(population_[,"F_ALLCAUSE"])
    #create matrix of ture = dead too for unit tests
    dead <- is.na(population_[,"F_ALLCAUSE"])==F
    
    #Estimate Diabetes Related complications and all cause deaths for this year
    population_ <- update_events_UKPDS82(population_,parameters_, treatment_, year, alive, random_numbs, LifeTables_)
    #Estimate Depression
    population_ <- update_events_SPHR_depression(population_,parameters_,year,alive,random_numbs)
    #Estimate Osteoarthritis
    population_ <- update_events_SPHR_osteoarthritis(population_,parameters_,year,alive,random_numbs)
    #Estimate Cancer incidence
    population_ <- update_events_SPHR_cancer(population_, parameters_,year, alive,random_numbs)
    
    #Stop the model if dead people have events
    if(sum(population_[,"MI_E"][dead])  !=0|
       sum(population_[,"MI2_E"][dead]) != 0|
       sum(population_[,"AMP_E"][dead]) != 0| 
       sum(population_[,"AMP2_E"][dead]) != 0|
       sum(population_[,"BLIND_E"][dead]) != 0|
       sum(population_[,"ULCER_E"][dead]) != 0|
       sum(population_[,"RENAL_E"][dead]) != 0|
       sum(population_[,"CHF_E"][dead]) != 0|
       sum(population_[,"IHD_E"][dead]) != 0|
       sum(population_[,"STRO_E"][dead]) != 0|
       sum(population_[,"STRO2_E"][dead]) != 0|
       sum(population_[,"ATFIB_E"][dead]) != 0|
       sum(population_[,"PVD_E"][dead]) != 0|
       sum(population_[,"MMALB_E"][dead]) != 0|
       sum(population_[,"CANB_E"][dead]) != 0|
       sum(population_[,"CANC_E"][dead]) != 0|
       sum(population_[,"DEP_E"][dead]) != 0|
       sum(population_[,"OST_E"][dead]) != 0){
      #stop the model if the number of people with dead people are recorded as
      #having an event
      #push everything to the global enviroment for bug checking
      for (variable in ls()) {
        assign(variable, get(variable), envir = .GlobalEnv)
      }
      
      stop("dead people are getting events")
    }
    
    ##QALYs
    population_ <- calculate_QALYs_MtHood2025_C2(population_,  year, alive, GlobalVars_)
    ##Costs
    population_ <- calculate_costs_MtHood2025_C2(population_, year, alive, GlobalVars_)
    
    #Record results
    results <- GenerateDetailedresults(results,population_, year, alive, GlobalVars_)
    
    
    #update histories
    population_ <- update_history_MtHood2025C2(population_,
                                               alive,
                                               A1c,
                                               BMI,
                                               SBP,
                                               SMO,       
                                               HDL,
                                               LDL,       
                                               ATF,
                                               MMALB,
                                               HEAM,
                                               WBC,
                                               HR,
                                               PVD,
                                               eGFR,
                                               treatment_,
                                               year,
                                               GlobalVars_)
    
    population_ <- update_patchars(population_, parameters_, alive)
    
    #Unit test
    #Stop the model if there are events still in the population matrix
    if(sum(population_[,"MI_E"])  !=0|
       sum(population_[,"MI2_E"]) != 0|
       sum(population_[,"AMP_E"]) != 0| 
       sum(population_[,"AMP2_E"]) != 0|
       sum(population_[,"BLIND_E"]) != 0|
       sum(population_[,"ULCER_E"]) != 0|
       sum(population_[,"RENAL_E"]) != 0|
       sum(population_[,"CHF_E"]) != 0|
       sum(population_[,"IHD_E"]) != 0|
       sum(population_[,"STRO_E"]) != 0|
       sum(population_[,"STRO2_E"]) != 0){
      
      #push everything to the global enviroment for bug checking
      for (variable in ls()) {
        assign(variable, get(variable), envir = .GlobalEnv)
      }
      #stop the model if the events are not reset to 0
      stop("events are not set to zero before reset")
    }
    
    #In years 3 and 4, record Life Years, QALYs and Costs for the challenge
    #Also, collect these data for everyone, not just the dead people
    if(year == 2){
      population_[,"LY_Y3"]   <- population_[,"YearsLived"]
      population_[,"QALY_Y3"] <- population_[,"QALY"]
      population_[,"Cost_Y3"] <- population_[,"COST"]
    }else if (year == 3){
      population_[,"LY_Y4"]   <- population_[,"YearsLived"]
      population_[,"QALY_Y4"] <- population_[,"QALY"]
      population_[,"Cost_Y4"] <- population_[,"COST"]
    }
    
    year <- year + 1
    
  }
  
  
  #For now return the Detailed results table or the population matrix if the run is deterministic
  if(GlobalVars_["Results_output", "Value"] == "Summary"){
    psaresults <- matrix(data=NA,nrow=1,ncol=24)
    #Life Years
    psaresults[,1] <- sum(results["Undiscounted life years accrued",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,2] <- sum(results["Discounted life years accrued",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,3] <- sum(results["Undiscounted QALYs",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,4] <- sum(results["Discounted QALYs",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,5] <- sum(results["Undiscounted Costs",],na.rm=TRUE)/length(population_[,"ID"])
    psaresults[,6] <- sum(results["Discounted Costs",],na.rm=TRUE)/length(population_[,"ID"])
    #10 year histories
    psaresults[,7] <- results["1st MI Hist",11]
    psaresults[,8] <- results["2nd MI Hist",11]
    psaresults[,9] <- results["1st Stroke Hist",11]
    psaresults[,10] <- results["2nd Stroke Hist",11]
    psaresults[,11] <- results["CHF Hist",11]
    psaresults[,12] <- results["IHD Hist",11]
    psaresults[,13] <- results["Blindness Hist",11]
    psaresults[,14] <- results["Ulcer Hist",11]
    psaresults[,15] <- results["1st Amputation Hist",11]
    psaresults[,16] <- results["2nd Amputation Hist",11]
    psaresults[,17] <- results["Renal Failure Hist",11]
    psaresults[,18] <- results["PVD Hist",11]
    psaresults[,19] <- results["MMALB Hist",11]
    psaresults[,20] <- results["ATFIB Hist",11]
    psaresults[,21] <- results["Breast Cancer Hist",11]
    psaresults[,22] <- results["Colorectal Cancer Hist",11]
    psaresults[,23] <- results["Depression Hist",11]
    psaresults[,24] <- results["Osteoarthritis Hist",11]
    #delete the original results matrix
    rm(results)
    return(psaresults)
  }else if(GlobalVars_["Results_output","Value"]=="Patient Level"&
           GlobalVars_["run_psa","Value"]==F){#option to produce the patient 
    #charateristics matrix for checking results stability for number of patients
    return(population_)
  }else if (GlobalVars_["Results_output", "Value"] == "MtHood2025C2"){
    MtHoodresults <- matrix(data = NA, nrow = length(population_[,"ID"]), ncol = 6)
    names <- c("3 year life expectancy",
               "3 year Quality Adjusted Life Years",
               "3 year costs",
               "4 year life expectancy",
               "4 year Quality Adjusted Life Years",
               "4 year costs")
    colnames(MtHoodresults) <- names
    rm(names)
    MtHoodresults[,"3 year life expectancy"] <- population_[,"LY_Y3"]
    MtHoodresults[,"3 year Quality Adjusted Life Years"] <- population_[,"QALY_Y3"]
    MtHoodresults[,"3 year costs"] <- population_[,"Cost_Y3"]
    MtHoodresults[,"4 year life expectancy"] <- population_[,"LY_Y4"]
    MtHoodresults[,"4 year Quality Adjusted Life Years"] <- population_[,"QALY_Y4"]
    MtHoodresults[,"4 year costs"] <- population_[,"Cost_Y4"]
    
    return(MtHoodresults)
  }
  else{#default is the 
    return(results)
  }
}



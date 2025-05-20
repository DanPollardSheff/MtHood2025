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



#' Generate the baseline population
#' @param diab_diab_population_ is a matrix of baseline characteristics read into 
#' this function
#' @param PopulationVariables_ is a list of all the variables in the population of this 
#' model
#' @param GlobalVars_ is the global variables matrix
#' @return population is a matrix of patient characteristics for use in the model
build_population <- function(diag_diab_population_, PopulationVariables_, GlobalVars_) {
  #limit the raw population characteristics to be the same as the number se
  #in the global varaibles
  diag_diab_population_ <- diag_diab_population_[1:as.numeric(GlobalVars_["n","Value"]),]
  #stop the simulation if there are too few patients in the csv files
  if(length(diag_diab_population_[,"CURR_AGE"]) > as.numeric(GlobalVars_["n","Value"])){
    stop("there are too few patients in the csv files, some patients will get 0 for
         their characteristics in the simulation")
  }
  
  
  n <- nrow(diag_diab_population_)  
  year <- 0
  population <- matrix(0, nrow = n, ncol = length(PopulationVariables_[,"Variable"]))
  colnames(population) <- PopulationVariables_[,"Variable"]
  population[, "ID"] <- 1:n
  population[, "CONS"] <- 1
  population[, "T"] <- 0
  population[, "AGE_0"] <- floor(diag_diab_population_[, "CURR_AGE"])
  population[, "MALE"] <- (-1*diag_diab_population_[, "Female"])+1
  population[, "FEMALE"] <- diag_diab_population_[, "Female"]
  
  #Read in a dummy variable for Afro-caribean descent (1= afro-caribean, 0=otherwise)
  population[, "AFRO"] <- diag_diab_population_[,"AFRO"]
  population[,"INDIAN"] <- diag_diab_population_[,"INDIAN"]
  #Set the binary variable for smoking status
  population[, "SMO"] <- ifelse(diag_diab_population_[, "Smoking"]>0,1,0)
  #record the history of diabetes (1 = yes, 0 = no)
  
  population[, "AGE"] <- floor(diag_diab_population_[, "CURR_AGE"])
  population[, "MEN"] <- replace(population[, "MEN"], population[, "AGE"] > 51 & population[, "FEMALE"] == 1, 1)
  population[, "HBA"] <- round(diag_diab_population_[, "HbA1c"],1)
  #This has been kept the same as in the SPHR diabetes model
  #Record the BMI of the population
  population[, "BMI"] <- diag_diab_population_[, "BMI"]
  #Apply a logical constraint to the BMI so that it cannot be less than 5kg/m2
  population[, "BMI"] <- replace(diag_diab_population_[, "BMI"], population[, "BMI"] < 5, 5)
  #Record the HDL cholesterol (add units)
  population[, "HDL"] <- diag_diab_population_[, "HDL"]
  #Record systolic blood pressure(add units)
  population[, "SBP"] <- diag_diab_population_[, "SBP"]
  
  
  #Use the Ara formula to estimate baseline QALYs
  #as there is a parameter determining baseline utility, set these two values to 0 for the time being
  population[, "QALY"] <- 0 
  population[, "EQ5D"] <- 0
  
  population[, "DEP_H"] <- runif(n) < (435/4781) #Source: Ali et al 2009. Prevalence of diagnosed depression in South Asian
  #and white European people with type 1 and type 2 diabetes mellitus in a UK secondary care population
  
  #Record the baseline BMI of the population in the HSE data
  population[,"BMI_0"] <- diag_diab_population_[,"Baseline_BMI"]
  
  #Record whether someone is on 1st line therapy for T2DM
  population[, "MET"]<-diag_diab_population_[,"MET"]
  #Record whether someone is on 2nd line therapy for T2DM
  population[,"MET2"]<- diag_diab_population_[,"MET2"]
  #Record whether someone is on 3rd line therapy for T2DM
  population[,"INSU"]<- diag_diab_population_[,"INSU"]
  
  ##parameters for the treatment version of the SPHR diabetes model 
  population[,"AMP_E"] <- diag_diab_population_[,"AMP_E"]
  population[,"AMP_H"] <- diag_diab_population_[,"AMP_H"]
  population[,"AMP2_H"] <- diag_diab_population_[,"AMP2_H"]
  population[,"IHD_E"] <- diag_diab_population_[,"IHD_E"]
  population[,"IHD_H"] <- diag_diab_population_[,"IHD_H"]
  population[,"CHF_E"] <- diag_diab_population_[,"CHF_E"]
  population[,"CHF_H"] <- diag_diab_population_[,"CHF_H"]
  population[,"RENAL_E"] <- diag_diab_population_[,"RENAL_E"]
  population[,"RENAL_H"] <- diag_diab_population_[,"RENAL_H"]
  population[,"STRO_E"] <- diag_diab_population_[,"STROKE_E"]
  population[,"STRO_H"] <- diag_diab_population_[,"STROKE_H"]
  population[,"MI_E"] <- diag_diab_population_[,"MI_E"]
  population[,"MI_H"] <- diag_diab_population_[,"MI_H"]
  population[,"MI2_E"] <- diag_diab_population_[,"MI2_E"]
  population[,"MI2_H"] <- diag_diab_population_[,"MI2_H"]
  population[,"BLIND_E"] <- diag_diab_population_[,"BLIND_E"]
  population[,"BLIND_H"] <- diag_diab_population_[,"BLIND_H"]
  population[,"MMALB_H"] <- diag_diab_population_[,"MIC_ALB"]
  population[,"ATFIB_H"] <- diag_diab_population_[,"ATFIB"]
  population[,"PVD_H"] <- diag_diab_population_[,"PVD"]
  
  #Add in diabetes treatment specific biomarkers
  population[,"HAEM"] <- diag_diab_population_[,"HAEM"]
  population[,"WBC"] <- diag_diab_population_[,"WBC"]
  population[,"eGFR"] <- diag_diab_population_[,"eGFR"]
  population[,"eGFR_U_60"] <- ifelse(diag_diab_population_[,"eGFR"]<60,diag_diab_population_[,"eGFR"],60)
  population[,"eGFR_O_60"] <- ifelse(diag_diab_population_[,"eGFR"]>60,diag_diab_population_[,"eGFR"]-60,0)
  population[,"HEART_R"] <- diag_diab_population_[,"Heart.rate"]
  population[,"BMI_U_18_5"] <- ifelse(population[,"BMI"]<18.5,1,0)
  population[,"BMI_O_E_25"] <- ifelse(population[,"BMI"]>=25,1,0)
  population[,"LDL"] <- diag_diab_population_[,"LDL"]
  population[,"LDL_O_35"] <- ifelse(population[,"LDL"]>3.5, population[,"LDL"]-3.5,0)
  #add in diabetes treatment specific chars
  population[,"DIAB_DUR"] <- diag_diab_population_[,"DIAB_DUR"]
  
  #add in values for first observations for each risk factor needed
  #in the absence of information this will be the baseline value
  population[,"HBA_0"] <- diag_diab_population_[,"Baseline_HbA1c"]
  population[,"LDL_0"] <- diag_diab_population_[,"Baseline_LDL"]
  population[,"HDL_0"] <- diag_diab_population_[,"Baseline_HDL"]
  population[,"HEART_R_0"] <- population[,"HEART_R"]
  population[,"HAEM_0"] <- population[,"HAEM"]
  population[,"WBC_0"] <- population[,"WBC"]
  population[,"eGFR_0"] <- population[,"eGFR"]
  population[,"SMO_0"] <- population[,"SMO"]
  population[,"SBP_0"] <- diag_diab_population_[,"Baseline_SBP"]
  
  #Give noone a history of Ulcers
  population[,"ULCER_H"] <-  0
  
  #Make all cause death missing, as missing indicates someone is alive in the code
  population[, "F_ALLCAUSE"] <- NA
  #Make smoking missing, so you can easily see whether smoking is assessed
  population[,"p_SMO"]<- NA
  
  return(population)
}

#' Generate the baseline population for MtHood 2025 challenge instructions
#' On the 16th May 2025, the instructions are at  https://www.mthooddiabeteschallenge.com/_files/ugd/4e5824_14eed2de8aa94f48b894eef0f5460d83.pdf
#' @param n_ is the number of replications
#' @param Female_ is a logical value. TRUE = the whole population is Female, FALSE = the whole population is male
#' @param PopulationVariables_ is a list of all the population varaibles used in this model along with a description of units and meanings
#' @return population is a matrix of patient characteristics for use in the model
#' 
build_population_MtHood2025 <- function(n_, Female_, PopulationVariables_) {
  
  #limit the raw population characteristics to be the same as the number se
  #in the global varaibles
    
  year <- 0
  population <- matrix(0, nrow = n_, ncol = length(PopulationVariables_[,"Variable"]))
  colnames(population) <- PopulationVariables_[,"Variable"]
  population[, "ID"] <- 1:n_
  population[, "CONS"] <- 1
  population[, "T"] <- 0
  population[, "AGE_0"] <- 66 #Mt Hood 2025 Reference Simulation instructions
  population[, "MALE"] <- ifelse(Female_==T,0,1)
  population[, "FEMALE"] <- ifelse(Female_==T,1,0)
  
  #Read in a dummy variable for Afro-caribean descent (1= afro-caribean, 0=otherwise)
  population[, "AFRO"] <- 0  #Mt Hood 2025 Reference Simulation instructions, person has white ehtnicity.
  population[,"INDIAN"] <- 0 #Mt Hood 2025 Reference Simulation instructions, person has white ehtnicity.
  #Set the binary variable for smoking status
  population[, "SMO"] <- 0   #Mt Hood 2025 Reference Simulation instructions, person is a non-smoker.
  #record the history of diabetes (1 = yes, 0 = no)
  
  population[, "AGE"] <- 66
  population[, "MEN"] <- replace(population[, "MEN"], population[, "AGE"] > 51 & population[, "FEMALE"] == 1, 1)
  population[, "HBA"] <- 7.5 #Mt Hood 2025 Reference Simulation instructions.
  
  #Record the BMI of the population
  population[, "BMI"] <- 28
  #Record the HDL cholesterol (add units)
  population[, "HDL"] <- 1.3
  #Record systolic blood pressure(add units)
  population[, "SBP"] <- 145
  
  population[, "QALY"] <- 0 
  population[, "EQ5D"] <- 0
  
  population[, "DEP_H"] <- runif(n_) < (435/4781) #Source: Ali et al 2009. Prevalence of diagnosed depression in South Asian
  #and white European people with type 1 and type 2 diabetes mellitus in a UK secondary care population
  
  #Record the baseline BMI of the population in the HSE data
  population[,"BMI_0"] <- 28
  
  #Instructions do not specify, so set everyone to be on 1st line therapy
  population[, "MET"] <- 1
  #Record whether someone is on 2nd line therapy for T2DM
  population[,"MET2"] <- 0
  #Record whether someone is on 3rd line therapy for T2DM
  population[,"INSU"]<- 0
  
  #Mt Hood 2025 Reference Simulation instructions, no history at baseline. 
  population[,"AMP_E"] <- 0
  population[,"AMP_H"] <- 0
  population[,"AMP2_H"] <- 0
  population[,"IHD_E"] <- 0
  population[,"IHD_H"] <- 0
  population[,"CHF_E"] <- 0
  population[,"CHF_H"] <- 0
  population[,"RENAL_E"] <- 0
  population[,"RENAL_H"] <- 0
  population[,"STRO_E"] <- 0
  population[,"STRO_H"] <- 0
  population[,"MI_E"] <- 0
  population[,"MI_H"] <- 0
  population[,"MI2_E"] <- 0
  population[,"MI2_H"] <- 0
  population[,"BLIND_E"] <- 0
  population[,"BLIND_H"] <- 0
  population[,"MMALB_H"] <- 0
  population[,"ATFIB_H"] <- 0
  population[,"PVD_H"] <- 0
  
  #Add in diabetes treatment specific biomarkers
  population[,"HAEM"] <- 14 
  population[,"WBC"] <- 7
  population[,"eGFR"] <- 70
  population[,"eGFR_U_60"] <- ifelse(population[,"eGFR"]<60,population[,"eGFR"],60)
  population[,"eGFR_O_60"] <- ifelse(population[,"eGFR"]>60,population[,"eGFR"]-60,0)
  population[,"HEART_R"] <- 79
  population[,"BMI_U_18_5"] <- ifelse(population[,"BMI"]<18.5,1,0)
  population[,"BMI_O_E_25"] <- ifelse(population[,"BMI"]>=25,1,0)
  population[,"LDL"] <- 3.0
  population[,"LDL_O_35"] <- ifelse(population[,"LDL"]>3.5, population[,"LDL"]-3.5,0)
  #add in diabetes treatment specific chars
  population[,"DIAB_DUR"] <- 8
  
  #add in values for first observations for each risk factor needed
  #in the absence of information this will be the baseline value
  population[,"HBA_0"] <- population[,"HBA"]
  population[,"LDL_0"] <- population[,"LDL"]
  population[,"HDL_0"] <- population[, "HDL"]
  population[,"HEART_R_0"] <- population[,"HEART_R"]
  population[,"HAEM_0"] <- population[,"HAEM"]
  population[,"WBC_0"] <- population[,"WBC"]
  population[,"eGFR_0"] <- population[,"eGFR"]
  population[,"SMO_0"] <- population[,"SMO"]
  population[,"SBP_0"] <- population[,"SBP"]
  
  #Give noone a history of Ulcers
  population[,"ULCER_H"] <-  0
  
  #Make all cause death missing, as missing indicates someone is alive in the code
  population[, "F_ALLCAUSE"] <- NA
  #Make smoking missing, so you can easily see whether smoking is assessed
  population[,"p_SMO"]<- NA
  
  return(population)
}

#' @param data_ is a list containing all baseline characteristcs and risk factor time paths
#' for Mt Hoot Diabetes Challenge, Challenge 2
#' @param treatment_ is a text varaible determining the treatment - this shouldn't 
#' be necessary as the bl characteristics should be the same
#' @param PopulationVariables_ is a list of all varaibles required in the microsimulation
#' and their meaning (although the meaning is not used in the code)
#' @return Population is the baseline characteristics of the population to use in the microsimulation

build_population_MtHood2025_C2 <- function(data_, treatment_, PopulationVariables_) {
  
  #extract only the baseline characteristics
  bldata <- data_[["BaselineCharacteristics.csv"]]
  #subset by Group
  if (treatment_=="Control"){
    bldata <- subset(bldata, Group==1)
  }else{
    bldata <- subset(bldata, Group==2)
  }
  #Remove Group, as it is no longer needed
  bldata <- bldata[,!(names(bldata)%in%"Group")]
    
  #limit the raw population characteristics to be the same as the number se
  #in the global varaibles
  
  year <- 0
  population <- matrix(0, nrow = length(bldata[,"ID"]), ncol = length(PopulationVariables_[,"Variable"]))
  colnames(population) <- PopulationVariables_[,"Variable"]
  population[, "ID"] <- 1:length(bldata[,"ID"])
  population[, "CONS"] <- 1
  population[, "T"] <- 0
  population[, "AGE_0"] <- floor(bldata$Age.now)
  population[, "MALE"] <- ifelse(bldata$Gender=="M",0,1)
  population[, "FEMALE"] <- ifelse(bldata$Gender=="F",1,0)
  
  #Read in a dummy variable for Afro-caribean descent (1= afro-caribean, 0=otherwise)
  population[, "AFRO"] <- ifelse(bldata$Ethnicity==2,1,0)  
  population[,"INDIAN"] <- ifelse(bldata$Ethnicity==3,1,0)   
  #Set the binary variable for smoking status
  population[, "SMO"] <- ifelse(bldata$Current.smoker=="Y",1,0)  
  #record the history of diabetes (1 = yes, 0 = no)
  
  population[, "AGE"] <- population[, "AGE_0"]
  population[, "MEN"] <- replace(population[, "MEN"], population[, "AGE"] > 51 & population[, "FEMALE"] == 1, 1)
  population[, "HBA"] <- bldata$HbA1c
  
  #Record the BMI of the population
  #No BMI, so do weight/height^2. Apply ^2 as a multiple of itself, as this is faster in R
  population[, "BMI"] <- bldata$Weight/(bldata$Height*bldata$Height)
  #Record the HDL cholesterol (add units)
  population[, "HDL"] <- bldata$HDL
  #Record systolic blood pressure(add units)
  population[, "SBP"] <- bldata$Systolic.BP
  
  population[, "QALY"] <- 0 
  population[, "EQ5D"] <- 0
  
  population[, "DEP_H"] <- runif(n_) < (435/4781) #Source: Ali et al 2009. Prevalence of diagnosed depression in South Asian
  #and white European people with type 1 and type 2 diabetes mellitus in a UK secondary care population
  
  #Record the baseline BMI of the population in the HSE data
  population[,"BMI_0"] <- population[, "BMI"]
  
  #Instructions do not specify, so set everyone to be on 1st line therapy
  population[, "MET"] <- 1
  #Record whether someone is on 2nd line therapy for T2DM
  population[,"MET2"] <- 0
  #Record whether someone is on 3rd line therapy for T2DM
  population[,"INSU"]<- 0
  
  #Mt Hood 2025 Reference Simulation instructions, no history at baseline. 
  #If the vlaues are missing there is no event of history of an event
  #If the value is not missing, check if the value is 1. If it is it is an event year, 
  #otherwise it is a historical event
  population[,"AMP_E"] <- ifelse(is.na(bldata$Amputation), 0,
                                 ifelse(bldata$Amputation==1,1,0))
  population[,"AMP_H"] <- ifelse(is.na(bldata$Amputation), 0,
                                 ifelse(bldata$Amputation==1,0,1))
  #No data on 2nd amputation, so set to 0
  population[,"AMP2_E"] <- 0
  population[,"AMP2_H"] <- 0
  
  population[,"IHD_E"] <- ifelse(is.na(bldata$IHD), 0,
                                 ifelse(bldata$IHD==1,1,0))
  population[,"IHD_H"] <- ifelse(is.na(bldata$IHD), 0,
                                 ifelse(bldata$IHD==1,0,1))
  population[,"CHF_E"] <- ifelse(is.na(bldata$Heart.failure), 0,
                                 ifelse(bldata$Heart.failure==1,1,0))
  population[,"CHF_H"] <- ifelse(is.na(bldata$Heart.failure), 0,
                                 ifelse(bldata$Heart.failure==1,0,1))
  population[,"RENAL_E"] <- ifelse(is.na(bldata$Renal.failure), 0,
                                   ifelse(bldata$Renal.failure==1,1,0))
  population[,"RENAL_H"] <- ifelse(is.na(bldata$Renal.failure), 0,
                                   ifelse(bldata$Renal.failure==1,0,1))
  population[,"STRO_E"] <- ifelse(is.na(bldata$Stroke), 0,
                                  ifelse(bldata$Stroke==1,1,0))
  population[,"STRO_H"] <- ifelse(is.na(bldata$Stroke), 0,
                                  ifelse(bldata$Stroke==1,0,1))
  population[,"MI_E"] <- ifelse(is.na(bldata$MI), 0,
                                ifelse(bldata$MI==1,1,0))
  population[,"MI_H"] <- ifelse(is.na(bldata$MI), 0,
                                ifelse(bldata$MI==1,0,1))
  population[,"MI2_E"] <- 0
  population[,"MI2_H"] <- 0
  
  population[,"BLIND_E"] <- ifelse(is.na(bldata$Blindness), 0,
                                   ifelse(bldata$Blindness==1,1,0))
  population[,"BLIND_H"] <- ifelse(is.na(bldata$Blindness), 0,
                                   ifelse(bldata$Blindness==1,0,1))
  population[,"ULCER_E"] <- ifelse(is.na(bldata$Ulcer), 0,
                                   ifelse(bldata$Ulcer==1,1,0))
  population[,"ULCER_H"] <- ifelse(is.na(bldata$Ulcer), 0,
                                   ifelse(bldata$Ulcer==1,0,1))
  
  
  population[,"MMALB_H"] <- ifelse(bldata$Albuminuria=="Y",1,0)
  population[,"ATFIB_H"] <- ifelse(bldata$AF=="Y",1,0)
  population[,"PVD_H"] <- ifelse(bldata$PVD=="Y",1,0)
  
  #Add in diabetes treatment specific biomarkers
  population[,"HAEM"] <- bldata$Haemoglobin 
  population[,"WBC"] <- bldata$WBC
  population[,"eGFR"] <- bldata$eGFR
  population[,"eGFR_U_60"] <- ifelse(population[,"eGFR"]<60,population[,"eGFR"],60)
  population[,"eGFR_O_60"] <- ifelse(population[,"eGFR"]>60,population[,"eGFR"]-60,0)
  population[,"HEART_R"] <- bldata$Heart.rate
  population[,"BMI_U_18_5"] <- ifelse(population[,"BMI"]<18.5,1,0)
  population[,"BMI_O_E_25"] <- ifelse(population[,"BMI"]>=25,1,0)
  population[,"LDL"] <- bldata$LDL
  population[,"LDL_O_35"] <- ifelse(population[,"LDL"]>3.5, population[,"LDL"]-3.5,0)
  #add in diabetes treatment specific chars
  population[,"DIAB_DUR"] <- bldata$Duration.of.diabetes
  
  #add in values for first observations for each risk factor needed
  #in the absence of information this will be the baseline value
  population[,"HBA_0"] <- population[,"HBA"]
  population[,"LDL_0"] <- population[,"LDL"]
  population[,"HDL_0"] <- population[, "HDL"]
  population[,"HEART_R_0"] <- population[,"HEART_R"]
  population[,"HAEM_0"] <- population[,"HAEM"]
  population[,"WBC_0"] <- population[,"WBC"]
  population[,"eGFR_0"] <- population[,"eGFR"]
  population[,"SMO_0"] <- population[,"SMO"]
  population[,"SBP_0"] <- population[,"SBP"]
  
  #Make all cause death missing, as missing indicates someone is alive in the code
  population[, "F_ALLCAUSE"] <- NA
  #Make smoking missing, so you can easily see whether smoking is assessed
  population[,"p_SMO"]<- NA
  
  #Remove the temporary bldata object
  rm(bldata)
  
  return(population)
}


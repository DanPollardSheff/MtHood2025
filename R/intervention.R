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


##'@param n_ is the number of patients in the model
##'@param treatment_ is a text term indicating the current treatment option
##'@param parameter_ is the row of the parameter matrix
##'@param endtime_ is a number indicating how many years the simulation is being run for
##'@param GlobalVars_ is the Global Variables matrix, this allows the duration of treatment 
##' effect to be controlled
##'@param attend_se_ is a TRUE/FALSE vector that indicates whether a patient attends
##' a structured education course in each year 
##'@return INTE_A1c is a matrix that gives the reduction in A1c for each patient 
##'in each year compared to a patient on a "normal" trajectory

initialise_intervention_dt_HbA1c <- function(n_,
                                             treatment_, 
                                             parameter_,
                                             endtime_, 
                                             GlobalVars_, 
                                             attend_se_) {
  INTE_A1c <- matrix(data=0, nrow = n_, ncol =endtime_+2)
  INTE_A1c[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  if(treatment_ == "Mt_HOOD_RS_A1c" | treatment_ == "Mt_HOOD_RS_ALL"){
    INTE_A1c[,2:(endtime_+2)] <- -0.5
  }else{ #if no treatment option is selected, leave them at baseline values
    INTE_A1c[,2:(endtime_+2)] <- 0
  }
  
  return(INTE_A1c)
}

##'@param n_ is the number of patients in the model
##'@param treatment_ is a text term indicating the current treatment option
##'@param parameter_ is the row of the parameter matrix
##'@param endtime_ is a number indicating how many years the simulation is being run for
##'@param GlobalVars_
##'@param attend_se_ is three column matrix. Column 1 is patient ID, column 2 is 
##'a 0,1 vector with 1 indicating attendance at an SE course in the first model year,
##'column 3 is a 0,1 vector with 1 indicating attendance at an SE course in the second
##'model year 
##'@return INTE_BMI is a matrix that gives the reduction in BMI for each patient 
##'in each year

initialise_intervention_dt_BMI <- function(n_,
                                           treatment_, 
                                           parameter_,
                                           endtime_,
                                           GlobalVars_, 
                                           attend_se_) {
  INTE_BMI <- matrix(data=0, nrow = n_, ncol =endtime_+2)
  INTE_BMI[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  if(treatment_ == "Mt_HOOD_RS_BMI"|treatment_ == "Mt_HOOD_RS_ALL"){
    INTE_BMI[,2:(endtime_+2)] <- -1
  }else { #If a valid option is not selected, do not use a treatment effect
    INTE_BMI[,2:(endtime_+2)] <- 0
  }
  
  return(INTE_BMI)
}

##'@param n_ is the number of patients in the model
##'@param treatment_ is a text term indicating the current treatment option
##'@param parameter_ is the row of the parameter matrix
##'@param endtime_ is a number indicating how many years the simulation is being run for
##'@param GlobalVars_
##'@param attend_se_ is three column matrix. Column 1 is patient ID, column 2 is 
##'a 0,1 vector with 1 indicating attendance at an SE course in the first model year,
##'column 3 is a 0,1 vector with 1 indicating attendance at an SE course in the second
##'model year 
##'@return INTE_SBP is a matrix that gives the reduction in SBP(mmHg) for each patient 
##'in each year

initialise_intervention_dt_SBP <- function(n_,
                                           treatment_, 
                                           parameter_,
                                           endtime_,
                                           GlobalVars_, 
                                           attend_se_) {
  INTE_SBP <- matrix(data=0, nrow = n_, ncol =endtime_+2)
  INTE_SBP[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  if(treatment_ == "Mt_HOOD_RS_SBP" | treatment_ == "Mt_HOOD_RS_ALL"){
    INTE_SBP[,2:(endtime_+2)] <- -10
  }else{
    #If no option is specified, set treatment effects to 0
    INTE_SBP[,2:(endtime_+2)] <-0 
  }
  
  return(INTE_SBP)
}

##'@param n_ is the number of patients in the model
##'@param treatment_ is a text term indicating the current treatment option
##'@param parameter_ is the row of the parameter matrix
##'@param endtime_ is a number indicating how many years the simulation is being run for
##'@param GlobalVars_
##'@param attend_se_ is three column matrix. Column 1 is patient ID, column 2 is 
##'a 0,1 vector with 1 indicating attendance at an SE course in the first model year,
##'column 3 is a 0,1 vector with 1 indicating attendance at an SE course in the second
##'model year 
##'@return INTE_HDL is a matrix that gives the change in high density lipoprotein cholesterol 
##'for each patient in each year
initialise_intervention_dt_HDL <- function(n_,
                                           treatment_, 
                                           parameter_,
                                           endtime_,
                                           GlobalVars_, 
                                           attend_se_) {
  INTE_HDL <- matrix(data=0, nrow = n_, ncol =(endtime_+2))
  INTE_HDL[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  
  if(treatment_ == "Mt_HOOD_RS_HDL" | treatment_ == "Mt_HOOD_RS_ALL"){
    INTE_HDL[,2:(endtime_+2)] <- 0
  }else{
    INTE_HDL[,2:(endtime_+2)] <- 0
  }
  
  return(INTE_HDL)
}

##'@param n_ is the number of patients in the model
##'@param treatment_ is a text term indicating the current treatment option
##'@param parameter_ is the row of the parameter matrix
##'@param endtime_ is a number indicating how many years the simulation is being run for
##'@return INTE_LDL is a matrix that gives the reduction in HDL cholesterol 
##'for each patient in each year

initialise_intervention_dt_LDL <- function(n_,
                                           treatment_, 
                                           parameter_,
                                           endtime_,
                                           GlobalVars_, 
                                           attend_se_) {
  INTE_LDL <- matrix(data=0, nrow = n_, ncol = (endtime_+2))
  INTE_LDL[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  
  if(treatment_ == "Mt_HOOD_RS_LDL" | treatment_ == "Mt_HOOD_RS_ALL"){
    INTE_LDL[,2:(endtime_+2)] <- -0.5
  }else{#If no valid option is selected, model no intervention effect
    INTE_LDL[,2:(endtime_+2)] <- 0
  }
  
  return(INTE_LDL)
}

##'@param Input_Prob_ is a numeric number between 0 and 1
##'@param HR_ is a hazard ratio 
##'@return output_prob is the new probability after applying the hazard ratio

Hazard_Ratio_Intervention <- function(Input_Prob_,
                                      HR_){

  input_rate <- -log(1-Input_Prob_)/1
  output_rate <- input_rate*HR_
  output_prob <- 1-exp(-output_rate*1)
  
  return(output_prob)
  
}
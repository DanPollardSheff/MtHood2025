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




##'@param population_, is the population matrix
##'@param endtime_, is the number of years to run this model for
##'@return A1c is the matrix giving the HbA1c trajectory for each individual 
##'over time

constant_A1c <- function(population_, endtime_){
  
  #set up a matrix to store the results
  A1c <- matrix(data=NA, nrow = length(population_[,"ID"]), ncol = endtime_ + 2)
  #make the first column the ID's of the patients
  A1c[,1] <- population_[,"ID"]
  #Make the second column the baseline HbA1c
  A1c[,2:(endtime_ + 2)] <- population_[,"HBA"]
  #Return the A1c matrix
  return(A1c)
  
}

##'@param population_, is the population matrix
##'@param endtime_, is the number of years to run this model for
##'@return BMI is the matrix giving the BMI trajectory for each individual 
##'over time

constant_BMI <- function(population_, endtime_){
  
  #set up a matrix to store the results
  BMI <- matrix(data=NA, nrow = length(population_[,"ID"]), ncol = endtime_ + 2)
  #make the first column the ID's of the patients
  BMI[,1] <- population_[,"ID"]
  #Make the second column the baseline BMI
  BMI[,2:(endtime_ + 2)] <- population_[,"BMI"]
  #Return the BMI matrix
  return(BMI)
  
}

##'@param population_, is the population matrix
##'@param endtime_, is the number of years to run this model for
##'@return LDL is the matrix giving the LDL trajectory for each individual 
##'over time

constant_LDL <- function(population_, endtime_){
  
  #set up a matrix to store the results
  LDL <- matrix(data=NA, nrow = length(population_[,"ID"]), ncol = endtime_ + 2)
  #make the first column the ID's of the patients
  LDL[,1] <- population_[,"ID"]
  #Make the second column the baseline LDL
  LDL[,2:(endtime_ + 2)] <- population_[,"LDL"]
  #Return the LDL matrix
  return(LDL)
  
}

##'@param population_, is the population matrix
##'@param endtime_, is the number of years to run this model for
##'@return HDL is the matrix giving the HDL trajectory for each individual 
##'over time

constant_HDL <- function(population_, endtime_){
  
  #set up a matrix to store the results
  HDL <- matrix(data=NA, nrow = length(population_[,"ID"]), ncol = endtime_ + 2)
  #make the first column the ID's of the patients
  HDL[,1] <- population_[,"ID"]
  #Make the second column the baseline HDL
  HDL[,2:(endtime_ + 2)] <- population_[,"HDL"]
  #Return the HDL matrix
  return(HDL)
  
}

##'@param population_, is the population matrix
##'@param endtime_, is the number of years to run this model for
##'@return SBP is the matrix giving the SBP trajectory for each individual 
##'over time

constant_SBP <- function(population_, endtime_){
  
  #set up a matrix to store the results
  SBP <- matrix(data=NA, nrow = length(population_[,"ID"]), ncol = endtime_ + 2)
  #make the first column the ID's of the patients
  SBP[,1] <- population_[,"ID"]
  #Make the second column the baseline Systolic Blood Pressure
  SBP[,2:(endtime_ + 2)] <- population_[,"SBP"]
  #Return the A1c matrix
  return(SBP)
  
}

##'@param population_, is the population matrix
##'@param endtime_, is the number of years to run this model for
##'@return HR is the matrix giving the Heart Rate trajectory for each individual 
##'over time


constant_HEARTRATE <- function(population_, endtime_){
  
  #set up a matrix to store the results
  HR <- matrix(data=NA, nrow = length(population_[,"ID"]), ncol = endtime_ + 2)
  #make the first column the ID's of the patients
  HR[,1] <- population_[,"ID"]
  #Make the second column the baseline Heart Rate
  HR[,2:(endtime_ + 2)] <- population_[,"HEART_R"]
  #Return the HR matrix
  return(HR)
  
}

##'@param population_, is the population matrix
##'@param endtime_, is the number of years to run this model for
##'@return WBC is the matrix giving the White Blood Cell Count trajectory for each individual 
##'over time

constant_WBC <- function(population_, endtime_){
  
  #set up a matrix to store the results
  WBC <- matrix(data=NA, nrow = length(population_[,"ID"]), ncol = endtime_ + 2)
  #make the first column the ID's of the patients
  WBC[,1] <- population_[,"ID"]
  #Make the second column the baseline White Blood Cell count
  WBC[,2:(endtime_ + 2)] <- population_[,"WBC"]
  #Return the WBC matrix
  return(WBC)
  
}

##'@param population_, is the population matrix
##'@param endtime_, is the number of years to run this model for
##'@return HAEM is the matrix giving the Hemoglobin Count trajectory for each individual 
##'over time

constant_HAEM <- function(population_, endtime_){
  
  #set up a matrix to store the results
  HAEM <- matrix(data=NA, nrow = length(population_[,"ID"]), ncol = endtime_ + 2)
  #make the first column the ID's of the patients
  HAEM[,1] <- population_[,"ID"]
  #Make the column baseline Hemoglobin
  HAEM[,2:(endtime_ + 2)] <- population_[,"HAEM"]
  #Return the HAEM matrix
  return(HAEM)
  
}

##'@param population_, is the population matrix
##'@param endtime_, is the number of years to run this model for
##'@return eGFR is the matrix giving the estimated glomerular filtration rate trajectory for each individual 
##'over time

constant_eGFR <- function(population_, endtime_){
  
  #set up a matrix to store the results
  eGFR <- matrix(data=NA, nrow = length(population_[,"ID"]), ncol = endtime_ + 2)
  #make the first column the ID's of the patients
  eGFR[,1] <- population_[,"ID"]
  #Make the column baseline eGFR
  eGFR[,2:(endtime_ + 2)] <- population_[,"eGFR"]
  #Return the eGFR matrix
  return(eGFR)
  
}
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

#    Mt Hood 2025 analyses using School for Public Health Research Diabetes Treatment Model version 3.1
#    Copyright (C) 2025   Pollard, Brennan

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

#Load everything needed to run the model
library(MASS)
library(VGAM)
library(parallel)
library(dplyr)
library(ggplot2)
set.seed(429)

#Create the Global options matrix
source("Global Options.R")
#Load in population and data files
source("all_model_files.R")

################################################################################
#Check if there is a folder called Results, and if not make it (it will be 
#ignored by git version control and anything saved here will not be on Github)
folder_check <- file.exists("Results/Challenge 2")
if(folder_check==F){
  dir.create("Results/Challenge 2")
}
################################################################################
pop_cont <- build_population_MtHood2025_C2(MtHood2025C2Data,"BC_Control",PopulationVariables)
#Produce standard results for error checking
GlobalVars["Results_output", "Value"] <- "Patient Level"

test_pl_bc_cont <- run_simulation_MtHood2025_C2(pop_cont,
                                     parameter,
                                     5,
                                     "BC_Control",
                                     GlobalVars,
                                     LifeTables,
                                     1,#SOUR, so 1 means deterministic model run
                                     MtHood2025C2Data
)

#get a true false vector for being in the control arm
controlarm <- MtHood2025C2Data$HDL.csv$Group == 1
interventionarm <- MtHood2025C2Data$HDL.csv$Group == 2
alive <- is.na(test_pl_bc_cont[,"F_ALLCAUSE"])
#Check the people alive for the 1st 5 years of simulation have their risk factor values at
#the start of year 4
sum(ifelse(test_pl_bc_cont[,"HDL"][alive]==MtHood2025C2Data$HDL.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_cont[,"HBA"][alive]==MtHood2025C2Data$A1cBC.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_cont[,"BMI"][alive]==MtHood2025C2Data$BMIBC.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_cont[,"SBP"][alive]==MtHood2025C2Data$SBPBC.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_cont[,"LDL"][alive]==MtHood2025C2Data$LDL.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_cont[,"WBC"][alive]==MtHood2025C2Data$WBC.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_cont[,"HAEM"][alive]==MtHood2025C2Data$HEAM.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_cont[,"eGFR"][alive]==MtHood2025C2Data$eGFR.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_cont[,"HEART_R"][alive]==MtHood2025C2Data$HR.csv[controlarm,"Year.4"][alive],0,1))

#Make a plan to check categorical variables, PL data frame is based on 0 (no event), 1 event. But data is Y, N

test_pl_bc_intv <- run_simulation_MtHood2025_C2(pop_cont,
                                                parameter,
                                                5,
                                                "BC_INTV",
                                                GlobalVars,
                                                LifeTables,
                                                1,#SOUR, so 1 means deterministic model run
                                                MtHood2025C2Data
)


alive <- is.na(test_pl_bc_intv[,"F_ALLCAUSE"])
#Check the people alive for the 1st 5 years of simulation have their risk factor values at
#the start of year 4
sum(ifelse(test_pl_bc_intv[,"HDL"][alive]==MtHood2025C2Data$HDL.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_intv[,"HBA"][alive]==MtHood2025C2Data$A1cBC.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_intv[,"BMI"][alive]==MtHood2025C2Data$BMIBC.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_intv[,"SBP"][alive]==MtHood2025C2Data$SBPBC.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_intv[,"LDL"][alive]==MtHood2025C2Data$LDL.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_intv[,"WBC"][alive]==MtHood2025C2Data$WBC.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_intv[,"HAEM"][alive]==MtHood2025C2Data$HEAM.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_intv[,"eGFR"][alive]==MtHood2025C2Data$eGFR.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_bc_intv[,"HEART_R"][alive]==MtHood2025C2Data$HR.csv[interventionarm,"Year.4"][alive],0,1))

#Make a plan to check categorical variables, PL data frame is based on 0 (no event), 1 event. But data is Y, N

test_pl_s2_cont <- run_simulation_MtHood2025_C2(pop_cont,
                                                parameter,
                                                5,
                                                "S2_Control",
                                                GlobalVars,
                                                LifeTables,
                                                1,#SOUR, so 1 means deterministic model run
                                                MtHood2025C2Data
)

#get a true false vector for being in the control arm
controlarm <- MtHood2025C2Data$HDL.csv$Group == 1
interventionarm <- MtHood2025C2Data$HDL.csv$Group == 2
alive <- is.na(test_pl_s2_cont[,"F_ALLCAUSE"])
#Check the people alive for the 1st 5 years of simulation have their risk factor values at
#the start of year 4
sum(ifelse(test_pl_s2_cont[,"HDL"][alive]==MtHood2025C2Data$HDL.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_cont[,"HBA"][alive]==MtHood2025C2Data$A1cS2.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_cont[,"BMI"][alive]==MtHood2025C2Data$BMIS2.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_cont[,"SBP"][alive]==MtHood2025C2Data$SBPS2.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_cont[,"LDL"][alive]==MtHood2025C2Data$LDL.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_cont[,"WBC"][alive]==MtHood2025C2Data$WBC.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_cont[,"HAEM"][alive]==MtHood2025C2Data$HEAM.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_cont[,"eGFR"][alive]==MtHood2025C2Data$eGFR.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_cont[,"HEART_R"][alive]==MtHood2025C2Data$HR.csv[controlarm,"Year.4"][alive],0,1))

#Make a plan to check categorical variables, PL data frame is based on 0 (no event), 1 event. But data is Y, N

test_pl_s2_intv <- run_simulation_MtHood2025_C2(pop_cont,
                                                parameter,
                                                5,
                                                "S2_INTV",
                                                GlobalVars,
                                                LifeTables,
                                                1,#SOUR, so 1 means deterministic model run
                                                MtHood2025C2Data
)


alive <- is.na(test_pl_s2_intv[,"F_ALLCAUSE"])
#Check the people alive for the 1st 5 years of simulation have their risk factor values at
#the start of year 4
sum(ifelse(test_pl_s2_intv[,"HDL"][alive]==MtHood2025C2Data$HDL.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_intv[,"HBA"][alive]==MtHood2025C2Data$A1cS2.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_intv[,"BMI"][alive]==MtHood2025C2Data$BMIS2.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_intv[,"SBP"][alive]==MtHood2025C2Data$SBPS2.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_intv[,"LDL"][alive]==MtHood2025C2Data$LDL.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_intv[,"WBC"][alive]==MtHood2025C2Data$WBC.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_intv[,"HAEM"][alive]==MtHood2025C2Data$HEAM.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_intv[,"eGFR"][alive]==MtHood2025C2Data$eGFR.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_s2_intv[,"HEART_R"][alive]==MtHood2025C2Data$HR.csv[interventionarm,"Year.4"][alive],0,1))

#Make a plan to check categorical variables, PL data frame is based on 0 (no event), 1 event. But data is Y, N

test_pl_ne_cont <- run_simulation_MtHood2025_C2(pop_cont,
                                                parameter,
                                                5,
                                                "NoEffect_Control",
                                                GlobalVars,
                                                LifeTables,
                                                1,#SOUR, so 1 means deterministic model run
                                                MtHood2025C2Data
)

#get a true false vector for being in the control arm
controlarm <- MtHood2025C2Data$HDL.csv$Group == 1
interventionarm <- MtHood2025C2Data$HDL.csv$Group == 2
alive <- is.na(test_pl_ne_cont[,"F_ALLCAUSE"])
#Check the people alive for the 1st 5 years of simulation have their risk factor values at
#the start of year 4
sum(ifelse(test_pl_ne_cont[,"HDL"][alive]==MtHood2025C2Data$HDL.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_cont[,"HBA"][alive]==MtHood2025C2Data$A1cNoEffect.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_cont[,"BMI"][alive]==MtHood2025C2Data$BMINoEffect.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_cont[,"SBP"][alive]==MtHood2025C2Data$SBPNoEffect.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_cont[,"LDL"][alive]==MtHood2025C2Data$LDL.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_cont[,"WBC"][alive]==MtHood2025C2Data$WBC.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_cont[,"HAEM"][alive]==MtHood2025C2Data$HEAM.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_cont[,"eGFR"][alive]==MtHood2025C2Data$eGFR.csv[controlarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_cont[,"HEART_R"][alive]==MtHood2025C2Data$HR.csv[controlarm,"Year.4"][alive],0,1))

#Make a plan to check categorical variables, PL data frame is based on 0 (no event), 1 event. But data is Y, N

test_pl_ne_intv <- run_simulation_MtHood2025_C2(pop_cont,
                                                parameter,
                                                5,
                                                "NoEffect_INTV",
                                                GlobalVars,
                                                LifeTables,
                                                1,#SOUR, so 1 means deterministic model run
                                                MtHood2025C2Data
)


alive <- is.na(test_pl_ne_intv[,"F_ALLCAUSE"])
#Check the people alive for the 1st 5 years of simulation have their risk factor values at
#the start of year 4
sum(ifelse(test_pl_ne_intv[,"HDL"][alive]==MtHood2025C2Data$HDL.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_intv[,"HBA"][alive]==MtHood2025C2Data$A1cNoEffect.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_intv[,"BMI"][alive]==MtHood2025C2Data$BMINoEffect.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_intv[,"SBP"][alive]==MtHood2025C2Data$SBPNoEffect.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_intv[,"LDL"][alive]==MtHood2025C2Data$LDL.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_intv[,"WBC"][alive]==MtHood2025C2Data$WBC.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_intv[,"HAEM"][alive]==MtHood2025C2Data$HEAM.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_intv[,"eGFR"][alive]==MtHood2025C2Data$eGFR.csv[interventionarm,"Year.4"][alive],0,1))
sum(ifelse(test_pl_ne_intv[,"HEART_R"][alive]==MtHood2025C2Data$HR.csv[interventionarm,"Year.4"][alive],0,1))

#Make a plan to check categorical variables, PL data frame is based on 0 (no event), 1 event. But data is Y, N

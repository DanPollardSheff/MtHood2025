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
folder_check <- file.exists("Results")
if(folder_check==F){
  dir.create("Results")
}
folder_check <- file.exists("Results/Challenge 1/Constant")
if(folder_check==F){
  dir.create("Results/Challenge 1/Constant")
}
folder_check <- file.exists("Results/Challenge 1/Trajectories")
if(folder_check==F){
  dir.create("Results/Challenge 1/Trajectories")
}
folder_check <- file.exists("Results/Challenge 2")
if(folder_check==F){
  dir.create("Results/Challenge 2")
}
################################################################################

##Model run - Men
set.seed(123)
popM  <- build_population_MtHood2025(as.numeric(GlobalVars["n", "Value"]), F, PopulationVariables)
rands <- generate_random(length(popM[,"ID"]))



controlM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Baseline", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)
write.csv(controlM, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, no intv, Male.csv")

A1cINTVM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_A1c", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(A1cINTVM, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, A1c only, Male.csv")


BMIINTVM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_BMI", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(BMIINTVM, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, BMI only, Male.csv")


SBPINTVM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_SBP", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(SBPINTVM, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, SBP only, Male.csv")

LDLINTVM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_LDL", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(LDLINTVM, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, LDL only, Male.csv")


allINTVM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_ALL", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(allINTVM, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, all interventions, Male.csv")


#Women
popF  <- build_population_MtHood2025(as.numeric(GlobalVars["n", "Value"]), T, PopulationVariables)


controlF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Baseline", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(controlF, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, no intv, Female.csv")


A1cINTVF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_A1c", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(A1cINTVF, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, A1c only, Female.csv")


BMIINTVF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_BMI", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(BMIINTVF, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, BMI only, Female.csv")


SBPINTVF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_SBP", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(SBPINTVF, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, SBP only, Female.csv")


LDLINTVF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_LDL", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(LDLINTVF, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, LDL only, Female.csv")


allINTVF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_ALL", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(allINTVF, "Results/Challenge 1/Constant/MtHoodReferenceSimulation, all interventions, Female.csv")

#Repeat with Trajectories
#Men
controlMtraj <- run_simulation_MtHood2025_C1_v2(popM,
                                             parameter, 
                                             40, 
                                             "Baseline", 
                                             GlobalVars,
                                             rands,
                                             LifeTables,
                                             1, #Probably not needed, but set SOUR to 1 so determinsitic will be set regardless
                                             MtHood2025C1Data)
write.csv(controlMtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, no intv, Male.csv")

A1cINTVMtraj <- run_simulation_MtHood2025_C1_v2(popM,
                           parameter, 
                           40, 
                           "Mt_HOOD_RS_A1c", 
                           GlobalVars,
                           rands,
                           LifeTables,
                           1, #Probably not needed, but set SOUR to 1 so determinsitic will be set regardless
                           MtHood2025C1Data)

write.csv(A1cINTVMtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, A1c only, Male.csv")


BMIINTVMtraj <- run_simulation_MtHood2025_C1_v2(popM,
                           parameter, 
                           40, 
                           "Mt_HOOD_RS_BMI", 
                           GlobalVars,
                           rands,
                           LifeTables,
                           1, #Probably not needed, but set SOUR to 1 so deterministic will be set regardless
                           MtHood2025C1Data)

write.csv(BMIINTVMtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, BMI only, Male.csv")


SBPINTVMtraj <- run_simulation_MtHood2025_C1_v2(popM,
                           parameter, 
                           40, 
                           "Mt_HOOD_RS_SBP", 
                           GlobalVars,
                           rands,
                           LifeTables,
                           1, #Probably not needed, but set SOUR to 1 so deterministic will be set regardless
                           MtHood2025C1Data)

write.csv(SBPINTVMtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, SBP only, Male.csv")

LDLINTVMtraj <- run_simulation_MtHood2025_C1_v2(popM,
                           parameter, 
                           40, 
                           "Mt_HOOD_RS_LDL", 
                           GlobalVars,
                           rands,
                           LifeTables,
                           1, #Probably not needed, but set SOUR to 1 so deterministic will be set regardless
                           MtHood2025C1Data)

write.csv(LDLINTVMtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, LDL only, Male.csv")


allINTVMtraj <- run_simulation_MtHood2025_C1_v2(popM,
                           parameter, 
                           40, 
                           "Mt_HOOD_RS_ALL", 
                           GlobalVars,
                           rands,
                           LifeTables,
                           1, #Probably not needed, but set SOUR to 1 so deterministic will be set regardless
                           MtHood2025C1Data)

write.csv(allINTVMtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, all interventions, Male.csv")
#Women
controlFtraj <- run_simulation_MtHood2025_C1_v2(popF,
                                                parameter, 
                                                40, 
                                                "Baseline", 
                                                GlobalVars,
                                                rands,
                                                LifeTables,
                                                1, #Probably not needed, but set SOUR to 1 so determinsitic will be set regardless
                                                MtHood2025C1Data)
write.csv(controlFtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, no intv, Female.csv")

A1cINTVFtraj <- run_simulation_MtHood2025_C1_v2(popF,
                                             parameter, 
                                             40, 
                                             "Mt_HOOD_RS_A1c", 
                                             GlobalVars,
                                             rands,
                                             LifeTables,
                                             1, #Probably not needed, but set SOUR to 1 so determinsitic will be set regardless
                                             MtHood2025C1Data)

write.csv(A1cINTVFtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, A1c only, Female.csv")


BMIINTVFtraj <- run_simulation_MtHood2025_C1_v2(popF,
                                             parameter, 
                                             40, 
                                             "Mt_HOOD_RS_BMI", 
                                             GlobalVars,
                                             rands,
                                             LifeTables,
                                             1, #Probably not needed, but set SOUR to 1 so deterministic will be set regardless
                                             MtHood2025C1Data)

write.csv(BMIINTVFtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, BMI only, Female.csv")


SBPINTVFtraj <- run_simulation_MtHood2025_C1_v2(popF,
                                             parameter, 
                                             40, 
                                             "Mt_HOOD_RS_SBP", 
                                             GlobalVars,
                                             rands,
                                             LifeTables,
                                             1, #Probably not needed, but set SOUR to 1 so deterministic will be set regardless
                                             MtHood2025C1Data)

write.csv(SBPINTVFtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, SBP only, Female.csv")

LDLINTVFtraj <- run_simulation_MtHood2025_C1_v2(popF,
                                             parameter, 
                                             40, 
                                             "Mt_HOOD_RS_LDL", 
                                             GlobalVars,
                                             rands,
                                             LifeTables,
                                             1, #Probably not needed, but set SOUR to 1 so deterministic will be set regardless
                                             MtHood2025C1Data)

write.csv(LDLINTVFtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, LDL only, Female.csv")


allINTVFtraj <- run_simulation_MtHood2025_C1_v2(popF,
                                             parameter, 
                                             40, 
                                             "Mt_HOOD_RS_ALL", 
                                             GlobalVars,
                                             rands,
                                             LifeTables,
                                             1, #Probably not needed, but set SOUR to 1 so deterministic will be set regardless
                                             MtHood2025C1Data)

write.csv(allINTVFtraj, "Results/Challenge 1/Trajectories/MtHoodReferenceSimulation, all interventions, Female.csv")

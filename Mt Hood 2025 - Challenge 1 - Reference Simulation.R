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
write.csv(controlM, "Results/MtHoodReferenceSimulation, no intv, Male.csv")

A1cINTVM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_A1c", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(A1cINTVM, "Results/MtHoodReferenceSimulation, A1c only, Male.csv")


BMIINTVM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_BMI", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(BMIINTVM, "Results/MtHoodReferenceSimulation, BMI only, Male.csv")


SBPINTVM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_SBP", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(SBPINTVM, "Results/MtHoodReferenceSimulation, SBP only, Male.csv")

LDLINTVM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_LDL", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(LDLINTVM, "Results/MtHoodReferenceSimulation, LDL only, Male.csv")


allINTVM <- run_simulation(popM,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_ALL", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(allINTVM, "Results/MtHoodReferenceSimulation, all interventions, Male.csv")


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

write.csv(controlF, "Results/MtHoodReferenceSimulation, no intv, Female.csv")


A1cINTVF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_A1c", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(A1cINTVF, "Results/MtHoodReferenceSimulation, A1c only, Female.csv")


BMIINTVF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_BMI", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(BMIINTVF, "Results/MtHoodReferenceSimulation, BMI only, Female.csv")


SBPINTVF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_SBP", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(SBPINTVF, "Results/MtHoodReferenceSimulation, SBP only, Female.csv")


LDLINTVF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_LDL", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(LDLINTVF, "Results/MtHoodReferenceSimulation, LDL only, Female.csv")


allINTVF <- run_simulation(popF,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_ALL", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

write.csv(allINTVF, "Results/MtHoodReferenceSimulation, all interventions, Female.csv")


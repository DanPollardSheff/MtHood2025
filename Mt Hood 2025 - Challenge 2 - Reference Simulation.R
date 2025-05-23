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
GlobalVars_["Mt Hood Utilities", "Value"] <- "C2"

test_pl <- run_simulation_MtHood2025_C2(pop_cont,
                                     parameter,
                                     5,
                                     "BC_Control",
                                     GlobalVars,
                                     LifeTables,
                                     1,#SOUR, so 1 means deterministic model run
                                     MtHood2025C2Data
)

GlobalVars["Results_output", "Value"] <- "MtHood2025C2"


test <- run_simulation_MtHood2025_C2(pop_cont,
                            parameter,
                            5,
                            "BC_Control",
                            GlobalVars,
                            LifeTables,
                            1,#SOUR, so 1 means determinsitic model run
                            MtHood2025C2Data
)

#Return to default results
GlobalVars["Results_output", "Value"] <- ""

test_boot_cont <- run_model_bootstrap_MtHood2025_C2(pop_cont,
                                               parameter,
                                               50,
                                               "BC_Control",
                                               GlobalVars,
                                               LifeTables,
                                               MtHood2025C2Data,
                                               20)
#Note SOUR not set, as it is fixed at 1 within the bootstrapping function for Mt Hood problem

test_boot_intv <- run_model_bootstrap_MtHood2025_C2(pop_cont,
                                               parameter,
                                               50,
                                               "BC_INTV",
                                               GlobalVars,
                                               LifeTables,
                                               MtHood2025C2Data,
                                               20)
#Note SOUR not set, as it is fixed at 1 within the bootstrapping function for Mt Hood problem
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


library(MASS)
library(VGAM)
library(parallel)
library(dplyr)
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
GlobalVars["Trajectory", "Value"] <- "Constant"

#Make population
population <- build_population_MtHood2025(as.numeric(GlobalVars["n", "Value"]), F, PopulationVariables)
#Make the random numbers
rands <- generate_random(length(population[,"ID"]))
#Set Global options
GlobalVars["Results_output", "Value"] <- "Detailed"
GlobalVars["run_psa", "Value"] <- F
GlobalVars["Mt Hood Utilities", "Value"] <- "TRUE"

MtHoodUtil <- run_simulation(population,
                             parameter, 
                             9, 
                             "Baseline", 
                             GlobalVars,
                             rands,
                             LifeTables,
                             1)

MtHoodUtil

GlobalVars["Mt Hood Utilities", "Value"] <- "FALSE"

DefaultUtil <- run_simulation(population,
                              parameter, 
                              9, 
                              "Baseline", 
                              GlobalVars,
                              rands,
                              LifeTables,
                              1)
DefaultUtil 
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

#Change the Global options to have the right values
#I want the patient level results - i.e. their values at the time of their death
GlobalVars["Results_output", "Value"] <- "Patient Level"
#Do this in 10,000 person groups to assess satbility
GlobalVars["n", "Value"] <- 10000

##Assess Stability
###Make 50,000 men to assess stability 
###Target: Incremental QALYs of all interventions v control
set.seed(123)
popM1  <- build_population_MtHood2025(as.numeric(GlobalVars["n", "Value"]), F, PopulationVariables)
rands <- generate_random(length(popM1[,"ID"]))



control <- run_simulation(popM1,
                          parameter, 
                          40, 
                          "Baseline", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

set.seed(123)
write.csv(control,"Results/ControlMen.csv")

allINTV <- run_simulation(popM1,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_ALL", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)


control <- control[,c("ID","QALY")]
allINTV <- allINTV[,c("ID", "QALY")]

#2nd set of 10,000
set.seed(63166)
popM2  <- build_population_MtHood2025(as.numeric(GlobalVars["n", "Value"]), F, PopulationVariables)
rands <- generate_random(length(popM1[,"ID"]))

control2 <- run_simulation(popM2,
                          parameter, 
                          40, 
                          "Baseline", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)

allINTV2 <- run_simulation(popM2,
                          parameter, 
                          40, 
                          "Mt_HOOD_RS_ALL", 
                          GlobalVars,
                          rands,
                          LifeTables,
                          1)


control2 <- control2[,c("ID","QALY")]
allINTV2 <- allINTV2[,c("ID", "QALY")]

control <- rbind(control,control2)
allINTV <- rbind(allINTV,allINTV2)

set.seed(617)
popM3  <- build_population_MtHood2025(as.numeric(GlobalVars["n", "Value"]), F, PopulationVariables)
rands <- generate_random(length(popM1[,"ID"]))

control3 <- run_simulation(popM3,
                           parameter, 
                           40, 
                           "Baseline", 
                           GlobalVars,
                           rands,
                           LifeTables,
                           1)

allINTV3 <- run_simulation(popM3,
                           parameter, 
                           40, 
                           "Mt_HOOD_RS_ALL", 
                           GlobalVars,
                           rands,
                           LifeTables,
                           1)


control3 <- control3[,c("ID","QALY")]
allINTV3 <- allINTV3[,c("ID", "QALY")]

control <- rbind(control,control3)
allINTV <- rbind(allINTV,allINTV3)


incremental <- control
incremental[,"QALY"] <- allINTV[,"QALY"] - control[,"QALY"]
incremental[,"ID"] <- 1:length(incremental[,"ID"]) #reset the ID variable, so that each person has a unique ID, e.g. 1st patient in simulation 2, is now 10001 not 1. 
running_mean_QALYs <- cumsum(incremental[,"QALY"])/seq_along(incremental[,"QALY"])
incremental <- cbind(incremental,running_mean_QALYs)

ggplot(incremental,
       aes(x=ID,
           y=running_mean_QALYs))+
  geom_line()+
  ylim(0.6,0.65)

##Looks stable at ~10k onwards. Run with 20k
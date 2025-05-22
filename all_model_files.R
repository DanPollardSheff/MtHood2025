#    Embedding RCT Health Economic Analysis using the Sheffield Type 2 Diabetes Treatment Model - version 3
#    Copyright (C) 2023  Pollard, Pidd, Breeze, Brennan, Thomas

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



sapply(list.files("R",full.names=T), source)

#Load parameters
load("data/parameter.rda")
#Load population variables
PopulationVariables <- read.csv("data/PopulationVariables.csv")
#Load life tables
LifeTables <- as.matrix(read.csv("data/LifeTables.csv"))
#Makes the Lifetables numeric, like the population matrix
LifeTables[,"AGE"] <- as.numeric(LifeTables[,"AGE"] )
LifeTables[,"FEMALE"] <- as.numeric(LifeTables[,"FEMALE"])

##Mt Hood 2025 Challenge 2, load everything in the data MtHood2025 challenge2 into the Global Envrioment as 1 list
#Get the names of all files for Mt Hood 2025 Challenge 2
allfiles <- list.files(path = "dataMtHood2025challenge2/", pattern = "\\.csv$", full.names=TRUE)
#Read all files into the global enviroment, making them a list
MtHood2025C2Data <- lapply(allfiles, read.csv)
#Make each object in the list have the same name as the .csv file
names(MtHood2025C2Data) <- basename(allfiles)


##Mt Hood 2025 Challenge 1, load everything in the data MtHood2025 challenge2 into the Global Envrioment as 1 list
#Get the names of all files for Mt Hood 2025 Challenge 2
allfiles <- list.files(path = "dataMtHood2025challenge1/", pattern = "\\.csv$", full.names=TRUE)
#Read all files into the global enviroment, making them a list
MtHood2025C1Data <- lapply(allfiles, read.csv)
#Make each object in the list have the same name as the .csv file
names(MtHood2025C1Data) <- basename(allfiles)
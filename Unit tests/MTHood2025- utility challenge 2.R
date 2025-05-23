#Create the Global options matrix
source("Global Options.R")
#Load in population and data files
source("all_model_files.R")

#Make a population that is one person and Male
popM <- build_population_MtHood2025 (1,F, PopulationVariables)

#Default in MTHood for BMI is 28Kg/m2, change to 25 for testing
popM[,"BMI"] <- 25
popM[,"BMI_0"] <- 25

#see if the baseline value is correct using the same function used in the model

nocomps_notoverweight <- calculate_QALYs_MtHood2025_C2(popM,
                                                    0,
                                                    T,
                                                    GlobalVars)
if(nocomps_notoverweight[,"EQ5D"]!= 0.807){
  stop("Error in the baseline probability, for the mean values in  Mt Hood challenge 2 simulation  2025")
}

rm(nocomps_notoverweight)

#Change BMI to 26
popM[,"BMI"] <- 26
popM[,"BMI_0"] <- 26

nocomps_BMI26 <- calculate_QALYs_MtHood2025_C2(popM,
                                            0,
                                            T,
                                            GlobalVars)

if(nocomps_BMI26[,"EQ5D"]!= 0.807){
  stop("Error in the application of BMI decrements, for the mean values in Mt Hood challenge 2 simulation 2025")
}
rm(nocomps_BMI26)

#Change BMI to 27 & check for errors
popM[,"BMI"] <- 27
popM[,"BMI_0"] <- 27

nocomps_BMI27 <- calculate_QALYs_MtHood2025_C2(popM,
                                            0,
                                            T,
                                            GlobalVars)

if(nocomps_BMI27[,"EQ5D"]!= 0.807){
  stop("Error in the application of BMI decrements, for the mean values in Mt Hood reference simulation 2025")
}
rm(nocomps_BMI27)

#Reset BMI back to 25, makes checking of events much easier
popM[,"BMI"] <- 25
popM[,"BMI_0"] <- 25


##Set BMI back to 25, and check for Blindness
popM[,"BLIND_E"] <- 1
popM[,"BLIND_H"] <- 0

blindthisyear <- calculate_QALYs_MtHood2025_C2(popM,
                                            0,
                                            T,
                                            GlobalVars)

if(blindthisyear[,"EQ5D"]!= 0.807){
  stop("Error in the application of blindness decrement, for the mean values in Mt Hood reference simulation 2025") 
}
rm(blindthisyear)

#Change to a blindness history
popM[,"BLIND_E"] <- 0
popM[,"BLIND_H"] <- 1

blindlastyear <- calculate_QALYs_MtHood2025_C2(popM,
                                            0,
                                            T,
                                            GlobalVars)

if(blindlastyear[,"EQ5D"]!= 0.807){
  stop("Error in the application of blindness decrement, for the mean values in Mt Hood reference simulation 2025") 
}
rm(blindlastyear)

#Change to a blindness history
popM[,"BLIND_E"] <- 0
popM[,"BLIND_H"] <- 0

###MMALB
#This year
popM[,"MMALB_E"] <- 1
popM[,"MMALB_H"] <- 0

MMALBthisyear <- calculate_QALYs_MtHood2025_C2(popM,
                                            0,
                                            T,
                                            GlobalVars)

if(MMALBthisyear[,"EQ5D"]!= 0.807){
  stop("Error in the application of MMALB decrement, for the mean values in Mt Hood reference simulation 2025") 
}
rm(MMALBthisyear)

#Last year
popM[,"MMALB_H"] <- 0
popM[,"MMALB_E"] <- 1

MMALBlastyear <- calculate_QALYs_MtHood2025_C2(popM,
                                            0,
                                            T,
                                            GlobalVars)

if(MMALBlastyear[,"EQ5D"]!= 0.807){
  stop("Error in the application of MMALB decrement, for the mean values in Mt Hood reference simulation 2025") 
}
rm(MMALBlastyear)



##
popM[,"RENAL_E"] <- 1
popM[,"RENAL_H"] <- 0

RENALthisyear <- calculate_QALYs_MtHood2025_C2(popM,
                                            0,
                                            T, 
                                            GlobalVars)


if(RENALthisyear[,"EQ5D"]!= 0.807 - 0.330){
  stop("Error in the application of renal decrement, year 1, for the mean values in Mt Hood reference simulation 2025") 
}

rm(RENALthisyear)

popM[,"RENAL_E"] <- 0
popM[,"RENAL_H"] <- 1

RENALnextyear <- calculate_QALYs_MtHood2025_C2(popM,
                                            0,
                                            T, 
                                            GlobalVars)

if(RENALnextyear[,"EQ5D"]!= 0.807 - 0.330){
  stop("Error in the application of renal decrement, subsequent years, for the mean values in Mt Hood reference simulation 2025") 
}

rm(RENALnextyear)

#Reset renal failures / MMALB
popM[,"RENAL_E"] <- 0
popM[,"RENAL_H"] <- 0
popM[,"MMALB_H"] <- 0
popM[,"MMALB_E"] <- 0

##PVD
popM[,"PVD_E"] <- 1
popM[,"PVD_H"] <- 0


PVDthisyear <- calculate_QALYs_MtHood2025_C2(popM,
                                            0,
                                            T, 
                                            GlobalVars)

if(PVDthisyear[,"EQ5D"]!= 0.807){
  stop("Error in the application of PVD decrement, for the mean values in Mt Hood reference simulation 2025")
}

rm(PVDthisyear)

popM[,"PVD_E"] <- 0
popM[,"PVD_H"] <- 1

PVDnextyear <- calculate_QALYs_MtHood2025_C2(popM,
                                          0,
                                          T, 
                                          GlobalVars)

if(PVDnextyear[,"EQ5D"]!= 0.807){
  stop("Error in the application of PVD decrement, for the mean values in Mt Hood reference simulation 2025")
}

rm(PVDnextyear)



###Active ulcer
popM[,"ULCER_E"] <- 1
popM[,"ULCER_H"] <- 0

ulceryear <- calculate_QALYs_MtHood2025_C2(popM,
                                        0,
                                        T, 
                                        GlobalVars)

if(ulceryear[,"EQ5D"]!= 0.807 -0.210){
  stop("Error in the application of Ulcer decrement, for the mean values in Mt Hood reference simulation 2025")
}

rm(ulceryear)

##History of ulcer
popM[,"ULCER_E"] <- 0
popM[,"ULCER_H"] <- 1

ulcerhist <- calculate_QALYs_MtHood2025_C2(popM,
                                        0,
                                        T, 
                                        GlobalVars)

#For ulcer history, the ulcer decrement should be ignored so the decrmenet for PVD should still apply
if(ulcerhist[,"EQ5D"]!= 0.807 -0.210){
  stop("Error in the application of Ulcer decrement, for ulcers over 1 year ago, for the mean values in Mt Hood reference simulation 2025")
}

rm(ulcerhist)

popM[,"ULCER_E"] <- 0
popM[,"ULCER_H"] <- 0

#Amputation event
popM[,"AMP_E"] <- 1
popM[,"AMP_H"] <- 0

ampevent <- calculate_QALYs_MtHood2025_C2(popM,
                                       0,
                                       T, 
                                       GlobalVars)

if(ampevent[,"EQ5D"]!= 0.807 -0.172){
  stop("Error in the application of amputation decrement, for amputations in the last year, for the mean values in Mt Hood reference simulation 2025")
}

rm(ampevent)

##amp hist
popM[,"AMP_E"] <- 0
popM[,"AMP_H"] <- 1

amphist <- calculate_QALYs_MtHood2025_C2(popM,
                                       0,
                                       T, 
                                       GlobalVars)
#If an amputation history, PVD decrements should still apply
if(amphist[,"EQ5D"]!= 0.807 -0.172){
  stop("Error in the application of amputation decrement, for amputations over 1 year ago, for the mean values in Mt Hood reference simulation 2025")
}

rm(amphist)
popM[,"AMP_E"] <- 0
popM[,"AMP_H"] <- 0

##2nd amputation event
popM[,"AMP2_E"] <- 1
popM[,"AMP2_H"] <- 0

amp2event <- calculate_QALYs_MtHood2025_C2(popM,
                                        0,
                                        T, 
                                        GlobalVars)

if(amp2event[,"EQ5D"]!= 0.807 -0.172){
  stop("Error in the application of amputation decrement, for a 2nd amputation this year, for the mean values in Mt Hood reference simulation 2025")
}

rm(amp2event)

popM[,"AMP2_E"] <- 0
popM[,"AMP2_H"] <- 1

amp2hist <- calculate_QALYs_MtHood2025_C2(popM,
                                        0,
                                        T, 
                                        GlobalVars)
#If an amputation history, PVD decrements should still apply
if(amp2hist[,"EQ5D"]!= 0.807 -0.172){
  stop("Error in the application of amputation decrement, for a 2nd amputation in previous years, for the mean values in Mt Hood reference simulation 2025")
}

rm(amp2hist)

#Reset Amputation
popM[,"AMP2_E"] <- 0
popM[,"AMP2_H"] <- 0
#Reset PVD
popM[,"PVD_E"] <- 0
popM[,"PVD_H"] <- 0

##Stroke
popM[,"STRO_E"] <- 1
popM[,"STRO_H"] <- 0

stroevent <- calculate_QALYs_MtHood2025_C2(popM,
                                        0,
                                        T, 
                                        GlobalVars)

if(stroevent[,"EQ5D"]!= 0.807 -0.165){
  stop("Error in the application of stroke decrement, for 1st stroke this year, for the mean values in Mt Hood reference simulation 2025")
}

rm(stroevent)

popM[,"STRO_E"] <- 0
popM[,"STRO_H"] <- 1

strohist <- calculate_QALYs_MtHood2025_C2(popM,
                                        0,
                                        T, 
                                        GlobalVars)

if(strohist[,"EQ5D"]!= 0.807 -0.165){
  stop("Error in the application of stroke decrement, for 1st stroke previous years, for the mean values in Mt Hood reference simulation 2025")
}

rm(strohist)

popM[,"STRO2_E"] <- 1
popM[,"STRO2_H"] <- 0

stro2event <- calculate_QALYs_MtHood2025_C2(popM,
                                       0,
                                       T, 
                                       GlobalVars)

if(stro2event[,"EQ5D"]!= 0.807 -0.165){
  stop("Error in the application of stroke decrement, for 2nd stroke this years, for the mean values in Mt Hood reference simulation 2025")
}

rm(stro2event)

popM[,"STRO2_E"] <- 0
popM[,"STRO2_H"] <- 1

stro2hist <- calculate_QALYs_MtHood2025_C2(popM,
                                         0,
                                         T, 
                                         GlobalVars)

if(stro2hist[,"EQ5D"]!= 0.807 -0.165){
  stop("Error in the application of stroke decrement, for 2nd stroke this years, for the mean values in Mt Hood reference simulation 2025")
}

rm(stro2hist)
#Reset Stroke History
popM[,"STRO2_E"] <- 0
popM[,"STRO2_H"] <- 0
popM[,"STRO_E"] <- 0
popM[,"STRO_H"] <- 0

##CHD

popM[,"MI_E"] <- 1
popM[,"MI_H"] <- 0

MI1event <- calculate_QALYs_MtHood2025_C2(popM,
                                       0,
                                       T, 
                                       GlobalVars)

if(MI1event[,"EQ5D"]!= 0.807 - 0.065){
  stop("Error in the application of 1st MI decrement, for the event year, for the mean values in Mt Hood reference simulation 2025")
}

rm(MI1event)

popM[,"MI_E"] <- 0
popM[,"MI_H"] <- 1

MI1hist <- calculate_QALYs_MtHood2025_C2(popM,
                                       0,
                                       T, 
                                       GlobalVars)

if(MI1hist[,"EQ5D"]!= 0.807){
  stop("Error in the application of 1st MI decrement, for the event year, for the mean values in Mt Hood reference simulation 2025")
}

rm(MI1hist)

popM[,"MI2_E"] <- 1
popM[,"MI2_H"] <- 0

MI2event <- calculate_QALYs_MtHood2025_C2(popM,
                                      0,
                                      T, 
                                      GlobalVars)

if(MI2event[,"EQ5D"]!= 0.807 - 0.065){
  stop("Error in the application of 2nd MI decrement, for the event year, for the mean values in Mt Hood reference simulation 2025")
}

rm(MI2event)

popM[,"MI2_E"] <- 0
popM[,"MI2_H"] <- 1

MI2hist <- calculate_QALYs_MtHood2025_C2(popM,
                                       0,
                                       T, 
                                       GlobalVars)

if(MI2hist[,"EQ5D"]!= 0.807){
  stop("Error in the application of 2nd MI decrement, for historical events, for the mean values in Mt Hood reference simulation 2025")
}

rm(MI2hist)

popM[,"IHD_E"] <- 1
popM[,"IHD_H"] <- 0

IHDevent <- calculate_QALYs_MtHood2025_C2(popM,
                                      0,
                                      T, 
                                      GlobalVars)

if(IHDevent[,"EQ5D"]!= 0.807){
  stop("Error in the application of IHD decrement, for event year, for the mean values in Mt Hood reference simulation 2025")
}

rm(IHDevent)

popM[,"IHD_E"] <- 0
popM[,"IHD_H"] <- 1

IHDhist <- calculate_QALYs_MtHood2025_C2(popM,
                                       0,
                                       T, 
                                       GlobalVars)

if(IHDhist[,"EQ5D"]!= 0.807){
  stop("Error in the application of IHD decrement, for previous years, for the mean values in Mt Hood reference simulation 2025")
}

rm(IHDhist)

popM[,"CHF_E"] <- 1
popM[,"CHF_H"] <- 0

CHFevent <- calculate_QALYs_MtHood2025_C2(popM,
                                      0,
                                      T, 
                                      GlobalVars)

if(CHFevent[,"EQ5D"]!= 0.807-0.101){
  stop("Error in the application of CHF decrement, for event year, for the mean values in Mt Hood reference simulation 2025")
}

rm(CHFevent)

CHFhist <- calculate_QALYs_MtHood2025_C2(popM,
                                       0,
                                       T, 
                                       GlobalVars)

if(CHFhist[,"EQ5D"]!= 0.807-0.101){
  stop("Error in the application of CHF decrement, for previous years, for the mean values in Mt Hood reference simulation 2025")
}

rm(CHFhist)

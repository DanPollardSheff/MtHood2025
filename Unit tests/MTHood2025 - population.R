#Create the Global options matrix
source("Global Options.R")
#Load in population and data files
source("all_model_files.R")

#Make a population that is one person and Male
popM <- build_population_MtHood2025 (1,F, PopulationVariables)

#Start checks

#Check they are male and not female
if(popM[,"MALE"]!=1 & popM[,"FEMALE"]!=0 &popM[,"MEN"]!=0){
  stop("Error in Gender/Menopause status, male population")
}

#Age
if(popM[,"AGE"]!=66|popM[,"AGE_0"]!=66){
  stop("Error in Age")
}

#Diabetes Duration
if(popM[,"DIAB_DUR"]!=8){
  stop("Error in Diabetes Duration")
}

#Smoking
if(popM[,"SMO"]!=0|popM[,"SMO_0"]!=0){
  stop("Error in Smoking Status")
}

#Ethnicity - check they are white ethnicity
if(popM[,"AFRO"]!=0|popM[,"INDIAN"]!=0){
  stop("Error in Ethnic Status")
}

##HbA1c - check it is 7.5%
if(popM[,"HBA_0"]!=7.5|popM[,"HBA"]!=7.5){
  stop("Error in HbA1c values")
}

##SBP - check it is 145
if(popM[,"SBP_0"]!=145|popM[,"SBP"]!=145){
  stop("Error in SBP")
}

##Diastolic blood pressure is not used in this model

##Total cholesterol is not used in this model

##HDL

if(popM[,"HDL_0"]!=1.3|popM[,"HDL"]!=1.3){
  stop("Error in HDL values")
}

#LDL
if(popM[,"LDL_0"]!=3.0|popM[,"LDL"]!=3.0){
  stop("Error in LDL values")
}

##Triglycerides are not used in this model

##BMI - check it is 28Kg/m2
#Continuous
if(popM[,"BMI_0"]!=28|popM[,"BMI"]!=28){
  stop("Error in BMI values")
}
#Categorical
if(popM[,"BMI_U_18_5"]!=0|popM[,"BMI_O_E_25"]!=1){
  stop("Error in categorical BMI values")
}

#Albumin:creatine ratio is not used in this model

#PVD status
if(popM[,"PVD_E"]!=0|popM[,"PVD_H"]!=0){
  stop("Error in PVD status")
}

#MMALB status
if(popM[,"MMALB_E"]!=0|popM[,"MMALB_H"]!=0){
  stop("Error in MMALB status")
}

#ATFIB status
if(popM[,"ATFIB_E"]!=0|popM[,"ATFIB_H"]!=0){
  stop("Error in ATFIB status")
}

#eGFR
#Continuous
if(popM[,"eGFR"]!=70|popM[,"eGFR_0"]!=70){
  stop("Error in eGFR values")
}
#Knot
if(popM[,"eGFR_U_60"]!=60|popM[,"eGFR_O_60"]!=(70-60)){
  stop("Error in eGFR knot values")
}

##WBC 
if(popM[,"WBC"]!=7|popM[,"WBC_0"]!=7){
  stop("Error in WBC values")
}

##Haemoglobin
if(popM[,"HAEM"]!= 14 |popM[,"HAEM_0"]!=14){
  stop("Error in Haemoglobin values")
}

##Macrovascular disease
if(
  popM[,"MI_E"]!=0|
  popM[,"MI_H"]!=0|
  popM[,"MI2_E"]!=0|
  popM[,"MI2_H"]!=0|
  popM[,"STRO_E"]!=0|
  popM[,"STRO_H"]!=0|
  popM[,"STRO2_E"]!=0|
  popM[,"STRO2_H"]!=0|
  popM[,"IHD_E"]!=0|
  popM[,"IHD_H"]!=0|
  popM[,"CHF_E"]!=0|
  popM[,"CHF_H"]!=0
){
stop("Error in Macrovascular disease variables")  
}

##Microvascular disease
if(
  popM[,"BLIND_E"]!=0|
  popM[,"BLIND_H"]!=0|
  popM[,"ULCER_E"]!=0|
  popM[,"ULCER_H"]!=0|
  popM[,"RENAL_E"]!=0|
  popM[,"RENAL_H"]!=0|
  popM[,"AMP_E"]!=0|
  popM[,"AMP_H"]!=0|
  popM[,"AMP2_E"]!=0|
  popM[,"AMP2_H"]!=0
){
  stop("Error in Microvascular disease variables")
}

##Other health states
if(
  popM[,"CANB_E"]!=0|
  popM[,"CANB_H"]!=0|
  popM[,"CANC_E"]!=0|
  popM[,"CANC_H"]!=0|
  popM[,"OST_E"]!=0|
  popM[,"OST_H"]!=0|
  popM[,"DEP_E"]!=0
){
  stop("Error in other health state variables")
}

#Make a population that is one person and Female
popF <- build_population_MtHood2025 (1,T, PopulationVariables)

##Get the indices for which the columns are different
diffcols <- which(apply(popF!=popM, 2, all))

##check only these columns from popF
popF <- popF[,diffcols]
#Note this turns the object into a named number and not a matrix so the notation for looking up values does not need a preceding, to tell R to look at all rows, as there are no rows in this object now

if(popF["MALE"]!=0 & popF["FEMALE"]!=1 &popF["MEN"]!=1){
  stop("Error in Gender/Menopause status, female population")
}

##Repeat for 10 people
#Remove the original objects
rm(popM, popF)

#Make a population that is one person and Male
popM <- build_population_MtHood2025 (10,F, PopulationVariables)

#Start checks

#Check they are male and not female
if(mean(popM[,"MALE"])!=1 & mean(popM[,"FEMALE"])!=0 &mean(popM[,"MEN"])!=0){
  stop("Error in Gender/Menopause status, male population")
}

#Age
if(mean(popM[,"AGE"])!=66|mean(popM[,"AGE_0"])!=66){
  stop("Error in Age")
}

#Diabetes Duration
if(mean(popM[,"DIAB_DUR"])!=8){
  stop("Error in Diabetes Duration")
}

#Smoking
if(mean(popM[,"SMO"])!=0|mean(popM[,"SMO_0"])!=0){
  stop("Error in Smoking Status")
}

#Ethnicity - check they are white ethnicity
if(mean(popM[,"AFRO"])!=0|mean(popM[,"INDIAN"])!=0){
  stop("Error in Ethnic Status")
}

##HbA1c - check it is 7.5%
if(mean(popM[,"HBA_0"])!=7.5|mean(popM[,"HBA"])!=7.5){
  stop("Error in HbA1c values")
}

##SBP - check it is 145
if(mean(popM[,"SBP_0"])!=145|mean(popM[,"SBP"])!=145){
  stop("Error in SBP")
}

##Diastolic blood pressure is not used in this model

##Total cholesterol is not used in this model

##HDL

if(mean(popM[,"HDL_0"])!=1.3|mean(popM[,"HDL"])!=1.3){
  stop("Error in HDL values")
}

#LDL
if(mean(popM[,"LDL_0"])!=3.0|mean(popM[,"LDL"])!=3.0){
  stop("Error in LDL values")
}

##Triglycerides are not used in this model

##BMI - check it is 28Kg/m2
#Continuous
if(mean(popM[,"BMI_0"])!=28|mean(popM[,"BMI"])!=28){
  stop("Error in BMI values")
}
#Categorical
if(mean(popM[,"BMI_U_18_5"])!=0|mean(popM[,"BMI_O_E_25"])!=1){
  stop("Error in categorical BMI values")
}

#Albumin:creatine ratio is not used in this model

#PVD status
if(mean(popM[,"PVD_E"])!=0|mean(popM[,"PVD_H"])!=0){
  stop("Error in PVD status")
}

#MMALB status
if(mean(popM[,"MMALB_E"])!=0|mean(popM[,"MMALB_H"])!=0){
  stop("Error in MMALB status")
}

#ATFIB status
if(mean(popM[,"ATFIB_E"])!=0|mean(popM[,"ATFIB_H"])!=0){
  stop("Error in ATFIB status")
}

#eGFR
#Continuous
if(mean(popM[,"eGFR"])!=70|mean(popM[,"eGFR_0"])!=70){
  stop("Error in eGFR values")
}
#Knot
if(mean(popM[,"eGFR_U_60"])!=60|mean(popM[,"eGFR_O_60"])!=(70-60)){
  stop("Error in eGFR knot values")
}

##WBC 
if(mean(popM[,"WBC"])!=7|mean(popM[,"WBC_0"])!=7){
  stop("Error in WBC values")
}

##Haemoglobin
if(mean(popM[,"HAEM"])!= 14 | mean(popM[,"HAEM_0"])!=14){
  stop("Error in Haemoglobin values")
}

##Macrovascular disease
if(
  mean(popM[,"MI_E"])!=0|
  mean(popM[,"MI_H"])!=0|
  mean(popM[,"MI2_E"])!=0|
  mean(popM[,"MI2_H"])!=0|
  mean(popM[,"STRO_E"])!=0|
  mean(popM[,"STRO_H"])!=0|
  mean(popM[,"STRO2_E"])!=0|
  mean(popM[,"STRO2_H"])!=0|
  mean(popM[,"IHD_E"])!=0|
  mean(popM[,"IHD_H"])!=0|
  mean(popM[,"CHF_E"])!=0|
  mean(popM[,"CHF_H"])!=0
){
  stop("Error in Macrovascular disease variables")  
}

##Microvascular disease
if(
  mean(popM[,"BLIND_E"])!=0|
  mean(popM[,"BLIND_H"])!=0|
  mean(popM[,"ULCER_E"])!=0|
  mean(popM[,"ULCER_H"])!=0|
  mean(popM[,"RENAL_E"])!=0|
  mean(popM[,"RENAL_H"])!=0|
  mean(popM[,"AMP_E"])!=0|
  mean(popM[,"AMP_H"])!=0|
  mean(popM[,"AMP2_E"])!=0|
  mean(popM[,"AMP2_H"])!=0
){
  stop("Error in Microvascular disease variables")
}

##Other health states
if(
  mean(popM[,"CANB_E"])!=0|
  mean(popM[,"CANB_H"])!=0|
  mean(popM[,"CANC_E"])!=0|
  mean(popM[,"CANC_H"])!=0|
  mean(popM[,"OST_E"])!=0|
  mean(popM[,"OST_H"])!=0|
  mean(popM[,"DEP_E"])!=0
){
  stop("Error in other health state variables")
}

#Make a population that is one person and Female
popF <- build_population_MtHood2025 (10,T, PopulationVariables)

##Get the indices for which the columns are different
diffcols <- which(apply(popF!=popM, 2, all))

##check only these columns from popF
popF <- popF[,diffcols]
#Note this turns the object into a named number and not a matrix so the notation for looking up values does not need a preceding, to tell R to look at all rows, as there are no rows in this object now

if(mean(popF[,"MALE"])!=0 & mean(popF[,"FEMALE"])!=1 &mean(popF[,"MEN"])!=1){
  stop("Error in Gender/Menopause status, female population")

}

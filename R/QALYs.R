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


##'@param population_ is the population matrix
##'@param parameters_ is a single row of the parameters matrix
##'@param year_ model year
##'@param alive_ is a vector of TRUEs and FALSEs to indicate whether or not people
##'are alive
##'@param GlobalVars_ is the global parameters matrix
##'@return population_ is the revised population matrix
calculate_QALYs_default <- function(population_, parameters_,  year_, alive_, GlobalVars_) {
  #Calculate multiplier to adjust for a population with T2D without any complications
  
  #Declare mean BMI, age and proportion female from Hayes et al
  Mean_BMI_Hayes <-28.4 #Mean BMI in the source of our baseline utilities 
  Mean_Age_Hayes <- 65.8
  Mean_pFemale_Hayes <- 4729/(6401+4729)
  
  Util_bl_mult <- parameters_[,"UTIL_BL"] / 
    (0.9454933+0.0256466*(1-Mean_pFemale_Hayes)+
       -0.0002213*Mean_Age_Hayes+
       -0.0000294*(Mean_Age_Hayes^2))
  
  #Calculate utility prior to adjusting for BMI and events 
  population_[,"EQ5D"][alive_] <- (0.9454933 + 
                             0.0256466*population_[,"MALE"][alive_]+
                             -0.0002213 * population_[,"AGE"][alive_] + 
                             -0.0000294 * (population_[,"AGE"][alive_]^2))*Util_bl_mult
  
  #apply the BMI decrements to this
  population_[,"EQ5D"][alive_] <- population_[,"EQ5D"][alive_] + 
    parameters_[,"UTIL_BMI"]*(population_[,"BMI"][alive_]-Mean_BMI_Hayes)
  
  #Ensure that QALYs are not greater than 1 after adjusting for BMI, age and gender
  population_[,"EQ5D"][alive_] <- ifelse(population_[,"EQ5D"][alive_]>1,1,population_[,"EQ5D"][alive_])
  
  #Adjust for model events
  #This is multiplicative, with no other option in line with ISPOR good practice 
  #guidance 
  #If you want decrements you must rewrite this code
  #MI 
  #event years
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"MI_E"][alive_]==1|population_[,"MI2_E"][alive_]==1,
                                                     parameters_[,"UTIL_MI_E"],1)
  #History of a first MI
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"MI_H"][alive_]==1,
                                                     parameters_[,"UTIL_MI_H"],1)
  #History of a second MI
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"MI2_H"][alive_]==1,
                                                     parameters_[,"UTIL_MI_H"],1)
  
  #Stroke
  #event years
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"STRO_E"][alive_]==1|population_[,"STRO2_E"][alive_]==1,
                                                     parameters_[,"UTIL_STRO_E"],1)
  #History of a first MI
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"STRO_H"][alive_]==1,
                                                     parameters_[,"UTIL_STRO_H"],1)
  #History of a second MI
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"STRO2_H"][alive_]==1,
                                                     parameters_[,"UTIL_STRO_H"],1)
  
  #CHF
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"CHF_E"][alive_]==1,
                                                     parameters_[,"UTIL_CHF"],1)
  #History of CHF
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"CHF_H"][alive_]==1,
                                                     parameters_[,"UTIL_CHF"],1)
  
  #IHD
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"IHD_E"][alive_]==1,
                                                     parameters_[,"UTIL_IHD"],1)
  #History of IHD
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"IHD_H"][alive_]==1,
                                                     parameters_[,"UTIL_IHD"],1)
  
  #Blind
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"BLIND_E"][alive_]==1,
                                                                     parameters_[,"UTIL_BLIND"],1)
  #History of blindness
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"BLIND_H"][alive_]==1,
                                                                     parameters_[,"UTIL_BLIND"],1)
  
  #Ulcers
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"ULCER_E"][alive_]==1,
                                                                     parameters_[,"UTIL_ULCER"],1)
  #History of Ulcer
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"ULCER_H"][alive_]==1,
                                                                     parameters_[,"UTIL_ULCER"],1)
  
  #Amputations
  #event year
  #event years
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"AMP_E"][alive_]==1|population_[,"AMP2_E"][alive_]==1,
                                                                     parameters_[,"UTIL_AMP"],1)
  #History of a first MI
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"AMP_H"][alive_]==1,
                                                                     parameters_[,"UTIL_AMP"],1)
  #History of a second MI
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"AMP2_H"][alive_]==1,
                                                                     parameters_[,"UTIL_AMP"],1)
  
  #Renal Failure
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"RENAL_E"][alive_]==1,
                                                                     parameters_[,"UTIL_RENAL"],1)
  #History of Renal Failure
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"RENAL_H"][alive_]==1,
                                                                     parameters_[,"UTIL_RENAL"],1)
  
  #Micro or Macro albuminurea
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"MMALB_E"][alive_]==1,
                                                                     parameters_[,"UTIL_MMALB"],1)
  #History of MMALB
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"MMALB_H"][alive_]==1,
                                                                     parameters_[,"UTIL_MMALB"],1)
  
  #ATFIB
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"ATFIB_E"][alive_]==1,
                                                                     parameters_[,"UTIL_ATFIB"],1)
  #History of ATFIB
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"ATFIB_H"][alive_]==1,
                                                                     parameters_[,"UTIL_ATFIB"],1)
  
  #PVD
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"PVD_E"][alive_]==1,
                                                                     parameters_[,"UTIL_PVD"],1)
  #History of MMALB
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"PVD_H"][alive_]==1,
                                                                     parameters_[,"UTIL_PVD"],1)
  
  #Breast Cancer
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"CANB_E"][alive_]==1,
                                                                     parameters_[,"UTIL_CANB"],1)
  #History of MMALB
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"CANB_H"][alive_]==1,
                                                                     parameters_[,"UTIL_CANB"],1)
  
  #Colorectal Cancer
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"CANC_E"][alive_]==1,
                                                                     parameters_[,"UTIL_CANC"],1)
  #History of MMALB
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"CANC_H"][alive_]==1,
                                                                     parameters_[,"UTIL_CANC"],1)
  
  #Depression
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"DEP_E"][alive_]==1,
                                                                     parameters_[,"UTIL_DEP"],1)
  #History of MMALB
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"DEP_H"][alive_]==1,
                                                                     parameters_[,"UTIL_DEP"],1)
  
  #Osteoarthritis
  #event year
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"OST_E"][alive_]==1,
                                                                     parameters_[,"UTIL_CANC"],1)
  #History of MMALB
  population_[,"EQ5D"][alive_] <-population_[,"EQ5D"][alive_]*ifelse(population_[,"OST_H"][alive_]==1,
                                                                     parameters_[,"UTIL_CANC"],1)
  
  #Half anyone who died this year QALYs, this assumes they died halfway through the year
  dead_thisyear <- is.na(population_[,"F_ALLCAUSE"][alive_])==F
  population_[,"EQ5D"][alive_] <- population_[,"EQ5D"][alive_]*ifelse(dead_thisyear,
                                                                      0.5,
                                                                      1)
  #record QALYs and discounted QALYs cumulatively
  population_[,"QALY"][alive_] <- population_[,"QALY"][alive_] + 
                                  population_[,"EQ5D"][alive_]
  population_[,"DiscQALY"][alive_] <- population_[,"DiscQALY"][alive_] + 
                                      (population_[,"EQ5D"][alive_]/
                                         ((1+as.numeric(GlobalVars_["disc_rate_QALYs", "Value"]))^year_))
  return(population_)
}



##'@param population_ is the population matrix
##'@param parameters_ is a single row of the parameters matrix
##'@param year_ model year
##'@param alive_ is a vector of TRUEs and FALSEs to indicate whether or not people
##'are alive
##'@param GlobalVars_ is the global parameters matrix
##'@return population_ is the revised population matrix
calculate_QALYs_MtHood2025_C1 <- function(population_,  year_, alive_, GlobalVars_) {
  #Warning, so this code is not accidentitally used in an model run for an actual project.
  warning("This model run uses the function for MtHood 2025 to calculate utilities/QALYs. This uses fixed values & uses utility decrements against good practice recommendations from ISPOR 2019 (10.1016/j.jval.2019.01.004) and against the prefered approach NICE methods guide (point 4.3.7, https://www.nice.org.uk/process/pmg36/resources/nice-health-technology-evaluations-the-manual-pdf-72286779244741, accessed 16th May 2025). This should only be used when running MtHood challenge simulations")
  #Calculate each person's utility this year
  
  if(GlobalVars_["Mt Hood Utility Values", "Value"]== "95% CI low"){
    #Constant
    population_[,"EQ5D"][alive_] <- 0.681
    #BMI, per unit of BMI above 25 Kg/m2
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + -0.008*(population_[,"BMI"][alive_]-25)
    
    #Retinopathy
    dec_ret <- min(-0.124*(population_[,"BLIND_E"][alive_]),
                   -0.124*(population_[,"BLIND_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_ret
    
    #Nephropathy
    #Renal Failure event year weights for hemodialysis, peritoneal dialysis and transplant come from https://www.ukkidney.org/sites/default/files/UK%20Renal%20Registry%20Annual%20Report%202022%20Patient%20Summary.pdf (accessed 16th May 2025), page 3
    #Renal Failure historical event weights  hemodialysis, peritoneal dialysis and transplant come from https://www.ukkidney.org/sites/default/files/UK%20Renal%20Registry%20Annual%20Report%202022%20Patient%20Summary.pdf (accessed 16th May 2025), page 4
    dec_neph <- min(-0.091*(population_[,"MMALB_E"][alive_]),
                    -0.091*(population_[,"MMALB_H"][alive_]),
                    (-0.137*(540/(540+39+6127+1548))+
                       -0.274*((39+6127)/(540+39+6127+1548))+
                       -0.342*(1548/(540+39+6127+1548)))*(population_[,"RENAL_E"][alive_]),
                    (-0.137*(39874/(39874+1452+25825+3800))+
                       -0.274*((1452+25825)/(39874+1452+25825+3800))+
                       -0.342*(3800/(39874+1452+25825+3800)))*(population_[,"RENAL_H"][alive_]))
    #Apply the biggest reduction in utility for nephropathy
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_neph
    
    rm(dec_neph)
    
    #Neuropathy
    #This model only tracks an ulcer state, so only apply in year 1
    dec_neuro <- min(-0.090*(population_[,"PVD_E"][alive_]),
                     -0.090*(population_[,"PVD_H"][alive_]),
                     -0.207*(population_[,"ULCER_E"][alive_]),
                     -0.389*(population_[,"AMP_E"][alive_]),
                     -0.389*(population_[,"AMP2_E"][alive_]))
    
    #Apply the biggest reduction in utility for nephropathy
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_neuro
    rm(dec_neuro)
    
    
    #Stroke
    dec_cebro <- min(-0.222*(population_[,"STRO_E"][alive_]),
                     -0.222*(population_[,"STRO_H"][alive_]),
                     -0.222*(population_[,"STRO2_E"][alive_]),
                     -0.222*(population_[,"STRO2_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_cebro
    rm(dec_cebro)
    
    #Coronary heart disease
    
    dec_CHD <- min( -0.067*(population_[,"MI_E"][alive_]),
                    -0.067*(population_[,"MI_H"][alive_]),
                    -0.067*(population_[,"MI2_E"][alive_]),
                    -0.067*(population_[,"MI2_H"][alive_]),
                    -0.169*(population_[,"CHF_E"][alive_]),
                    -0.169*(population_[,"CHF_H"][alive_]),
                    -0.126*(population_[,"IHD_E"][alive_]),
                    -0.126*(population_[,"IHD_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_CHD
    rm(dec_CHD) 
  }
  else if (GlobalVars_["Mt Hood Utility Values", "Value"]== "95% CI high"){
    #Constant
    population_[,"EQ5D"][alive_] <- 0.889
    #BMI, per unit of BMI above 25 Kg/m2
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + -0.004*(population_[,"BMI"][alive_]-25)
    
    #Retinopathy
    dec_ret <- min(-0.025*(population_[,"BLIND_E"][alive_]),
                   -0.025*(population_[,"BLIND_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_ret
    
    #Nephropathy
    #Renal Failure event year weights for hemodialysis, peritoneal dialysis and transplant come from https://www.ukkidney.org/sites/default/files/UK%20Renal%20Registry%20Annual%20Report%202022%20Patient%20Summary.pdf (accessed 16th May 2025), page 3
    #Renal Failure historical event weights  hemodialysis, peritoneal dialysis and transplant come from https://www.ukkidney.org/sites/default/files/UK%20Renal%20Registry%20Annual%20Report%202022%20Patient%20Summary.pdf (accessed 16th May 2025), page 4
    dec_neph <- min(-0.005*(population_[,"MMALB_E"][alive_]),
                    -0.005*(population_[,"MMALB_H"][alive_]),
                    (-0.027*(540/(540+39+6127+1548))+
                       -0.054*((39+6127)/(540+39+6127+1548))+
                       -0.066*(1548/(540+39+6127+1548)))*(population_[,"RENAL_E"][alive_]),
                    (-0.027*(39874/(39874+1452+25825+3800))+
                       -0.054*((1452+25825)/(39874+1452+25825+3800))+
                       -0.066*(3800/(39874+1452+25825+3800)))*(population_[,"RENAL_H"][alive_]))
    #Apply the biggest reduction in utility for nephropathy
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_neph
    
    rm(dec_neph)
    
    #Neuropathy
    #This model only tracks an ulcer state, so only apply in year 1
    dec_neuro <- min(-0.032*(population_[,"PVD_E"][alive_]),
                     -0.032*(population_[,"PVD_H"][alive_]),
                     -0.133*(population_[,"ULCER_E"][alive_]),
                     -0.170*(population_[,"AMP_E"][alive_]),
                     -0.170*(population_[,"AMP2_E"][alive_]))
    
    #Apply the biggest reduction in utility for nephropathy
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_neuro
    rm(dec_neuro)
    
    
    #Stroke
    dec_cebro <- min(-0.105*(population_[,"STRO_E"][alive_]),
                     -0.105*(population_[,"STRO_H"][alive_]),
                     -0.105*(population_[,"STRO2_E"][alive_]),
                     -0.105*(population_[,"STRO2_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_cebro
    rm(dec_cebro)
    
    #Coronary heart disease
    
    dec_CHD <- min( -0.042*(population_[,"MI_E"][alive_]),
                    -0.042*(population_[,"MI_H"][alive_]),
                    -0.042*(population_[,"MI2_E"][alive_]),
                    -0.042*(population_[,"MI2_H"][alive_]),
                    -0.048*(population_[,"CHF_E"][alive_]),
                    -0.048*(population_[,"CHF_H"][alive_]),
                    -0.054*(population_[,"IHD_E"][alive_]),
                    -0.054*(population_[,"IHD_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_CHD
    rm(dec_CHD) 
  }
  #Mean Values - default if either the 95% CI low or high are not specified
  else{
    #Constant
    population_[,"EQ5D"][alive_] <- 0.785
    #BMI, per unit of BMI above 25 Kg/m2
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + -0.006*(population_[,"BMI"][alive_]-25)
    
    #Retinopathy
    dec_ret <- min(-0.074*(population_[,"BLIND_E"][alive_]),
                   -0.074*(population_[,"BLIND_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_ret

    #Nephropathy
    #Renal Failure event year weights for hemodialysis, peritoneal dialysis and transplant come from https://www.ukkidney.org/sites/default/files/UK%20Renal%20Registry%20Annual%20Report%202022%20Patient%20Summary.pdf (accessed 16th May 2025), page 3
    #Renal Failure historical event weights  hemodialysis, peritoneal dialysis and transplant come from https://www.ukkidney.org/sites/default/files/UK%20Renal%20Registry%20Annual%20Report%202022%20Patient%20Summary.pdf (accessed 16th May 2025), page 4
    dec_neph <- min(-0.048*(population_[,"MMALB_E"][alive_]),
               -0.048*(population_[,"MMALB_H"][alive_]),
               (-0.082*(540/(540+39+6127+1548))+
                  -0.164*((39+6127)/(540+39+6127+1548))+
                  -0.204*(1548/(540+39+6127+1548)))*(population_[,"RENAL_E"][alive_]),
               (-0.082*(39874/(39874+1452+25825+3800))+
                  -0.164*((1452+25825)/(39874+1452+25825+3800))+
                  -0.204*(3800/(39874+1452+25825+3800)))*(population_[,"RENAL_H"][alive_]))
    #Apply the biggest reduction in utility for nephropathy
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_neph
    
    rm(dec_neph)
    
    #Neuropathy
    #This model only tracks an ulcer state, so only apply in year 1
    dec_neuro <- min(-0.061*(population_[,"PVD_E"][alive_]),
               -0.061*(population_[,"PVD_H"][alive_]),
               -0.170*(population_[,"ULCER_E"][alive_]),
               -0.280*(population_[,"AMP_E"][alive_]),
               -0.280*(population_[,"AMP2_E"][alive_]))
    
    #Apply the biggest reduction in utility for nephropathy
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_neuro
    rm(dec_neuro)
    
    
    #Stroke
    dec_cebro <- min(-0.164*(population_[,"STRO_E"][alive_]),
                     -0.164*(population_[,"STRO_H"][alive_]),
                     -0.164*(population_[,"STRO2_E"][alive_]),
                     -0.164*(population_[,"STRO2_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_cebro
    rm(dec_cebro)
    
    #Coronary heart disease
    
    dec_CHD <- min( -0.055*(population_[,"MI_E"][alive_]),
                    -0.055*(population_[,"MI_H"][alive_]),
                    -0.055*(population_[,"MI2_E"][alive_]),
                    -0.055*(population_[,"MI2_H"][alive_]),
                    -0.108*(population_[,"CHF_E"][alive_]),
                    -0.108*(population_[,"CHF_H"][alive_]),
                    -0.090*(population_[,"IHD_E"][alive_]),
                    -0.090*(population_[,"IHD_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_CHD
    rm(dec_CHD) 
  }
  
  #Calculate QALYs
  population_[,"QALY"][alive_] <- population_[,"QALY"][alive_] + 
    population_[,"EQ5D"][alive_]
  population_[,"DiscQALY"][alive_] <- population_[,"DiscQALY"][alive_] + 
    (population_[,"EQ5D"][alive_]/
       ((1+as.numeric(GlobalVars_["disc_rate_QALYs", "Value"]))^year_))
  
  
  return(population_)
  }

calculate_QALYs_MtHood2025_C2 <- function(population_,  year_, alive_, GlobalVars_) {
  #Warning, so this code is not accidentitally used in an model run for an actual project.
  warning("This model run uses the function for MtHood 2025 to calculate utilities/QALYs. This uses fixed values & uses utility decrements against good practice recommendations from ISPOR 2019 (10.1016/j.jval.2019.01.004) and against the prefered approach NICE methods guide (point 4.3.7, https://www.nice.org.uk/process/pmg36/resources/nice-health-technology-evaluations-the-manual-pdf-72286779244741, accessed 16th May 2025). This should only be used when running MtHood challenge simulations")
  #Calculate each person's utility this year
 
    #Constant
    population_[,"EQ5D"][alive_] <- 0.807
    
    #Retinopathy
    dec_ret <- min(-0.000*(population_[,"BLIND_E"][alive_]),
                   -0.000*(population_[,"BLIND_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_ret
    
    #Nephropathy
    #Renal Failure event year weights for hemodialysis, peritoneal dialysis and transplant come from https://www.ukkidney.org/sites/default/files/UK%20Renal%20Registry%20Annual%20Report%202022%20Patient%20Summary.pdf (accessed 16th May 2025), page 3
    #Renal Failure historical event weights  hemodialysis, peritoneal dialysis and transplant come from https://www.ukkidney.org/sites/default/files/UK%20Renal%20Registry%20Annual%20Report%202022%20Patient%20Summary.pdf (accessed 16th May 2025), page 4
    dec_neph <- min(-0.330*(population_[,"RENAL_E"][alive_]),
                    -0.330*(population_[,"RENAL_H"][alive_]))
    #Apply the biggest reduction in utility for nephropathy
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_neph
    
    rm(dec_neph)
    
    #Neuropathy
    #This model only tracks an ulcer state, so only apply in year 1
    dec_neuro <- min(-0.210*(population_[,"ULCER_E"][alive_]),
                     -0.210*(population_[,"ULCER_H"][alive_])
                     -0.172*(population_[,"AMP_E"][alive_]),
                     -0.172*(population_[,"AMP2_E"][alive_]))
    
    #Apply the biggest reduction in utility for nephropathy
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_neuro
    rm(dec_neuro)
    
    #Coronary heart disease
    
    dec_CHD <- min( -0.065*(population_[,"MI_E"][alive_]),
                    -0.000*(population_[,"MI_H"][alive_]),
                    -0.065*(population_[,"MI2_E"][alive_]),
                    -0.000*(population_[,"MI2_H"][alive_]),
                    -0.101*(population_[,"CHF_E"][alive_]),
                    -0.101*(population_[,"CHF_H"][alive_]),
                    -0.000*(population_[,"IHD_E"][alive_]),
                    -0.000*(population_[,"IHD_H"][alive_]),
                    -0.165*(population_[,"STRO_E"][alive_]),
                    -0.165*(population_[,"STRO_H"][alive_]),
                    -0.165*(population_[,"STRO2_E"][alive_]),
                    -0.165*(population_[,"STRO2_H"][alive_]))
    
    population_[,"EQ5D"][alive_] <-  population_[,"EQ5D"][alive_] + dec_CHD
    rm(dec_CHD) 
  
  #Calculate QALYs
  population_[,"QALY"][alive_] <- population_[,"QALY"][alive_] + 
    population_[,"EQ5D"][alive_]
  population_[,"DiscQALY"][alive_] <- population_[,"DiscQALY"][alive_] + 
    (population_[,"EQ5D"][alive_]/
       ((1+as.numeric(GlobalVars_["disc_rate_QALYs", "Value"]))^year_))
  
  
  return(population_)
}





##'@param population_ is the population matrix
##'@param parameters_ is a single row of the parameters matrix
##'@param year_ model year
##'@param alive_ is a vector of TRUEs and FALSEs to indicate whether or not people
##'are alive
##'@param GlobalVars_ is the global parameters matrix
##'@return population_ is the revised population matrix
calculate_QALYs <- function(population_, parameters_,  year_, alive_, GlobalVars_){
  #If the option is selected use the Mt Hood Challenge 1 utility function
  if(GlobalVars_["Mt Hood Utilities", "Value"]=="C1"){
  population_ <- calculate_QALYs_MtHood2025C1(population_,  year_, alive_, GlobalVars_)
  }else if (GlobalVars_["Mt Hood Utilities", "Value"]=="C2"){  #If the option is selected use the Mt Hood Challenge 2 utility function
  population_ <- calculate_QALYs_MtHood2025C2(population_,  year_, alive_, GlobalVars_)
  } else{#Otherwise use the default utilities function
  population_ <- calculate_QALYs_default(population_, parameters_,  year_, alive_, GlobalVars_)
  }
  return(population_)
}
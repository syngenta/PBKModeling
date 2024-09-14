library(httk)
###############################################################################
#########                  calculation of Kw and Ky                 ###########
###############################################################################
# Kw #
#### function calculating partitioning coefficien; Schmitt 2008
calc_partition <- function(pKa_Donor, pKa_Accept, LogP,
                           Fcell, Fint, FW, FL, FP, Fn_L, Fn_PL, Fa_PL, 
                           Funbound.plasma, tissue.pH, plasma.pH){
  
  if (pKa_Donor != 'None' | pKa_Accept != 'None'){
    ionization        <- calc_ionization(pKa_Donor = pKa_Donor,pKa_Accept = pKa_Accept, pH = tissue.pH)
    fraction_neutral  <- ionization[["fraction_neutral"]]
    fraction_charged  <- ionization[["fraction_charged"]]
    fraction_negative <- ionization[["fraction_negative"]]
    fraction_positive <- ionization[["fraction_positive"]]
    fraction_zwitter  <- ionization[["fraction_zwitter"]]
    
    plasma                  <- calc_ionization(pKa_Donor = pKa_Donor,pKa_Accept =pKa_Accept, pH = plasma.pH)
    fraction_neutral_plasma <- plasma[["fraction_neutral"]]
    fraction_zwitter_plasma <- plasma[["fraction_zwitter"]]
    fraction_charged_plasma <- plasma[["fraction_charged"]]
    
  }else{
    fraction_neutral  <- 1
    fraction_charged  <- 0
    fraction_negative <- 0
    fraction_positive <- 0
    fraction_zwitter  <- 0
    
    fraction_neutral_plasma <- 1
    fraction_zwitter_plasma <- 0
    fraction_charged_plasma <- 0
  }
  
  Pow                 <- 10^(LogP)
  
  alpha               <- 0.001                                                          # Schmitt 2008
  Fprotein.plasma     <- 0.0383                                                         # Table 1 Scanes et al 2022 Quantiative comparison
  SD_Fprotein.plasma  <- 0.052/100                                                      # Table 1 Scanes et al 2022 Quantiative comparison
  FWpl                <- 1 - Fprotein.plasma                                            # the water fractions in plasma
  FWint               <- FWpl                                                           # the water fractions in interstitium
  Kn_L                <- Pow * (fraction_neutral + fraction_zwitter + 
                                  alpha * fraction_charged)                             # K_(n_L)=D_(o:w)=P_(o:w) eq 13
  Kn_PL               <- 10^(1.294 + 0.304 * log10(Pow))                                # membrane affinity
  Ka_PL               <- Kn_PL * (fraction_neutral + fraction_zwitter + 
                                    20 * fraction_positive + 0.05 * fraction_negative)  # eq 17 & 18
  KP                  <- 0.163 + 0.0221 * Kn_PL                                         # Proteins, eq 19 & 9
  KAPPAcell2pu        <- (fraction_neutral_plasma + fraction_zwitter_plasma + 
                            alpha * fraction_charged_plasma)/(fraction_neutral + 
                                                                fraction_zwitter + alpha * fraction_charged)
  
  Kint       <- (FWint + 0.37 * (1/Funbound.plasma - FWpl))                             # eq 4  
  Kcell      <- (FW + Kn_L * Fn_L + Kn_PL * Fn_PL + Ka_PL * Fa_PL + KP * FP)            # eq 7
  
  Ktissue2pu <- Fint * Kint + KAPPAcell2pu * Fcell * Kcell                              # eq 12
  
  return(Ktissue2pu)
}


#################################################################################################################
calc_unionized_fraction  <- function(pKa_Donor, pKa_Accept, pH){
  
  ion               <- calc_ionization(pKa_Donor = pKa_Donor, pKa_Accept = pKa_Accept, pH = pH)
  fraction_neutral  <- ion[["fraction_neutral"]]
  fraction_charged  <- ion[["fraction_charged"]]
  
  return(fraction_neutral)
}


#################################################################################################################
calc_ky      <- function(pKa, pKb, fub){
  
  alpha_6             <- calc_unionized_fraction(pKa_Donor = pKa, pKa_Accept = pKb, pH = 6)
  alpha_7.4           <- calc_unionized_fraction(pKa_Donor = pKa, pKa_Accept = pKb, pH = 7.4)
  
  # less than 10% CV for protein fraction in plasma
  Fprotein.plasma     <- 0.0383                         # fraction of protein in hen plasma; Table 1 Scanes et al, 2022
  Funbound.plasma     <- fub
  FWpl                <- 1 - Fprotein.plasma            # the water fractions in plasma
  #FWint               <- FWpl                           # the water fractions in interstitium (egg white)
  FWint               <- FW_y                           # fraction of water within egg yolk
  
  ky                  <- fub * (FWint + Fpr_y/Fprotein.plasma  * (1/Funbound.plasma - FWpl)) * alpha_7.4 / alpha_6    # eq4 in scheferlie, hekman 2016, eq4 in Schimitt 2008
} # eq 10 in Hekman 2016


calc_kw      <- function(pKa, pKb, fub){
  
  alpha_7.6           <- calc_unionized_fraction(pKa_Donor = pKa, pKa_Accept = pKb, pH = 7.6)
  alpha_7.4           <- calc_unionized_fraction(pKa_Donor = pKa, pKa_Accept = pKb, pH = 7.4)
  
  # less than 10% CV for protein fraction in plasma
  Fprotein.plasma     <- 0.0383                         # fraction of protein in hen plasma; Table 1 Scanes et al, 2022
  Funbound.plasma     <- fub
  FWpl                <- 1 - Fprotein.plasma            # the water fractions in plasma
  #FWint               <- FWpl                          # the water fractions in interstitium (egg white)
  FWint               <- FW_w                           # fraction of water within egg white
  
  kw                  <- fub * (FWint + Fpr_w/Fprotein.plasma  * (1/Funbound.plasma - FWpl)) * alpha_7.4 / alpha_7.6    # eq4 in scheferlie, hekman 2016, eq4 in Schimitt 2008
}# eq 7 in Hekman 2016


calc_kw_Kdw      <- function(pKa, pKb, fub){
  
  alpha_7.6           <- calc_unionized_fraction(pKa_Donor = pKa, pKa_Accept = pKb, pH = 7.6)
  alpha_7.4           <- calc_unionized_fraction(pKa_Donor = pKa, pKa_Accept = pKb, pH = 7.4)
  
  # less than 10% CV for protein fraction in plasma
  Fprotein.plasma     <- 0.0383                         # fraction of protein in hen plasma; Table 1 Scanes et al, 2022
  Funbound.plasma     <- fub
  FWpl                <- 1 - Fprotein.plasma            # the water fractions in plasma
  #FWint               <- FWpl                          # the water fractions in interstitium (egg white)
  FWint               <- FW_w                           # fraction of water within egg white
  
  fuw                 <- 1/(FWint + Fpr_w/Fprotein.plasma  * (1/Funbound.plasma - FWpl))
  kw                  <- fub / fuw * alpha_7.4 / alpha_7.6    # eq4 in scheferlie, hekman 2016, eq4 in Schimitt 2008
  kw                  <- kw * (1 - exp(-1 * 0.26 * (fub * alpha_7.4 + fuw * alpha_7.6)))
}


##################################################################################
########                calculation of Rblood2plasma                      ########
##################################################################################
calc_Rblood2plasma_chicken  <- function(LogP, Compound_type, fub, pKa, Hematocrit_fraction){
  
  # Rodgers and Rowland 2006
  pHIW        <- 7.22             # pH of intracellular of blood cells
  pHEW        <- 7.4              # pH of extracellular water 
  plasma.pH   <- 7.49             # pH of plasma for chicken Hugo Chiodi, and James W. Terman; Arterial blood gases of the domestic hen
  Hematocrit_fraction   <- Hematocrit_fraction 
  
  if(Compound_type == 'acid'){
    X  <- 1 + 10^(pHIW - pKa)
    Y  <- 1 + 10^(plasma.pH - pKa)
  }else if(Compound_type == 'base'){
    X  <- 1 + 10^(pKa - pHIW) 
    Y  <- 1 + 10^(pKa - plasma.pH)
  }else if(Compound_type == 'neutral'){
    X  <- 1
    Y  <- 1
  }else if(Compound_type == 'zwitterions'){
    X  <- 1 + 10^(pKa_BASE - pHIW) + 10^(pHIW - pKa_ACID)
    Y  <- 1 + 10^(pKa_BASE - plasma.pH)  + 10^(plasma.pH  - pKa_ACID)
  }
  
  fracIWBC      <- 0.603             # Table 1, Rodgers and Rowland 2005 using rat for laying hen, blood cells
  fracNLBC      <- 0.0017            # Table 1, Rodgers and Rowland 2005
  fracNPBC      <- 0.0029            # Table 1, Rodgers and Rowland 2005
  ChemicalP     <- 10^(LogP)
  Pow           <- 10^(LogP)
  
  Krbc2pu       <- X / Y * fracIWBC + ((ChemicalP * fracNLBC + (0.3 * ChemicalP + 0.7) * fracNPBC) / Y)   # eq B2; Rodgers and Rowland 2006
  Rblood2plasma <- 1 - Hematocrit_fraction/100 + Hematocrit_fraction/100 * Krbc2pu * fub                  # eq B1; Rodgers and Rowland 2006
  
  return(Rblood2plasma)
}


#######################################################################
###                     partition coefficients                      ###  to unbound plasma; 12 for rat, 14 for human
#######################################################################
####### calculation of tissue partition coefficient (Parent)  ##################
PC_parent            <- predict_partitioning_schmitt(chem.name  = Compound_name_parent,
                                                     chem.cas= CAS_parent,
                                                     species= 'Rat',
                                                     regression = F, 
                                                     adjusted.Funbound.plasma = F) 
Kgut2pu_parent        <- PC_parent$Kgut2pu      
Kkidney2pu_parent     <- PC_parent$Kkidney2pu  
Kliver2pu_parent      <- PC_parent$Kliver2pu     
Klung2pu_parent       <- PC_parent$Klung2pu   
Kmuscle2pu_parent     <- PC_parent$Kmuscle2pu 
Kadipose2pu_parent    <- PC_parent$Kadipose2pu 
Krep2pu_parent        <- Kgut2pu_parent 
Krest2pu_parent       <- (PC_parent$Kbone2pu * Vbone + PC_parent$Kbrain2pu * Vbrain + PC_parent$Kheart2pu * Vheart +
                            PC_parent$Kspleen2pu * Vspleen + PC_parent$Krest2pu * Vrest_a) / (Vbone + Vbrain + Vheart + Vspleen + Vrest_a)

####### calculation of tissue partition coefficient (daughter)  ##################
PC_daughter            <- predict_partitioning_schmitt(chem.name= Compound_name_daughter,
                                                       chem.cas= CAS_daughter,
                                                       species= 'Rat',
                                                       regression = F, 
                                                       adjusted.Funbound.plasma = F) 
Kgut2pu_daughter        <- PC_daughter$Kgut2pu      
Kkidney2pu_daughter     <- PC_daughter$Kkidney2pu  
Kliver2pu_daughter      <- PC_daughter$Kliver2pu      
Klung2pu_daughter       <- PC_daughter$Klung2pu   
Kmuscle2pu_daughter     <- PC_daughter$Kmuscle2pu
Kadipose2pu_daughter    <- PC_daughter$Kadipose2pu
Krep2pu_daughter        <- Kgut2pu_daughter 
Krest2pu_daughter       <- (PC_daughter$Kbone2pu * Vbone + PC_daughter$Kbrain2pu * Vbrain + PC_daughter$Kheart2pu * Vheart +
                              PC_daughter$Kspleen2pu * Vspleen + PC_daughter$Krest2pu * Vrest_a) / (Vbone + Vbrain + Vheart + Vspleen + Vrest_a)

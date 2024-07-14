
###################################################################################################################
#             Metabolism: non-specific binding in metabolism  (microsomal binding & hepatocyte)                   #
###################################################################################################################

calc_hep_mic_fu        <- function(pH, media_conc_metabolic, media_conc_binding,fuinc_binding,
                                   Compound_type, LogP, pKa){
  # nonspecific microsomal binding at 0.5 mg/mL based on regression model (only used if original fuinc is not known)
  # Kmic            <- 10^(0.072*(LogP)^2+0.067*(LogP)-1.126)
  # fuinc_1         <- 1 / (Kmic+1)'
  # liver density: 1.05 g/mL (Pearce et al. 2017)
  # Vmax is in the unit of umol/h/mg protein; Km is in the unit of um (umol/L); C_tissue is in the unit of umol/L.
  
  if(tolower(Compound_type) == 'acid'){
    LogP       <- LogP + log10(1/(1+10^(pH-pKa)))           # logD calculation for acid
  }else{ 
    LogP       <- LogP
  }                                     # logP for bases and neutrals 
  
  if(tolower(fuinc_binding) != 'none'){
    fuinc           <- 1 / ( media_conc_metabolic / media_conc_binding * ((1 - fuinc_binding) / fuinc_binding) + 1)     
  }else{
    Kmic            <- 10^(0.072*(LogP)^2+0.067*(LogP)-1.126)
    fuinc           <- 1 / ( 125*0.005*Kmic+1)              # calc_hep_fu(); Ratio of cell volume to incubation volume. Default (0.005) 
  }
  
  return(fuinc)
}

####################################################################
###                    clearance (L/d/kg BW)                     ###  
####################################################################
# hepatic and intestinal clearance for parent; hepatic clearance for daughter (rat only)
calc_metabolic_clearance  <- function(type_clearance, incubation, fuinc, Vmax_unit,       # if media conc not available, use 1 instead for both
                                      Vmax, Km, Clint_ori, 
                                      C_tissue, Ktissue2pu, tissue_specific_volume, 
                                      MPPG, fub, BW_scaledfrom){
  
  if(tolower(type_clearance) == 'liver'){
    
    if(is.numeric(Vmax)){
      if(Vmax_unit == 'umol/h/kg bw'){
        if(incubation == 'scaled'){
          Vmax_scaledfrom   <- Vmax * BW_scaledfrom * 24                                                           # umol/h; volume per unit time
          Vmax_scaled       <- Vmax_scaledfrom * ((BW / BW_scaledfrom)^0.75) /BW                               # umol/h/kg bw; 0.67 or 0.75
          Clint             <- Vmax_scaled /(Km * fuinc + C_tissue / Ktissue2pu) 
        }else{
          Clint        <- Vmax /(Km * fuinc + C_tissue / Ktissue2pu) * 24
        }
      }else{
        # Vmax is in the unit of umol/h/mg protein or umol/h/million cells; Km is in the unit of um (umol/L); C_tissue is in the unit of umol/L.
        Clint          <- Vmax /(Km * fuinc + C_tissue / Ktissue2pu) * MPPG * (tissue_specific_volume * 1.05 * 1000) * 24  # unit: L/h/kg BW
      }
    }else if(is.numeric(Clint_ori)){
      if(incubation == 'scaled'){
        # # Clint_scaled from is in the unit of L/d/kg BW
        Clint_scaledfrom   <- Clint_ori * BW_scaledfrom
        Clint              <- Clint_scaledfrom * ((BW / BW_scaledfrom)^0.75) /BW       
      }else{
        # Clint_ori is in the unit of uL/h/million cell
        Clint          <- Clint_ori / fuinc * 1e-6 * MPPG * tissue_specific_volume * 1.05 * 1000 * 24               # L/h/million cells * million cells/g tissue * g/kg  -->  L/h/kg BW
      }
    }
    
  }else if(tolower(type_clearance) == 'intestine'){
    if(is.numeric(Vmax)){
      if(Vmax_unit == 'umol/h/kg bw'){
        if(incubation == 'scaled'){
          Vmax_scaledfrom   <- Vmax * BW_scaledfrom * 24                                                           # umol/h; volume per unit time
          Vmax_scaled       <- Vmax_scaledfrom * ((BW / BW_scaledfrom)^0.75) /BW                               # umol/h/kg bw; 0.67 or 0.75
          Clint             <- Vmax_scaled /(Km * fuinc + C_tissue / (Ktissue2pu * fub))
        }else{
          Clint        <- Vmax /(Km * fuinc + C_tissue / (Ktissue2pu * fub)) * 24
        }
      }else{
        # Vmax is in the unit of umol/h/mg protein or umol/h/million cells; Km is in the unit of um (umol/L); C_tissue is in the unit of umol/L.
        Clint          <- Vmax /(Km * fuinc + C_tissue / (Ktissue2pu * fub)) * MPPG * (tissue_specific_volume * 1.05 * 1000)  * 24               # unit: L/d/kg BW
      }
    }else if(is.numeric(Clint_ori)){
      if(incubation == 'scaled'){
        # # Clint_scaled from is in the unit of L/h/kg BW
        Clint_scaledfrom   <- Clint_ori * BW_scaledfrom
        Clint              <- Clint_scaledfrom * ((BW / BW_scaledfrom)^0.75) /BW   
      }else if(incubation == 'half-life') {
        # Clint_ori is in the unit of h (half-life)
        # Clint           <- log(2) * 0.2 / 1000 / Clint_ori                                                   # case by case; check the unit; equation A14 in appenxix
        Clint           <-  log(2) / Clint_ori                                                                 # h-1         
      }else{
        # Clint_ori is in the unit of uL/h/million cell or uL/h/mg protein
        Clint          <- Clint_ori / fuinc * 1e-6 * MPPG * tissue_specific_volume * 1.05 * 1000 * 24              # L/h/million cells * million cells/g tissue * g/kg  -->  L/h/kg BW
      }
    }
  }
  
  return(Clint)                                                                                          # unit: (L/h/kg BW)
}

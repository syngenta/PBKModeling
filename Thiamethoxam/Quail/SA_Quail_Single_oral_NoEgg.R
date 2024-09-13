library(stationaRy)
library(isdparser)
library(tidyr)
library(lubridate)
library(WriteXLS)
library(readxl)
library(writexl)
library("tidyr")
library(dplyr)
library(Metrics)
library(PKNCA)
library(cowplot)
library(knitr)
library(reshape2)
library(httk)
library(deSolve)
library(pksensi)
library(ggplot2)
library(httk)
library(reshape2)


rm(list = ls())
SPECIES  <- 'Quail'

current_dir <- getwd()
print(current_dir)
source(file.path("HTTK", "Bird", "General code", "Avian_Noegg.R"))
source(file.path("HTTK", "Bird", "General code", paste(SPECIES, "_PhyData.R", sep = '')))
source(file.path("HTTK", "Bird", "General code", 'Metabolism.R'))

#======================================================================================================
#                          Compound physico-chemical & in vitro properties                            #
#======================================================================================================
####################### Parent (Table 1)    TMX ######################
# physico-chemical properties  
species        <- 'Rat'                  # (either "Rat", "Rabbit", "Dog", "Mouse", or "Human").
# Table 1, Gasthuys, Elke 2019 Comparative physiology of glomerular filtration rate by plasma clearance of exogenous creatinine and exo-iohexol in six different avian species
BW             <- 0.186                  # kg, body weight of 18 weeks 
BW_hen         <- 1.45                   # kg
SDBW           <- 0.08

# TMX
MW_parent             <- 291.72
CAS_parent            <- '153719-23-4a'
S_parent              <- 4100                    # mg/L; http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/631.htm
pKa_a_parent          <- 'NA'                    # No dissociation; http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/631.htm
pKa_b_parent          <- 'NA'                    # pka accept (base)
fub_parent            <- 0.55                    # fraction unbound in plasma of chicken; from SED
LogP_parent           <- -0.13                   # LogP FROM SED
Peff_parent           <- mean(c(8.71E-6, 19.21E-6, 16.4E-6))   # Cacao permeability (cm/s) from SED for human
Compound_type_parent  <- 'neutral'
Compound_name_parent  <- 'Thiamethoxam (TMX)'    # has to be different with HTTK
Density.parent        <- 1.57                    # g/cm3 # http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/631.htm
Rblood2plasma_parent  <- 1.342 
type_hep_clearance    <- 'liver'

Vmax_unit_parent          <- 'none'               # 'umol/h/kg bw' or 'none'
incubation_hep_parent     <- 'scaled'             
Vmax_hep_parent           <- "None"
Km_hep_parent             <- 'None'
Clint_ori_hep_parent      <- 259.3772 
parent_hep_conc_metabolic        <- 0.5                                   # 0.5 mg/ml protein; Shen et al. 2013
fuinc_hep_parent                 <- 9999                                  # binding was not measured in in vitro metabolsim assays
parent_hep_conc_binding          <- 9999


######################        Clothianidin        ##########################
MW_daughter             <- 249.7
CAS_daughter            <- '210880-92-5a'
S_daughter              <- 340                    # mg/L; http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/631.htm
pKa_a_daughter          <- 11.1                   # No dissociation; http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/171.htm; pKa_Donor="pka.a") Compound H dissociation equilibirum constant
pKa_b_daughter          <- 14-11.1                # pka accept (base) Compound H association equilibirum constant
fub_daughter            <- 0.71                   # fraction unbound in plasma of rat; from SED
LogP_daughter           <- 0.702                  # LogP httk
Peff_daughter           <- 19.8E-6                # Cacao permeability (cm/s) from SED for human
Compound_type_daughter  <- 'base'
Compound_name_daughter  <- 'Clothianidin (CLO)'
Density.daughter        <- 1.61                   # g/cm3 # http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/631.htm
Rblood2plasma_daughter  <-'None' 

Vmax_unit_daughter          <- 'none'                 # 'umol/h/kg bw' or 'none'
incubation_hep_daughter     <- 'scaled'        
Vmax_hep_daughter           <- "None"
Km_hep_daughter             <- 'None'
Clint_ori_hep_daughter      <- 41.00266 
daughter_hep_conc_metabolic        <- 0.5                                   # 0.5 mg/ml protein; Shen et al. 2013
fuinc_hep_daughter                 <- 9999                                  # binding was not measured in in vitro metabolsim assays
daughter_hep_conc_binding          <- 9999


##########################################################
###                    absorption                      ###  
##########################################################
# Parent
Peff_human        <- 0.4926 * log10(Peff_parent * 1e6) - 0.1454         # Peff is in the unit of e-6 cm/s; 1e-4 cm/s
Peff_rat          <- (Peff_human - 0.03)/3.6                            # 1e-4 cm/s
ka                <-  2 * Peff_rat / R_rat / 1e4 * 3600 * 24            # d-1
fa                <-  0.9                                               # absorption fraction;

####################################################################
###                    clearance (L/d/kg BW)                     ###  
####################################################################
# hepatic and intestinal clearance for parent; hepatic clearance for daughter (rat only)
liver.density           <-  1.05

##################################################################
###                    excretion    L/d/kg BW                  ###  
##################################################################
# Table 3, Scanes el al 2022 Quantitative morphometric 
# Abstract, Gasthuys, Elke 2019 Comparative physiology of glomerular filtration rate by plasma clearance of exogenous creatinine and exo-iohexol in six different avian species
Qgfr_parent       <-  2.57 /1000 * 60 * 24                      
Qgfr_daughter     <-  2.57 /1000 * 60 * 24  

#=====================     End of user input    =======================#


#===============================================================================
#                      Input to HTTK for parameterization                      #
#===============================================================================
# Add to Httk database
# Reference: add_chemtable (Httk manual dated Sep. 22, 2022, page 7)
# Add acibenzolar and acibenzolar acid data to HTTK chemical table for further compound-specific parameterization using HTTK
# Add to Httk database
# Reference: add_chemtable (Httk manual dated Sep. 22, 2022, page 7)
# Add acibenzolar and acibenzolar acid data to HTTK chemical table for further compound-specific parameterization using HTTK
my.new.data           <- as.data.frame(cbind(c(Compound_name_parent, Compound_name_daughter), 
                                             c(CAS_parent, CAS_daughter), 
                                             c("DTX_parent","DTX_daughter"), 
                                             c(MW_parent, MW_daughter), 
                                             c(LogP_parent, LogP_daughter), 
                                             c(fub_parent, fub_daughter), 
                                             c(0, 0), 
                                             c(pKa_a_parent, pKa_a_daughter), 
                                             c(pKa_b_parent, pKa_b_daughter)))

colnames(my.new.data)     <- c("Name","CASRN","DTXSID","MW","LogP","Fup","CLint", "pKa.a", 'pKa.b')
chem.physical_and_invitro.data <- add_chemtable(my.new.data,
                                                current.table=
                                                  chem.physical_and_invitro.data,
                                                data.list=list(
                                                  Compound  ="Name",
                                                  CAS       ="CASRN",
                                                  DTXSID    ="DTXSID",
                                                  MW        ="MW",
                                                  logP      ="LogP",
                                                  Funbound.plasma="Fup",
                                                  Clint     ="CLint",
                                                  pKa_Donor = 'pKa.a',
                                                  pKa_Accept= 'pKa.b'),
                                                species     = species,
                                                reference   ="SED")

source("C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/General Code/Partition.R")



#######################################################################
###                     partition coefficients                      ###  to unbound plasma; 12 for rat, 14 for human
#######################################################################
## predicting the tissue to unbound plasma partition coefficients for the tissues contained in the tissue.data table 
# predicting partitioning cofficient for different tissues, rightnow assume using adjusted Funbound 
pHIW        <- 7.22             # pH of intracellular of blood cells
pHEW        <- 7.4              # pH of extracellular water 
plasma.pH   <- 7.49             # pH of plasma for chicken Hugo Chiodi, and James W. Terman; Arterial blood gases of the domestic hen

Rblood2plasma_parent     <- if(is.numeric(Rblood2plasma_parent)){Rblood2plasma_parent}else{calc_Rblood2plasma_chicken(LogP_parent, Compound_type_parent, fub_parent, pKa_a_parent, Hematocrit_fraction )}
Rblood2plasma_daughter   <- if(is.numeric(Rblood2plasma_daughter)){Rblood2plasma_daughter}else{calc_Rblood2plasma_chicken(LogP_daughter, Compound_type_daughter, fub_daughter, pKa_a_daughter, Hematocrit_fraction )}

ky_parent            <-  calc_ky(pKa_a_parent,   pKa_b_parent,   fub_parent) 
ky_daughter          <-  calc_ky(pKa_a_daughter, pKa_b_daughter, fub_daughter) 
kw_parent            <-  calc_kw_Kdw(pKa_a_parent,   pKa_b_parent,   fub_parent) 
kw_daughter          <-  calc_kw_Kdw(pKa_a_daughter, pKa_b_daughter, fub_daughter) 

#=========================================================================================
#                                 PBTK model equations                                   #
#=========================================================================================

########################################################################################################################
#                                                  PBPK MODEL RUN                                                      #
########################################################################################################################
initState <- c(Agutlumen_parent  = 0, 
               Agut_parent       = 0,  
               Aliver_parent     = 0,
               Akidney_parent    = 0,
               Aadipose_parent = 0,
               Amuscle_parent    = 0,
               Arest_parent      = 0,
               Arep_parent       = 0,
               Aven_parent        = 0,
               Alung_parent       = 0,
               Aart_parent        = 0,
               Aurine_parent      = 0,
               AUC_Cplasma_parent = 0,
               AUC_Cblood_parent  = 0,
               
               Agut_daughter       = 0,  
               Aliver_daughter     = 0,
               Akidney_daughter    = 0,
               Aadipose_daughter   = 0,
               Amuscle_daughter    = 0,
               Arest_daughter      = 0,
               Arep_daughter       = 0,
               Aven_daughter       = 0,
               Alung_daughter      = 0,
               Aart_daughter       = 0,
               Aurine_daughter     = 0,
               AUC_Cplasma_daughter = 0,
               AUC_Cblood_daughter  = 0,
               
               Agut_parent_in = 0, 
               Aliver_parent_out = 0,
               Aven_parent_out   = 0,
               Aint_parent_out   = 0,
               Aart_parent_out   = 0,
               
               Aurine_daughter_out = 0,
               Aliver_daughter_out = 0)



#########################    parameters for exposure scenario ###############
factor               <- 1
fraction_daughter    <- 1

parms <- c( ka                   = 12.1847922, 
            Clint_parent         = 692.6002544, 
            Clint_daughter       = 79.2898978, 
            
            Rblood2plasma_parent   = 0.9027282,  
            Rblood2plasma_daughter = 1.3597193,   

            fub_parent           = fub_parent,   
            fub_daughter         = fub_daughter, 

            # ratio
            Kgut2pu_parent       = Kgut2pu_parent   ,      
            Kkidney2pu_parent    = Kkidney2pu_parent ,   
            Kliver2pu_parent     = Kliver2pu_parent,     
            Klung2pu_parent      = Klung2pu_parent  ,      
            Krest2pu_parent      = Krest2pu_parent  , 
            Kadipose2pu_parent   = Kadipose2pu_parent,
            Kmuscle2pu_parent    = Kmuscle2pu_parent, 
            Krep2pu_parent       = Krep2pu_parent   ,
            
            Kgut2pu_daughter       = Kgut2pu_daughter   ,      
            Kkidney2pu_daughter    = Kkidney2pu_daughter,   
            Kliver2pu_daughter     = Kliver2pu_daughter,    
            Klung2pu_daughter      = Klung2pu_daughter  ,      
            Krest2pu_daughter      = Krest2pu_daughter  , 
            Kadipose2pu_daughter   = Kadipose2pu_daughter,
            Kmuscle2pu_daughter    = Kmuscle2pu_daughter, 
            Krep2pu_daughter       = Krep2pu_daughter   ,
            
            # Parametes for flow 
            Qart          = Qart, 
            Qgfr_parent   = 6.7386862,  
            Qgfr_daughter = 15.3186752 , 
            
            Qgut          = Qgut,
            Qkidney       = Qkidney,
            Qliver        = Qliver,
            Qadipose      = Qadipose,
            Qmuscle       = Qmuscle,
            Qrep          = Qrep,
            
            # volume 
            Vart          = Vart,
            Vgut          = Vgut,
            Vkidney       = Vkidney,
            Vliver        = Vliver,
            Vadipose      = Vadipose,
            Vmuscle       = Vmuscle,
            Vrep          = Vrep,
            Vlung         = Vlung,
            Vrest         = Vrest,
            Vven          = Vven)



Pred  <- function(pars, Qrest.new){
  
  parms_rest <- c( fa                   = fa,  
                   fraction_daughter    = 1,
                   Qrest                = Qrest.new)
  parms      <- c(pars,parms_rest)
  
  # Modeling duration
  days          <- 2
  Time_max      <- days                                   
  StartTime     <- 0                                        # Time_d
  StopTime      <- Time_max                                 
  dt            <- 0.01
  Times         <- seq(StartTime, StopTime, dt)
  # dose
  dose          <- 100                                      # mg/kg bw
  Oral_input    <- dose * 1000 / MW_parent                  # umol/kg/d, STUDY 2

  Dose_events    <- data.frame(var    = "Agutlumen_parent",
                               time   =  0,
                               value  = Oral_input,
                               method = "add")   
  
  df  <- ode(y       = initState, 
             times   = Times, 
             func    = pbtk5cpt_repeated, 
             parms   = parms,
             atol    = 1e-6, rtol = 1e-8,
             events  = list(data = Dose_events),
             method  = 'lsoda')
  
  df               <- as.data.frame(df)
  colnames(df)[1]  <- "Time"
  outdf  <- cbind.data.frame(Time = df$Time,
                             AUC  = df$AUC_Cblood_parent + df$AUC_Cblood_daughter,
                             Cmax = max(df$C_blood_parent + df$C_blood_daughter))
  
  
  return (outdf)
}

NSC.AUC              <- matrix(nrow=42,ncol=1)
NSC.Cmax             <- matrix(nrow=42,ncol=1)
rownames(NSC.AUC)    <- names(parms)
rownames(NSC.Cmax)   <- names(parms)

fold                 <- 1.1

for (i in 1: 42){
  
  Qrest            <- parms[['Qart']] - (parms[['Qliver']] + parms[['Qkidney']] + parms[['Qmuscle']] + parms[['Qadipose']] + parms[['Qgut']] + parms[['Qrep']])
  pars.changed     <- parms[i]  * fold 
  pars.rest        <- parms[-i]
  pars             <- c(pars.changed, pars.rest)    
  Qrest.new        <- pars[['Qart']] - (pars[['Qliver']] + pars[['Qkidney']] + pars[['Qmuscle']] + pars[['Qadipose']] + pars[['Qgut']] + pars[['Qrep']])
  
  cat('i = ', i,  ', Pars: ',  names(parms[i]), ', Original value: ', parms[i], ', Changed: ',  pars.changed, '\n')
  
  delta.AUC        <- (Pred(pars, Qrest.new)$AUC  - Pred(parms, Qrest)$AUC)  / Pred(parms, Qrest)$AUC  / (fold-1) * 100
  delta.Cmax       <- (Pred(pars, Qrest.new)$Cmax - Pred(parms, Qrest)$Cmax) / Pred(parms, Qrest)$Cmax / (fold-1) * 100

  NSC.AUC[i,1]        <- tail(delta.AUC, n = 1)
  NSC.Cmax[i,1]       <- tail(delta.Cmax, n = 1)
  
  colnames(NSC.AUC)[1]   <- 'AUC'
  colnames(NSC.Cmax)[1]  <- 'Cmax'
  
}

NSC               <- cbind(NSC.AUC, NSC.Cmax)
file_path         <- file.path(current_dir, "HTTK", "Bird", 'TMX', 'Quail', 'Sensitivity.csv')
write.csv(NSC, file_path , row.names = TRUE)
#=============================================================================
############                End of PBPK model                    #############
#=============================================================================

##########################      pbpk output     #############################  
# Mass balance
# it is difficult to evaluate becase state variable related to egg white and yolk were set to zero routinely.
ggplot() +
  geom_line (data = df, aes(Time, AUC), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (d)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14)

ggplot() +
  geom_line (data = df, aes(Time, Cmax), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (d)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14)

ggplot() +
  geom_line (data = df, aes(Time, C_blood_parent+C_blood_daughter), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (d)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14)



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
#dev.off()

#==================================================================================================
#                                      Species:Mallard duck                                       # 
#                                     No of compartment: 5                                        #
#==================================================================================================
# Reference
# PKNCA: https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk: https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk

#=============================================================================================================
#                                   Generic Maternal PBPK Model Parameters                                   #
#                                  Tissue and physiological data for species                                 #
#=============================================================================================================

#======================================================================================================
#                          Compound physico-chemical & in vitro properties                            #
#======================================================================================================
####################### Parent (Table 1)    TMX ######################
# physico-chemical properties  
SPECIES        <- 'Duck'
species        <- 'Rat'                 # (either "Rat", "Rabbit", "Dog", "Mouse", or "Human").
# Table 1, Gasthuys, Elke 2019 Comparative physiology of glomerular filtration rate by plasma clearance of exogenous creatinine and exo-iohexol in six different avian species
BW             <- 1.05                     # 
BW_hen         <- 1.45                   # kg
BW_human       <- 70                    # kg
SDBW           <- 0.08

source("C:/XXX/HTTK/Bird/General Code/Avian_Noegg.R")
source(paste("C:/XXX/OneDrive - Syngenta/HTTK/Bird/General Code/", SPECIES, "_PhyData.R", sep = ''))
source("C:/XXX/OneDrive - Syngenta/HTTK/Bird/General Code/Metabolism.R")


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
BW_scaledfrom         <- BW_hen
type_hep_clearance    <- 'liver'

Vmax_unit_parent          <- 'none'                   # 'umol/h/kg bw' or 'none'
incubation_hep_parent     <- 'scaled'             # human
Vmax_hep_parent           <- "None"
Km_hep_parent             <- 'None'
Clint_ori_hep_parent      <- 258.3149 
parent_hep_conc_metabolic        <- 0.5                                   # 0.5 mg/ml protein; Shen et al. 2013
fuinc_hep_parent                 <- 9999                                  # binding was not measured in in vitro metabolsim assays
parent_hep_conc_binding          <- 9999





######################        Clothianidin        ##########################
MW_daughter             <- 249.7
CAS_daughter            <- '210880-92-5a'
S_daughter              <- 340                    # mg/L; http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/631.htm
pKa_a_daughter          <- 11.1                   # No dissociation; http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/171.htm; pKa_Donor="pka.a") Compound H dissociation equilibirum constant
pKa_b_daughter          <- 14-11.1                  # pka accept (base) Compound H association equilibirum constant
fub_daughter            <- 0.71                   # fraction unbound in plasma of rat; from SED
LogP_daughter           <- 0.70                   # LogP httk
Peff_daughter           <- 19.8E-6                # Cacao permeability (cm/s) from SED for human
Compound_type_daughter  <- 'base'
Compound_name_daughter  <- 'Clothianidin (CLO)'
Density.daughter        <- 1.61                   # g/cm3 # http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/631.htm
Rblood2plasma_daughter  <-'None' 

Vmax_unit_daughter          <- 'none'                 # 'umol/h/kg bw' or 'none'
incubation_hep_daughter     <- 'scaled'           # human
Vmax_hep_daughter           <- "None"
Km_hep_daughter             <- 'None'
Clint_ori_hep_daughter      <- 40.83473 
daughter_hep_conc_metabolic        <- 0.5                                   # 0.5 mg/ml protein; Shen et al. 2013
fuinc_hep_daughter                 <- 9999                                  # binding was not measured in in vitro metabolsim assays
daughter_hep_conc_binding          <- 9999


#===============================================================================
#                      Input to HTTK for parameterization                      #
#===============================================================================
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

source("C:/xxx/OneDrive - Syngenta/HTTK/Bird/General Code/Partition.R")

## Permeability-surface area product L/d kg
PSc_adipose_parent      <- Peff_parent * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000            # Peff * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000  # cm/s -> cm/d * g/kg * cm2/g -> cm3/d kg -> L/d kg
PSc_adipose_daughter    <- Peff_daughter * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000          # Peff * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000  # cm/s -> cm/d * g/kg * cm2/g -> cm3/d kg -> L/d kg


##########################################################
###                    absorption                      ###  
##########################################################
# Parent
# calculate of ka; use rat ka first 
# R_rat             <- 0.2             # 0.2cm # radius of rat jejunum
Peff_human        <- 10^(0.4926 * log10(Peff_parent * 1e6) - 0.1454)         # Peff is in the unit of e-6 cm/s; 1e-4 cm/s
#Peff_human        <- (0.4926 * log10(Peff_parent * 1e6) - 0.1454)  
Peff_rat          <- (Peff_human - 0.03)/3.6                            # 1e-4 cm/s
ka                <-  12.1847922 #2 * Peff_rat / R_rat / 1e4 * 3600 * 24            # d-1
#ka                <- 2 * Peff_rat / R_rat / 1e4 * 3600 * 24  

#R_hen         <-  10/10/2                                     # radius of chicken jejunum,cm; Table 6 Kamil, Dariusz  2020
#ka_broiler    <-  22                                          # 13.03
#ka            <-  ka_broiler * (BW / BW_broiler)^(-0.25)      # 2 * Peff / R_hen * 3600 * 24                # h-1   ---> d-1
# Mean retention time (MRT) in the small intestine (3h), Birger Svihus KhaledItani 2019
# Table 10, Wang et al, 2020 chicken and turkey; not used
#factor        <- 1
#kt_broiler    <- log(0.98/0.02)/3 * 24 * factor                      # Transit rate constant (1/h) 1.3  ---> 1/d
#kt            <-  kt_broiler * ((BW / BW_broiler)^(-0.25))

fa        <- 0.9 #249.7 / 291.72                                   # absorption fraction;

####################################################################
###                    clearance (L/d/kg BW)                     ###  
####################################################################
# hepatic and intestinal clearance for parent; hepatic clearance for daughter (rat only)
## Clearance
# hepatic, intestinal and plasma clearance rate L/h/BW for parent compound
fu_hep_parent                    <- calc_hep_mic_fu(pH                   = pH,
                                                    media_conc_metabolic = if(parent_hep_conc_metabolic == 9999){parent_hep_conc_metabolic ='none'}else{parent_hep_conc_metabolic},  #parent_hep_conc_metabolic,
                                                    media_conc_binding   = if(parent_hep_conc_binding == 9999){parent_hep_conc_binding ='none'}else{parent_hep_conc_binding},  #parent_hep_conc_binding, 
                                                    fuinc_binding        = if(fuinc_hep_parent == 9999){fuinc_hep_parent ='none'}else{fuinc_hep_parent},  #fuinc_hep_parent, 
                                                    Compound_type        = Compound_type_parent, 
                                                    LogP                 = LogP_parent, 
                                                    pKa                  = pKa_a_parent)

fu_hep_daughter                  <- calc_hep_mic_fu(pH                   = pH,
                                                    media_conc_metabolic = if(daughter_hep_conc_metabolic == 9999){daughter_hep_conc_metabolic ='none'}else{daughter_hep_conc_metabolic},  #daughter_hep_conc_metabolic,
                                                    media_conc_binding   = if(daughter_hep_conc_binding == 9999){daughter_hep_conc_binding ='none'}else{daughter_hep_conc_binding},  #daughter_hep_conc_binding, 
                                                    fuinc_binding        = if(fuinc_hep_daughter == 9999){fuinc_hep_daughter ='none'}else{fuinc_hep_daughter},  #fuinc_hep_daughter, 
                                                    Compound_type        = Compound_type_daughter, 
                                                    LogP                 = LogP_daughter, 
                                                    pKa                  = pKa_a_daughter)


# clearance
Clint_hep_parent                 <- calc_metabolic_clearance(type_clearance       = type_hep_clearance,
                                                             incubation           = incubation_hep_parent, 
                                                             fuinc                = 1, #fu_hep_parent,
                                                             Vmax_unit            = Vmax_unit_parent,
                                                             Vmax =  Vmax_hep_parent , Km = Km_hep_parent, 
                                                             Clint_ori            = if(Clint_ori_hep_parent==9999){Clint_ori_hep_parent = 'none'}else{Clint_ori_hep_parent},
                                                             C_tissue             = 1, Ktissue2pu = 1,
                                                             tissue_specific_volume = Vliver, 
                                                             MPPG = MPPGL_chicken , fub = fub_parent, BW_scaledfrom = BW_scaledfrom) 


Clint_hep_daughter               <- calc_metabolic_clearance(type_clearance       = type_hep_clearance,
                                                             incubation           = incubation_hep_daughter, 
                                                             fuinc                = 1, #fu_hep_daughter,
                                                             Vmax_unit            = Vmax_unit_daughter,
                                                             Vmax =  Vmax_hep_daughter , Km = Km_hep_daughter, 
                                                             Clint_ori            = if(Clint_ori_hep_daughter==9999){Clint_ori_hep_daughter = 'none'}else{Clint_ori_hep_daughter},
                                                             C_tissue             = 1, Ktissue2pu = 1,
                                                             tissue_specific_volume = Vliver, 
                                                             MPPG = MPPGL_chicken , fub = fub_daughter, BW_scaledfrom = BW_scaledfrom) 


##################################################################
###                    excretion    L/d/kg BW                  ###  
##################################################################
# Table 3, Scanes el al 2022 Quantitative morphometric 
# Abstract, Gasthuys, Elke 2019 Comparative physiology of glomerular filtration rate by plasma clearance of exogenous creatinine and exo-iohexol in six different avian species
factor              <-  1
Qgfr_hen_parent     <-  5.18112  * BW_hen                                                 ## L/d/kg BW
Qgfr_hen_daughter   <-  21.83472 * BW_hen                                                 # 
Qgfr_parent         <-  Qgfr_hen_parent   * ((BW / BW_hen)^0.75) /BW                       ## 2.6 / 1000 * 60  * 24 * factor   L/d/kg BW L for laying hen and duck
Qgfr_daughter       <-  Qgfr_hen_daughter * ((BW / BW_hen)^0.75) /BW  

#=====================     End of user input    =======================#


#######################################################################
###                     partition coefficients                      ###  to unbound plasma; 12 for rat, 14 for human
#######################################################################
## predicting the tissue to unbound plasma partition coefficients for the tissues contained in the tissue.data table 
# predicting partitioning cofficient for different tissues, rightnow assume using adjusted Funbound 
Rblood2plasma_parent     <- if(is.numeric(Rblood2plasma_parent)){Rblood2plasma_parent}else{calc_Rblood2plasma_chicken(LogP_parent, Compound_type_parent, fub_parent, pKa_a_parent, Hematocrit_fraction )}
Rblood2plasma_daughter   <- if(is.numeric(Rblood2plasma_daughter)){Rblood2plasma_daughter}else{calc_Rblood2plasma_chicken(LogP_daughter, Compound_type_daughter, fub_daughter, pKa_a_daughter, Hematocrit_fraction )}

######            calculation of Kw and Ky                 ######
# Yolk & white
# ky_parent            <-  calc_partition(pKa_Donor = pKa_a_parent, pKa_Accept = pKa_b_parent, LogP_parent,
#                                         Fcell_y, Fint_y, FW_y, FL_y, Fpr_y, Fn_L_y, Fn_PL_y, Fa_PL_y, 
#                                         fub_parent, tissue.pH = pH_y, plasma.pH = plasma.pH ) * fub_parent          # rate constant (unitless) for egg yolk
# 
# ky_daughter          <-  calc_partition(pKa_Donor = pKa_a_daughter, pKa_Accept = pKa_b_daughter, LogP_daughter,
#                                         Fcell_y, Fint_y, FW_y, FL_y, Fpr_y, Fn_L_y, Fn_PL_y, Fa_PL_y, 
#                                         fub_daughter, tissue.pH = pH_y, plasma.pH = plasma.pH ) * fub_daughter      # rate constant (unitless) for egg yolk

#kw_parent            <-  calc_kw(pKa_a_parent,   pKa_b_parent,   fub_parent)                                        # partition to egg white
#kw_daughter          <-  calc_kw(pKa_a_daughter, pKa_b_daughter, fub_daughter)
ky_parent            <-  calc_ky(pKa_a_parent,   pKa_b_parent,   fub_parent) 
ky_daughter          <-  calc_ky(pKa_a_daughter, pKa_b_daughter, fub_daughter) 
kw_parent            <-  calc_kw_Kdw(pKa_a_parent,   pKa_b_parent,   fub_parent) 
kw_daughter          <-  calc_kw_Kdw(pKa_a_daughter, pKa_b_daughter, fub_daughter) 
  
#=========================================================================================
#                                 PBTK model equations                                   #
#=========================================================================================

########################################################################
###                    define modeling variables                     ###  
########################################################################
days          <- 28
Time_max      <- days                                     # 5 d  

StartTime     <- 0                                        # Time_d
StopTime      <- Time_max                                 
dt            <- 0.01
Times         <- seq(StartTime, StopTime, dt)
#==============================================
dose_ppm      <- 5200   # 5200; 300
#==============================================                               # ug/kg/d

var                <- rep("Agutlumen_parent", 24 * days/1)                                   # dose pattern every one hour
time               <- seq(0, (days), by = 1/24)
time               <- time[-length(time)]
value_1            <-  c(0,
                         0.002022362,
                         0.004044725,
                         0.002022362,
                         0.026560358,
                         0.047265253,
                         0.048256832,
                         0.048752621,
                         0.053875778,
                         0.059164198,
                         0.063461039,
                         0.062799987,
                         0.059494724,
                         0.059494724,
                         0.064948407,
                         0.063791565,
                         0.073376826,
                         0.089076823,
                         0.087919981,
                         0.069245248,
                         0.009302867,
                         0.005123318,
                         0,
                         0)   

dose_5200ppm      <- c(23, 43, 32, 34, 38)                                    # group mean food consumption: g food/bird/d ; ppm: ug/g food; for sub-acute study
dose_300ppm       <- c(2197/7/2, 2197/7/2, 2197/7/2, 2197/7/2, 2197/7/2)      # g food/cage/week  -> g/bird/d * ug/g -> ug/bird/d

# for tox studies
var                <- rep("Agutlumen_parent", 24 * 5/1)                       # dose pattern every one hour
time               <- seq(0, (5), by = 1/24)
time               <- time[-length(time)]
dose_5d            <- matrix(nrow = 5, ncol = length(value_1))

for (i in 1:5){
  cat("day = ", i, '\n')
  if(dose_ppm == 300) {dose   <- dose_300ppm}
  if(dose_ppm == 5200){dose   <- dose_5200ppm}
  
  dose_5d[i,]   <- value_1 * dose_ppm * dose[i] / MW_parent / 0.15           # FOR 5200 ppm subacute study; 0.15 0.188 is the average body weight of the sub-acute study for the 5 days
  #dose_5d[i,]   <- value_1 * dose_ppm * dose[i] / MW_parent / BW              # umol/kg bw; for chronic study
}
df_dose_5d  <- as.vector(t(dose_5d))



Dose_events    <- data.frame(var    = rep("Agutlumen_parent", 5 * 24),
                             time   =  time,
                             value  = df_dose_5d,
                             method = "add")
                    

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
               Aadipose_daughter = 0,
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
factor <- 1

parms <- c( ka            = ka,                 # h-1
            fa            = fa,  
            Clint_parent         = Clint_hep_parent ,
            Clint_daughter       = Clint_hep_daughter ,
            
            Rblood2plasma_parent   = Rblood2plasma_parent,  
            Rblood2plasma_daughter = Rblood2plasma_daughter,   
            BW = BW,              
            
            fub_parent           = fub_parent,   
            fub_daughter         = fub_daughter, 
            fraction_daughter    = 1,
            
            # ratio
            Kgut2pu_parent       = Kgut2pu_parent    * factor,      
            Kkidney2pu_parent    = Kkidney2pu_parent * factor,   
            Kliver2pu_parent     = Kliver2pu_parent  * factor,     
            Klung2pu_parent      = Klung2pu_parent   * factor,      
            Krest2pu_parent      = Krest2pu_parent   * factor, 
            Kadipose2pu_parent   = Kadipose2pu_parent * factor,
            Kmuscle2pu_parent    = Kmuscle2pu_parent * factor,
            Krep2pu_parent       = Krep2pu_parent    * factor,
            
            Kgut2pu_daughter       = Kgut2pu_daughter    * factor,      
            Kkidney2pu_daughter    = Kkidney2pu_daughter * factor,   
            Kliver2pu_daughter     = Kliver2pu_daughter  * factor,     
            Klung2pu_daughter      = Klung2pu_daughter   * factor,      
            Krest2pu_daughter      = Krest2pu_daughter   * factor, 
            Kadipose2pu_daughter   = Kadipose2pu_daughter * factor,
            Kmuscle2pu_daughter    = Kmuscle2pu_daughter * factor,
            Krep2pu_daughter       = Krep2pu_daughter    * factor,

            # Parametes for flow 
            Qart          = Qart, 
            Qgfr_parent   = Qgfr_parent ,   
            Qgfr_daughter = Qgfr_daughter ,
            
            Qgut          = Qgut,
            Qkidney       = Qkidney,
            Qliver        = Qliver,
            Qrest         = Qrest, 
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
            Vven          = Vven,
            
            V_Vadipose    = V_Vadipose ,
            V_Tadipose    = V_Tadipose ,
            PSc_adipose_parent   = PSc_adipose_parent,
            PSc_adipose_daughter = PSc_adipose_daughter)



df  <- ode(y       = initState, 
           times   = Times, 
           func    = pbtk5cpt_repeated, 
           parms   = parms,
           atol    = 1e-6, rtol = 1e-8,
           events  = list(data = Dose_events),
           method  = 'lsoda')

df               <- as.data.frame(df)
colnames(df)[1]  <- "Time"



#=============================================================================
############                End of PBPK model                    #############
#=============================================================================

##########################      pbpk output     #############################  
######           Plot           ###########
# Mass balance
ggplot() +
  geom_line (data = df, aes(Time, Mass_parent_bal), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14) + ylim(-0.01, 0.01)

ggplot() +
  geom_line (data = df, aes(Time, Mass_daughter_bal), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14) + ylim(-0.01, 0.01)

ggplot() +
  geom_line (data = df, aes(Time, Mass_bal), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14) + ylim(-0.01, 0.01)

ggplot() +
  geom_line(data = df, aes(Time, C_blood_parent), col="#00AFBB", lwd=1) +
  theme(text = element_text(size = 20)) + xlim(0, 6)

ggplot() +
  geom_line(data = df, aes(Time, C_blood_daughter), col="#00AFBB", lwd=1) +
  theme(text = element_text(size = 20)) + xlim(0, 6) 

#==========================================================================================================
###                                   Across-species extrapolation                                      ###  
#==========================================================================================================

##################################        Toxicity endpoint derivation       ################################
# AUC of blood (acute, single oral dose)
AUC24_parent     <- df[df$Time == 2,]$AUC_Cblood_parent 
AUC24_parent                                                   # 10.83266 umol*L/d
AUC24_daughter   <- df[df$Time == 2,]$AUC_Cblood_daughter 
AUC24_daughter                                                 # 94.12444 umol*L/d
AUC24_total      <- AUC24_parent + AUC24_daughter              # 104.9571 umol*L/d

Cmax_parent      <- max(df$C_blood_parent)            # 96.76751
Cmax_daughter    <- max(df$C_blood_daughter)          # 678.6988
Cmax_total       <- max(df$C_blood_total)             # 765.3447

stat  <- data.frame('Name'  = c('AUC24_parent', 'AUC24_daughter', 'AUC24_total', 'Cmax_parent', 'Cmax_daughter', 'Cmax_total' ),
                    'Value' = c(AUC24_parent, AUC24_daughter, AUC24_total, Cmax_parent, Cmax_daughter, Cmax_total))
stat
'           Name     Value
1   AUC24_parent  47.20760
2 AUC24_daughter 117.54598
3    AUC24_total 164.75358
4    Cmax_parent  69.63179
5  Cmax_daughter 152.39238
6     Cmax_total 219.18220'

# chronic
'           Name    Value
1   AUC24_parent 2.320941
2 AUC24_daughter 5.784932
3    AUC24_total 8.105873
4    Cmax_parent 2.624991
5  Cmax_daughter 5.744910
6     Cmax_total 8.262767'

# AUC TMX
AUC24_TMX_blood   <- df[df$Time == 1,]$AUC_Cblood_parent 
AUC24_TMX_blood
# 1.435053 ug*L/h
AUC24_TMX_plasma   <- df[df$Time == 1,]$AUC_Cplasma_parent 
AUC24_TMX_plasma
# 1.847343 ug*L/h
AUC24_TMX_plasma_unbound   <- df[df$Time == 1,]$AUC_plasma_parent_unbound 
AUC24_TMX_plasma_unbound
# 1.016039 ug*L/h

'# AUC TMX
> AUC24_TMX_blood   <- df[df$Time == 1,]$AUC_Cblood_parent 
> AUC24_TMX_blood
[1] 2.052923
> # 1.435053 ug*L/h
> AUC24_TMX_plasma   <- df[df$Time == 1,]$AUC_Cplasma_parent 
> AUC24_TMX_plasma
[1] 2.642727
> # 1.847343 ug*L/h
> AUC24_TMX_plasma_unbound   <- df[df$Time == 1,]$AUC_plasma_parent_unbound 
> AUC24_TMX_plasma_unbound
[1] 1.4535
> # 1.016039 ug*L/h'

# AUC CLO
AUC24_CLO_blood   <- df[df$Time == 1,]$AUC_Cblood_daughter
AUC24_CLO_blood 
# 14.24027 ug*L/h
AUC24_CLO_plasma  <- df[df$Time == 1,]$AUC_Cplasma_daughter
AUC24_CLO_plasma 
# 15.27864 ug*L/h
AUC24_CLO_plasma_unbound   <- df[df$Time == 1,]$AUC_plasma_daughter_unbound 
AUC24_CLO_plasma_unbound
# 10.84784 ug*L/h


Cmax_Tmx_blood    <- max(df$C_blood_parent)
Cmax_Tmx_blood  
# 12.05475
Cmax_Tmx_plasma    <- max(df$C_plasma_parent)
Cmax_Tmx_plasma  
# 15.51807
Cmax_Tmx_plasma_unbound    <- max(df$C_plasma_parent_unbound)
Cmax_Tmx_plasma_unbound  
# 8.534939

Cmax_CLO_blood    <- max(df$C_blood_daughter)
Cmax_CLO_blood 
# 80.59109
Cmax_CLO_plasma    <- max(df$C_plasma_daughter)
Cmax_CLO_plasma
# 86.46763
Cmax_CLO_plasma_unbound    <- max(df$C_plasma_daughter_unbound)
Cmax_CLO_plasma_unbound  
# 61.39202

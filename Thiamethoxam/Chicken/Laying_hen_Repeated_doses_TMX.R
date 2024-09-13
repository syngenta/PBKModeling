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

rm(list = ls())

#==================================================================================================
#                                      Species: Laying Hen                                        # 
#==================================================================================================
# Referencedata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
# PKNCA: https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk: https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk

SPECIES  <- 'Laying_hen'
BW             <- 1.45                  # kg, body weight of 18 weeks female laying hen

current_dir <- getwd()
print(current_dir)

source(file.path("HTTK", "Bird", "General code", "Avian_Wegg.R"))
source(file.path("HTTK", "Bird", "General code", paste(SPECIES, "_PhyData.R", sep = '')))
source(file.path("HTTK", "Bird", "General code", 'Metabolism.R'))

df_tissue_human                     <- tissue.data[which(tissue.data$Species == 'Human'),]
Vliver_human                        <- subset(df_tissue_human , variable == "Vol (L/kg)" &
                                              tolower(Tissue) == 'liver')$value
#======================================================================================================
#                          Compound physico-chemical & in vitro properties                            #
#======================================================================================================
####################### Parent (Table 1)    TMX ######################
# physico-chemical properties  
species        <- 'Rat'                 # (either "Rat", "Rabbit", "Dog", "Mouse", or "Human").
# Table 1, Gasthuys, Elke 2019 Comparative physiology of glomerular filtration rate by plasma clearance of exogenous creatinine and exo-iohexol in six different avian species

SDBW           <- 0.08
BW_broiler     <- 2.5                   # kg
BW_mouse       <- 0.02
BW_rat         <- 0.25
BW_human       <- 70

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
BW_scaledfrom         <- BW_human
type_hep_clearance    <- 'liver'

Vmax_unit_parent          <- 'none'                   # 'umol/h/kg bw' or 'none'
incubation_hep_parent     <- 'scaled'             # human
Vmax_hep_parent           <- "None"
Km_hep_parent             <- 'None'
Clint_human_parent        <- 0.6115 / 1E6 * 60 * HPGL_human * (Vliver_human * 1.05 * 1000) * 24     # ul/min/million hepatocytes; human; HTTK database--> ul/h/million hepatocytes -> L/h/kg
Clint_ori_hep_parent      <- Clint_human_parent 
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
LogP_daughter           <- 0.70                   # LogP httk
Peff_daughter           <- 19.8E-6                # Cacao permeability (cm/s) from SED for human
Compound_type_daughter  <- 'base'
Compound_name_daughter  <- 'Clothianidin (CLO)'
Density.daughter        <- 1.61                   # g/cm3 # http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/631.htm
Rblood2plasma_daughter  <-'None' 

Vmax_unit_daughter          <- 'none'             # 'umol/h/kg bw' or 'none'
incubation_hep_daughter     <- 'scaled'           # human
Vmax_hep_daughter           <- "None"
Km_hep_daughter             <- 'None'
Clint_human_daughter        <- 10.73 / 1E6 * 60 * HPGL_human * (Vliver_human * 1.05 * 1000)  * 24             # ul/min/million hepatocytes; human; HTTK database--> ul/h/million hepatocytes;
Clint_ori_hep_daughter      <- Clint_human_daughter 
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

source(file.path("HTTK", "Bird", "General code", "Partition.R"))

## Permeability-surface area product L/d kg; not used
PSc_adipose_parent      <- Peff_parent * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000            # Peff * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000  # cm/s -> cm/d * g/kg * cm2/g -> cm3/d kg -> L/d kg
PSc_adipose_daughter    <- Peff_daughter * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000          # Peff * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000  # cm/s -> cm/d * g/kg * cm2/g -> cm3/d kg -> L/d kg

##########################################################
###                    absorption                      ###  
##########################################################
# Parent
# calculate of ka; use rat ka first 
# R_rat             <- 0.2             # 0.2cm # radius of rat jejunum
Peff_human        <- 10^(0.4926 * log10(Peff_parent * 1e6) - 0.1454)         # Peff is in the unit of e-6 cm/s; 1e-4 cm/s
Peff_rat          <- (Peff_human - 0.03)/3.6                                 # 1e-4 cm/s
ka                <-  12.1847922 
fa                <- 0.9                                                     # absorption fraction;

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
factor            <-  1
Qgfr_broiler      <-  3.4  * BW_broiler                                                 ## L/d

# laying hen Abstract, Gasthuys, Elke 2019 Comparative physiology of
Qgfr_parent       <-  2.57 /1000 * 60 * 24                       ## 2.6 / 1000 * 60  * 24 * factor   L/d/kg BW L for laying hen and duck
Qgfr_daughter     <-  2.57 /1000 * 60 * 24  
#=====================     End of user input    =======================#



#######################################################################
###                     partition coefficients                      ###  to unbound plasma; 
#######################################################################
## predicting the tissue to unbound plasma partition coefficients for the tissues contained in the tissue.data table 
# predicting partitioning cofficient for different tissues, rightnow assume using adjusted Funbound 
Rblood2plasma_parent     <- if(is.numeric(Rblood2plasma_parent)){Rblood2plasma_parent}else{calc_Rblood2plasma_chicken(LogP_parent, Compound_type_parent, fub_parent, pKa_a_parent, Hematocrit_fraction )}
Rblood2plasma_daughter   <- if(is.numeric(Rblood2plasma_daughter)){Rblood2plasma_daughter}else{calc_Rblood2plasma_chicken(LogP_daughter, Compound_type_daughter, fub_daughter, pKa_a_daughter, Hematocrit_fraction )}

######            calculation of Kw and Ky                 ######
# Yolk & white
ky_parent            <-  calc_ky(pKa_a_parent,   pKa_b_parent,   fub_parent) 
ky_daughter          <-  calc_ky(pKa_a_daughter, pKa_b_daughter, fub_daughter) 
kw_parent            <-  calc_kw_Kdw(pKa_a_parent,   pKa_b_parent,   fub_parent) 
kw_daughter          <-  calc_kw_Kdw(pKa_a_daughter, pKa_b_daughter, fub_daughter) 

#=================================================================================
#                                    PBPK                                        #
#=================================================================================
########################################################################
###                    define modeling variables                     ###  
########################################################################
days          <- 5
Time_max      <- days                                     # 5 d  

talbumen      <- 10 /24                                   # 10 h
tlag          <- 1                                        # time interval between ovulation and egg lay (approximately one day or 24 h)
tsig          <- 2                                        # pre-ovulation time of maximum follicle growth rate (48h, time for maximum growth achieved)

StartTime     <- 0                                        # Time_d
StopTime      <- Time_max                                 
dt            <- 0.01
Times         <- seq(StartTime, StopTime, dt)

dose          <- 7.66
#dose          <- 7.92
#Oral_input    <- 7.66 * 1000 / MW_parent                  # umol/kg/d, STUDY 1
Oral_input    <- dose * 1000 / MW_parent                  # umol/kg/d, STUDY 2

# dose pattern \HTTK\Bird\Intake pattern\pattern.xlsx (Clark et al. 2019; Choi et al. 2004)
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


value_percent  <- c(rep(value_1, days))
value          <- value_percent  *  Oral_input 
time           <- seq(0, (days), by = 1/24)

source(file.path("HTTK", "Bird", "TMX", 'Laying hen', 'event_repeated.R'))

########################################################################################################################
#                                                  PBPK MODEL RUN                                                      #
########################################################################################################################
initState <- c(Agutlumen_parent  = 0, 
               Agut_parent       = 0,  
               Aliver_parent     = 0,
               Akidney_parent    = 0,
               Aadipose_parent   = 0,
               Amuscle_parent    = 0,
               Arest_parent      = 0,
               Ww                = egg_white/1000,
               Arep_parent       = 0,
               Ayolk1_parent     = 0,
               Ayolk2_parent     = 0,
               Ayolk3_parent     = 0,
               Ayolk4_parent     = 0,
               Ayolk5_parent     = 0,
               Ayolk6_parent     = 0,
               Ayolk7_parent     = 0,
               Ayolk8_parent     = 0,
               Ayolk9_parent     = 0,
               Awhite_parent     = 0,  # no 19
               Wyolk1     = egg_yolk / (1+exp(-(0-(0.38-tsig-tlag)))) /1000,
               Wyolk2     = egg_yolk / (1+exp(-(0-(1.38-tsig-tlag)))) /1000,
               Wyolk3     = egg_yolk / (1+exp(-(0-(2.38-tsig-tlag)))) /1000,
               Wyolk4     = egg_yolk / (1+exp(-(0-(3.38-tsig-tlag)))) /1000,
               Wyolk5     = egg_yolk / (1+exp(-(0-(4.38-tsig-tlag)))) /1000,
               Wyolk6     = egg_yolk / (1+exp(-(0-(5.38-tsig-tlag)))) /1000,
               Wyolk7     = egg_yolk / (1+exp(-(0-(6.38-tsig-tlag)))) /1000,
               Wyolk8     = egg_yolk / (1+exp(-(0-(7.38-tsig-tlag)))) /1000,
               Wyolk9     = egg_yolk / (1+exp(-(0-(8.38-tsig-tlag)))) /1000,
               Aven_parent        = 0,
               Alung_parent       = 0,
               Aart_parent        = 0,
               Aurine_parent      = 0,
               AUC_Cplasma_parent = 0,
               AUC_Cblood_parent  = 0,  #34
               
               Agut_daughter       = 0,  
               Aliver_daughter     = 0,
               Akidney_daughter    = 0,
               Aadipose_daughter   = 0,
               Amuscle_daughter    = 0,
               Arest_daughter      = 0,
               Arep_daughter       = 0, #41
               Ayolk1_daughter     = 0,
               Ayolk2_daughter     = 0,
               Ayolk3_daughter     = 0,
               Ayolk4_daughter     = 0,
               Ayolk5_daughter     = 0,
               Ayolk6_daughter     = 0,
               Ayolk7_daughter     = 0,
               Ayolk8_daughter     = 0,
               Ayolk9_daughter     = 0, #50
               Awhite_daughter     = 0, 
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
parms <- c( dt            = dt,
            ka            = ka,                 
            fa            = fa,  
            Clint_parent         = Clint_hep_parent   * 37, 
            Clint_daughter       = Clint_hep_daughter / 3,   
            
            Rblood2plasma_parent   = Rblood2plasma_parent,  
            Rblood2plasma_daughter = Rblood2plasma_daughter,   
            BW = BW,              
            
            fub_parent           = fub_parent,   
            fub_daughter         = fub_daughter, 
            fraction_daughter    = 1,
            
            # ratio
            Kgut2pu_parent       = Kgut2pu_parent    ,      
            Kkidney2pu_parent    = Kkidney2pu_parent ,   
            Kliver2pu_parent     = Kliver2pu_parent  ,     
            Klung2pu_parent      = Klung2pu_parent   ,      
            Krest2pu_parent      = Krest2pu_parent   , 
            Kadipose2pu_parent   = Kadipose2pu_parent ,
            Kmuscle2pu_parent    = Kmuscle2pu_parent ,
            Krep2pu_parent       = Krep2pu_parent    ,
            
            Kgut2pu_daughter       = Kgut2pu_daughter    ,      
            Kkidney2pu_daughter    = Kkidney2pu_daughter ,   
            Kliver2pu_daughter     = Kliver2pu_daughter  ,     
            Klung2pu_daughter      = Klung2pu_daughter   ,      
            Krest2pu_daughter      = Krest2pu_daughter   , 
            Kadipose2pu_daughter   = Kadipose2pu_daughter ,
            Kmuscle2pu_daughter    = Kmuscle2pu_daughter ,
            Krep2pu_daughter       = Krep2pu_daughter    ,

            tsig          = tsig,
            tlag          = tlag,
            ky_parent     = ky_parent, 
            ky_daughter   = ky_daughter,
            kw_parent     = kw_parent,
            kw_daughter   = kw_daughter,
            
            # Parametes for flow 
            Qart          = Qart, 
            Qgfr_parent   = Qgfr_parent * 1.4 ,   
            Qgfr_daughter = Qgfr_daughter * 5.9,
            
            Qgut          = Qgut,
            Qkidney       = Qkidney,
            Qliver        = Qliver,
            Qrest         = Qart - (Qliver + Qkidney + Qmuscle + Qadipose + Qgut + Qrep), 
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
            Vadipose    = Vadipose )



df  <- ode(y       = initState, 
           times   = Times, 
           func    = pbtk5cpt_repeated, 
           parms   = parms,
           atol    = 1e-6, rtol = 1e-8,
           events  = list(func = event_repeated, root = TRUE), 
           rootfun = root_repeated ,
           method  = 'lsoda')

df               <- as.data.frame(df)
colnames(df)[1]  <- "Time"

file_path     <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", "Hen")
write.csv(df, paste(file_path, '/hen.oral.', dose, 'mgkg.csv', sep = ''), row.names = FALSE)

#=============================================================================
############                End of PBPK model                    #############
#=============================================================================
study1 <- "Metabolism_O_1998_hen"  # 7.66
study2 <- "Metabolism_T_1998_hen"  # 7.92

######           Plot           ###########
file_path                     <- file.path(current_dir, "HTTK", "Bird", "TMX", "InvivoData.xlsx")
Metabolism_1998_TMX           <- read_excel(file_path, sheet = study1)
#Metabolism_1998_TMX           <- read_excel(file_path, sheet = study2)

# urine
Metabolism_1998_excreta       <- Metabolism_1998_TMX[Metabolism_1998_TMX$Matrix == 'Excreta',]
Metabolism_1998_excreta$TMX_excreta  <- Metabolism_1998_excreta$'%TMX' /100 *  Metabolism_1998_excreta$Concentration        # Unit: % of dose
Metabolism_1998_excreta$CLO_excreta  <- Metabolism_1998_excreta$'%CGA322704' /100 *  Metabolism_1998_excreta$Concentration  # Unit: % of dose
df_Metabolism_1998_excreta    <- subset(Metabolism_1998_excreta, select = c('Time_h', "TMX_excreta", "CLO_excreta"))

# egg white
Metabolism_1998_egg_white               <- Metabolism_1998_TMX[Metabolism_1998_TMX$Matrix == 'Egg_white',]
Metabolism_1998_egg_white$TMXConc_umolL   <- Metabolism_1998_egg_white$Concentration * 1000 * Metabolism_1998_egg_white$'%TMX' / 100 / MW_parent     # [umol/L]
Metabolism_1998_egg_white$CLOConc_umolL   <- Metabolism_1998_egg_white$Concentration * 1000 * Metabolism_1998_egg_white$'%CGA322704' / 100 / MW_parent
df_Metabolism_1998_eggwhite             <- subset(Metabolism_1998_egg_white, select = c('Time_h_range', 'Time_h', "TMXConc_umolL", "CLOConc_umolL", "Weight_g", "Matrix"))

# egg yolk
Metabolism_1998_egg_yolk               <- Metabolism_1998_TMX[Metabolism_1998_TMX$Matrix == 'Egg_yolk',] 
Metabolism_1998_egg_yolk$TMXConc_umolL   <- Metabolism_1998_egg_yolk$Concentration * 1000 * Metabolism_1998_egg_yolk$'%TMX' / 100 / MW_parent
Metabolism_1998_egg_yolk$CLOConc_umolL   <- Metabolism_1998_egg_yolk$Concentration * 1000 * Metabolism_1998_egg_yolk$'%CGA322704' / 100 / MW_parent
df_Metabolism_1998_yolk                <- subset(Metabolism_1998_egg_yolk, select = c('Time_h_range', 'Time_h', "TMXConc_umolL", "CLOConc_umolL", "Weight_g", "Matrix"))

# blood
Metabolism_1998_blood                  <- Metabolism_1998_TMX[Metabolism_1998_TMX$Matrix == 'Blood' ,]
Metabolism_1998_blood$TMXBloodConc_umolL <- Metabolism_1998_blood$Concentration * 1000 * Metabolism_1998_blood$'%TMX' / 100/ MW_parent
Metabolism_1998_blood$CLOBloodConc_umolL <- Metabolism_1998_blood$Concentration * 1000 * Metabolism_1998_blood$'%CGA322704'/100 / MW_parent
df_Metabolism_1998_blood               <- subset(Metabolism_1998_blood, select = c('Time_h', "Matrix", "TMXBloodConc_umolL", "CLOBloodConc_umolL")) 

# liver
Metabolism_1998_liver                  <- Metabolism_1998_TMX[Metabolism_1998_TMX$Matrix == 'Liver' ,]
Metabolism_1998_liver$TMXliverConc_umolL <- Metabolism_1998_liver$Concentration * 1000 * Metabolism_1998_liver$'%TMX' / 100/ MW_parent
Metabolism_1998_liver$CLOliverConc_umolL <- Metabolism_1998_liver$Concentration * 1000 * Metabolism_1998_liver$'%CGA322704'/100 / MW_parent 
df_Metabolism_1998_liver               <- subset(Metabolism_1998_liver, select = c('Time_h', "Matrix", "TMXliverConc_umolL", "CLOliverConc_umolL")) 

# fat
Metabolism_1998_fat                  <- Metabolism_1998_TMX[Metabolism_1998_TMX$Matrix == 'Fat' ,]
Metabolism_1998_fat$TMXfatConc_umolL <- Metabolism_1998_fat$Concentration * 1000 * Metabolism_1998_fat$'%TMX' / 100/ MW_parent
Metabolism_1998_fat$CLOfatConc_umolL <- Metabolism_1998_fat$Concentration * 1000 * Metabolism_1998_fat$'%CGA322704'/100 / MW_parent 
df_Metabolism_1998_fat               <- subset(Metabolism_1998_fat, select = c('Time_h', "Matrix", "TMXfatConc_umolL", "CLOfatConc_umolL")) 



# Mass balance
# it is difficult to evaluate for laying hen because state variable related to egg white and yolk were set to zero routinly.
ggplot() +
  geom_line (data = df, aes(Time, Mass_parent_bal), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14)

ggplot() +
  geom_line (data = df, aes(Time, Mass_daughter_bal), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14)

ggplot() +
  geom_line (data = df, aes(Time, Mass_bal), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14)




ggplot() +
  geom_line(data = df, aes(Time-0.3, C_blood_parent, color = "Prediction"),  lwd=0.7) +
  geom_point(data = df_Metabolism_1998_blood  , aes(Time_h/24, TMXBloodConc_umolL, shape = 'Observation'),col="#ff5044",  size=1.8) +
  theme(text = element_text(size = 20)) + xlim(0, 4) + 
  xlab("Time (d)") + 
  ylab(expression("TMX Blood Concentration ("*mu*"mol/L)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#00AFBB', "#ff5044"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(17)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Thiamethoxam (TMX)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen.", dose, ".TMXBlood_0106.tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 16, height = 10,units = "cm", compression = 'lzw')


ggplot() +
  geom_line(data = df, aes(Time-0.3, C_blood_daughter, color = "Prediction"),  lwd=0.7) +
  geom_point(data = df_Metabolism_1998_blood  , aes(Time_h/24, CLOBloodConc_umolL, shape = 'Observation'),col="#d8923d",  size=1.8) +
  theme(text = element_text(size = 20)) + xlim(0, 4) + 
  xlab("Time (d)") + 
  ylab(expression("CTD Blood Concentration ("*mu*"mol/L)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#276DC2', "#d8923d"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(17)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Clothianidin (CTD)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen.", dose, ".CTDBlood_0106.tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 16, height = 10,units = "cm", compression = 'lzw')


df_blood_parent           <- as.data.frame(df %>% dplyr::select(Time, C_blood_parent))
df_blood_daughter         <- as.data.frame(df %>% dplyr::select(Time, C_blood_daughter))
colnames(df_blood_parent)[2]   <- 'Value'
colnames(df_blood_daughter)[2] <- 'Value'
df_blood_parent$Time      <- df_blood_parent$Time - 0.3
df_blood_daughter$Time    <- df_blood_daughter$Time - 0.3
df_blood_parent           <- df_blood_parent[df_blood_parent$Time >= 0,]
df_blood_daughter         <- df_blood_daughter[df_blood_daughter$Time >= 0,]
df_blood_parent$Class     <- 'TMX'
df_blood_daughter$Class   <- "CTD"
df_blood_parent$Matrix    <- 'Blood'
df_blood_daughter$Matrix  <- "Blood"

obs_1998_blood_parent   <- as.data.frame(df_Metabolism_1998_blood  %>% dplyr::select(Time_h, TMXBloodConc_umolL))
obs_1998_blood_daughter <- as.data.frame(df_Metabolism_1998_blood  %>% dplyr::select(Time_h, CLOBloodConc_umolL))
obs_1998_blood_parent$Time_h     <- obs_1998_blood_parent$Time_h/24
obs_1998_blood_daughter$Time_h   <- obs_1998_blood_daughter$Time_h/24
colnames(obs_1998_blood_parent)[2]    <- 'Value'
colnames(obs_1998_blood_daughter)[2]  <- 'Value'
obs_1998_blood_parent$Class      <- 'TMX'
obs_1998_blood_daughter$Class    <- "CTD"
obs_1998_blood_parent$Matrix     <- 'Blood'
obs_1998_blood_daughter$Matrix   <- "Blood"

### panel plot
# http://zevross.com/blog/2019/04/02/easy-multi-panel-plots-in-r-using-facet_wrap-and-facet_grid-from-ggplot2/
ggplot() +
  geom_line(data = df, aes(Time-0.3, C_liver_parent, color = "Prediction"), lwd=0.7) +
  geom_point(data = df_Metabolism_1998_liver  , aes(Time_h/24, TMXliverConc_umolL, shape = 'Observation'),col="#ff5044", size=1.8) + 
  theme(text = element_text(size = 20)) + xlim(0, 4) + 
  xlab("Time (d)") + 
  ylab(expression("TMX Liver Concentration ("*mu*"mol/kg)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#00AFBB', "#ff5044"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(17)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Thiamethoxam (TMX)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen.", dose, ".TMXLiver.tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 16, height = 10,units = "cm", compression = 'lzw')

  
ggplot() +
  geom_line(data = df, aes(Time-0.3, C_liver_daughter, color = "Prediction"), lwd=0.7) +
  geom_point(data = df_Metabolism_1998_liver  , aes(Time_h/24, CLOliverConc_umolL, shape = 'Observation'),col="#d8923d", size=1.8) + ylab("Concentration")+
  theme(text = element_text(size = 20)) + xlim(0, 4) + 
  xlab("Time (d)") + 
  ylab(expression("CTD Liver Concentration ("*mu*"mol/kg)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#276DC2', "#d8923d"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(17)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Clothianidin (CTD)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen.", dose, ".CTDLiver.tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 16, height = 10,units = "cm", compression = 'lzw')


ggplot() +
  geom_line (data = df, aes(Time-0.3, Aurine_parent / (Oral_input  * 4) * 100,  color = "Prediction"),  lwd=0.8) + 
  geom_point(data = df_Metabolism_1998_excreta, aes(Time_h/24, TMX_excreta, shape = 'Observation'),size=1.8, colour = "#FC4E07") + 
  theme(text = element_text(size = 20)) + xlim(0, 4) + 
  xlab("Time (d)") + 
  ylab(("Accumulated TMX mass in excreta (%)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#00AFBB', "#FC4E07"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(17)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Thiamethoxam (TMX)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen.", dose, ".TMXurine_0106.tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 16, height = 10,units = "cm", compression = 'lzw')


ggplot() +
  geom_line (data = df, aes(Time-0.3, Aurine_daughter/ (Oral_input  * 4)  * 100,  color = "Prediction"),  lwd=0.8) + 
  geom_point(data = df_Metabolism_1998_excreta, aes(Time_h/24, CLO_excreta, shape = 'Observation'),size=1.8, colour = "#d8923d") + 
  xlab("Time (d)") + 
  ylab(("Accumulated CTD mass in excreta (%)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#276DC2', "#d8923d"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(17)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Clothianidin (CTD)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen.", dose, ".CTDurine_0106.tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 16, height = 10,units = "cm", compression = 'lzw')


df_urine_parent           <- as.data.frame(df %>% dplyr::select(Time, Aurine_parent))
df_urine_daughter         <- as.data.frame(df %>% dplyr::select(Time, Aurine_daughter))
df_urine_parent$Aurine_parent      <- df_urine_parent$Aurine_parent/ (Oral_input  * 4)  * 100
df_urine_daughter$Aurine_daughter  <- df_urine_daughter$Aurine_daughter/ (Oral_input  * 4)  * 100
df_urine_parent$Time      <- df_urine_parent$Time - 0.3
df_urine_daughter$Time    <- df_urine_daughter$Time -0.3
df_urine_parent           <- df_urine_parent[df_urine_parent$Time >= 0,]
df_urine_daughter         <- df_urine_daughter[df_urine_daughter$Time >= 0,]
colnames(df_urine_parent)[2]   <- 'Value'
colnames(df_urine_daughter)[2] <- 'Value'
df_urine_parent$Class     <- 'TMX'
df_urine_daughter$Class    <- "CTD"
df_urine_parent$Matrix    <- 'Excreta'
df_urine_daughter$Matrix   <- "Excreta"

obs_1998_urine_parent   <- as.data.frame(df_Metabolism_1998_excreta  %>% dplyr::select(Time_h, TMX_excreta))
obs_1998_urine_daughter <- as.data.frame(df_Metabolism_1998_excreta  %>% dplyr::select(Time_h, CLO_excreta))
obs_1998_urine_parent$Time_h     <- obs_1998_urine_parent$Time_h/24
obs_1998_urine_daughter$Time_h   <- obs_1998_urine_daughter$Time_h/24
colnames(obs_1998_urine_parent)[2]    <- 'Value'
colnames(obs_1998_urine_daughter)[2]  <- 'Value'
obs_1998_urine_parent$Class      <- 'TMX'
obs_1998_urine_daughter$Class    <- "CTD"
obs_1998_urine_parent$Matrix     <- 'Excreta'
obs_1998_urine_daughter$Matrix   <- "Excreta"

#===========       Egg white     =============
ggplot() +
  geom_line (data = df, aes(Time, Awhite_parent), col="#00AFBB", lwd=2) + ylab("Concentration (ug/L)") +
  theme(text = element_text(size = 20))
ggplot() +
  geom_line (data = df, aes(Time, C_white_parent), col="#00AFBB", lwd=2) + ylab("Concentration (ug/L)") +
  theme(text = element_text(size = 20))

df_white_out                       <- df[(df$Time - 0.37)%%1 < 1e-5,]
file_path     <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", "Hen")
write.csv(df_white_out, paste(file_path, '/hen.eggwhite.oral.', dose, 'mgkg.csv', sep = ''), row.names = FALSE)

ggplot() +
  geom_line(data = df_white_out, aes(Time, C_white_parent,  color = "Prediction"),  lwd=1) + 
  geom_point(data = df_Metabolism_1998_eggwhite  , aes(Time_h/24 +0.3, TMXConc_umolL, shape = 'Observation' ),size=2, colour = "#FC4E07") + 
  theme(text = element_text(size = 20)) +
  xlab("Time (d)") + 
  ylab(expression("TMX Egg White Concentration ("*mu*"mol/L)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#00AFBB', "#FC4E07"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(17)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Thiamethoxam (TMX)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen.", dose, ".TMXEggWhite.tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 16, height = 10,units = "cm", compression = 'lzw')


ggplot() +
  geom_line(data = df_white_out, aes(Time, C_white_daughter,  color = "Prediction"), lwd=1) + 
  geom_point(data = df_Metabolism_1998_eggwhite  , aes(Time_h/24 + 0.3, CLOConc_umolL, shape = 'Observation' ), color = '#d8923d',size=2) + 
  theme(text = element_text(size = 20)) +
  xlab("Time (d)") + 
  ylab(expression("CTD Egg White Concentration ("*mu*"mol/L)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#276DC2', "#d8923d"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(17)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Clothianidin (CTD)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen.", dose, ".CTDEggWhite.tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 16, height = 10,units = "cm", compression = 'lzw')


df_white_parent           <- as.data.frame(df_white_out %>% dplyr::select(Time, C_white_parent))
df_white_daughter         <- as.data.frame(df_white_out %>% dplyr::select(Time, C_white_daughter))
colnames(df_white_parent)[2]   <- 'Value'
colnames(df_white_daughter)[2] <- 'Value'
df_white_parent$Class      <- 'TMX'
df_white_daughter$Class    <- "CTD"
df_white_parent$Matrix     <- 'Egg white'
df_white_daughter$Matrix   <- "Egg white"

obs_1998_white_parent   <- as.data.frame(df_Metabolism_1998_eggwhite  %>% dplyr::select(Time_h, TMXConc_umolL))
obs_1998_white_daughter <- as.data.frame(df_Metabolism_1998_eggwhite  %>% dplyr::select(Time_h, CLOConc_umolL))
obs_1998_white_parent$Time_h     <- obs_1998_white_parent$Time_h/24 + 0.3
obs_1998_white_daughter$Time_h   <- obs_1998_white_daughter$Time_h/24 + 0.3
colnames(obs_1998_white_parent)[2]    <- 'Value'
colnames(obs_1998_white_daughter)[2]  <- 'Value'
obs_1998_white_parent$Class      <- 'TMX'
obs_1998_white_daughter$Class    <- "CTD"
obs_1998_white_parent$Matrix     <- 'Egg white'
obs_1998_white_daughter$Matrix   <- "Egg white"


ggplot() +
  geom_line (data = df, aes(Time, Ayolk1_parent), col="#00AFBB", lwd=2) + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Ayolk2_parent), col="blue", lwd=2)    + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Ayolk3_parent), col="red", lwd=2)     + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Ayolk4_parent), col="black", lwd=2)   + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Ayolk5_parent), col="orange", lwd=2)  + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Ayolk6_parent), col="green", lwd=2)   + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Ayolk7_parent), col="yellow", lwd=2)  + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Ayolk8_parent), col="grey", lwd=2)    + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Ayolk9_parent), col="white", lwd=2)   + ylab("Mass (ug/kg)") +
  theme(text = element_text(size = 20))#+ xlim(0, 100)

######### get yolk concentration ###########
df_yolk_out         <- subset(df, (Time - 0.37)%%1 < 1e-4 | (Time - 0.38)%%1 < 1e-4)
df_yolk_out_parent   <- df_yolk_out[c('Time', 'Ayolk1_parent', 'Ayolk2_parent', 'Ayolk3_parent', 'Ayolk4_parent', 
                                            'Ayolk5_parent', 'Ayolk6_parent', 'Ayolk7_parent', 'Ayolk8_parent', 'Ayolk9_parent',
                                            'Cyolk1_parent', 'Cyolk2_parent', 'Cyolk3_parent', 'Cyolk4_parent', 'Cyolk5_parent', 
                                            'Cyolk6_parent', 'Cyolk7_parent', 'Cyolk8_parent', 'Cyolk9_parent')]
df_Cyolk_out_parent  <- df_yolk_out_parent[c('Time', 'Cyolk1_parent', 'Cyolk2_parent', 'Cyolk3_parent', 'Cyolk4_parent', 'Cyolk5_parent', 
                                             'Cyolk6_parent', 'Cyolk7_parent', 'Cyolk8_parent', 'Cyolk9_parent')]

rownames(df_Cyolk_out_parent) <- 1:nrow(df_Cyolk_out_parent)
mdata                         <- melt(df_Cyolk_out_parent, id=c("Time"))
mdata2                        <- mdata[which(mdata$value == 0)-1,] 

file_path     <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", "Hen")
write.csv(mdata2, paste(file_path, '/hen.eggyolk.TMX.oral.', dose, 'mgkg.csv', sep = ''), row.names = FALSE)


ggplot() +
  geom_line(data = mdata2, aes(Time, value, col="Prediction"), lwd=1,  size = 3) + 
  geom_point(data = df_Metabolism_1998_yolk , aes(Time_h/24+ 0.3, TMXConc_umolL, shape = 'Observation'),size=2, colour = "#FC4E07") + 
  theme(text = element_text(size = 20)) +
  xlab("Time (d)") + 
  ylab(expression("TMX Egg Yolk Concentration ("*mu*"mol/L)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#00AFBB', "#FC4E07"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(17)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Thiamethoxam (TMX)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen.", dose, ".TMXEggYolk.tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 16, height = 10,units = "cm", compression = 'lzw')


df_yolk_out_daughter   <- df_yolk_out[c('Time', 'Ayolk1_daughter', 'Ayolk2_daughter', 'Ayolk3_daughter', 'Ayolk4_daughter', 
                                      'Ayolk5_daughter', 'Ayolk6_daughter', 'Ayolk7_daughter', 'Ayolk8_daughter', 'Ayolk9_daughter',
                                      'Cyolk1_daughter', 'Cyolk2_daughter', 'Cyolk3_daughter', 'Cyolk4_daughter', 'Cyolk5_daughter', 
                                      'Cyolk6_daughter', 'Cyolk7_daughter', 'Cyolk8_daughter', 'Cyolk9_daughter')]
df_Cyolk_out_daughter  <- df_yolk_out_daughter[c('Time', 'Cyolk1_daughter', 'Cyolk2_daughter', 'Cyolk3_daughter', 'Cyolk4_daughter', 'Cyolk5_daughter', 
                                             'Cyolk6_daughter', 'Cyolk7_daughter', 'Cyolk8_daughter', 'Cyolk9_daughter')]

rownames(df_Cyolk_out_daughter) <- 1:nrow(df_Cyolk_out_daughter)
mdata_daughter                         <- melt(df_Cyolk_out_daughter, id=c("Time"))
mdata2_daughter                        <- mdata_daughter[which(mdata_daughter$value == 0)-1,] 

file_path     <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", "Hen")
write.csv(mdata2_daughter, paste(file_path, '/hen.eggyolk.CTD.oral.', dose, 'mgkg.csv', sep = ''), row.names = FALSE)

ggplot() +
  geom_line(data = mdata2_daughter, aes(Time, value, col="Prediction"), lwd=1, size=3) + 
  geom_point(data = df_Metabolism_1998_yolk , aes(Time_h/24 + 0.3, CLOConc_umolL, shape = 'Observation'),size=2.5, col = "#d8923d") + 
  theme(text = element_text(size = 20)) +
  xlab("Time (d)") + 
  ylab(expression("CTD Egg Yolk Concentration ("*mu*"mol/L)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#276DC2', "#d8923d"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(17)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Clothianidin (CTD)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen.", dose, ".CTDEggYolk.tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 16, height = 10,units = "cm", compression = 'lzw')


df_yolk_parent           <- as.data.frame(mdata2 %>% dplyr::select(Time, value))
df_yolk_daughter         <- as.data.frame(mdata2_daughter %>% dplyr::select(Time, value))
colnames(df_yolk_parent)[2]   <- 'Value'
colnames(df_yolk_daughter)[2] <- 'Value'
df_yolk_parent$Class      <- 'TMX'
df_yolk_daughter$Class    <- "CTD"
df_yolk_parent$Matrix     <- 'Egg yolk'
df_yolk_daughter$Matrix   <- "Egg yolk"

obs_1998_yolk_parent   <- as.data.frame(df_Metabolism_1998_yolk %>% dplyr::select(Time_h, TMXConc_umolL))
obs_1998_yolk_daughter <- as.data.frame(df_Metabolism_1998_yolk %>% dplyr::select(Time_h, CLOConc_umolL))
obs_1998_yolk_parent$Time_h     <- obs_1998_yolk_parent$Time_h/24 + 0.3
obs_1998_yolk_daughter$Time_h   <- obs_1998_yolk_daughter$Time_h/24 + 0.3

colnames(obs_1998_yolk_parent)[2]    <- 'Value'
colnames(obs_1998_yolk_daughter)[2]  <- 'Value'
obs_1998_yolk_parent$Class      <- 'TMX'
obs_1998_yolk_daughter$Class    <- "CTD"
obs_1998_yolk_parent$Matrix     <- 'Egg yolk'
obs_1998_yolk_daughter$Matrix   <- "Egg yolk"

df_byw                <- data.frame(rbind(df_blood_parent, df_blood_daughter, 
                                    df_white_parent, df_white_daughter,
                                    df_yolk_parent, df_yolk_daughter))
df_byw$Class          <- factor(df_byw$Class, levels =  c('TMX', 'CTD'))
df_byw$Matrix         <- factor(df_byw$Matrix, levels =  c('Blood', 'Egg white', "Egg yolk"))
obs_byw               <- data.frame(rbind(obs_1998_blood_parent, obs_1998_blood_daughter, 
                                    obs_1998_white_parent,obs_1998_white_daughter, 
                                    obs_1998_yolk_parent, obs_1998_yolk_daughter))
obs_byw$Class         <- factor(obs_byw$Class, levels =  c('TMX', 'CTD'))
obs_byw$Matrix        <- factor(obs_byw$Matrix, levels =  c('Blood', 'Egg white', "Egg yolk"))


# Panel plot
ggplot(obs_byw) +
  geom_point(aes(Time_h, Value, color = Class), size = 2, shape = 24, fill = 'white', stroke = 1) +
  geom_line(data = df_byw, aes(Time, Value, color = Class), lwd=0.6) +
  #facet_wrap(~Class) +  
  #facet_grid(Matrix~Class, space="free", scales="free")+
  facet_grid(Matrix~Class, space="fixed", scales="free")+
  theme_bw(base_size = 16) + theme(legend.position = 'none') + 
  theme(strip.text = element_text(size = 16)) +
  scale_color_brewer(palette = "Set1") + 
  labs(title = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''))+
  ylab(expression("Concentration ("*mu*"mol/L)")) + xlab('Time (d)')
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Hen', paste("plot.layinghen_Panel.", dose, ".tiff", sep = ''))
ggsave(file.path,
       dpi = 300, width = 26, height = 20,units = "cm",compression = 'lzw')


ggplot(df_byw) +
  geom_line(aes(Time, Value, color = Class), lwd=1) +
  geom_point(data = obs_byw, aes(Time_h, Value, color = Class), size = 2, shape = 24, fill = 'white', stroke = 1) +
   #facet_wrap(~Class) +  
  #facet_grid(Matrix~Class, space="free", scales="free")+
  facet_grid(Matrix~Class, space="free")+
  theme_bw(base_size = 14) + theme(legend.position = 'none') + 
  scale_color_brewer(palette = "Set1") + 
  labs(title = paste("Laying hen [Repeated oral dose for 4 consecutive days at ", dose, " mg/kg BW/d]", sep = ''))+
  ylab(expression("Concentration ("*mu*"mol/L)")) + xlab('Time (d)')


# urine
ggplot() +
  geom_line (data = df, aes(Time,Ww), col="#00AFBB", lwd=2) + ylab("Mass (kg)") +
  theme(text = element_text(size = 20))+ xlim(0, 10)

ggplot() +
  geom_line (data = df, aes(Time,Wyolk1), col="#00AFBB", lwd=2) + 
  geom_line (data = df, aes(Time,Wyolk2), col="blue", lwd=2) +
  geom_line (data = df, aes(Time,Wyolk3), col="red", lwd=2) +
  geom_line (data = df, aes(Time,Wyolk4), col="black", lwd=2) + 
  geom_line (data = df, aes(Time,Wyolk5), col="orange", lwd=2) +
  geom_line (data = df, aes(Time,Wyolk6), col="green", lwd=2) +
  geom_line (data = df, aes(Time,Wyolk7), col="yellow", lwd=2) +
  geom_line (data = df, aes(Time,Wyolk8), col="grey", lwd=2) + 
  geom_line (data = df, aes(Time,Wyolk9), col="white", lwd=2) +ylab("Mass (kg)") +
  theme(text = element_text(size = 20))# + xlim(120, 135)

##########     end   ############




ggplot() +
  geom_line (data = df, aes(Time, Wyolk1), col="#00AFBB", lwd=2) + ylab("Mass (kg)") +
  theme(text = element_text(size = 20))+ xlim(0, 10)

ggplot() +
  geom_line (data = df, aes(Time, Ayolk5), col="#00AFBB", lwd=2) + ylab("Concentration (ug/L)") +
  theme(text = element_text(size = 20))+ xlim(0, 10)
ggplot() +
  geom_line (data = df, aes(Time, Wyolk5), col="#00AFBB", lwd=2) + ylab("Concentration (ug/L)") +
  theme(text = element_text(size = 20))+ xlim(0, 10)

ggplot() +
  geom_line (data = df, aes(Time, Cyolk1), col="#00AFBB", lwd=2) + ylab("Concentration (ug/L)") +
  geom_line (data = df, aes(Time, Cyolk2), col="blue", lwd=2) + ylab("Concentration (ug/L)") +
  geom_line (data = df, aes(Time, Cyolk3), col="red", lwd=2) + ylab("Concentration (ug/L)") +
  geom_line (data = df, aes(Time, Cyolk4), col="black", lwd=2) + ylab("Concentration (ug/L)") +
  geom_line (data = df, aes(Time, Cyolk5), col="orange", lwd=2) + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Cyolk6), col="green", lwd=2) + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Cyolk7), col="yellow", lwd=2) + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Cyolk8), col="grey", lwd=2) + ylab("Mass (ug/kg)") +
  geom_line (data = df, aes(Time, Cyolk9), col="white", lwd=2) + ylab("Mass (ug/kg)") +
  theme(text = element_text(size = 20)) #+ xlim(0, 20)
  

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
#dev.off()

#==================================================================================================
#                                      Species: Laying Hen                                        # 
#                                     No of compartment: 5                                        #
#==================================================================================================
# Reference
# PKNCA: https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk: https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk

SPECIES  <- 'Laying_hen'
BW             <- 1.474                  # kg, body weight of 18 weeks female laying hen, Residue study 2000
BW_broiler     <- 2.5                   # kg
BW_human       <- 70                    # kg
SDBW           <- 0.08

source("C:/xxx/OneDrive - Syngenta/HTTK/Bird/General Code/Avian_Wegg.R")
source(paste("C:/xxx/OneDrive - Syngenta/HTTK/Bird/General Code/", SPECIES, "_PhyData.R", sep = ''))
source("C:/xxx/OneDrive - Syngenta/HTTK/Bird/General Code/Metabolism.R")

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

source("C:/xxx/OneDrive - Syngenta/HTTK/Bird/General Code/Partition.R")


## Permeability-surface area product L/d kg
PSc_adipose_parent      <- Peff_parent * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000            # Peff * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000  # cm/s -> cm/d * g/kg * cm2/g -> cm3/d kg -> L/d kg
PSc_adipose_daughter    <- Peff_daughter * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000          # Peff * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000  # cm/s -> cm/d * g/kg * cm2/g -> cm3/d kg -> L/d kg


##########################################################
###                    absorption                      ###  
##########################################################
# Parent
# calculate of ka; use rat ka first 
#R_rat             <- 0.2             # 0.2cm # radius of rat jejunum
Peff_human        <- 10^(0.4926 * log10(Peff_parent * 1e6) - 0.1454)         # Peff is in the unit of e-6 cm/s; 1e-4 cm/s
Peff_rat          <- (Peff_human - 0.03)/3.6                            # 1e-4 cm/s
ka                <-  12.1847922 #2 * Peff_rat / R_rat / 1e4 * 3600 * 24            # d-1

fa        <- 0.9 #249.7 / 291.72                                   # absorption fraction;

####################################################################
###                    clearance (L/d/kg BW)                     ###  
####################################################################
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
Qgfr_broiler      <-  3.4  * BW_broiler                                                 ## L/d/kg BW

# laying hen Abstract, Gasthuys, Elke 2019 Comparative physiology of
Qgfr_parent       <-  2.57 /1000 * 60 * 24                       ## 2.6 / 1000 * 60  * 24 * factor   L/d/kg BW L for laying hen and duck
Qgfr_daughter     <-  2.57 /1000 * 60 * 24  
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
days          <- 28#+9 # actual dose start at day 9 because egg start lay at day 9
Time_max      <- days                                     # 5 d  

talbumen      <- 10 /24                                   # 10 h
tlag          <- 1                                        # time interval between ovulation and egg lay (approximately one day or 24 h)
tsig          <- 2                                        # pre-ovulation time of maximum follicle growth rate (48h, time for maximum growth achieved)

StartTime     <- 0                                        # Time_d
StopTime      <- Time_max                                 
dt            <- 0.01
Times         <- seq(StartTime, StopTime, dt)

Oral_input    <- 0.901 * 1000 / MW_parent                 # ug/kg/d  -> umol/kg/d

# dose pattern 89/91 Sitovition report; C:\xxx\OneDrive - Syngenta\HTTK\Bird\Intake pattern\pattern.xlsx
var                <- rep("Agutlumen_parent", 24 * days/1)                                   # dose pattern every one hour
time               <- seq(0, (days), by = 1/24)
time               <- time[-length(time)]
#value_1        <-     c(2.7,4.8,8.1,12.4,16.5,22.4,33.2,47.4,62.9,78.7,91.8,100)            # Table A6 in appendix; not used
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

time               <- seq(0, (days), by = 1/24)

Dose_events    <- data.frame(var    = rep("Agutlumen_parent", 4),
                             time   = c(0.3, 1.3, 2.3, 3.3),
                             value  = rep(Oral_input, 4),
                             method = "add")

### Approxfun is only appropriate for external input use (not for state variable change!!!)
######## scenario 1: pulse input daily, lasting for a certain period
'signal        <- data.frame(times = dose_times , import = rep(0, length(dose_times)))
signal$import <- ifelse((signal$times) %% 1 == 0, Dose  , 0)  # 100 mg  ---> ug/kg bw
signal
                           

input         <- approxfun(signal,  method = "constant", yleft = 0, yright = 0, f = 0, rule = 2)
input(seq(from = 11, to = 15, by = 0.1))
plot(input(seq(1, 20, by = 0.1)), type="p", xlab = "Time (min)", ylab = "Butadiene inhaled concentration (ppm)")'

######## scenario 2: realistic feeding scenario, following feed pattern (Clark et al. 2019; Choi et al. 2004)
'Intake_pattern    <- read.csv("C:/xxx/OneDrive - Syngenta/HTTK/Bird/Intake pattern/Pattern.csv")
Intake_pattern    <- na.omit(Intake_pattern )
x                 <- Intake_pattern $Time/24
y                 <- Intake_pattern $g * Dose_conc / BW             # g * ug/g -> ug/bw

xx                <- seq(from = 0, to = (Time_max - 1/24), by = 1/24)
yy                <- rep(y, times = Time_max, each = 1)

input             <- approxfun(xx, yy) 
curve(input(x), 0, Time_max, col = "green2")
curve(input(x), 0, 5, col = "green2")
points(xx, yy)'


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
               Awhite_parent     = 0, 
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
               AUC_Cblood_parent  = 0,
               
               Agut_daughter       = 0,  
               Aliver_daughter     = 0,
               Akidney_daughter    = 0,
               Aadipose_daughter   = 0,
               Amuscle_daughter    = 0,
               Arest_daughter      = 0,
               Arep_daughter       = 0,
               Ayolk1_daughter     = 0,
               Ayolk2_daughter     = 0,
               Ayolk3_daughter     = 0,
               Ayolk4_daughter     = 0,
               Ayolk5_daughter     = 0,
               Ayolk6_daughter     = 0,
               Ayolk7_daughter     = 0,
               Ayolk8_daughter     = 0,
               Ayolk9_daughter     = 0,
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


##########################      pbpk output     #############################  
## event triggered if time %%24 = 0
## https://stat.ethz.ch/pipermail/r-sig-dynamic-models/2013q1/000160.html
## https://stat.ethz.ch/pipermail/r-sig-dynamic-models/2013q1/000158.html
# https://stat.ethz.ch/pipermail/r-sig-dynamic-models/2018/000606.html          continuous infusion
# https://cran.r-project.org/web/packages/pksensi/vignettes/pbtk1cpt.html
# https://rstudio-pubs-static.s3.amazonaws.com/557239_2e4541dc478b4ef69a2ff2664da352bd.html  (inhalation)
# https://tpetzoldt.github.io/deSolve-forcing/deSolve-forcing.html    approxfun (useful for infusion and inhalation)

# https://stackoverflow.com/questions/71009453/combining-root-and-non-root-event-functions-in-r-desolve

'root <- function(t, y, parms) {return(t%%1)}
## set state variable Ww = 0
event <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  return(y)
}'

# DOSING WAS FOR 5 CONSECUTIVE DAYS FOLLOWED BY A WITHDRAWAL PERIOD OF 7 DAYS
eventdose <- data.frame(var    = rep("Agutlumen_parent", 4),
                        time   = c(0.3, 1.3, 2.3, 3.3),
                        value  = rep(Oral_input, 4),
                        method = "add")





eventW <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white parent
  y[51]      <- 0                                             # chemical mass in egg white daughter
  return(y)
}

event1 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white parent
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[10]      <- 0                                             # chemical mass in egg yolk 1 parent
  y[42]      <- 0                                             # chemical mass in egg yolk 1 daughter
  y[20]      <- 0.1/1000                                      # mass in egg yolk 1
  return(y)
}

event2 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[11]      <- 0                                             # chemical mass in egg yolk 2
  y[43]      <- 0                                             # chemical mass in egg yolk 2 daughter
  y[21]      <- 0.1/1000                                      # mass in egg yolk 2
  return(y)
}

event3 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[12]      <- 0                                             # chemical mass in egg yolk 3
  y[44]      <- 0                                             # chemical mass in egg yolk 3 daughter
  y[22]      <- 0.1/1000                                      # mass in egg yolk 3
  return(y)
}

event4 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[13]      <- 0                                             # chemical mass in egg yolk 4
  y[45]      <- 0                                             # chemical mass in egg yolk 41 daughter
  y[23]      <- 0.1/1000                                      # mass in egg yolk 4
  return(y)
}

event5 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[14]      <- 0                                             # chemical mass in egg yolk 5
  y[46]      <- 0                                             # chemical mass in egg yolk 5 daughter
  y[24]      <- 0.1/1000                                      # mass in egg yolk 5
  return(y)
}

event6 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[15]      <- 0                                             # chemical mass in egg yolk 6
  y[47]      <- 0                                             # chemical mass in egg yolk 6 daughter
  y[25]      <- 0.1/1000                                      # mass in egg yolk 6
  return(y)
}

event7 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[16]      <- 0                                             # chemical mass in egg yolk 7
  y[48]      <- 0                                             # chemical mass in egg yolk 7 daughter
  y[26]      <- 0.1/1000                                      # mass in egg yolk 7
  return(y)
}

event8 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[17]      <- 0                                             # chemical mass in egg yolk 8
  y[49]      <- 0                                             # chemical mass in egg yolk 8 daughter
  y[27]      <- 0.1/1000                                      # mass in egg yolk 8
  return(y)
}

event9 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[18]      <- 0                                             # chemical mass in egg yolk 9
  y[50]      <- 0                                             # chemical mass in egg yolk 9 daughter
  y[28]      <- 0.1/1000                                      # mass in egg yolk 9
  return(y)
}




event_dose0 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0                               
  return(y)
}

event_dose1 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.002022362                              
  return(y)
}

event_dose2 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.004044725                               
  return(y)
}

event_dose3 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.002022362                                
  return(y)
}

event_dose4 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.026560358                                 
  return(y)
}

event_dose5 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.047265253                               
  return(y)
}

event_dose6 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input *  0.048256832                                
  return(y)
}

event_dose7 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.048752621                                
  return(y)
}

event_dose8 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.053875778                                
  return(y)
}

event_dose9 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.059164198                                 
  return(y)
}

event_dose10 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.063461039                                
  return(y)
}

event_dose11 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.062799987                                
  return(y)
}

event_dose12 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.059494724                                
  return(y)
}

event_dose13 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.059494724                               
  return(y)
}

event_dose14 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.064948407                                
  return(y)
}

event_dose15 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.063791565                               
  return(y)
}

event_dose16 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.073376826                               
  return(y)
}

event_dose17 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.089076823                                
  return(y)
}

event_dose18 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.087919981                                
  return(y)
}

event_dose19 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.069245248                                
  return(y)
}

event_dose20 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.009302867                                
  return(y)
}

event_dose21 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.005123318                              
  return(y)
}

event_dose22 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0                              
  return(y)
}

event_dose23 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0                               
  return(y)
}

# egg lay start at day 9
event_NoDose <- function(t, y, parms) {
  y[1]       <- 0                               
  return(y)
}

root_repeated <- function(t, y, parms) {
  yroot      <- c(((t-0.38)%%1),   ((t-0.38)/9)%%1, ((t-1.38)/9)%%1, ((t-2.38)/9)%%1, 
                  ((t-3.38)/9)%%1, ((t-4.38)/9)%%1, ((t-5.38)/9)%%1, ((t-6.38)/9)%%1, 
                  ((t-7.38)/9)%%1, ((t-8.38)/9)%%1, (t-17.38),
                  (t-64.38),        (t-65.38),      (t-66.38),        (t-67.38),
                  (t-68.38),        (t-69.38),      (t-70.38),        (t-71.38),
                  (t-126.38),       (t-127.38), 
                  (t-123.38),       (t-124.38),     (t-125.38),                       # egg laying events
                  
                  (t-dt)%%1,     (t-0.04)%%1,           (t-0.08)%%1,     (t-0.12)%%1, 
                  (t-0.16)%%1,   (t-0.21)%%1,           (t-0.25)%%1,     (t-0.29)%%1, 
                  (t-0.33)%%1,   (t-0.37)%%1,           (t-0.42)%%1,     (t-0.46)%%1, 
                  (t-0.5)%%1,    (t-0.54)%%1,           (t-0.58)%%1,     (t-0.62)%%1, 
                  (t-0.67)%%1,   (t-0.71)%%1,           (t-0.75)%%1,     (t-0.79)%%1, 
                  (t-0.83)%%1,   (t-0.87)%%1,           (t-0.92)%%1,     (t-0.96)%%1,  # feeding events
                  (t-1.67),      
                  (t-8.29),              (t-8.37),        (t-8.62),
                  (t-9.29),              (t-9.37),        (t-9.62),
                  (t-10.29),              (t-10.37),        (t-10.62),
                  (t-11.29),              (t-11.37),        (t-11.62),
                  (t-12.29),              (t-12.37),        (t-12.62),
                  (t-13.29),              (t-13.37),        (t-13.62),
                  (t-14.29),              (t-14.37),        (t-14.62),
                  (t-15.29),              (t-15.37),        (t-15.62),
                  (t-16.33),              (t-16.58),
                  (t-17.33),              (t-17.58),
                  (t-18.33),              (t-18.58),
                  (t-19.33),              (t-19.58),
                  (t-20.33),              (t-20.58),
                  (t-21.33),              (t-21.58),
                  (t-22.33),              (t-22.58),
                  (t-23.33),              (t-23.58),
                  (t-24.33),              (t-24.58),
                  (t-25.33),              (t-25.58),
                  (t-26.33),              (t-26.58),
                  (t-27.33),              (t-27.58),
                  (t-9)
                  
                  )
  return(yroot)
}

event_repeated <- function(t, y, parms) {
  ret     <- y
  # SEQUENCE MATTERS
  # dose event
  if(abs((t-dt)%%1) < 1e-6)          ret <- event_dose0(t, y, parms)
  if(abs((t-0.04)%%1) < 1e-6)        ret <- event_dose1(t, y, parms)
  if(abs((t-0.08)%%1) < 1e-6)        ret <- event_dose2(t, y, parms)
  if(abs((t-0.12)%%1) < 1e-6)        ret <- event_dose3(t, y, parms)
  if(abs((t-0.16)%%1) < 1e-6)        ret <- event_dose4(t, y, parms)
  if(abs((t-0.21)%%1) < 1e-6)        ret <- event_dose5(t, y, parms)
  if(abs((t-0.25)%%1) < 1e-6)        ret <- event_dose6(t, y, parms)
  if(abs((t-0.29)%%1) < 1e-6)        ret <- event_dose7(t, y, parms)
  if(abs((t-0.33)%%1) < 1e-6)        ret <- event_dose8(t, y, parms)
  if(abs((t-0.37)%%1) < 1e-6)        ret <- event_dose9(t, y, parms)
  if(abs((t-0.42)%%1) < 1e-6)        ret <- event_dose10(t, y, parms)
  if(abs((t-0.46)%%1) < 1e-6)        ret <- event_dose11(t, y, parms)
  if(abs((t-0.50)%%1) < 1e-6)        ret <- event_dose12(t, y, parms)
  if(abs((t-0.54)%%1) < 1e-6)        ret <- event_dose13(t, y, parms)
  if(abs((t-0.58)%%1) < 1e-6)        ret <- event_dose14(t, y, parms)
  if(abs((t-0.62)%%1) < 1e-6)        ret <- event_dose15(t, y, parms)
  if(abs((t-0.67)%%1) < 1e-6)        ret <- event_dose16(t, y, parms)
  if(abs((t-0.71)%%1) < 1e-6)        ret <- event_dose17(t, y, parms)
  if(abs((t-0.75)%%1) < 1e-6)        ret <- event_dose18(t, y, parms)
  if(abs((t-0.79)%%1) < 1e-6)        ret <- event_dose19(t, y, parms)
  if(abs((t-0.83)%%1) < 1e-6)        ret <- event_dose20(t, y, parms)
  if(abs((t-0.87)%%1) < 1e-6)        ret <- event_dose21(t, y, parms)
  if(abs((t-0.92)%%1) < 1e-6)        ret <- event_dose22(t, y, parms)
  if(abs((t-0.96)%%1) < 1e-6)        ret <- event_dose23(t, y, parms)
  
  if(abs(t-1.67) < 1e-6)             ret <- event_dose16(t, y, parms)
  if(abs(t-8.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-8.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-8.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-9.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-9.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-9.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-10.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-10.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-10.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-11.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-11.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-11.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-12.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-12.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-12.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-13.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-13.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-13.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-14.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-14.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-14.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-15.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-15.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-15.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-16.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-16.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-17.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-17.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-18.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-18.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-19.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-19.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-20.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-20.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-21.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-21.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-22.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-22.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-23.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-23.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-24.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-24.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-25.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-25.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-26.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-26.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-27.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-27.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if((t-9) < 1e-6)                    ret <- event_NoDose(t, y, parms)
  
  # egg laying event
  if((abs(t-0.38)%%1)     < 1e-6)  ret <- eventW(t, y, parms)
  if(abs(((t-0.38)/9)%%1) < 1e-6)  ret <- event1(t, y, parms)
  if(abs(((t-1.38)/9)%%1) < 1e-6)  ret <- event2(t, y, parms)
  if(abs(((t-2.38)/9)%%1) < 1e-6)  ret <- event3(t, y, parms)
  if(abs(((t-3.38)/9)%%1) < 1e-6)  ret <- event4(t, y, parms)
  if(abs(((t-4.38)/9)%%1) < 1e-6)  ret <- event5(t, y, parms)
  if(abs(((t-5.38)/9)%%1) < 1e-6)  ret <- event6(t, y, parms)
  if(abs(((t-6.38)/9)%%1) < 1e-6)  ret <- event7(t, y, parms)
  if(abs(((t-7.38)/9)%%1) < 1e-6)  ret <- event8(t, y, parms)
  if(abs(((t-8.38)/9)%%1) < 1e-6)  ret <- event9(t, y, parms)
  if(abs(t-17.38) < 1e-6)          ret <- event9(t, y, parms)
  if(abs(t-64.38) < 1e-6)          ret <- event2(t, y, parms)
  if(abs(t-65.38) < 1e-6)          ret <- event3(t, y, parms)
  if(abs(t-66.38) < 1e-6)          ret <- event4(t, y, parms)
  if(abs(t-67.38) < 1e-6)          ret <- event5(t, y, parms)
  if(abs(t-68.38) < 1e-6)          ret <- event6(t, y, parms)
  if(abs(t-69.38) < 1e-6)          ret <- event7(t, y, parms)
  if(abs(t-70.38) < 1e-6)          ret <- event8(t, y, parms)
  if(abs(t-71.38) < 1e-6)          ret <- event9(t, y, parms)
  
  if(abs(t-126.38) < 1e-6)         ret <- event1(t, y, parms)
  if(abs(t-127.38) < 1e-6)         ret <- event2(t, y, parms)
  if(abs(t-123.38) < 1e-6)         ret <- event7(t, y, parms)
  if(abs(t-124.38) < 1e-6)         ret <- event8(t, y, parms)
  if(abs(t-125.38) < 1e-6)         ret <- event9(t, y, parms)   # egg laying events above
  return((ret))
}


#########################    parameters for exposure scenario ###############
factor <- 1


parms <- c( dt            = dt,
            ka            = ka,                 # h-1
            fa            = fa,  
            Clint_parent         = Clint_hep_parent   * 37, # 50
            Clint_daughter       = Clint_hep_daughter / 3,   # /4
            
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
            
            tsig          = tsig,
            tlag          = tlag,
            ky_parent     = ky_parent, 
            ky_daughter   = ky_daughter,
            kw_parent     = kw_parent,
            kw_daughter   = kw_daughter,
            
            # Parametes for flow 
            Qart          = Qart, 
            Qgfr_parent   = Qgfr_parent * 1.4 ,# * 2.3,   
            Qgfr_daughter = Qgfr_daughter * 5.9,
            
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
            Vadipose      = Vadipose)



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

ggplot() +
  geom_line(data = df, aes(Time, Agutlumen_parent), col="#00AFBB", lwd=1) +
  xlim(0,28) + 
  theme(text = element_text(size = 20)) #+ xlim(3, 4)


ggplot() +
  geom_line(data = df, aes(Time, C_blood_parent), col="#00AFBB", lwd=1) +
  xlim(0, 28) + #geom_point(data = df_Metabolism_O_1998_blood  , aes(Time_h/24+0.3, TMXBloodConc_ugL),col="red",  size=2.5) + ylab("Concentration")+
  theme(text = element_text(size = 20)) #+ xlim(13, 14)

ggplot() +
  geom_line(data = df, aes(Time, C_blood_daughter), col="#00AFBB", lwd=2) +
  ylab("Concentration")+
  theme(text = element_text(size = 20)) #+ xlim(3, 4)

ggplot() +
  geom_line(data = df, aes(Time, C_liver_parent), col="#00AFBB", lwd=1) +
  xlim(0, 28) + #geom_point(data = df_Metabolism_O_1998_blood  , aes(Time_h/24+0.3, TMXBloodConc_ugL),col="red",  size=2.5) + ylab("Concentration")+
  theme(text = element_text(size = 20)) + xlim(27, 28)

ggplot() +
  geom_line(data = df, aes(Time, C_liver_daughter), col="#00AFBB", lwd=2) +
  ylab("Concentration")+
  theme(text = element_text(size = 20)) #+ xlim(3, 4)


#=============================================================================
############                End of PBPK model                    #############
#=============================================================================

######           Plot           ###########
Metabolism_O_1998_TMX           <- read_excel("C:xxx/OneDrive - Syngenta/HTTK/Bird/TMX/InvivoData.xlsx", sheet = "Metabolism_O_1998_hen")


# egg white
Metabolism_O_1998_egg_white               <- Metabolism_O_1998_TMX[Metabolism_O_1998_TMX$Matrix == 'Egg_white',]
Metabolism_O_1998_egg_white$TMXConc_ugL   <- Metabolism_O_1998_egg_white$Concentration * 1000 * Metabolism_O_1998_egg_white$'%TMX' / 100
Metabolism_O_1998_egg_white$CLOConc_ugL   <- Metabolism_O_1998_egg_white$Concentration * 1000 * Metabolism_O_1998_egg_white$'%CGA322704' / 100
df_Metabolism_O_1998_eggwhite             <- subset(Metabolism_O_1998_egg_white, select = c('Time_h_range', 'Time_h', "TMXConc_ugL", "CLOConc_ugL", "Weight_g", "Matrix"))



########### Egg yolk concentration ###########
df_yolk_out         <- subset(df, (Time - 0.37)%%1 < 1e-4 | (Time - 0.38)%%1 < 1e-4)
df_yolk_out2         <- subset(df_yolk_out, (Time - 9) > 1e-4)
df_yolk_out_parent   <- df_yolk_out2[c('Time', 'Ayolk1_parent', 'Ayolk2_parent', 'Ayolk3_parent', 'Ayolk4_parent', 
                                            'Ayolk5_parent', 'Ayolk6_parent', 'Ayolk7_parent', 'Ayolk8_parent', 'Ayolk9_parent',
                                            'Cyolk1_parent', 'Cyolk2_parent', 'Cyolk3_parent', 'Cyolk4_parent', 'Cyolk5_parent', 
                                            'Cyolk6_parent', 'Cyolk7_parent', 'Cyolk8_parent', 'Cyolk9_parent')]
df_Cyolk_out_parent  <- df_yolk_out_parent[c('Time', 'Cyolk1_parent', 'Cyolk2_parent', 'Cyolk3_parent', 'Cyolk4_parent', 'Cyolk5_parent', 
                                             'Cyolk6_parent', 'Cyolk7_parent', 'Cyolk8_parent', 'Cyolk9_parent')]

rownames(df_Cyolk_out_parent) <- 1:nrow(df_Cyolk_out_parent)
mdata                         <- melt(df_Cyolk_out_parent, id=c("Time"))
mdata2                        <- mdata[which(mdata$value == 0)-1,] 
ggplot() +
  geom_point(data = mdata2, aes(Time-9, value, col=variable), lwd=5) + ylab("Mass (kg)") +
  #geom_point(data = df_Metabolism_O_1998_yolk , aes(Time_h/24+ 0.3, TMXConc_ugL),size=2.5) + ylab("Concentration")+
  theme(text = element_text(size = 20))#+ xlim(0, 4)


df_yolk_out_daughter   <- df_yolk_out2[c('Time', 'Ayolk1_daughter', 'Ayolk2_daughter', 'Ayolk3_daughter', 'Ayolk4_daughter', 
                                      'Ayolk5_daughter', 'Ayolk6_daughter', 'Ayolk7_daughter', 'Ayolk8_daughter', 'Ayolk9_daughter',
                                      'Cyolk1_daughter', 'Cyolk2_daughter', 'Cyolk3_daughter', 'Cyolk4_daughter', 'Cyolk5_daughter', 
                                      'Cyolk6_daughter', 'Cyolk7_daughter', 'Cyolk8_daughter', 'Cyolk9_daughter')]
df_Cyolk_out_daughter  <- df_yolk_out_daughter[c('Time', 'Cyolk1_daughter', 'Cyolk2_daughter', 'Cyolk3_daughter', 'Cyolk4_daughter', 'Cyolk5_daughter', 
                                             'Cyolk6_daughter', 'Cyolk7_daughter', 'Cyolk8_daughter', 'Cyolk9_daughter')]

rownames(df_Cyolk_out_daughter) <- 1:nrow(df_Cyolk_out_daughter)
mdata_daughter                         <- melt(df_Cyolk_out_daughter, id=c("Time"))
mdata2_daughter                        <- mdata_daughter[which(mdata_daughter$value == 0)-1,] 
ggplot() +
  geom_point(data = mdata2_daughter, aes(Time-9, value, col=variable), lwd=5) + ylab("Mass (kg)") +
  #geom_point(data = df_Metabolism_O_1998_yolk , aes(Time_h/24 + 0.3, CLOConc_ugL),size=2.5) + ylab("Concentration")+
  theme(text = element_text(size = 20))#+ xlim(0, 4)



############  Egg concentration  ############
df_egg_out1         <- subset(df, (Time - 0.37)%%1 < 1e-4 | (Time - 0.38)%%1 < 1e-4)
df_egg_out          <- subset(df_egg_out1, (Time - 9) > 1e-4)
df_egg_out_parent   <- df_egg_out[c('Time', 'Cegg1_parent', 'Cegg2_parent', 'Cegg3_parent', 'Cegg4_parent', 'Cegg5_parent', 
                                      'Cegg6_parent', 'Cegg7_parent', 'Cegg8_parent', 'Cegg9_parent')]
df_Cegg_out_parent  <- df_egg_out_parent[c('Time', 'Cegg1_parent', 'Cegg2_parent', 'Cegg3_parent', 'Cegg4_parent', 'Cegg5_parent', 
                                             'Cegg6_parent', 'Cegg7_parent', 'Cegg8_parent', 'Cegg9_parent')]

rownames(df_Cegg_out_parent) <- 1:nrow(df_Cegg_out_parent)
mdata                         <- melt(df_Cegg_out_parent, id=c("Time"))
egg_8                         <- mdata[mdata$Time == 16.37 & mdata$variable == 'Cegg8_parent',]
mdata2                        <- mdata[which(mdata$value == 0)-1,] 
mdata2                        <- rbind(mdata2, egg_8)
ggplot() +
  geom_point(data = mdata2, aes(Time-9, value * MW_parent/1000, col=variable), lwd=5) + ylab("Concentration (ppm)") +
  #geom_point(data = df_Metabolism_O_1998_egg , aes(Time_h/24+ 0.3, TMXConc_ugL),size=2.5) + ylab("Concentration")+
  geom_hline(yintercept = 0.01, linetype="dashed", color="red", lwd = 2) +
  theme(text = element_text(size = 20))#+ xlim(0, 4)

## plot for parent
# Create a new data frame for the hline
hline_data <- data.frame(x = range(mdata2$Time-9), y = 0.01)
ggplot() +
  geom_point(data = mdata2, aes(Time-9, value * MW_parent/1000, shape = "Prediction"), colour = '#00AFBB', size=2) + 
  geom_line(data = hline_data, aes(x, y, color = "Observation"), linetype="dashed", lwd = 1.2) + 
  theme(text = element_text(size = 20)) +
  xlab("Time (d)") + 
  ylab(expression("TMX Egg (whites+yolks) Concentration ("*mu*"g/mL)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#00AFBB', "#FC4E07"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(16)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [ad libitum daily for 28 consecutive days at 10 ppm/d]", sep = ''),
    subtitle = paste("Thiamethoxam (TMX)"))
ggsave(filename = paste("C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/plot.layinghen2000.10ppm.TMXEgg.tiff", sep = ''),
       dpi = 300, width = 16, height = 10,units = "cm",compression = 'lzw')



df_egg_out_daughter   <- df_egg_out[c('Time', 'Cegg1_daughter', 'Cegg2_daughter', 'Cegg3_daughter', 'Cegg4_daughter', 'Cegg5_daughter', 
                                        'Cegg6_daughter', 'Cegg7_daughter', 'Cegg8_daughter', 'Cegg9_daughter')]
df_Cegg_out_daughter  <- df_egg_out_daughter[c('Time', 'Cegg1_daughter', 'Cegg2_daughter', 'Cegg3_daughter', 'Cegg4_daughter', 'Cegg5_daughter', 
                                                 'Cegg6_daughter', 'Cegg7_daughter', 'Cegg8_daughter', 'Cegg9_daughter')]

rownames(df_Cegg_out_daughter) <- 1:nrow(df_Cegg_out_daughter)
mdata_daughter                         <- melt(df_Cegg_out_daughter, id=c("Time"))
egg_8_daughter                         <- mdata_daughter[mdata_daughter$Time == 16.37 & mdata_daughter$variable == 'Cegg8_daughter',]
mdata2_daughter                        <- mdata_daughter[which(mdata_daughter$value == 0)-1,] 
mdata2_daughter                        <- rbind(mdata2_daughter, egg_8_daughter)
ggplot() +
  geom_point(data = mdata2_daughter, aes(Time-9, value * MW_daughter/1000, col=variable), lwd=5) + ylab("Concentration (ppm)") +
  #geom_point(data = df_Metabolism_O_1998_egg , aes(Time_h/24 + 0.3, CLOConc_ugL),size=2.5) + ylab("Concentration")+
  geom_hline(yintercept = 0.01, linetype="dashed", color="red", lwd = 2) +
  theme(text = element_text(size = 20))#+ xlim(0, 4)

## plot for daughter
# Create a new data frame for the hline
hline_data <- data.frame(x = range(mdata2$Time-9), y = 0.01)
ggplot() +
  geom_point(data = mdata2_daughter, aes(Time-9, value * MW_daughter/1000, shape = "Prediction"), colour = '#276DC2', size=2) + 
  geom_line(data = hline_data, aes(x, y, color = "Observation"), linetype="dashed", lwd = 1.2) + 
  theme(text = element_text(size = 20)) +
  xlab("Time (d)") + 
  ylab(expression("CTD Egg (whites+yolks) Concentration ("*mu*"g/mL)")) +
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 12))  + 
  # legend
  scale_color_manual(name = NULL,
                     values = c('#276DC2', "#d8923d"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(16)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = 'right',
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Laying hen [ad libitum daily for 28 consecutive days at 10 ppm/d]", sep = ''),
    subtitle = paste("Clothianidin (CTD)"))
ggsave(filename = paste("C:/Uxxx/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/plot.layinghen2000.10ppm.CTDEgg.tiff", sep = ''),
       dpi = 300, width = 16, height = 10,units = "cm",compression = 'lzw')




ggplot() +
  geom_line (data = df, aes(Time,Ww), col="#00AFBB", lwd=2) + ylab("Mass (kg)") +
  theme(text = element_text(size = 20))#+ xlim(0, 10)

ggplot() +
  geom_line (data = df, aes(Time,Agutlumen_parent), col="#00AFBB", lwd=2) + ylab("Mass (kg)") +
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



  

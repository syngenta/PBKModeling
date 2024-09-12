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
#==================================================================================================
#                                      Species: Quail                                             # 
#==================================================================================================
# Reference
# PKNCA: https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk: https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk

#=============================================================================================================
#                                   Generic Maternal PBPK Model Parameters                                   #
#                                  Tissue and physiological data for species                                 #
#=============================================================================================================
rm(list = ls())
SPECIES        <- 'Quail'
species        <- 'Rat'                  # (either "Rat", "Rabbit", "Dog", "Mouse", or "Human").
BW             <- 0.2                    # kg, body weight of 18 weeks 
BW_hen         <- 1.45                   # kg
SDBW           <- 0.08

current_dir <- getwd()
print(current_dir)
source(file.path("HTTK", "Bird", "General code", "Avian_Noegg.R"))
source(file.path("HTTK", "Bird", "General code", paste(SPECIES, "_PhyData.R", sep = '')))
source(file.path("HTTK", "Bird", "General code", 'Metabolism.R'))


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


######################        Clothianidin        ##########################
MW_daughter             <- 249.7
CAS_daughter            <- '210880-92-5a'
S_daughter              <- 340                    # mg/L; http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/631.htm
pKa_a_daughter          <- 11.1                   # No dissociation; http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/171.htm; pKa_Donor="pka.a") Compound H dissociation equilibirum constant
pKa_b_daughter          <- 14-11.1                  # pka accept (base) Compound H association equilibirum constant
fub_daughter            <- 0.71                   # fraction unbound in plasma of rat; from SED
LogP_daughter           <- 0.702                  # LogP httk
Peff_daughter           <- 19.8E-6                # Cacao permeability (cm/s) from SED for human
Compound_type_daughter  <- 'base'
Compound_name_daughter  <- 'Clothianidin (CLO)'
Density.daughter        <- 1.61                   # g/cm3 # http://sitem.herts.ac.uk/aeru/ppdb/en/Reports/631.htm
Rblood2plasma_daughter  <-'None' 


## Permeability-surface area product L/d kg; not used
PSc_adipose_parent      <- Peff_parent * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000            # Peff * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000  # cm/s -> cm/d * g/kg * cm2/g -> cm3/d kg -> L/d kg
PSc_adipose_daughter    <- Peff_daughter * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000          # Peff * 3600 * 24 * (Vadipose * 1000) * S_adipose /1000  # cm/s -> cm/d * g/kg * cm2/g -> cm3/d kg -> L/d kg

fa        <- 0.9 


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

source(file.path("HTTK", "Bird", "General code", "Partition.R"))


  
########################################################################
###                    define modeling variables                     ###  
########################################################################
days          <- 28
Time_max      <- days                                     # 5 d  

StartTime     <- 0                                        # Time_d
StopTime      <- Time_max+28                                 
dt            <- 0.01
Times         <- seq(StartTime, StopTime, dt)
#==============================================
dose_ppm      <- 320   # 320, 1000, 2000 or 3200
#==============================================
# dose pattern ; \HTTK\Bird\Intake pattern\pattern.xlsx
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

dose_320          <- c(30.9, 
                       21.9,
                       27.6,
                       17,
                       13.9,
                       22.3,
                       20.3,
                       
                       23.7,
                       23.7,
                       23.7,
                       23.7,
                       23.7,
                       23.7,
                       23.7,
                       
                       25.6,
                       25.6,
                       25.6,
                       25.6,
                       25.6,
                       25.6,
                       25.6,
                       
                       26.1,
                       26.1,
                       26.1,
                       26.1,
                       26.1,
                       26.1,
                       26.1) * 1000    # ug/kg/d for 28 days; Table 1 female

dose_1000         <- c(89.4, 
                       56.5,
                       53.4,
                       57.5,
                       63.7,
                       71.9,
                       77.1,
                       
                       88.2,
                       88.2,
                       88.2,
                       88.2,
                       88.2,
                       88.2,
                       88.2,
                       
                       87.5,
                       87.5,
                       87.5,
                       87.5,
                       87.5,
                       87.5,
                       87.5,
                       
                       72.4,
                       72.4,
                       72.4,
                       72.4,
                       72.4,
                       72.4,
                       72.4) * 1000   # ug/kg/d for 28 days

dose_2000         <- c(105, 
                       88.6,
                       97,
                       112,
                       108,
                       131,
                       112,
                       
                       144,
                       144,
                       144,
                       144,
                       144,
                       144,
                       144,
                       
                       160,
                       160,
                       160,
                       160,
                       160,
                       160,
                       160,
                       
                       133,
                       133,
                       133,
                       133,
                       133,
                       133,
                       133) * 1000   # ug/kg/d for 28 days

dose_3200         <- c(135, 
                       104,
                       188,
                       158,
                       96.7,
                       192,
                       141,
                       
                       243,
                       243,
                       243,
                       243,
                       243,
                       243,
                       243,
                       
                       257,
                       257,
                       257,
                       257,
                       257,
                       257,
                       257,
                       
                       256,
                       256,
                       256,
                       256,
                       256,
                       256,
                       256) * 1000   # ug/kg/d for 28 days


dose_28d    <- matrix(nrow = days, ncol = length(value_1))

for (i in 1:days){
  
  cat("day = ", i, '\n')
  if(dose_ppm == 320){dose  <- dose_320}
  if(dose_ppm == 1000){dose <- dose_1000}
  if(dose_ppm == 2000){dose <- dose_2000}
  if(dose_ppm == 3200){dose <- dose_3200}
  
  dose_28d[i,]   <- value_1 * dose[i] / MW_parent
}

df_dose_28d  <- as.vector(t(dose_28d))


Dose_events    <- data.frame(var    = rep("Agutlumen_parent", days * 24),
                             time   =  time,
                             value  = df_dose_28d,
                             method = "add")

                

########################################################################################################################
#                                                  PBPK MODEL RUN                                                      #
########################################################################################################################
initState <- c(Agutlumen_parent  = 0, 
               Agut_parent       = 0,  
               Aliver_parent     = 0,
               Akidney_parent    = 0,
               Aadipose_parent = 0,
               #A_V_adipose_parent = 0,
               #A_T_adipose_parent = 0,
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
               #A_V_adipose_daughter = 0,
               #A_T_adipose_daughter = 0,
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
parms <- c( ka                   = 12.1847922,                 # d-1
            fa                   = fa,  
            Clint_parent         = 692.6002544, 
            Clint_daughter       =  79.2898978, 
            
            Rblood2plasma_parent   = 0.9027282,   
            Rblood2plasma_daughter = 1.3597193,   
            BW = BW,              
            
            fub_parent           = fub_parent,   
            fub_daughter         = fub_daughter, 
            fraction_daughter    = 1,
            # ratio
            Kgut2pu_parent       = Kgut2pu_parent    ,      
            Kkidney2pu_parent    = Kkidney2pu_parent ,   
            Kliver2pu_parent     = Kliver2pu_parent,     
            Klung2pu_parent      = Klung2pu_parent   ,      
            Krest2pu_parent      = Krest2pu_parent   , 
            Kadipose2pu_parent   = Kadipose2pu_parent ,
            Kmuscle2pu_parent    = Kmuscle2pu_parent, 
            Krep2pu_parent       = Krep2pu_parent    ,
            
            Kgut2pu_daughter       = Kgut2pu_daughter    ,      
            Kkidney2pu_daughter    = Kkidney2pu_daughter, 
            Kliver2pu_daughter     = Kliver2pu_daughter,     
            Klung2pu_daughter      = Klung2pu_daughter   ,      
            Krest2pu_daughter      = Krest2pu_daughter   , 
            Kadipose2pu_daughter   = Kadipose2pu_daughter ,
            Kmuscle2pu_daughter    = Kmuscle2pu_daughter, 
            Krep2pu_daughter       = Krep2pu_daughter    ,

            # Parametes for flow 
            Qart          = Qart, 
            Qgfr_parent   = 6.7386862,    
            Qgfr_daughter = 15.3186752 , 
            
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

file_path    <- file.path(current_dir, "HTTK", "Bird",'TMX', 'Plots', 'Quail')
write.csv(df, paste(file_path, '/quail.adlib.',  dose_ppm, 'ppm.csv', sep = ''), row.names = FALSE)

#=============================================================================
############                End of PBPK model                    #############
#=============================================================================
# AUC of blood (chronic, adlib dose)
AUC24_parent     <- df[df$Time == 5,]$AUC_Cblood_parent     - df[df$Time == 4,]$AUC_Cblood_parent 
AUC24_daughter   <- df[df$Time == 5,]$AUC_Cblood_daughter   - df[df$Time == 4,]$AUC_Cblood_daughter 
AUC24_total      <- AUC24_parent + AUC24_daughter            

Cmax_parent      <- max(df$C_blood_parent)            
Cmax_daughter    <- max(df$C_blood_daughter)        
Cmax_total       <- max(df$C_blood_total)            

stat  <- data.frame('Name'  = c('AUC24_parent', 'AUC24_daughter', 'AUC24_total', 'Cmax_parent', 'Cmax_daughter', 'Cmax_total' ),
                    'Value' = c(AUC24_parent, AUC24_daughter, AUC24_total, Cmax_parent, Cmax_daughter, Cmax_total))
stat


##########################      pbpk output     #############################  
######           Plot           ###########
file_path                   <- file.path(current_dir, "HTTK", "Bird", "TMX", 'InvivoData.xlsx')
Residue_2016_quail          <- read_excel(file_path, sheet = "Residue_2016_quail")

# blood
Pre_egg_female                  <- Residue_2016_quail[Residue_2016_quail$Dose_ppm == dose_ppm & Residue_2016_quail$group == 'pre-egg' & Residue_2016_quail$Gender == 'female',]
Pre_egg_female$TMXBloodConc_umol <- Pre_egg_female$Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female$CLOBloodConc_umol <- Pre_egg_female$Blood_conc_ng.ml_CGA322704 /MW_daughter
Pre_egg_female$TMXBloodConc_umol_SD <- Pre_egg_female$SD.Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female$CLOBloodConc_umol_SD <- Pre_egg_female$SD.Blood_conc_ng.ml_CGA322704 /MW_daughter
df_Pre_egg_female               <- subset(Pre_egg_female, select = c('Time_d', "TMXBloodConc_umol", 'TMXBloodConc_umol_SD',
                                                                           'CLOBloodConc_umol', "CLOBloodConc_umol_SD")) 
df_Pre_egg_female   <- na.omit(df_Pre_egg_female)

##########   OVER VIREW PLOT  ################
ggplot() +
  geom_line(data = df, aes(Time, C_blood_parent, color = "Prediction"),  lwd=0.3) +
  geom_errorbar(data = df_Pre_egg_female , aes(x=(Time_d - 14/24), 
                                                  ymin= pmax(TMXBloodConc_umol  - TMXBloodConc_umol_SD, 0) , ymax = TMXBloodConc_umol  + TMXBloodConc_umol_SD) , 
                size = 0.6, width=0.52, position=position_dodge(0.05),colour="#ff5044") + 
  geom_point(data = df_Pre_egg_female , aes((Time_d - 14/24), TMXBloodConc_umol, shape = 'Observation'),col="#ff5044",  size=1.5) + 
  xlab("Time (d)") + 
  ylab(expression("TMX Blood Concentration ("*mu*"mol/L)")) + 
  theme_bw() + 
  theme(plot.title = element_text(size=12), plot.subtitle = element_text(size=10))+
  theme(text = element_text(size = 12))  + xlim(0,35) +
  # legend
  scale_color_manual(name = NULL,
                     values = c('#00AFBB', "#ff5044"),
                     breaks=c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(16)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = c(0.9, 0.88),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Northern bobwhite [ad libitum daily for 28 consecutive days at ", dose_ppm, " ppm/d]", sep = ''),
    subtitle = paste("Thiamethoxam (TMX)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Quail', paste("plot.bobwhite.", dose_ppm, ".TMXBlood_35d.tiff", sep = ''))
ggsave(file_path,
       dpi = 300, width = 18, height = 10,units = "cm",compression = 'lzw')

ggplot() +
  geom_line(data  = df, aes(Time, C_blood_daughter, color = "Prediction"),  lwd=0.3) +
  geom_errorbar(data = df_Pre_egg_female , aes(x=(Time_d - 14/24), 
                                                  ymin= pmax(CLOBloodConc_umol - CLOBloodConc_umol_SD, 0), ymax = CLOBloodConc_umol + CLOBloodConc_umol_SD), 
                size = 0.6,width=0.52, position=position_dodge(0.05),colour="#d8923d") + 
  geom_point(data = df_Pre_egg_female , aes((Time_d - 14/24), CLOBloodConc_umol, shape = 'Observation'),col="#d8923d", size=1.5) +
  xlab("Time (d)") + 
  ylab(expression("CTD Blood Concentration ("*mu*"mol/L)")) + 
  theme_bw() + 
  theme(plot.title = element_text(size=12), plot.subtitle = element_text(size=10))+
  theme(text = element_text(size = 12))  + xlim(0,35) +
  # legend
  scale_color_manual(name = NULL,
                     values = c('#276DC2', "#d8923d"),
                     breaks = c("Prediction", "Observation")) +
  scale_shape_manual(name=NULL, values=c(16)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = c(0.9, 0.88),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm'),
        legend.text=element_text(size=8)) +
  # Titles and Labels
  ggtitle(
    label = paste("Northern bobwhite [ad libitum daily for 28 consecutive days at ", dose_ppm, " ppm/d]", sep = ''),
    subtitle = paste("Clothianidin (CTD)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Quail', paste("plot.bobwhite.", dose_ppm, ".CTDBlood_35d.tiff", sep = ''))
ggsave(file_path,
       dpi = 300, width = 18, height = 10,units = "cm",compression = 'lzw')


##########   7d PLOT  ################
ggplot() +
  geom_line(data = df, aes(Time, C_blood_parent, color = "Prediction"), lwd=0.3) +
  geom_errorbar(data = df_Pre_egg_female , aes(x=(Time_d - 14/24), 
                                               ymin= pmax(TMXBloodConc_umol  - TMXBloodConc_umol_SD, 0), ymax = TMXBloodConc_umol  + TMXBloodConc_umol_SD) , 
                size = 0.6, width=0.2, position=position_dodge(0.05),colour="#ff5044") + 
  geom_point(data = df_Pre_egg_female , aes((Time_d - 14/24), TMXBloodConc_umol, shape = 'Observation'),col="#ff5044",  size=1.5) + 
  xlab("Time (d)") + 
  ylab(expression("TMX Blood Concentration ("*mu*"mol/L)")) + 
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 10))  + xlim(0,7)+
  # legend
  scale_color_manual(name = NULL,
                     values = c('#00AFBB', "#ff5044"),
                     breaks = c("Prediction", "Observation")) +
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
    label = paste("Northern bobwhite [ad libitum daily for 28 consecutive days at ", dose_ppm, " ppm/d]", sep = ''),
    subtitle = paste("Thiamethoxam (TMX)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Quail', paste("plot.bobwhite.", dose_ppm, ".TMXBlood_7d.tiff", sep = ''))
ggsave(file_path,
       dpi = 300, width = 16, height = 10,units = "cm",compression = 'lzw')

ggplot() +
  geom_line(data  = df, aes(Time, C_blood_daughter, color = "Prediction"), lwd=0.3) +
  geom_errorbar(data = df_Pre_egg_female , aes(x     = (Time_d - 14/24), 
                                               ymin  = pmax(CLOBloodConc_umol - CLOBloodConc_umol_SD, 0), ymax = CLOBloodConc_umol + CLOBloodConc_umol_SD), 
                size = 0.6,width=0.2, position=position_dodge(0.05),colour="#d8923d") + 
  geom_point(data = df_Pre_egg_female , aes((Time_d - 14/24), CLOBloodConc_umol, shape = 'Observation'),col="#d8923d", size=1.5) +
  xlab("Time (d)") + 
  ylab(expression("CTD Blood Concentration ("*mu*"mol/L)")) + 
  theme_bw() + 
  theme(plot.title = element_text(size=9), plot.subtitle = element_text(size=7))+
  theme(text = element_text(size = 10))  + xlim(0, 7)+
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
    label = paste("Northern bobwhite [ad libitum daily for 28 consecutive days at ", dose_ppm, " ppm/d]", sep = ''),
    subtitle = paste("Clothianidin (CTD)"))
file_path           <- file.path(current_dir, "HTTK", "Bird", "TMX", "Plots", 'Quail', paste("plot.bobwhite.", dose_ppm, ".CTDBlood_7d.tiff", sep = ''))  
ggsave(file_path,
       dpi = 300, width = 16, height = 10,units = "cm",compression = 'lzw')



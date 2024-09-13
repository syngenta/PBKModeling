library(stationaRy)
library(isdparser)
library(tidyr)
library(lubridate)
library(WriteXLS)
library(readxl)
library(writexl)
library(FME)
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
library(minpack.lm)

current_dir <- getwd()
print(current_dir)
source(file.path("HTTK", "Bird", "TMX", 'Quail', "Quail_Single_oral_NoEgg.R"))

#################################################################
##########            calibration dataset             ###########
#################################################################
file_path                   <- file.path("HTTK", "Bird", "TMX", 'InvivoData.xlsx')
Residue_2016_quail          <- read_excel(file_path, sheet = "Residue_2016_quail")

# blood
Pre_egg_female_320                  <- Residue_2016_quail[Residue_2016_quail$Dose_ppm == 320 & Residue_2016_quail$group == 'pre-egg' & Residue_2016_quail$Gender == 'female',]
Pre_egg_female_320$Time             <- round(Pre_egg_female_320$Time_d - 14/24, 2)
Pre_egg_female_320$C_blood_parent   <- round(Pre_egg_female_320$Blood_conc_ng.ml_TMX /MW_parent, 4)
Pre_egg_female_320$C_blood_daughter <- round(Pre_egg_female_320$Blood_conc_ng.ml_CGA322704 /MW_daughter, 4)
df_Pre_egg_female_320                   <-  as.data.frame(Pre_egg_female_320 %>% dplyr::select(Time, C_blood_parent, C_blood_daughter))
df_Pre_egg_female_320                   <- na.omit(df_Pre_egg_female_320)

Pre_egg_female_1000                  <- Residue_2016_quail[Residue_2016_quail$Dose_ppm == 1000 & Residue_2016_quail$group == 'pre-egg' & Residue_2016_quail$Gender == 'female',]
Pre_egg_female_1000$Time             <- round(Pre_egg_female_1000$Time_d - 14/24, 2)
Pre_egg_female_1000$C_blood_parent   <- round(Pre_egg_female_1000$Blood_conc_ng.ml_TMX /MW_parent, 4)
Pre_egg_female_1000$C_blood_daughter <- round(Pre_egg_female_1000$Blood_conc_ng.ml_CGA322704 /MW_daughter, 4)
df_Pre_egg_female_1000                   <-  as.data.frame(Pre_egg_female_1000 %>% dplyr::select(Time, C_blood_parent, C_blood_daughter))
df_Pre_egg_female_1000                   <- na.omit(df_Pre_egg_female_1000)



########################################################################
#                            Calibration                               #
########################################################################
days          <- 28
Time_max      <- days                                     # 5 d  

StartTime     <- 0                                        # Time_d
StopTime      <- Time_max + 28                                 
dt            <- 0.01
Times         <- seq(StartTime, StopTime, dt)
time          <- seq(0, (days), by = 1/24)
time          <- time[-length(time)]

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


dose_28d_320     <- matrix(nrow = days, ncol = length(value_1))
dose_28d_1000    <- matrix(nrow = days, ncol = length(value_1))

for (i in 1:days){
  
  cat("day = ", i, '\n')
  
  dose_28d_320[i,]    <- value_1 * dose_320[i]  / MW_parent
  dose_28d_1000[i,]   <- value_1 * dose_1000[i] / MW_parent
}

df_dose_28d_320   <- as.vector(t(dose_28d_320))
df_dose_28d_1000  <- as.vector(t(dose_28d_1000))


# function 
Pred_320mgkg.adlib.2016     <- function(pars){
  
  names(pars) <- names(pars)
  
  ########################################################################
  ###                    define modeling variables                     ###  
  ########################################################################
  Dose_events_320    <- data.frame(var    = rep("Agutlumen_parent", days * 24),
                                   time   =  time,
                                   value  = df_dose_28d_320,
                                   method = "add")                  
  
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
  
  
  out_320mgkg.adlib.2016       <- ode(y       = initState, 
                                       times   = Times, 
                                       func    = pbtk5cpt_repeated, 
                                       parms   = pars,                        # make sure use pars here
                                       atol    = 1e-6, rtol = 1e-8,
                                       events  = list(data = Dose_events_320),
                                       method  = 'lsoda')
  
  out_320mgkg.adlib.2016                  <- as.data.frame(out_320mgkg.adlib.2016)
  colnames(out_320mgkg.adlib.2016)[1]     <- "Time"
  out_320mgkg.adlib.2016       <- cbind.data.frame (Time              = out_320mgkg.adlib.2016$Time,
                                                     C_blood_parent   = out_320mgkg.adlib.2016$C_blood_parent,
                                                     C_blood_daughter = out_320mgkg.adlib.2016$C_blood_daughter)    # umol/L 
  
  return(list("out_320mgkg.adlib.2016"     = out_320mgkg.adlib.2016))
}


Pred_1000mgkg.adlib.2016     <- function(pars){
  
  names(pars) <- names(pars)
  
  ########################################################################
  ###                    define modeling variables                     ###  
  ########################################################################
  Dose_events_1000   <- data.frame(var    = rep("Agutlumen_parent", days * 24),
                                   time   =  time,
                                   value  = df_dose_28d_1000,
                                   method = "add")                  
  
  initState <- c(Agutlumen_parent  = 0, 
                 Agut_parent       = 0,  
                 Aliver_parent     = 0,
                 Akidney_parent    = 0,
                 Aadipose_parent   = 0,
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
  
  
  out_1000mgkg.adlib.2016      <- ode(y       = initState, 
                                      times   = Times, 
                                      func    = pbtk5cpt_repeated, 
                                      parms   = pars,                        # make sure use pars here
                                      atol    = 1e-6, rtol = 1e-8,
                                      events  = list(data = Dose_events_1000),
                                      method  = 'lsoda')
  
  out_1000mgkg.adlib.2016                  <- as.data.frame(out_1000mgkg.adlib.2016)
  colnames(out_1000mgkg.adlib.2016)[1]     <- "Time"
  out_1000mgkg.adlib.2016       <- cbind.data.frame (Time              = out_1000mgkg.adlib.2016$Time,
                                                    C_blood_parent     = out_1000mgkg.adlib.2016$C_blood_parent,
                                                    C_blood_daughter   = out_1000mgkg.adlib.2016$C_blood_daughter)    # umol/L 
  
  return(list("out_1000mgkg.adlib.2016"     = out_1000mgkg.adlib.2016))
}


## Cost fuction (from FME pckage) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  
  out_320mgkg.adlib.2016                  <-  Pred_320mgkg.adlib.2016(pars)
  out_1000mgkg.adlib.2016                 <-  Pred_1000mgkg.adlib.2016(pars)
  #print(out_repeated_30mgkg.oral.Afifi.1997)
  
  cost  <- modCost(model = out_320mgkg.adlib.2016$out_320mgkg.adlib.2016,      obs= df_Pre_egg_female_320,   weight='std', x="Time")               # weight std can only be used if there is more than 1 data point
  cost  <- modCost(model = out_1000mgkg.adlib.2016$out_1000mgkg.adlib.2016,    obs= df_Pre_egg_female_1000,  weight='std', x="Time", cost = cost)
 
  return(cost)
}

########################################################################
#                      Senstivity function (FME)                       #
########################################################################
pars <-  c( ka                   = 17.48444,                
            Clint_parent         = 425.6144, 
            Clint_daughter       = 67.28163,   
            
            Rblood2plasma_parent   = 1.342,  
            Rblood2plasma_daughter = 0.932,   
            
            Qgfr_parent   = 5.181,  
            Qgfr_daughter = 21.8)


system.time(Fit<- modFit(f=MCcost, p=pars, method = "Marq" ,     #"Nelder-Mead" or Marq
                         lower = c(12,  350,   40,  0.9,  0.8,    4,  15.3),    
                         upper = c(24,  700,  100,  1.4,  1.4,  6.8,    24),
                         control = nls.lm.control(nprint=1)))


summary(Fit)
Fit$par

#==================================================================================================
#                                        check the fits                                           # 
#==================================================================================================
# Reference
# PKNCA: https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk: https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk
Sim.fit.out_320mgkg.adlib.2016                  = Pred_320mgkg.adlib.2016(Fit$par)$out_320mgkg.adlib.2016  
df.Sim.fit.out_320mgkg.adlib.2016.blood.TMX     = cbind.data.frame (Time=Sim.fit.out_320mgkg.adlib.2016$Time, C_blood_parent   = Sim.fit.out_320mgkg.adlib.2016$C_blood_parent)
df.Sim.fit.out_320mgkg.adlib.2016.blood.CLO     = cbind.data.frame (Time=Sim.fit.out_320mgkg.adlib.2016$Time, C_blood_daughter = Sim.fit.out_320mgkg.adlib.2016$C_blood_daughter)


Sim.fit.out_1000mgkg.adlib.2016                  = Pred_1000mgkg.adlib.2016(Fit$par)$out_1000mgkg.adlib.2016  
df.Sim.fit.out_1000mgkg.adlib.2016.blood.TMX     = cbind.data.frame (Time=Sim.fit.out_1000mgkg.adlib.2016$Time, C_blood_parent   = Sim.fit.out_1000mgkg.adlib.2016$C_blood_parent)
df.Sim.fit.out_1000mgkg.adlib.2016.blood.CLO     = cbind.data.frame (Time=Sim.fit.out_1000mgkg.adlib.2016$Time, C_blood_daughter = Sim.fit.out_1000mgkg.adlib.2016$C_blood_daughter)


##########################      pbpk output     #############################  
######           Plot           ###########
file_path                   <- file.path("HTTK", "Bird", "TMX", 'InvivoData.xlsx')
Residue_2016_quail          <- read_excel(file_path, sheet = "Residue_2016_quail")
# 320 ppm
# blood
Pre_egg_female                  <- Residue_2016_quail[Residue_2016_quail$Dose_ppm == 320 & Residue_2016_quail$group == 'pre-egg' & Residue_2016_quail$Gender == 'female',]
Pre_egg_female$TMXBloodConc_umol <- Pre_egg_female$Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female$CLOBloodConc_umol <- Pre_egg_female$Blood_conc_ng.ml_CGA322704 /MW_daughter
Pre_egg_female$TMXBloodConc_umol_SD <- Pre_egg_female$SD.Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female$CLOBloodConc_umol_SD <- Pre_egg_female$SD.Blood_conc_ng.ml_CGA322704 /MW_daughter
df_Pre_egg_female               <- subset(Pre_egg_female, select = c('Time_d', "TMXBloodConc_umol", 'TMXBloodConc_umol_SD',
                                                                           'CLOBloodConc_umol', "CLOBloodConc_umol_SD")) 
df_Pre_egg_female   <- na.omit(df_Pre_egg_female)
# assume sampling happens at 10 am
ggplot() +
  geom_line(data = df.Sim.fit.out_320mgkg.adlib.2016.blood.TMX, aes(Time, C_blood_parent), col="#00AFBB", lwd=2) +
  geom_point(data = df_Pre_egg_female , aes((Time_d - 14/24), TMXBloodConc_umol),col="red",  size=2.5) + 
  geom_errorbar(data = df_Pre_egg_female , aes(x=(Time_d - 14/24), 
                                                  ymin= TMXBloodConc_umol  - TMXBloodConc_umol_SD , ymax = TMXBloodConc_umol  + TMXBloodConc_umol_SD) , 
                size = 0.8,width=0.5, position=position_dodge(0.05),colour="black") + ylab("Concentration")+
  theme(text = element_text(size = 20)) + xlim(0,4)

ggplot() +
  geom_line(data = df.Sim.fit.out_320mgkg.adlib.2016.blood.CLO , aes(Time, C_blood_daughter), col="#00AFBB", lwd=2) +
  geom_point(data = df_Pre_egg_female , aes((Time_d - 14/24), CLOBloodConc_umol),col="red", size=2.5) +
  geom_errorbar(data = df_Pre_egg_female , aes(x=(Time_d - 14/24), 
                                                  ymin= CLOBloodConc_umol - CLOBloodConc_umol_SD, ymax = CLOBloodConc_umol + CLOBloodConc_umol_SD), 
                size = 0.8,width=0.5, position=position_dodge(0.05),colour="black") + ylab("Concentration") +
  ylab("Concentration")+
  theme(text = element_text(size = 20)) + xlim(0, 4) 


# 1000 ppm
# blood
Pre_egg_female                  <- Residue_2016_quail[Residue_2016_quail$Dose_ppm == 1000 & Residue_2016_quail$group == 'pre-egg' & Residue_2016_quail$Gender == 'female',]
Pre_egg_female$TMXBloodConc_umol <- Pre_egg_female$Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female$CLOBloodConc_umol <- Pre_egg_female$Blood_conc_ng.ml_CGA322704 /MW_daughter
Pre_egg_female$TMXBloodConc_umol_SD <- Pre_egg_female$SD.Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female$CLOBloodConc_umol_SD <- Pre_egg_female$SD.Blood_conc_ng.ml_CGA322704 /MW_daughter
df_Pre_egg_female               <- subset(Pre_egg_female, select = c('Time_d', "TMXBloodConc_umol", 'TMXBloodConc_umol_SD',
                                                                     'CLOBloodConc_umol', "CLOBloodConc_umol_SD")) 
df_Pre_egg_female   <- na.omit(df_Pre_egg_female)
# assume sampling happens at 10 am
ggplot() +
  geom_line(data = df.Sim.fit.out_1000mgkg.adlib.2016.blood.TMX, aes(Time, C_blood_parent), col="#00AFBB", lwd=2) +
  geom_point(data = df_Pre_egg_female , aes((Time_d - 14/24), TMXBloodConc_umol),col="red",  size=2.5) + 
  geom_errorbar(data = df_Pre_egg_female , aes(x=(Time_d - 14/24), 
                                               ymin= TMXBloodConc_umol  - TMXBloodConc_umol_SD , ymax = TMXBloodConc_umol  + TMXBloodConc_umol_SD) , 
                size = 0.8,width=0.5, position=position_dodge(0.05),colour="black") + ylab("Concentration")+
  theme(text = element_text(size = 20)) + xlim(0,4)

ggplot() +
  geom_line(data = df.Sim.fit.out_1000mgkg.adlib.2016.blood.CLO , aes(Time, C_blood_daughter), col="#00AFBB", lwd=2) +
  geom_point(data = df_Pre_egg_female , aes((Time_d - 14/24), CLOBloodConc_umol),col="red", size=2.5) +
  geom_errorbar(data = df_Pre_egg_female , aes(x=(Time_d - 14/24), 
                                               ymin= CLOBloodConc_umol - CLOBloodConc_umol_SD, ymax = CLOBloodConc_umol + CLOBloodConc_umol_SD), 
                size = 0.8,width=0.5, position=position_dodge(0.05),colour="black") + ylab("Concentration") +
  ylab("Concentration")+
  theme(text = element_text(size = 20)) + xlim(0, 4) 

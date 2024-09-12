library(FME)
library(invgamma) 
library(ggplot2)
library(httk)
library(deSolve)
library(pksensi)
library(ggplot2)
library(PKNCA)
library(sensitivity)
library(ODEsensitivity)
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(tidyr)
library(tidyverse)
library(dplyr)
library(minpack.lm)  # Package for model fitting

current_dir <- getwd()
print(current_dir)
source(file.path("HTTK", "Bird", "TMX", 'Quail', "Quail_Single_oral_NoEgg.R"))

#################################################################################
##########            calibration dataset (Japanese quail)            ###########
####  This R script only served as a reference for coding for calibration    ####
#################################################################################
file_path               <- file.path("HTTK", "Bird", "TMX", 'InvivoData.xlsx')
Residue_2022_quail      <- read_excel(file_path , sheet = "Pan_2022_quail")
# TMX Plasma
Plasma                  <- Residue_2022_quail[Residue_2022_quail$Matrix ==  'plasma',]
Plasma_TMX              <- Plasma[, c('Time_h', 'TMX')]
Plasma_CLO              <- Plasma[, c('Time_h', 'CLO')]
Plasma_TMX$Time               <- round(Plasma_TMX$Time_h/24, 2)
Plasma_CLO$Time               <- round(Plasma_CLO$Time_h/24, 2)
Plasma_TMX$C_plasma_parent    <- Plasma_TMX$TMX / MW_parent       # umol/L
Plasma_CLO$C_plasma_daughter  <- Plasma_CLO$CLO / MW_daughter     # umol/L

Plasma_TMX_5mgkg.oral         <- as.data.frame(Plasma_TMX %>% dplyr::select(Time, C_plasma_parent))
Plasma_CLO_5mgkg.oral         <- as.data.frame(Plasma_CLO %>% dplyr::select(Time, C_plasma_daughter))
Plasma_TMX_5mgkg.oral         <- na.omit(Plasma_TMX_5mgkg.oral )   
Plasma_CLO_5mgkg.oral         <- na.omit(Plasma_CLO_5mgkg.oral ) 

# muscle
muscle                  <- Residue_2022_quail[Residue_2022_quail$Matrix ==  'muscle',]
muscle_TMX              <- muscle[, c('Time_h', 'TMX')]
muscle_CLO              <- muscle[, c('Time_h', 'CLO')]
muscle_TMX$Time               <- round(muscle_TMX$Time_h/24, 2)
muscle_CLO$Time               <- round(muscle_CLO$Time_h/24, 2)
muscle_TMX$C_muscle_parent    <- muscle_TMX$TMX / MW_parent       # umol/L
muscle_CLO$C_muscle_daughter  <- muscle_CLO$CLO / MW_daughter     # umol/L

muscle_TMX_5mgkg.oral         <- as.data.frame(muscle_TMX %>% dplyr::select(Time, C_muscle_parent))
muscle_CLO_5mgkg.oral         <- as.data.frame(muscle_CLO %>% dplyr::select(Time, C_muscle_daughter))
muscle_TMX_5mgkg.oral         <- na.omit(muscle_TMX_5mgkg.oral )   
muscle_CLO_5mgkg.oral         <- na.omit(muscle_CLO_5mgkg.oral ) 

# kidney
kidney                  <- Residue_2022_quail[Residue_2022_quail$Matrix ==  'kidney',]
kidney_TMX              <- kidney[, c('Time_h', 'TMX')]
kidney_CLO              <- kidney[, c('Time_h', 'CLO')]
kidney_TMX$Time               <- round(kidney_TMX$Time_h/24, 2)
kidney_CLO$Time               <- round(kidney_CLO$Time_h/24, 2)
kidney_TMX$C_kidney_parent    <- kidney_TMX$TMX / MW_parent       # umol/L
kidney_CLO$C_kidney_daughter  <- kidney_CLO$CLO / MW_daughter     # umol/L

kidney_TMX_5mgkg.oral         <- as.data.frame(kidney_TMX %>% dplyr::select(Time, C_kidney_parent))
kidney_CLO_5mgkg.oral         <- as.data.frame(kidney_CLO %>% dplyr::select(Time, C_kidney_daughter))
kidney_TMX_5mgkg.oral         <- na.omit(kidney_TMX_5mgkg.oral )   
kidney_CLO_5mgkg.oral         <- na.omit(kidney_CLO_5mgkg.oral ) 

# liver
liver                  <- Residue_2022_quail[Residue_2022_quail$Matrix ==  'liver',]
liver_TMX              <- liver[, c('Time_h', 'TMX')]
liver_CLO              <- liver[, c('Time_h', 'CLO')]
liver_TMX$Time               <- round(liver_TMX$Time_h/24, 2)
liver_CLO$Time               <- round(liver_CLO$Time_h/24, 2)
liver_TMX$C_liver_parent    <- liver_TMX$TMX / MW_parent       # umol/L
liver_CLO$C_liver_daughter  <- liver_CLO$CLO / MW_daughter     # umol/L

liver_TMX_5mgkg.oral         <- as.data.frame(liver_TMX %>% dplyr::select(Time, C_liver_parent))
liver_CLO_5mgkg.oral         <- as.data.frame(liver_CLO %>% dplyr::select(Time, C_liver_daughter))
liver_TMX_5mgkg.oral         <- na.omit(liver_TMX_5mgkg.oral )   
liver_CLO_5mgkg.oral         <- na.omit(liver_CLO_5mgkg.oral ) 


########################################################################
#                            Calibration                               #
########################################################################
Pred_5mgkg.oral.Pan.2022     <- function(pars){
  
  names(pars) <- names(pars)
  
  ########################################################################
  ###                    define modeling variables                     ###  
  ########################################################################
  days          <- 5
  Time_max      <- days                                     # 5 d  
  
  StartTime     <- 0                                        # Time_d
  StopTime      <- Time_max                                 
  dt            <- 0.01
  Times         <- seq(StartTime, StopTime, dt)
  
  Oral_input    <- 5 * 1000 /MW_parent                              # ug/kg/d
  
  Dose_events    <- data.frame(var    = "Agutlumen_parent",
                               time   =  0,
                               value  = Oral_input,
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
  
  
  out_5mgkg.oral.Pan.2022       <- ode(y       = initState, 
                                       times   = Times, 
                                       func    = pbtk5cpt_repeated, 
                                       parms   = pars,                        # make sure use pars here
                                       atol    = 1e-6, rtol = 1e-8,
                                       events  = list(data = Dose_events),
                                       method  = 'lsoda')
  
  out_5mgkg.oral.Pan.2022                  <- as.data.frame(out_5mgkg.oral.Pan.2022)
  colnames(out_5mgkg.oral.Pan.2022)[1]     <- "Time"
  out_5mgkg.oral.Pan.2022       <- cbind.data.frame (Time              = out_5mgkg.oral.Pan.2022$Time,
                                                     C_plasma_parent   = out_5mgkg.oral.Pan.2022$C_plasma_parent,
                                                     C_plasma_daughter = out_5mgkg.oral.Pan.2022$C_plasma_daughter,
                                                     C_liver_parent    = out_5mgkg.oral.Pan.2022$C_liver_parent,
                                                     C_liver_daughter  = out_5mgkg.oral.Pan.2022$C_liver_daughter,
                                                     C_muscle_parent    = out_5mgkg.oral.Pan.2022$C_muscle_parent,
                                                     C_muscle_daughter  = out_5mgkg.oral.Pan.2022$C_muscle_daughter,
                                                     C_kidney_parent    = out_5mgkg.oral.Pan.2022$C_kidney_parent,
                                                     C_kidney_daughter  = out_5mgkg.oral.Pan.2022$C_kidney_daughter)    # ug/L 
  
  return(list("out_5mgkg.oral.Pan.2022"     = out_5mgkg.oral.Pan.2022))
}



## Cost fuction (from FME pckage) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  
  out_5mgkg.oral.Pan.2022                  <-  Pred_5mgkg.oral.Pan.2022(pars)
  #print(out_repeated_30mgkg.oral.Afifi.1997)
  
  cost  <- modCost(model = out_5mgkg.oral.Pan.2022$out_5mgkg.oral.Pan.2022,    obs= Plasma_TMX_5mgkg.oral,  weight='std', x="Time")               # weight std can only be used if there is more than 1 data point
  cost  <- modCost(model = out_5mgkg.oral.Pan.2022$out_5mgkg.oral.Pan.2022,    obs= Plasma_CLO_5mgkg.oral,  weight='std', x="Time", cost = cost)
  cost  <- modCost(model = out_5mgkg.oral.Pan.2022$out_5mgkg.oral.Pan.2022,    obs= muscle_TMX_5mgkg.oral,  weight='std', x="Time", cost = cost)
  cost  <- modCost(model = out_5mgkg.oral.Pan.2022$out_5mgkg.oral.Pan.2022,    obs= muscle_CLO_5mgkg.oral,  weight='std', x="Time", cost = cost)
  cost  <- modCost(model = out_5mgkg.oral.Pan.2022$out_5mgkg.oral.Pan.2022,    obs= kidney_TMX_5mgkg.oral,  weight='std', x="Time", cost = cost)
  cost  <- modCost(model = out_5mgkg.oral.Pan.2022$out_5mgkg.oral.Pan.2022,    obs= kidney_CLO_5mgkg.oral,  weight='std', x="Time", cost = cost)
  cost  <- modCost(model = out_5mgkg.oral.Pan.2022$out_5mgkg.oral.Pan.2022,    obs= liver_TMX_5mgkg.oral,   weight='std', x="Time", cost = cost)
  cost  <- modCost(model = out_5mgkg.oral.Pan.2022$out_5mgkg.oral.Pan.2022,    obs= liver_CLO_5mgkg.oral,   weight='std', x="Time", cost = cost)
  
  
 
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
            Qgfr_daughter = 21.8,
            
            Kmuscle2pu_parent      = 1.114,
            Kmuscle2pu_daughter    = 1.331,
            Kkidney2pu_parent      = 1.503,
            Kkidney2pu_daughter    = 2.041,
            Kliver2pu_parent       = 1.548,
            Kliver2pu_daughter     = 2.192)



system.time(Fit<- modFit(f=MCcost, p=pars, method = "Nelder-Mead" ,     #"Nelder-Mead" or Marq
                         lower = c( 6,  300,  40,  0.8,  0.8,   4,  15.3,  1,  1,  1,  1,  1,  1),    
                         upper = c(24,  800,  120,  1.4,  1.4, 6.8,    25,  5,  5,  5,  5,  5,  5),
                         control = nls.lm.control(nprint=1)))


summary(Fit)
Fit$par

#======================================================================================================
############                Prediction using the calibrated PBPK model                    #############
#======================================================================================================
## Model calibration plot using ggplot2 
Sim.fit.out_5mgkg.oral.Pan.2022               = Pred_5mgkg.oral.Pan.2022(Fit$par)$out_5mgkg.oral.Pan.2022  
df.sim.out_5mgkg.oral.Pan.2022.plasma.TMX     = cbind.data.frame (Time=Sim.fit.out_5mgkg.oral.Pan.2022$Time, C_plasma_parent   = Sim.fit.out_5mgkg.oral.Pan.2022$C_plasma_parent)
df.sim.out_5mgkg.oral.Pan.2022.plasma.CLO     = cbind.data.frame (Time=Sim.fit.out_5mgkg.oral.Pan.2022$Time, C_plasma_daughter = Sim.fit.out_5mgkg.oral.Pan.2022$C_plasma_daughter)

# ##########################      pbpk output     #############################  
# ######           Plot           ###########
# plasma
ggplot() +
   geom_line(data = df.sim.out_5mgkg.oral.Pan.2022.plasma.TMX, aes(Time, C_plasma_parent), col="#00AFBB", lwd=2) +
   geom_point(data = Plasma_TMX_5mgkg.oral  , aes((Time), C_plasma_parent),col="red",  size=2.5) + xlim(0,1)#

ggplot() +
   geom_line(data = df.sim.out_5mgkg.oral.Pan.2022.plasma.CLO, aes(Time, C_plasma_daughter), col="#00AFBB", lwd=2) +
   geom_point(data = Plasma_CLO_5mgkg.oral  , aes((Time), C_plasma_daughter),col="red",  size=2.5) + xlim(0,1)

 
# liver
df.sim.out_5mgkg.oral.Pan.2022.liver.TMX     = cbind.data.frame (Time=Sim.fit.out_5mgkg.oral.Pan.2022$Time, C_liver_parent   = Sim.fit.out_5mgkg.oral.Pan.2022$C_liver_parent)
df.sim.out_5mgkg.oral.Pan.2022.liver.CLO     = cbind.data.frame (Time=Sim.fit.out_5mgkg.oral.Pan.2022$Time, C_liver_daughter = Sim.fit.out_5mgkg.oral.Pan.2022$C_liver_daughter)

liver                  <- Residue_2022_quail[Residue_2022_quail$Matrix ==  'liver',]
liver_TMX              <- liver[, c('Time_h', 'TMX')]
liver_CLO              <- liver[, c('Time_h', 'CLO')]
liver_TMX$Time               <- round(liver_TMX$Time_h/24, 2)
liver_CLO$Time               <- round(liver_CLO$Time_h/24, 2)
liver_TMX$C_liver_parent    <- liver_TMX$TMX / MW_parent       # umol/L
liver_CLO$C_liver_daughter  <- liver_CLO$CLO / MW_daughter     # umol/L

liver_TMX_5mgkg.oral         <- as.data.frame(liver_TMX %>% dplyr::select(Time, C_liver_parent))
liver_CLO_5mgkg.oral         <- as.data.frame(liver_CLO %>% dplyr::select(Time, C_liver_daughter))
liver_TMX_5mgkg.oral         <- na.omit(liver_TMX_5mgkg.oral )   
liver_CLO_5mgkg.oral         <- na.omit(liver_CLO_5mgkg.oral ) 

ggplot() +
  geom_line(data = df.sim.out_5mgkg.oral.Pan.2022.liver.TMX, aes(Time, C_liver_parent), col="#00AFBB", lwd=2) +
  geom_point(data = liver_TMX_5mgkg.oral  , aes((Time), C_liver_parent),col="red",  size=2.5) + xlim(0,1)#

ggplot() +
  geom_line(data = df.sim.out_5mgkg.oral.Pan.2022.liver.CLO, aes(Time, C_liver_daughter), col="#00AFBB", lwd=2) +
  geom_point(data = liver_CLO_5mgkg.oral  , aes((Time), C_liver_daughter),col="red",  size=2.5) + xlim(0,1)


# muscle
df.sim.out_5mgkg.oral.Pan.2022.muscle.TMX     = cbind.data.frame (Time=Sim.fit.out_5mgkg.oral.Pan.2022$Time, C_muscle_parent   = Sim.fit.out_5mgkg.oral.Pan.2022$C_muscle_parent)
df.sim.out_5mgkg.oral.Pan.2022.muscle.CLO     = cbind.data.frame (Time=Sim.fit.out_5mgkg.oral.Pan.2022$Time, C_muscle_daughter = Sim.fit.out_5mgkg.oral.Pan.2022$C_muscle_daughter)

muscle                  <- Residue_2022_quail[Residue_2022_quail$Matrix ==  'muscle',]
muscle_TMX              <- muscle[, c('Time_h', 'TMX')]
muscle_CLO              <- muscle[, c('Time_h', 'CLO')]
muscle_TMX$Time               <- round(muscle_TMX$Time_h/24, 2)
muscle_CLO$Time               <- round(muscle_CLO$Time_h/24, 2)
muscle_TMX$C_muscle_parent    <- muscle_TMX$TMX / MW_parent       # umol/L
muscle_CLO$C_muscle_daughter  <- muscle_CLO$CLO / MW_daughter     # umol/L

muscle_TMX_5mgkg.oral         <- as.data.frame(muscle_TMX %>% dplyr::select(Time, C_muscle_parent))
muscle_CLO_5mgkg.oral         <- as.data.frame(muscle_CLO %>% dplyr::select(Time, C_muscle_daughter))
muscle_TMX_5mgkg.oral         <- na.omit(muscle_TMX_5mgkg.oral )   
muscle_CLO_5mgkg.oral         <- na.omit(muscle_CLO_5mgkg.oral ) 

ggplot() +
  geom_line(data = df.sim.out_5mgkg.oral.Pan.2022.muscle.TMX, aes(Time, C_muscle_parent), col="#00AFBB", lwd=2) +
  geom_point(data = muscle_TMX_5mgkg.oral  , aes((Time), C_muscle_parent),col="red",  size=2.5) + xlim(0,1)#

ggplot() +
  geom_line(data = df.sim.out_5mgkg.oral.Pan.2022.muscle.CLO, aes(Time, C_muscle_daughter), col="#00AFBB", lwd=2) +
  geom_point(data = muscle_CLO_5mgkg.oral  , aes((Time), C_muscle_daughter),col="red",  size=2.5) + xlim(0,1)


# kidney
df.sim.out_5mgkg.oral.Pan.2022.kidney.TMX     = cbind.data.frame (Time=Sim.fit.out_5mgkg.oral.Pan.2022$Time, C_kidney_parent   = Sim.fit.out_5mgkg.oral.Pan.2022$C_kidney_parent)
df.sim.out_5mgkg.oral.Pan.2022.kidney.CLO     = cbind.data.frame (Time=Sim.fit.out_5mgkg.oral.Pan.2022$Time, C_kidney_daughter = Sim.fit.out_5mgkg.oral.Pan.2022$C_kidney_daughter)

kidney                  <- Residue_2022_quail[Residue_2022_quail$Matrix ==  'kidney',]
kidney_TMX              <- kidney[, c('Time_h', 'TMX')]
kidney_CLO              <- kidney[, c('Time_h', 'CLO')]
kidney_TMX$Time               <- round(kidney_TMX$Time_h/24, 2)
kidney_CLO$Time               <- round(kidney_CLO$Time_h/24, 2)
kidney_TMX$C_kidney_parent    <- kidney_TMX$TMX / MW_parent       # umol/L
kidney_CLO$C_kidney_daughter  <- kidney_CLO$CLO / MW_daughter     # umol/L

kidney_TMX_5mgkg.oral         <- as.data.frame(kidney_TMX %>% dplyr::select(Time, C_kidney_parent))
kidney_CLO_5mgkg.oral         <- as.data.frame(kidney_CLO %>% dplyr::select(Time, C_kidney_daughter))
kidney_TMX_5mgkg.oral         <- na.omit(kidney_TMX_5mgkg.oral )   
kidney_CLO_5mgkg.oral         <- na.omit(kidney_CLO_5mgkg.oral ) 

ggplot() +
  geom_line(data = df.sim.out_5mgkg.oral.Pan.2022.kidney.TMX, aes(Time, C_kidney_parent), col="#00AFBB", lwd=2) +
  geom_point(data = kidney_TMX_5mgkg.oral  , aes((Time), C_kidney_parent),col="red",  size=2.5) + xlim(0,1)#

ggplot() +
  geom_line(data = df.sim.out_5mgkg.oral.Pan.2022.kidney.CLO, aes(Time, C_kidney_daughter), col="#00AFBB", lwd=2) +
  geom_point(data = kidney_CLO_5mgkg.oral  , aes((Time), C_kidney_daughter),col="red",  size=2.5) + xlim(0,1)


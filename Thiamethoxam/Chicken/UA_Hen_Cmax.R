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


#======================================================================================================
##########                            Uncertainty analysis (Cmax)                             ##########
#======================================================================================================

SA_hen       <- read.csv('C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Laying Hen/Sensitivity.csv')
SA_duck      <- read.csv('C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Duck/Sensitivity.csv')
SA_quail     <- read.csv('C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Quail/Sensitivity.csv')

SA_hen$Species        <- 'Chicken'
SA_duck$Species       <- 'Duck'
SA_quail$Species      <- 'Quail'

SA                    <- as.data.frame(rbind(SA_hen, SA_duck, SA_quail))
SA                    <- unique(SA)
SA$Species            <- factor(SA$Species, levels=c('Chicken', 'Duck', 'Quail'))
colnames(SA)[1]       <- 'parameter'
SA$parameter          <- gsub("fub", "fup", SA$parameter)

SA$AUC                <- (SA$AUC)/100
SA$Cmax               <- (SA$Cmax)/100

SA_AUC                <- subset(SA, select = - abs(Cmax))
SA_Cmax               <- subset(SA, select = - abs(AUC))
SA_AUC                <- SA_AUC[abs(SA_AUC$AUC) >= 0.1, ]
SA_Cmax               <- SA_Cmax[abs(SA_Cmax$Cmax) >= 0.1, ]

SA_Cmax_chicken          <- SA_Cmax[SA_Cmax$Species == 'Chicken',]
SA_Cmax_chicken   <- SA_Cmax_chicken   %>%
  mutate(group = case_when(
    abs(Cmax) >=0.5 ~ 'high',
    abs(Cmax) < 0.5 & abs(Cmax) >= 0.2 ~ 'medium',
    abs(Cmax) < 0.2 ~ 'low'
  ))

print(SA_Cmax_chicken$parameter)

library(truncnorm) 
library(ggridges)     # Used to create Figure
library(ggplot2)   # ggplot is the basic package for creating plots.
#install.packages("tkrplot")
library(tkrplot)
#install.packages("rriskDistributions")
library(rriskDistributions)
library(EnvStats)
library(httk)
library(mc2d)
library(fitdistrplus)

#rm(list = ls())
set.seed(10) 

# Normal distributions were implemented for physiological parameters; 
# chemical-specific parameters were assumed to be log-normally distributed (Henri et al., 2017; Li et al., 2017; Yang et al., 2015).

######################################################################
#####################  NOW run local sensitivity file (chicken) until line 330
#dev.off()

SPECIES  <- 'Laying_hen'


#### mean value of parameters
Rblood2plasma_parent.mean          <- 1.342
Rblood2plasma_daughter.mean        <- 0.9450215
fub_parent.mean                    <- 0.55
fub_daughter.mean                  <- 0.71
Clint_parent.mean                  <- 259.3772
Clint_daughter.mean                <- 41.00266
ka.mean                            <- 12.1847922
Kmuscle2pu_daughter.mean           <- 1.33
Vmuscle.mean                       <- 0.398
                 

#### std value of parameters

Rblood2plasma_parent.sd       <- 0.3 * Rblood2plasma_parent.mean
Rblood2plasma_daughter.sd     <- 0.3 * Rblood2plasma_daughter.mean
fub_parent.sd                 <- 0.3 * fub_parent.mean
fub_daughter.sd               <- 0.3 * fub_daughter.mean
Clint_parent.sd               <- 0.5 * Clint_parent.mean
Clint_daughter.sd             <- 0.5 * Clint_daughter.mean 
ka.sd                         <- 0.3 * ka.mean
Kmuscle2pu_daughter.sd        <- 0.3 * Kmuscle2pu_daughter.mean
Vmuscle.sd                    <- 0.12 * Vmuscle.mean                      



#############     define distribution
# Normal distributions were implemented for physiological parameters; 
m.log.Rblood2plasma_parent             <- log(Rblood2plasma_parent.mean^2/(Rblood2plasma_parent.sd^2+Rblood2plasma_parent.mean^2)^0.5) 
sd.log.Rblood2plasma_parent            <- (log(1+Rblood2plasma_parent.sd^2/Rblood2plasma_parent.mean^2))^0.5 

m.log.Rblood2plasma_daughter           <- log(Rblood2plasma_daughter.mean^2/(Rblood2plasma_daughter.sd^2+Rblood2plasma_daughter.mean^2)^0.5) 
sd.log.Rblood2plasma_daughter          <- (log(1+Rblood2plasma_daughter.sd^2/Rblood2plasma_daughter.mean^2))^0.5 

m.log.fub_parent             <- log(fub_parent.mean^2/(fub_parent.sd^2+fub_parent.mean^2)^0.5)
sd.log.fub_parent            <- (log(1+fub_parent.sd^2/fub_parent.mean^2))^0.5

m.log.fub_daughter           <- log(fub_daughter.mean^2/(fub_daughter.sd^2+fub_daughter.mean^2)^0.5) 
sd.log.fub_daughter          <- (log(1+fub_daughter.sd^2/fub_daughter.mean^2))^0.5 

m.log.Clint_parent         <- log(Clint_parent.mean^2/(Clint_parent.sd^2+Clint_parent.mean^2)^0.5)
sd.log.Clint_parent        <- (log(1+Clint_parent.sd^2/Clint_parent.mean^2))^0.5

m.log.Clint_daughter         <- log(Clint_daughter.mean^2/(Clint_daughter.sd^2+Clint_daughter.mean^2)^0.5)
sd.log.Clint_daughter        <- (log(1+Clint_daughter.sd^2/Clint_daughter.mean^2))^0.5

m.log.Kmuscle2pu_daughter    <- log(Kmuscle2pu_daughter.mean^2/(Kmuscle2pu_daughter.sd^2+Kmuscle2pu_daughter.mean^2)^0.5)
sd.log.Kmuscle2pu_daughter   <- (log(1+Kmuscle2pu_daughter.sd^2/Kmuscle2pu_daughter.mean^2))^0.5

m.log.ka                     <- log(ka.mean^2/(ka.sd^2+ka.mean^2)^0.5)
sd.log.ka                    <- (log(1+ka.sd^2/ka.mean^2))^0.5



####################
N=1000
set.seed(10) 

idata <- 
  tibble(ID=1:N) %>% 
  mutate(
    
    Vmuscle = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = Vmuscle.mean, sd = Vmuscle.sd),
      b = qnorm(0.975, mean = Vmuscle.mean, sd = Vmuscle.sd),
      mean = Vmuscle.mean,
      sd = Vmuscle.sd
    ),
    
    Clint_parent = rlnormTrunc(
      N,
      meanlog = m.log.Clint_parent,
      sdlog = sd.log.Clint_parent,
      min = qlnorm(0.025, meanlog = m.log.Clint_parent, sdlog = sd.log.Clint_parent),
      max = qlnorm(0.975, meanlog = m.log.Clint_parent, sdlog = sd.log.Clint_parent)
    ),
    
    ## Chemical properties
    Clint_daughter = rlnormTrunc(
      N,
      meanlog = m.log.Clint_daughter,
      sdlog = sd.log.Clint_daughter,
      min = qlnorm(0.025, meanlog = m.log.Clint_daughter, sdlog = sd.log.Clint_daughter),
      max = qlnorm(0.975, meanlog = m.log.Clint_daughter, sdlog = sd.log.Clint_daughter)
    ),
    
    fub_parent = rlnormTrunc(
      N,
      meanlog = m.log.fub_parent,
      sdlog = sd.log.fub_parent,
      min = qlnorm(0.025, meanlog = m.log.fub_parent, sdlog = sd.log.fub_parent),
      max = qlnorm(0.975, meanlog = m.log.fub_parent, sdlog = sd.log.fub_parent)
    ),
    
    Kmuscle2pu_daughter = rlnormTrunc(
      N,
      meanlog = m.log.Kmuscle2pu_daughter,
      sdlog = sd.log.Kmuscle2pu_daughter,
      min = qlnorm(0.025, meanlog = m.log.Kmuscle2pu_daughter, sdlog = sd.log.Kmuscle2pu_daughter),
      max = qlnorm(0.975, meanlog = m.log.Kmuscle2pu_daughter, sdlog = sd.log.Kmuscle2pu_daughter)
    ),
    
    ka = rlnormTrunc(
      N,
      meanlog = m.log.ka,
      sdlog = sd.log.ka,
      min = qlnorm(0.025, meanlog = m.log.ka, sdlog = sd.log.ka),
      max = qlnorm(0.975, meanlog = m.log.ka, sdlog = sd.log.ka)
    ),

    Rblood2plasma_parent = rlnormTrunc(
      N,
      meanlog = m.log.Rblood2plasma_parent,
      sdlog = sd.log.Rblood2plasma_parent,
      min = qlnorm(0.025, meanlog = m.log.Rblood2plasma_parent, sdlog = sd.log.Rblood2plasma_parent),
      max = qlnorm(0.975, meanlog = m.log.Rblood2plasma_parent, sdlog = sd.log.Rblood2plasma_parent)
    ),
    
    Rblood2plasma_daughter = rlnormTrunc(
      N,
      meanlog = m.log.Rblood2plasma_daughter,
      sdlog = sd.log.Rblood2plasma_daughter,
      min = qlnorm(0.025, meanlog = m.log.Rblood2plasma_daughter, sdlog = sd.log.Rblood2plasma_daughter),
      max = qlnorm(0.975, meanlog = m.log.Rblood2plasma_daughter, sdlog = sd.log.Rblood2plasma_daughter)
    ),
    
    fub_daughter = rlnormTrunc(
      N,
      meanlog = m.log.fub_daughter,
      sdlog = sd.log.fub_daughter,
      min = qlnorm(0.025, meanlog = m.log.fub_daughter, sdlog = sd.log.fub_daughter),
      max = qlnorm(0.975, meanlog = m.log.fub_daughter, sdlog = sd.log.fub_daughter)
    )

  )


N_refined  <- nrow(idata)
N_refined 



########################################################################
###                    define modeling variables                     ###  
########################################################################
tlag          <- 1                                        # time interval between ovulation and egg lay (approximately one day or 24 h)
tsig          <- 2                                        # pre-ovulation time of maximum follicle growth rate (48h, time for maximum growth achieved)
talbumen      <- 10 /24                                   # 10 h


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
               Awhite_parent     = 0,  #19
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
               AUC_Cblood_parent  = 0, # 34
               
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



#============================================================
######             Single dose population               #####
#============================================================
fraction_daughter            <-  1
pars  <-  c( ka                   = ka,                 # h-1
             Clint_parent         = 259.3772, # 50
             Clint_daughter       = 41.00266,   # /4
             
             Rblood2plasma_parent   = Rblood2plasma_parent,  
             Rblood2plasma_daughter = Rblood2plasma_daughter,   
             
             fub_parent           = fub_parent,   
             fub_daughter         = fub_daughter, 
             
             # ratio
             Kgut2pu_parent       = Kgut2pu_parent    ,      
             Kkidney2pu_parent    = Kkidney2pu_parent ,   
             Kliver2pu_parent     = Kliver2pu_parent  ,     
             Klung2pu_parent      = Klung2pu_parent   ,      
             Krest2pu_parent      = Krest2pu_parent   , 
             Kadipose2pu_parent   = Kadipose2pu_parent ,
             Kmuscle2pu_parent    = Kmuscle2pu_parent ,
             Krep2pu_parent       = Krep2pu_parent    ,
             
             Kgut2pu_daughter       = Kgut2pu_daughter  ,      
             Kkidney2pu_daughter    = Kkidney2pu_daughter ,   
             Kliver2pu_daughter     = Kliver2pu_daughter  ,     
             Klung2pu_daughter      = Klung2pu_daughter   ,      
             Krest2pu_daughter      = Krest2pu_daughter   , 
             Kadipose2pu_daughter   = Kadipose2pu_daughter ,
             Kmuscle2pu_daughter    = Kmuscle2pu_daughter ,
             Krep2pu_daughter       = Krep2pu_daughter    ,
             
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
             Vven          = Vven )



Pred     <- function(pars){
  
  names(pars) <- names(pars)
  
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
  
  # event
  event_dose <- function(t, y, parms) {
    y[1]       <- y[1] + Oral_input                                  
    return(y)
  }
  
  event_repeated <- function(t, y, parms) {
    ret     <- y
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
    if(abs(t-125.38) < 1e-6)         ret <- event9(t, y, parms)
    
    if(abs(t-0.3) < 1e-6)              ret <- event_dose(t, y, parms)
    #if(abs(t-1.3) < 1e-6)              ret <- event_dose(t, y, parms)
    
    return((ret))
  }
  
  
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
  out              <-        cbind.data.frame(Time = df$Time,
                                              AUC  = df$AUC_Cblood_parent + df$AUC_Cblood_daughter,
                                              Cmax = max(df$C_blood_parent + df$C_blood_daughter))
  
  return(list("out"  = out))
}



Newtime.r         = Pred(pars)$out$Time
nrwo.r            = length (Newtime.r)

# Create the matrix 


print(SA_Cmax_chicken$parameter)

UA_list  <- c('Clint_parent',
              'Clint_daughter',
              'Rblood2plasma_parent', 
              'Rblood2plasma_daughter',
              'fub_parent',
              'fub_daughter',
              'ka',
              'Kmuscle2pu_daughter',
              'Vmuscle')



MC.Cmax     = matrix(nrow = N_refined , ncol = length(UA_list))
colnames(MC.Cmax)   <- UA_list

for (j in seq(UA_list)){
  
  p = UA_list[j]
  
  for (i in 1:N_refined){
    
    cat("Parameter = ", p, "; iteration = ", i , "\n")
    
    new_pars       <- pars
    pars_UA        <- idata[i,]%>% dplyr :: select(-ID)
    
    if( p == 'Qadipose'){
      Qrest                 <-  Qcardiac - Qliver - Qgut - pars_UA[[p]] - Qbrain
      new_pars[['Qrest']]   <- Qrest
    }
    if( p == 'Qgut'){
      Qrest                 <-  Qcardiac - Qliver - Qadipose - pars_UA[[p]] - Qbrain
      new_pars[['Qrest']]   <- Qrest
    }
    
    new_pars[[p]]  <- pars_UA[[p]]
    
    
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
                   Awhite_parent     = 0,  #19
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
                   AUC_Cblood_parent  = 0, # 34
                   
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
    
    # event
    event_dose <- function(t, y, parms) {
      y[1]       <- y[1] + Oral_input                                  
      return(y)
    }
    
    event_repeated <- function(t, y, parms) {
      ret     <- y
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
      if(abs(t-125.38) < 1e-6)         ret <- event9(t, y, parms)
      
      if(abs(t-0.3) < 1e-6)              ret <- event_dose(t, y, parms)
      #if(abs(t-1.3) < 1e-6)              ret <- event_dose(t, y, parms)
      
      return((ret))
    }
    
    
    df  <- ode(y       = initState, 
               times   = Times, 
               func    = pbtk5cpt_repeated, 
               parms   = parms,
               atol    = 1e-6, rtol = 1e-8,
               events  = list(func = event_repeated, root = TRUE), 
               rootfun = root_repeated ,
               method  = 'lsoda')
    
    
    df                   <- as.data.frame(df)
    colnames(df)[1]      <- "Time"
    
    MC.Cmax[i,j]     <- max(df$C_blood_parent + df$C_blood_daughter)

    
  }
}

# Calculate the 95th percentile for each column
percentile_95 <- apply(MC.Cmax, 2, function(x) quantile(x, probs = 0.95))

# Calculate the median for each column
medians <- apply(MC.Cmax, 2, median)

# Create a new data frame with the 95th percentile and median
new_df <- rbind(percentile_95, medians)

stat  <- t(new_df)
stat   <- as.data.frame(stat)
stat$ratio  <- round((stat$percentile_95 /stat$medians), 2)

stat <- stat %>%
  mutate(group = case_when(
    ratio >=2 ~ 'high',
    ratio < 2 & ratio >= 1.3 ~ 'medium',
    ratio < 1.3 ~ 'low'
  ))

# end



write.csv(SA_Cmax_chicken, "C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Laying hen/SA_Cmax_chicken.csv")
write.csv(stat, "C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Laying hen/UA_Cmax_chicken.csv")









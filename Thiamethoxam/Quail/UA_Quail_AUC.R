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
##########                            Uncertainty analysis (AUC)                             ##########
#======================================================================================================
current_dir   <- getwd()
print(current_dir)
file_path     <- file.path(current_dir, "HTTK", "Bird", "TMX")
SA_hen        <- read.csv(paste(file_path, '/Laying Hen/Sensitivity.csv', sep = ''))
SA_duck       <- read.csv(paste(file_path, '/Duck/Sensitivity.csv', sep = ''))
SA_quail      <- read.csv(paste(file_path, '/Quail/Sensitivity.csv', sep = ''))

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

SA_AUC_quail          <- SA_AUC[SA_AUC$Species == 'Quail',]
SA_AUC_quail   <- SA_AUC_quail   %>%
  mutate(group = case_when(
    abs(AUC) >=0.5 ~ 'high',
    abs(AUC) < 0.5 & abs(AUC) >= 0.2 ~ 'medium',
    abs(AUC) < 0.2 ~ 'low'
  ))

print(SA_AUC_quail$parameter)
file_path     <- file.path(current_dir, "HTTK", "Bird", "TMX", "Quail")
write.csv(SA_AUC_quail, paste(file_path, '/SA_AUC_quail.csv', sep = ''), row.names = FALSE)


library(truncnorm) 
library(ggridges)     
library(ggplot2)  
library(tkrplot)
library(rriskDistributions)
library(EnvStats)
library(httk)
library(mc2d)
library(fitdistrplus)

source(file.path("HTTK", "Bird", "TMX", 'Quail', 'UASupporting_quail.R'))
set.seed(10) 

BW             <- 0.186                  # kg, body weight of 18 weeks 
BW_hen         <- 1.45                   # kg
BW_human       <- 70                     # kg
SDBW           <- 0.08

BW_mouse    <- 0.02
BW_rat      <- 0.25
BW_human    <- 70


#### mean value of parameters
Rblood2plasma_parent.mean          <- 0.9027282
Rblood2plasma_daughter.mean        <- 1.3597193
fub_daughter.mean                  <- 0.71
Clint_daughter.mean                <- 79.2898978
Qgfr_daughter.mean                 <- 15.3186752


#### std value of parameters
Rblood2plasma_parent.sd       <- 0.3 * Rblood2plasma_parent.mean
Rblood2plasma_daughter.sd     <- 0.3 * Rblood2plasma_daughter.mean
fub_daughter.sd               <- 0.3 * fub_daughter.mean
Clint_daughter.sd             <- 0.5 * Clint_daughter.mean 
Qgfr_daughter.sd              <- 0.3 * Qgfr_daughter.mean


#############     define distribution
# Normal distributions were implemented for physiological parameters; 
m.log.Rblood2plasma_parent             <- log(Rblood2plasma_parent.mean^2/(Rblood2plasma_parent.sd^2+Rblood2plasma_parent.mean^2)^0.5) 
sd.log.Rblood2plasma_parent            <- (log(1+Rblood2plasma_parent.sd^2/Rblood2plasma_parent.mean^2))^0.5 

m.log.Rblood2plasma_daughter           <- log(Rblood2plasma_daughter.mean^2/(Rblood2plasma_daughter.sd^2+Rblood2plasma_daughter.mean^2)^0.5) 
sd.log.Rblood2plasma_daughter          <- (log(1+Rblood2plasma_daughter.sd^2/Rblood2plasma_daughter.mean^2))^0.5 

m.log.fub_daughter             <- log(fub_daughter.mean^2/(fub_daughter.sd^2+fub_daughter.mean^2)^0.5) 
sd.log.fub_daughter            <- (log(1+fub_daughter.sd^2/fub_daughter.mean^2))^0.5 

m.log.Clint_daughter        <- log(Clint_daughter.mean^2/(Clint_daughter.sd^2+Clint_daughter.mean^2)^0.5)
sd.log.Clint_daughter       <- (log(1+Clint_daughter.sd^2/Clint_daughter.mean^2))^0.5

m.log.Qgfr_daughter          <- log(Qgfr_daughter.mean^2/(Qgfr_daughter.sd^2+Qgfr_daughter.mean^2)^0.5)
sd.log.Qgfr_daughter         <- (log(1+Qgfr_daughter.sd^2/Qgfr_daughter.mean^2))^0.5



####################
N=1000
set.seed(10) 

idata <- 
  tibble(ID=1:N) %>% 
  mutate(
    
    ## Chemical properties
    Clint_daughter = rlnormTrunc(
      N,
      meanlog = m.log.Clint_daughter,
      sdlog = sd.log.Clint_daughter,
      min = qlnorm(0.025, meanlog = m.log.Clint_daughter, sdlog = sd.log.Clint_daughter),
      max = qlnorm(0.975, meanlog = m.log.Clint_daughter, sdlog = sd.log.Clint_daughter)
    ),
    
    Qgfr_daughter = rlnormTrunc(
      N,
      meanlog = m.log.Qgfr_daughter,
      sdlog = sd.log.Qgfr_daughter,
      min = qlnorm(0.025, meanlog = m.log.Qgfr_daughter, sdlog = sd.log.Qgfr_daughter),
      max = qlnorm(0.975, meanlog = m.log.Qgfr_daughter, sdlog = sd.log.Qgfr_daughter)
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



###########   Single dose of 66 mg/kg BW  #########
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

#============================================================
######             Single dose population               #####
#============================================================
pars  <- c( ka                   = 12.1847922,          # d-1
            Clint_parent         = 692.6002544, 
            Clint_daughter       = 79.2898978, 
            
            Rblood2plasma_parent   = 0.9027282,   
            Rblood2plasma_daughter = 1.3597193,    
            
            fub_parent           = fub_parent,   
            fub_daughter         = fub_daughter, 
            
            # ratio
            Kgut2pu_parent       = Kgut2pu_parent ,      
            Kkidney2pu_parent    = Kkidney2pu_parent ,   
            Kliver2pu_parent     = Kliver2pu_parent,  
            Klung2pu_parent      = Klung2pu_parent ,      
            Krest2pu_parent      = Krest2pu_parent , 
            Kadipose2pu_parent   = Kadipose2pu_parent ,
            Kmuscle2pu_parent    = Kmuscle2pu_parent, 
            Krep2pu_parent       = Krep2pu_parent ,
            
            Kgut2pu_daughter       = Kgut2pu_daughter  ,      
            Kkidney2pu_daughter    = Kkidney2pu_daughter,   
            Kliver2pu_daughter     = Kliver2pu_daughter,    
            Klung2pu_daughter      = Klung2pu_daughter ,      
            Krest2pu_daughter      = Krest2pu_daughter , 
            Kadipose2pu_daughter   = Kadipose2pu_daughter,
            Kmuscle2pu_daughter    = Kmuscle2pu_daughter, 
            Krep2pu_daughter       = Krep2pu_daughter ,
            fraction_daughter      = 1,
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


UA_list  <- c('Clint_daughter',
              'Rblood2plasma_parent', 
              'Rblood2plasma_daughter',
              'fub_daughter',
              'Qgfr_daughter')



MC.AUC     = matrix(nrow = N_refined , ncol = length(UA_list))
colnames(MC.AUC)   <- UA_list

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
    
    # Modeling duration
    days          <- 2
    Time_max      <- days                                   
    StartTime     <- 0                                        # Time_d
    StopTime      <- Time_max                                 
    dt            <- 0.01
    Times         <- seq(StartTime, StopTime, dt)
    
    # dose
    dose          <- 100                                      # mg/kg bw
    Oral_input    <- dose * 1000 / MW_parent                  # umol/kg
    
    Dose_events    <- data.frame(var    = "Agutlumen_parent",
                                 time   =  0,
                                 value  = Oral_input,
                                 method = "add")    
    
    
    df  <- ode(y       = initState, 
               times   = Times, 
               func    = pbtk5cpt_repeated, 
               parms   = new_pars,
               atol    = 1e-6, rtol = 1e-8,
               events  = list(data = Dose_events),
               method  = 'lsoda')
    
    df                   <- as.data.frame(df)
    colnames(df)[1]      <- "Time"
    
    MC.AUC[i,j]     <- max(df$AUC_Cblood_parent + df$AUC_Cblood_daughter)

    
  }
}

# Calculate the 95th percentile for each column
percentile_95 <- apply(MC.AUC, 2, function(x) quantile(x, probs = 0.95))
# Calculate the median for each column
medians       <- apply(MC.AUC, 2, median)
# Create a new data frame with the 95th percentile and median
new_df        <- rbind(percentile_95, medians)

stat   <- t(new_df)
stat   <- as.data.frame(stat)
stat$ratio  <- round((stat$percentile_95 /stat$medians), 2)

stat <- stat %>%
  mutate(group = case_when(
    ratio >=2 ~ 'high',
    ratio < 2 & ratio >= 1.3 ~ 'medium',
    ratio < 1.3 ~ 'low'
  ))


file_path     <- file.path(current_dir, "HTTK", "Bird", "TMX", "Quail")
write.csv(stat, paste(file_path, '/UA_AUC_quail.csv', sep = ''), row.names = FALSE)

# end







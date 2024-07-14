library(stationaRy)
library(isdparser)
library(tidyr)
library(lubridate)
library(WriteXLS)
#install.packages(c("readxl","writexl")) 
library(readxl)
library(writexl)
library("tidyr")
library(dplyr)
library(Metrics)
library(PKNCA)
library(cowplot)
library(knitr)
library(reshape2)
rm(list = ls())
# PKNCA: https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
#dev.off()

MW_parent       <-   291.72        # for difen             <-  406.26
MW_daughter     <-   249.7         # for difen             <-  406.26

#####################             QUAIL                 #################    
#######################    Experimental data    #########################
Residue_2016_quail          <- read_excel("C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/InvivoData.xlsx", sheet = "Residue_2016_quail")

# blood
Pre_egg_female_320                   <- Residue_2016_quail[Residue_2016_quail$Dose_ppm == 320 & Residue_2016_quail$group == 'pre-egg' & Residue_2016_quail$Gender == 'female',]
Pre_egg_female_320$TMXBloodConc_umol <- Pre_egg_female_320$Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female_320$CLOBloodConc_umol <- Pre_egg_female_320$Blood_conc_ng.ml_CGA322704 /MW_daughter
df_Pre_egg_female_320                <- subset(Pre_egg_female_320, select = c('Time_d', "TMXBloodConc_umol", 'CLOBloodConc_umol')) 
df_Pre_egg_female_320                <- na.omit(df_Pre_egg_female_320)

Pre_egg_female_1000                   <- Residue_2016_quail[Residue_2016_quail$Dose_ppm == 1000 & Residue_2016_quail$group == 'pre-egg' & Residue_2016_quail$Gender == 'female',]
Pre_egg_female_1000$TMXBloodConc_umol <- Pre_egg_female_1000$Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female_1000$CLOBloodConc_umol <- Pre_egg_female_1000$Blood_conc_ng.ml_CGA322704 /MW_daughter
df_Pre_egg_female_1000                <- subset(Pre_egg_female_1000, select = c('Time_d', "TMXBloodConc_umol", 'CLOBloodConc_umol')) 
df_Pre_egg_female_1000                <- na.omit(df_Pre_egg_female_1000)

Pre_egg_female_2000                   <- Residue_2016_quail[Residue_2016_quail$Dose_ppm == 2000 & Residue_2016_quail$group == 'pre-egg' & Residue_2016_quail$Gender == 'female',]
Pre_egg_female_2000$TMXBloodConc_umol <- Pre_egg_female_2000$Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female_2000$CLOBloodConc_umol <- Pre_egg_female_2000$Blood_conc_ng.ml_CGA322704 /MW_daughter
df_Pre_egg_female_2000                <- subset(Pre_egg_female_2000, select = c('Time_d', "TMXBloodConc_umol", 'CLOBloodConc_umol')) 
df_Pre_egg_female_2000                <- na.omit(df_Pre_egg_female_2000)

Pre_egg_female_3200                   <- Residue_2016_quail[Residue_2016_quail$Dose_ppm == 3200 & Residue_2016_quail$group == 'pre-egg' & Residue_2016_quail$Gender == 'female',]
Pre_egg_female_3200$TMXBloodConc_umol <- Pre_egg_female_3200$Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female_3200$CLOBloodConc_umol <- Pre_egg_female_3200$Blood_conc_ng.ml_CGA322704 /MW_daughter
df_Pre_egg_female_3200                <- subset(Pre_egg_female_3200, select = c('Time_d', "TMXBloodConc_umol", 'CLOBloodConc_umol')) 
df_Pre_egg_female_3200                <- na.omit(df_Pre_egg_female_3200)

df_Pre_egg_female_320$Dose        <- '320 ppm'
df_Pre_egg_female_1000$Dose       <- '1000 ppm'
df_Pre_egg_female_2000$Dose       <- '2000 ppm'
df_Pre_egg_female_3200$Dose       <- '3200 ppm'

df_obs                     <- as.data.frame(rbind(df_Pre_egg_female_320, df_Pre_egg_female_1000, df_Pre_egg_female_2000, df_Pre_egg_female_3200))
colnames(df_obs)[1]        <- 'Time'
df_obs$Time                <- df_obs$Time- 14/24

# Assuming df is your data frame
df_obs_TMX   <- df_obs %>% dplyr::select(1, 2, 4)
df_obs_CTD   <- df_obs %>% dplyr::select(1, 3, 4)
df_obs_TMX$Compound  <- 'TMX'
df_obs_CTD$Compound  <- 'CTD'

# Rename the columns of df2 to match df1
names(df_obs_CTD) <- names(df_obs_TMX)

# Stack df1 and df2
df_obs                     <- bind_rows(df_obs_TMX, df_obs_CTD)
colnames(df_obs)[2]        <- 'BloodConc'

#######################   Modeling results   #######################
df_320ppm      <- read.csv('C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Quail/quail.adlib.320ppm.csv')
df_1000ppm     <- read.csv('C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Quail/quail.adlib.1000ppm.csv')
df_2000ppm     <- read.csv('C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Quail/quail.adlib.2000ppm.csv')
df_3200ppm     <- read.csv('C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Quail/quail.adlib.3200ppm.csv')

df_320ppm$Dose      <- '320 ppm'
df_1000ppm$Dose     <- '1000 ppm'
df_2000ppm$Dose     <- '2000 ppm'
df_3200ppm$Dose     <- '3200 ppm'

dat_320          <- as.data.frame(df_320ppm  %>% dplyr::select(Time, C_blood_parent,  C_blood_daughter, Dose))
dat_1000         <- as.data.frame(df_1000ppm %>% dplyr::select(Time, C_blood_parent,  C_blood_daughter, Dose))
dat_2000         <- as.data.frame(df_2000ppm %>% dplyr::select(Time, C_blood_parent,  C_blood_daughter, Dose))
dat_3200         <- as.data.frame(df_3200ppm %>% dplyr::select(Time, C_blood_parent,  C_blood_daughter, Dose))
dat              <- as.data.frame(rbind(dat_320, dat_1000, dat_2000, dat_3200))


# Assuming df is your data frame
dat_TMX   <- dat %>% dplyr::select(1, 2, 4)
dat_CTD   <- dat %>% dplyr::select(1, 3, 4)
dat_TMX$Compound  <- 'TMX'
dat_CTD$Compound  <- 'CTD'

# Rename the columns of df2 to match df1
names(dat_CTD) <- names(dat_TMX)

# Stack df1 and df2
dat                     <- bind_rows(dat_TMX, dat_CTD)
colnames(dat)[2]        <- 'C_Blood'



df                <- merge(df_obs, dat , by = c('Time', 'Compound', 'Dose'), all = FALSE)
df$log.obs        <- log10(df$BloodConc)
df$log.pred       <- log10(df$C_Blood)
df$Dose           <- factor(df$Dose, levels = c("320 ppm", "1000 ppm", "2000 ppm", "3200 ppm"))

#############################################################
#=======        plot  AUC diagnoal plot ================
# ## All data
# AUC$log.obs = log(AUC$obs,10)
# AUC$log.pre = log(AUC$pred,10)
# Cmax$log.obs = log(Cmax$obs,10)
# Cmax$log.pre = log(Cmax$pred,10)
# 
# AUC <- AUC[is.finite(AUC$log.obs) & is.finite(AUC$log.pre), ]
# 
# fit <- lm(log.obs ~ log.pre, data=AUC)                   # Fit the model
# AUC$residuals = residuals(fit)                           # Save the residual values
# AUC$predicted = predict(fit)                             # Save the predicted values
# AUC$OPratio = AUC$pred/AUC$obs          # Estimated the predicted-to-observed ratio


############## plot AUC #####################
p <- 
  ggplot(df, aes(x = log.obs, y = log.pred, color = Dose)) +
  geom_point(shape = 4, size = 3, stroke = 1.3) +
  #geom_point   ( shape=4, size = 3, color = '#00AFBB' ,stroke = 0.8)  +
  #scale_colour_manual(labels = c(expression(paste("Plasma (", mu,"g/L)")), "Urine (%)"),values=c('#2E9FDF', '#E7B800')) +
  #scale_shape_manual(labels = c(expression(paste("Plasma (", mu,"g/L)")), "Urine (%)"), values = c(19, 17)) +
  #scale_fill_discrete(labels = c(expression(paste("Plasma (", mu,"g/L)")),"Urine (%)")) +
  #scale_fill_manual(values=c('blue', 'green')) + # legend name
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", size = 0.6) +
  geom_abline (intercept = log10(3), linetype = "dashed",
               slope     = 1,
               color     ="black", size = 0.6) +
  geom_abline (intercept = log10(1/3), linetype = "dashed", 
               slope     = 1,
               color     ="black", size = 0.6) +
  geom_abline (intercept = log10(5), linetype = "dotted",
               slope     = 1,
               color     ="black", size = 0.8) +
  geom_abline (intercept = log10(0.2), linetype = "dotted", 
               slope     = 1,
               color     ="black", size = 0.8) +
  xlim(-1.2,1.5) + 
  ylim(-1.2,1.5)+
  #labs(y = expression(paste(Log[10] * ' [Predicted blood concentration (\u03BCmol/L)]'))) +
  #labs(x = expression(paste(Log[10] * ' [Observed blood concentration (\u03BCmol/L)]'))) +
  labs(y = expression(paste('Predicted blood concentration (\u03BCmol/L)'))) +
  labs(x = expression(paste('Observed blood concentration (\u03BCmol/L)'))) +
  theme_bw(base_size = 16)+
  theme(legend.position = c(0.17, 0.86),
        legend.title            = element_blank(),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm')) +
  annotate("text", x=Inf, y=-Inf, label="Northern Bobwhite", hjust=1.05, vjust=-1.5, size=5)
p

p1<-
  p +
  annotation_logticks() +
  scale_y_continuous(limits = c(-1.2,1.5), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-1.2,1.5),labels = scales::math_format(10^.x))
p1

ggsave("Evaluation_quail.tiff",scale = 1,
       plot = p1,
       path = "C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Quail/",
       width = 15, height = 15, units = "cm",dpi=600, compression = 'lzw')



#########################################################
###########  analysis of Cmax outlier ###################
#########################################################
value      <- 3 # one order of magnititude
df$comp    <- ifelse(df$BloodConc > df$C_Blood*value, '1',
                     ifelse(df$BloodConc < df$C_Blood/value, '1', '0'))
df_outlier <- df[df$comp == '1',]

print(paste0("Percentage of outlier (3 Fold, C): ", round(nrow(df_outlier)/nrow(df)*100, 1)))
print(paste0("Percentage of Non-outlier (3 Fold, C): ", 100 - round(nrow(df_outlier)/nrow(df)*100, 1)))
                                                           
   
value      <- 5 # one order of magnititude
df$comp    <- ifelse(df$BloodConc > df$C_Blood*value, '1',
                     ifelse(df$BloodConc < df$C_Blood/value, '1', '0'))
df_outlier <- df[df$comp == '1',]

print(paste0("Percentage of outlier (5 Fold, C): ", round(nrow(df_outlier)/nrow(df)*100, 1)))
print(paste0("Percentage of Non-outlier (5 Fold, C): ", 100 - round(nrow(df_outlier)/nrow(df)*100, 1)))                    

#####################              END  pf Quail              ########################


#####################               laying hen                 #######################

rm(list = ls())

MW_parent       <-   291.72        # for difen             <-  406.26
MW_daughter     <-   249.7         # for difen             <-  406.26
#######################    Experimental data    #########################
study1 <- "Metabolism_O_1998_hen"
study2 <- "Metabolism_T_1998_hen"
Metabolism_1998_7.66           <- read_excel("C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/InvivoData.xlsx", sheet = study1)
Metabolism_1998_7.92           <- read_excel("C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/InvivoData.xlsx", sheet = study2)

######### 7.66
# urine
Metabolism_1998_7.66_excreta              <- Metabolism_1998_7.66[Metabolism_1998_7.66$Matrix == 'Excreta',]
Metabolism_1998_7.66_excreta$TMX_excreta  <- Metabolism_1998_7.66_excreta$'%TMX' /100 *  Metabolism_1998_7.66_excreta$Concentration        # Unit: % of dose
Metabolism_1998_7.66_excreta$CLO_excreta  <- Metabolism_1998_7.66_excreta$'%CGA322704' /100 *  Metabolism_1998_7.66_excreta$Concentration  # Unit: % of dose
Metabolism_1998_7.66_excreta$Matrx        <- 'Excreta'
Metabolism_1998_7.66_excreta$Time         <- Metabolism_1998_7.66_excreta$Time_h/24 
df_Metabolism_1998_7.66_excreta           <- subset(Metabolism_1998_7.66_excreta, select = c('Time', "TMX_excreta", "CLO_excreta", "Matrix"))
colnames(df_Metabolism_1998_7.66_excreta) <- c('Time', "TMX", "CLO", "Matrix") 

# egg white
Metabolism_1998_7.66_egg_white                 <- Metabolism_1998_7.66[Metabolism_1998_7.66$Matrix == 'Egg_white',]
Metabolism_1998_7.66_egg_white$TMXConc_umolL   <- Metabolism_1998_7.66_egg_white$Concentration * 1000 * Metabolism_1998_7.66_egg_white$'%TMX' / 100 / MW_parent     # [umol/L]
Metabolism_1998_7.66_egg_white$CLOConc_umolL   <- Metabolism_1998_7.66_egg_white$Concentration * 1000 * Metabolism_1998_7.66_egg_white$'%CGA322704' / 100 / MW_parent
Metabolism_1998_7.66_egg_white$Time            <- Metabolism_1998_7.66_egg_white$Time_h/24 
df_Metabolism_1998_7.66_eggwhite               <- subset(Metabolism_1998_7.66_egg_white, select = c('Time', "TMXConc_umolL", "CLOConc_umolL", "Matrix"))
colnames(df_Metabolism_1998_7.66_eggwhite)     <- c('Time', "TMX", "CLO", "Matrix") 

# egg yolk
Metabolism_1998_7.66_egg_yolk                 <- Metabolism_1998_7.66[Metabolism_1998_7.66$Matrix == 'Egg_yolk',] 
Metabolism_1998_7.66_egg_yolk$TMXConc_umolL   <- Metabolism_1998_7.66_egg_yolk$Concentration * 1000 * Metabolism_1998_7.66_egg_yolk$'%TMX' / 100 / MW_parent
Metabolism_1998_7.66_egg_yolk$CLOConc_umolL   <- Metabolism_1998_7.66_egg_yolk$Concentration * 1000 * Metabolism_1998_7.66_egg_yolk$'%CGA322704' / 100 / MW_parent
Metabolism_1998_7.66_egg_yolk$Time            <- Metabolism_1998_7.66_egg_yolk$Time_h/24 
df_Metabolism_1998_7.66_yolk                  <- subset(Metabolism_1998_7.66_egg_yolk, select = c('Time', "TMXConc_umolL", "CLOConc_umolL", "Matrix"))
colnames(df_Metabolism_1998_7.66_yolk)        <- c('Time', "TMX", "CLO", "Matrix") 

# blood
Metabolism_1998_7.66_blood                    <- Metabolism_1998_7.66[Metabolism_1998_7.66$Matrix == 'Blood' ,]
Metabolism_1998_7.66_blood$TMXBloodConc_umolL <- Metabolism_1998_7.66_blood$Concentration * 1000 * Metabolism_1998_7.66_blood$'%TMX' / 100/ MW_parent
Metabolism_1998_7.66_blood$CLOBloodConc_umolL <- Metabolism_1998_7.66_blood$Concentration * 1000 * Metabolism_1998_7.66_blood$'%CGA322704'/100 / MW_parent
Metabolism_1998_7.66_blood$Time               <- Metabolism_1998_7.66_blood$Time_h/24 
df_Metabolism_1998_7.66_blood                 <- subset(Metabolism_1998_7.66_blood, select = c('Time',  "TMXBloodConc_umolL", "CLOBloodConc_umolL", "Matrix")) 
colnames(df_Metabolism_1998_7.66_blood)       <- c('Time', "TMX", "CLO", "Matrix") 

# liver
Metabolism_1998_7.66_liver                    <- Metabolism_1998_7.66[Metabolism_1998_7.66$Matrix == 'Liver' ,]
Metabolism_1998_7.66_liver$TMXliverConc_umolL <- Metabolism_1998_7.66_liver$Concentration * 1000 * Metabolism_1998_7.66_liver$'%TMX' / 100/ MW_parent
Metabolism_1998_7.66_liver$CLOliverConc_umolL <- Metabolism_1998_7.66_liver$Concentration * 1000 * Metabolism_1998_7.66_liver$'%CGA322704'/100 / MW_parent 
Metabolism_1998_7.66_liver$Time               <- Metabolism_1998_7.66_liver$Time_h/24 
df_Metabolism_1998_7.66_liver                 <- subset(Metabolism_1998_7.66_liver, select = c('Time', "TMXliverConc_umolL", "CLOliverConc_umolL", "Matrix")) 
colnames(df_Metabolism_1998_7.66_liver)       <- c('Time', "TMX", "CLO", "Matrix") 

# Format
df_obs_7.66                     <- as.data.frame(rbind(df_Metabolism_1998_7.66_excreta, df_Metabolism_1998_7.66_eggwhite, df_Metabolism_1998_7.66_yolk, 
                                                       df_Metabolism_1998_7.66_blood, df_Metabolism_1998_7.66_liver))
df_obs_7.66$Dose                <- '7.66 mg/kg/d'
df_obs_7.66                     <- na.omit(df_obs_7.66)


### 7.92
# urine
Metabolism_1998_7.92_excreta              <- Metabolism_1998_7.92[Metabolism_1998_7.92$Matrix == 'Excreta',]
Metabolism_1998_7.92_excreta$TMX_excreta  <- Metabolism_1998_7.92_excreta$'%TMX' /100 *  Metabolism_1998_7.92_excreta$Concentration        # Unit: % of dose
Metabolism_1998_7.92_excreta$CLO_excreta  <- Metabolism_1998_7.92_excreta$'%CGA322704' /100 *  Metabolism_1998_7.92_excreta$Concentration  # Unit: % of dose
Metabolism_1998_7.92_excreta$Matrx        <- 'Excreta'
Metabolism_1998_7.92_excreta$Time         <- Metabolism_1998_7.92_excreta$Time_h/24 
df_Metabolism_1998_7.92_excreta           <- subset(Metabolism_1998_7.92_excreta, select = c('Time', "TMX_excreta", "CLO_excreta", "Matrix"))
colnames(df_Metabolism_1998_7.92_excreta) <- c('Time', "TMX", "CLO", "Matrix") 

# egg white
Metabolism_1998_7.92_egg_white                 <- Metabolism_1998_7.92[Metabolism_1998_7.92$Matrix == 'Egg_white',]
Metabolism_1998_7.92_egg_white$TMXConc_umolL   <- Metabolism_1998_7.92_egg_white$Concentration * 1000 * Metabolism_1998_7.92_egg_white$'%TMX' / 100 / MW_parent     # [umol/L]
Metabolism_1998_7.92_egg_white$CLOConc_umolL   <- Metabolism_1998_7.92_egg_white$Concentration * 1000 * Metabolism_1998_7.92_egg_white$'%CGA322704' / 100 / MW_parent
Metabolism_1998_7.92_egg_white$Time            <- Metabolism_1998_7.92_egg_white$Time_h/24 
df_Metabolism_1998_7.92_eggwhite               <- subset(Metabolism_1998_7.92_egg_white, select = c('Time', "TMXConc_umolL", "CLOConc_umolL", "Matrix"))
colnames(df_Metabolism_1998_7.92_eggwhite)     <- c('Time', "TMX", "CLO", "Matrix") 

# egg yolk
Metabolism_1998_7.92_egg_yolk                 <- Metabolism_1998_7.92[Metabolism_1998_7.92$Matrix == 'Egg_yolk',] 
Metabolism_1998_7.92_egg_yolk$TMXConc_umolL   <- Metabolism_1998_7.92_egg_yolk$Concentration * 1000 * Metabolism_1998_7.92_egg_yolk$'%TMX' / 100 / MW_parent
Metabolism_1998_7.92_egg_yolk$CLOConc_umolL   <- Metabolism_1998_7.92_egg_yolk$Concentration * 1000 * Metabolism_1998_7.92_egg_yolk$'%CGA322704' / 100 / MW_parent
Metabolism_1998_7.92_egg_yolk$Time            <- Metabolism_1998_7.92_egg_yolk$Time_h/24 
df_Metabolism_1998_7.92_yolk                  <- subset(Metabolism_1998_7.92_egg_yolk, select = c('Time', "TMXConc_umolL", "CLOConc_umolL", "Matrix"))
colnames(df_Metabolism_1998_7.92_yolk)        <- c('Time', "TMX", "CLO", "Matrix") 

# blood
Metabolism_1998_7.92_blood                    <- Metabolism_1998_7.92[Metabolism_1998_7.92$Matrix == 'Blood' ,]
Metabolism_1998_7.92_blood$TMXBloodConc_umolL <- Metabolism_1998_7.92_blood$Concentration * 1000 * Metabolism_1998_7.92_blood$'%TMX' / 100/ MW_parent
Metabolism_1998_7.92_blood$CLOBloodConc_umolL <- Metabolism_1998_7.92_blood$Concentration * 1000 * Metabolism_1998_7.92_blood$'%CGA322704'/100 / MW_parent
Metabolism_1998_7.92_blood$Time               <- Metabolism_1998_7.92_blood$Time_h/24 
df_Metabolism_1998_7.92_blood                 <- subset(Metabolism_1998_7.92_blood, select = c('Time',  "TMXBloodConc_umolL", "CLOBloodConc_umolL", "Matrix")) 
colnames(df_Metabolism_1998_7.92_blood)       <- c('Time', "TMX", "CLO", "Matrix") 

# liver
Metabolism_1998_7.92_liver                    <- Metabolism_1998_7.92[Metabolism_1998_7.92$Matrix == 'Liver' ,]
Metabolism_1998_7.92_liver$TMXliverConc_umolL <- Metabolism_1998_7.92_liver$Concentration * 1000 * Metabolism_1998_7.92_liver$'%TMX' / 100/ MW_parent
Metabolism_1998_7.92_liver$CLOliverConc_umolL <- Metabolism_1998_7.92_liver$Concentration * 1000 * Metabolism_1998_7.92_liver$'%CGA322704'/100 / MW_parent 
Metabolism_1998_7.92_liver$Time               <- Metabolism_1998_7.92_liver$Time_h/24 
df_Metabolism_1998_7.92_liver                 <- subset(Metabolism_1998_7.92_liver, select = c('Time', "TMXliverConc_umolL", "CLOliverConc_umolL", "Matrix")) 
colnames(df_Metabolism_1998_7.92_liver)       <- c('Time', "TMX", "CLO", "Matrix") 

# Format
df_obs_7.92                     <- as.data.frame(rbind(df_Metabolism_1998_7.92_excreta, df_Metabolism_1998_7.92_eggwhite, df_Metabolism_1998_7.92_yolk, 
                                                       df_Metabolism_1998_7.92_blood, df_Metabolism_1998_7.92_liver))
df_obs_7.92$Dose                <- '7.92 mg/kg/d'
df_obs_7.92                     <- na.omit(df_obs_7.92)

df_obs                          <- as.data.frame(rbind(df_obs_7.66, df_obs_7.92))


# Assuming df is your data frame
df_obs_TMX           <- df_obs %>% dplyr::select(1, 2, 4, 5)
df_obs_CTD           <- df_obs %>% dplyr::select(1, 3, 4, 5)
df_obs_TMX$Compound  <- 'TMX'
df_obs_CTD$Compound  <- 'CTD'

# Rename the columns of df2 to match df1
names(df_obs_CTD)    <- names(df_obs_TMX)

# Stack df1 and df2
df_obs                     <- bind_rows(df_obs_TMX, df_obs_CTD)
colnames(df_obs)[2]        <- 'Value'


#######################   Modeling results   #######################
# 7.66
df_7.66_eggwhite          <- read.csv('C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/hen.eggwhite.oral.7.66mgkg.csv')
df_7.66_eggyolk.CTD       <- read.csv('C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/hen.eggyolk.CTD.oral.7.66mgkg.csv')
df_7.66_eggyolk.TMX       <- read.csv('C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/hen.eggyolk.TMX.oral.7.66mgkg.csv')
df_7.66_rest              <- read.csv('C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/hen.oral.7.66mgkg.csv')


df_7.66_eggwhite$Dose        <- '7.66 mg/kg/d'
df_7.66_eggyolk.CTD$Dose     <- '7.66 mg/kg/d'
df_7.66_eggyolk.TMX$Dose     <- '7.66 mg/kg/d'
df_7.66_rest$Dose            <- '7.66 mg/kg/d'

df_7.66_eggwhite$Time        <- df_7.66_eggwhite$Time - 0.3
df_7.66_eggyolk.CTD$Time     <- df_7.66_eggyolk.CTD$Time  - 0.3
df_7.66_eggyolk.TMX$Time     <- df_7.66_eggyolk.TMX$Time  - 0.3
df_7.66_rest$Time            <- df_7.66_rest$Time  - 0.3

Oral_input                     <- 7.66 * 1000 / MW_parent                  # umol/kg/d, STUDY 2
df_7.66_rest$excreta_parent    <-  df_7.66_rest$Aurine_parent / (Oral_input  * 4) * 100
df_7.66_rest$excreta_daughter  <-  df_7.66_rest$Aurine_daughter / (Oral_input  * 4) * 100

df_7.66_rest         <- as.data.frame(df_7.66_rest %>% dplyr::select(Time, C_blood_parent,  C_blood_daughter, C_liver_parent, C_liver_daughter, 
                                                                     excreta_parent, excreta_daughter, Dose))
df_7.66_eggwhite     <- as.data.frame(df_7.66_eggwhite  %>% dplyr::select(Time, C_white_parent,  C_white_daughter, Dose))
df_7.66_eggyolk.CTD  <- as.data.frame(df_7.66_eggyolk.CTD %>% dplyr::select(Time, value, Dose))
df_7.66_eggyolk.TMX  <- as.data.frame(df_7.66_eggyolk.TMX %>% dplyr::select(Time, value, Dose))
df_7.66_eggyolk.CTD$Matrix    <- 'Egg_yolk'
df_7.66_eggyolk.TMX$Matrix    <- 'Egg_yolk' 
df_7.66_eggyolk.CTD$Compound  <- 'CTD'
df_7.66_eggyolk.TMX$Compound  <- 'TMX' 
colnames(df_7.66_eggyolk.TMX)[2]  <- 'Value_pred'
colnames(df_7.66_eggyolk.CTD)[2]  <- 'Value_pred'

library(reshape2)
dfm_7.66_rest <- melt(df_7.66_rest, id.vars = c("Time", "Dose"))
dfm_7.66_rest$Matrix  <- ifelse(grepl("C_blood_parent", dfm_7.66_rest$variable), "Blood", 
                             ifelse(grepl("C_blood_daughter", dfm_7.66_rest$variable), "Blood", 
                                    ifelse(grepl("C_liver_parent", dfm_7.66_rest$variable), "Liver",
                                           ifelse(grepl("C_liver_daughter", dfm_7.66_rest$variable), "Liver", 
                                                  ifelse(grepl("excreta_parent", dfm_7.66_rest$variable), "Excreta", 
                                                         ifelse(grepl("excreta_daughter", dfm_7.66_rest$variable), "Excreta",dfm_7.66_rest$variable))))))

dfm_7.66_rest$Compound <- ifelse(grepl("C_blood_parent", dfm_7.66_rest$variable), "TMX", 
                               ifelse(grepl("C_blood_daughter", dfm_7.66_rest$variable), "CTD", 
                                      ifelse(grepl("C_liver_parent", dfm_7.66_rest$variable), "TMX",
                                             ifelse(grepl("C_liver_daughter", dfm_7.66_rest$variable), "CTD", 
                                                    ifelse(grepl("excreta_parent", dfm_7.66_rest$variable), "TMX", 
                                                           ifelse(grepl("excreta_daughter", dfm_7.66_rest$variable), "CTD",dfm_7.66_rest$variable))))))

colnames(dfm_7.66_rest)[4]  <- 'Value_pred'


dfm_7.66_eggwhite            <- melt(df_7.66_eggwhite , id.vars = c("Time", "Dose"))
dfm_7.66_eggwhite$Compound   <- ifelse(grepl("C_white_parent", dfm_7.66_eggwhite$variable), "TMX", 
                                 ifelse(grepl("C_white_daughter", dfm_7.66_eggwhite$variable), "CTD", dfm_7.66_eggwhite$variable))
colnames(dfm_7.66_eggwhite)[4]  <- 'Value_pred'
dfm_7.66_eggwhite$Matrix        <- 'Egg_white'


df_7.66_eggyolk.TMX  <- as.data.frame(df_7.66_eggyolk.TMX %>% dplyr::select(Time, Value_pred, Dose, Matrix, Compound))
df_7.66_eggyolk.CTD  <- as.data.frame(df_7.66_eggyolk.CTD %>% dplyr::select(Time, Value_pred, Dose, Matrix, Compound))
dfm_7.66_eggwhite    <- as.data.frame(dfm_7.66_eggwhite   %>% dplyr::select(Time, Value_pred, Dose, Matrix, Compound))
dfm_7.66_rest        <- as.data.frame(dfm_7.66_rest       %>% dplyr::select(Time, Value_pred, Dose, Matrix, Compound))


# 7.92
df_7.92_eggwhite          <- read.csv('C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/hen.eggwhite.oral.7.92mgkg.csv')
df_7.92_eggyolk.CTD       <- read.csv('C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/hen.eggyolk.CTD.oral.7.92mgkg.csv')
df_7.92_eggyolk.TMX       <- read.csv('C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/hen.eggyolk.TMX.oral.7.92mgkg.csv')
df_7.92_rest              <- read.csv('C:/xxx/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/hen.oral.7.92mgkg.csv')


df_7.92_eggwhite$Dose        <- '7.92 mg/kg/d'
df_7.92_eggyolk.CTD$Dose     <- '7.92 mg/kg/d'
df_7.92_eggyolk.TMX$Dose     <- '7.92 mg/kg/d'
df_7.92_rest$Dose            <- '7.92 mg/kg/d'

df_7.92_eggwhite$Time        <- df_7.92_eggwhite$Time - 0.3
df_7.92_eggyolk.CTD$Time     <- df_7.92_eggyolk.CTD$Time  - 0.3
df_7.92_eggyolk.TMX$Time     <- df_7.92_eggyolk.TMX$Time  - 0.3
df_7.92_rest$Time            <- df_7.92_rest$Time  - 0.3

Oral_input                     <- 7.92 * 1000 / MW_parent                  # umol/kg/d, STUDY 2
df_7.92_rest$excreta_parent    <-  df_7.92_rest$Aurine_parent / (Oral_input  * 4) * 100
df_7.92_rest$excreta_daughter  <-  df_7.92_rest$Aurine_daughter / (Oral_input  * 4) * 100

df_7.92_rest         <- as.data.frame(df_7.92_rest %>% dplyr::select(Time, C_blood_parent,  C_blood_daughter, C_liver_parent, C_liver_daughter, 
                                                                     excreta_parent, excreta_daughter, Dose))
df_7.92_eggwhite     <- as.data.frame(df_7.92_eggwhite  %>% dplyr::select(Time, C_white_parent,  C_white_daughter, Dose))
df_7.92_eggyolk.CTD  <- as.data.frame(df_7.92_eggyolk.CTD %>% dplyr::select(Time, value, Dose))
df_7.92_eggyolk.TMX  <- as.data.frame(df_7.92_eggyolk.TMX %>% dplyr::select(Time, value, Dose))
df_7.92_eggyolk.CTD$Matrix    <- 'Egg_yolk'
df_7.92_eggyolk.TMX$Matrix    <- 'Egg_yolk' 
df_7.92_eggyolk.CTD$Compound  <- 'CTD'
df_7.92_eggyolk.TMX$Compound  <- 'TMX' 
colnames(df_7.92_eggyolk.TMX)[2]  <- 'Value_pred'
colnames(df_7.92_eggyolk.CTD)[2]  <- 'Value_pred'

library(reshape2)
dfm_7.92_rest <- melt(df_7.92_rest, id.vars = c("Time", "Dose"))
dfm_7.92_rest$Matrix  <- ifelse(grepl("C_blood_parent", dfm_7.92_rest$variable), "Blood", 
                                ifelse(grepl("C_blood_daughter", dfm_7.92_rest$variable), "Blood", 
                                       ifelse(grepl("C_liver_parent", dfm_7.92_rest$variable), "Liver",
                                              ifelse(grepl("C_liver_daughter", dfm_7.92_rest$variable), "Liver", 
                                                     ifelse(grepl("excreta_parent", dfm_7.92_rest$variable), "Excreta", 
                                                            ifelse(grepl("excreta_daughter", dfm_7.92_rest$variable), "Excreta",dfm_7.92_rest$variable))))))

dfm_7.92_rest$Compound <- ifelse(grepl("C_blood_parent", dfm_7.92_rest$variable), "TMX", 
                                 ifelse(grepl("C_blood_daughter", dfm_7.92_rest$variable), "CTD", 
                                        ifelse(grepl("C_liver_parent", dfm_7.92_rest$variable), "TMX",
                                               ifelse(grepl("C_liver_daughter", dfm_7.92_rest$variable), "CTD", 
                                                      ifelse(grepl("excreta_parent", dfm_7.92_rest$variable), "TMX", 
                                                             ifelse(grepl("excreta_daughter", dfm_7.92_rest$variable), "CTD",dfm_7.92_rest$variable))))))

colnames(dfm_7.92_rest)[4]  <- 'Value_pred'


dfm_7.92_eggwhite            <- melt(df_7.92_eggwhite , id.vars = c("Time", "Dose"))
dfm_7.92_eggwhite$Compound   <- ifelse(grepl("C_white_parent", dfm_7.92_eggwhite$variable), "TMX", 
                                       ifelse(grepl("C_white_daughter", dfm_7.92_eggwhite$variable), "CTD", dfm_7.92_eggwhite$variable))
colnames(dfm_7.92_eggwhite)[4]  <- 'Value_pred'
dfm_7.92_eggwhite$Matrix        <- 'Egg_white'


df_7.92_eggyolk.TMX  <- as.data.frame(df_7.92_eggyolk.TMX %>% dplyr::select(Time, Value_pred, Dose, Matrix, Compound))
df_7.92_eggyolk.CTD  <- as.data.frame(df_7.92_eggyolk.CTD %>% dplyr::select(Time, Value_pred, Dose, Matrix, Compound))
dfm_7.92_eggwhite    <- as.data.frame(dfm_7.92_eggwhite   %>% dplyr::select(Time, Value_pred, Dose, Matrix, Compound))
dfm_7.92_rest        <- as.data.frame(dfm_7.92_rest       %>% dplyr::select(Time, Value_pred, Dose, Matrix, Compound))



# merge
dat                  <- as.data.frame(rbind(df_7.66_eggyolk.TMX, df_7.66_eggyolk.CTD, dfm_7.66_eggwhite, dfm_7.66_rest,
                                            df_7.92_eggyolk.TMX, df_7.92_eggyolk.CTD, dfm_7.92_eggwhite, dfm_7.92_rest))

df                <- merge(df_obs, dat , by = c('Time', 'Compound', 'Matrix', 'Dose'), all = FALSE)
df                <- df[df$Value != 0,]
df$log.obs        <- log10(df$Value)
df$log.pred       <- log10(df$Value_pred)
df$Dose           <- factor(df$Dose, levels = c("7.66 mg/kg/d", "7.92 mg/kg/d"))

############## plot AUC #####################
p <- 
  ggplot(df, aes(x = log.obs, y = log.pred, color = Dose)) +
  geom_point(shape = 4, size = 3, stroke = 1.3) +
  #geom_point   ( shape=4, size = 3, color = '#00AFBB' ,stroke = 0.8)  +
  #scale_colour_manual(labels = c(expression(paste("Plasma (", mu,"g/L)")), "Urine (%)"),values=c('#2E9FDF', '#E7B800')) +
  #scale_shape_manual(labels = c(expression(paste("Plasma (", mu,"g/L)")), "Urine (%)"), values = c(19, 17)) +
  #scale_fill_discrete(labels = c(expression(paste("Plasma (", mu,"g/L)")),"Urine (%)")) +
  #scale_fill_manual(values=c('blue', 'green')) + # legend name
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", size = 0.6) +
  geom_abline (intercept = log10(3), linetype = "dashed",
               slope     = 1,
               color     ="black", size = 0.6) +
  geom_abline (intercept = log10(1/3), linetype = "dashed", 
               slope     = 1,
               color     ="black", size = 0.6) +
  geom_abline (intercept = log10(5), linetype = "dotted",
               slope     = 1,
               color     ="black", size = 0.8) +
  geom_abline (intercept = log10(0.2), linetype = "dotted", 
               slope     = 1,
               color     ="black", size = 0.8) +
  xlim(-1.5,1.5) + 
  ylim(-1.5,1.5)+
  #labs(y = expression(paste(Log[10] * ' [Predicted blood concentration (\u03BCmol/L)]'))) +
  #labs(x = expression(paste(Log[10] * ' [Observed blood concentration (\u03BCmol/L)]'))) +
  labs(y = expression(paste('Predicted value (\u03BCmol/L or %)'))) +
  labs(x = expression(paste('Observed value (\u03BCmol/L or %)'))) +
  theme_bw(base_size = 16)+
  theme(legend.position = c(0.17, 0.91),
        legend.title            = element_blank(),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm')) +
  annotate("text", x=Inf, y=-Inf, label="Domestic Chicken", hjust=1.12, vjust=-1.5, size=5)
p

p1<-
  p +
  annotation_logticks() +
  scale_y_continuous(limits = c(-1.5,1.5), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-1.5,1.5),labels = scales::math_format(10^.x))
p1

ggsave("Evaluation_hen.tiff",scale = 1,
       plot = p1,
       path = "C:/Users/s1036120/OneDrive - Syngenta/HTTK/Bird/TMX/Plots/Hen/",
       width = 15, height = 15, units = "cm", dpi=600, compression = 'lzw')

max(df$log.obs)
min(df$log.obs)
max(df$log.pred)
min(df$log.pred)
#########################################################
###########  analysis of Cmax outlier ###################
#########################################################
value      <- 3 # one order of magnititude
df$comp    <- ifelse(df$Value > df$Value_pred*value, '1',
                     ifelse(df$Value< df$Value_pred/value, '1', '0'))
df_outlier <- df[df$comp == '1',]

print(paste0("Percentage of outlier (3 Fold, C): ", round(nrow(df_outlier)/nrow(df)*100, 1)))
print(paste0("Percentage of Non-outlier (3 Fold, C): ", 100 - round(nrow(df_outlier)/nrow(df)*100, 1)))


value      <- 5 # one order of magnititude
df$comp    <- ifelse(df$Value > df$Value_pred*value, '1',
                     ifelse(df$Value< df$Value_pred/value, '1', '0'))
df_outlier <- df[df$comp == '1',]

print(paste0("Percentage of outlier (5 Fold, C): ", round(nrow(df_outlier)/nrow(df)*100, 1)))
print(paste0("Percentage of Non-outlier (5 Fold, C): ", 100 - round(nrow(df_outlier)/nrow(df)*100, 1)))                    












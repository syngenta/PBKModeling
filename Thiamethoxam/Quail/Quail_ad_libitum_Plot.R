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
SPECIES                 <- 'Quail'
MW_parent               <- 291.72
MW_daughter             <- 249.7

current_dir   <- getwd()
print(current_dir)
file_path     <- file.path(current_dir, "HTTK", "Bird", "TMX", 'Plots', 'Quail')

df_320     <- read.csv(paste(file_path, '/quail.adlib.320ppm.csv', sep =''))
df_1000    <- read.csv(paste(file_path, '/quail.adlib.1000ppm.csv', sep =''))
df_320     <- subset(df_320, select = c('Time' , 'C_blood_parent', 'C_blood_daughter'))
df_1000    <- subset(df_1000, select = c('Time' , 'C_blood_parent', 'C_blood_daughter'))

df_320_TMX   <-  subset(df_320, select = c('Time' , 'C_blood_parent'))
df_320_CTD   <-  subset(df_320, select = c('Time' , 'C_blood_daughter'))
df_1000_TMX  <-  subset(df_1000, select = c('Time' , 'C_blood_parent'))
df_1000_CTD  <-  subset(df_1000, select = c('Time' , 'C_blood_daughter'))

colnames(df_320_TMX)[2]   <- 'Conc_umol'
colnames(df_320_CTD)[2]   <- 'Conc_umol'
colnames(df_1000_TMX)[2]  <- 'Conc_umol'
colnames(df_1000_CTD)[2]  <- 'Conc_umol'

df_320_TMX$Compound    <-  'TMX'
df_320_CTD$Compound    <-  'CTD'
df_1000_TMX$Compound   <-  'TMX'
df_1000_CTD$Compound   <-  'CTD'

df_320    <- as.data.frame(rbind(df_320_TMX, df_320_CTD))
df_1000   <- as.data.frame(rbind(df_1000_TMX, df_1000_CTD))
df_320$Dose   <- '320 ppm/d'
df_1000$Dose  <- '1000 ppm/d'
df            <- as.data.frame(rbind(df_320, df_1000))

df$Dose       <- factor(df$Dose, level = c('320 ppm/d', '1000 ppm/d'))
df$Compound   <- factor(df$Compound, level = c('TMX', 'CTD'))

file_path                   <- file.path(current_dir, "HTTK", "Bird", "TMX", 'InvivoData.xlsx')
Residue_2016_quail          <- read_excel(file_path , sheet = "Residue_2016_quail")

# blood
#dose_ppm      <- 1000  
Pre_egg_female                      <- Residue_2016_quail[Residue_2016_quail$group == 'pre-egg' & Residue_2016_quail$Gender == 'female',]
Pre_egg_female                      <- Pre_egg_female[Pre_egg_female$Dose_ppm == 1000 |Pre_egg_female$Dose_ppm == 320,]
Pre_egg_female$TMXBloodConc_umol    <- Pre_egg_female$Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female$CLOBloodConc_umol    <- Pre_egg_female$Blood_conc_ng.ml_CGA322704 /MW_daughter
Pre_egg_female$TMXBloodConc_umol_SD <- Pre_egg_female$SD.Blood_conc_ng.ml_TMX /MW_parent
Pre_egg_female$CLOBloodConc_umol_SD <- Pre_egg_female$SD.Blood_conc_ng.ml_CGA322704 /MW_daughter
df_Pre_egg_female                   <- subset(Pre_egg_female, select = c('Time_d','Dose_ppm',  "TMXBloodConc_umol", 'TMXBloodConc_umol_SD',
                                                                           'CLOBloodConc_umol', "CLOBloodConc_umol_SD")) 
df_Pre_egg_female      <- na.omit(df_Pre_egg_female)
df_Pre_egg_female_TMX  <- subset(df_Pre_egg_female, select = c("Time_d", "Dose_ppm", "TMXBloodConc_umol", "TMXBloodConc_umol_SD" ))
df_Pre_egg_female_CTD  <- subset(df_Pre_egg_female, select = c("Time_d", "Dose_ppm", "CLOBloodConc_umol", "CLOBloodConc_umol_SD" ))

df_Pre_egg_female_TMX$Compound      <- 'TMX'
df_Pre_egg_female_CTD$Compound      <- 'CTD'
colnames(df_Pre_egg_female_TMX)[3]  <- 'Conc_umol'
colnames(df_Pre_egg_female_TMX)[4]  <- 'Conc_umol_SD'
colnames(df_Pre_egg_female_CTD)[3]  <- 'Conc_umol'
colnames(df_Pre_egg_female_CTD)[4]  <- 'Conc_umol_SD'

df_Pre_egg        <- as.data.frame(rbind(df_Pre_egg_female_TMX, df_Pre_egg_female_CTD))
df_Pre_egg$Dose   <- paste0(df_Pre_egg$Dose_ppm, ' ppm/d')

df_Pre_egg$Dose       <- factor(df_Pre_egg$Dose, level = c('320 ppm/d', '1000 ppm/d'))
df_Pre_egg$Compound   <- factor(df_Pre_egg$Compound, level = c('TMX', 'CTD'))


##########   OVER VIREW PLOT  ################

ggplot(df_Pre_egg) +
  geom_point(aes((Time_d - 14/24), Conc_umol, color = Compound), size = 2, shape = 21, fill = 'white', stroke = 0.8) +
  geom_errorbar(aes((Time_d - 14/24),
                    ymin = pmax(Conc_umol - Conc_umol_SD, 0),
                    ymax = Conc_umol + Conc_umol_SD, color = Compound), width=0.12, linewidth = 0.6) +
  geom_line(data = df, aes(Time, Conc_umol, color = Compound), lwd=0.3) +
  #facet_wrap(~Class) +  
  #facet_grid(Matrix~Class, space="free", scales="free")+
  facet_grid(Dose~Compound, space="fixed", scales="free")+
  theme_bw(base_size = 16) + theme(legend.position = 'none') + 
  theme(strip.text = element_text(size = 16)) +
  scale_color_brewer(palette = "Set1") + xlim(0,7)+
  labs(title = paste("Northern bobwhite [daily dietary exposure for 28 consecutive days]"))+
  ylab(expression("Blood concentration ("*mu*"mol/L)")) + xlab('Time (d)')

file_path  <- file.path(current_dir, "HTTK", "Bird", "TMX", "Quail", 'plot.bobwhite.Blood_35d.tiff')
ggsave(file_path,
       dpi = 300, width = 26, height = 18,units = "cm",compression = 'lzw')


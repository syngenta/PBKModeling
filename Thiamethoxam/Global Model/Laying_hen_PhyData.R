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

#==================================================================================================
#                                      Species: Laying Hen                                        # 
#                                     No of compartment: 5                                        #
#==================================================================================================
# Reference
# PKNCA: https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk: https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk

#=============================================================================================================
#                                   Generic Maternal PBPK Model Parameters                                   #
#                                  Tissue and physiological data for species                                 #
#=============================================================================================================
R_rat                   <- 0.11             # 0.2cm # radius of rat jejunum; https://github.com/wfsrqivive/rat_PBK/blob/main/model/generate_input_data.r
R_human                 <- 1.25       
#######################################################################
###               tissue weight for laying hen  (% bw)              ### 
###               tissue volume for laying hen  (L/kg)              ###  
#######################################################################
egg_white             <- 33.6            # g internal study
egg_yolk              <- 12.4            # g
# Table 7, Scanes el al 2022 Quantitative morphometric,  CV information provided 
Hematocrit_fraction     <- 27.1  
SD_Hematocrit_fraction  <- 0.86/27.1

# L/kg (organ weight as %bw)
# Table 7, Scanes el al 2022 Quantitative morphometric (L/kg)
Vblood      <- 6.3   /100   # Table 7, Scanes el al 2022 Quantitative morphometric %
Vart        <- Vblood * 0.16              # L/kg
Vven        <- Vblood * 0.595 
Vplasma     <- Vblood * (1-Hematocrit_fraction/100)

Vheart      <- 0.45  /100   # Table 7, Scanes el al 2022 Quantitative morphometric
Vspleen     <- 0.16  /100   # Table 7, Scanes el al 2022 Quantitative morphometric

Vliver      <- 2.33  /100   # Table 7, Scanes el al 2022 Quantitative morphometric
Vkidney     <- 0.8   /100   # Table 7, Scanes el al 2022 Quantitative morphometric
Vovary      <- 2.38  /100   # Table 7, Scanes el al 2022 Quantitative morphometric
Voviduct    <- 2.84  /100   # Table 7, Scanes el al 2022 Quantitative morphometric
Vgizzard    <- 3.21  /100   # Table 7, Scanes el al 2022 Quantitative morphometric
Vcrop       <- (0.055 + 0.19)/100      # Table 3 Wang 2021 for laying hen crop + proventriculus

#Vskin       <- 12    /100   # Fereidoun et al. (2007) Mean Percentage of Skin and Visible Fat in 10 Chicken Carcass Weight 8-20%;  Table 1. Wang et al 2020 chicken and turkey
Vbrain      <- 0.2   /100   # # Table 3 L.S. Lautz et al. 2020 (Environmental International)
Vlung       <- 0.8   /100   # # Table 3 L.S. Lautz et al. 2020 (Environmental International)
Vmuscle     <- 39.8  /100   # # Table 3 L.S. Lautz et al. 2020 (Environmental International)
Vadipose    <- 13.6  /100   # # Average of Table 3 L.S. Lautz et al. 2020 (Environmental International) might include skin weight (12.1); Wang 2021: 15.1 (female)
Vint        <- 5.8   /100   # # intestine, Table 3 L.S. Lautz et al. 2020 (Environmental International)
Vbone       <- 20.3  /100   # # Table 3 L.S. Lautz et al. 2020 (Environmental International) Carcass

Vrep        <- Vovary + Voviduct 
Vgut        <- Vcrop + Vint + Vgizzard

Vrest       <- 1 - (Vgut + Vkidney + Vliver + Vadipose + Vmuscle + Vrep + Vlung + Vart + Vven)   # L/kg BW
Vrest_a     <- Vrest - (Vbone + Vbrain + Vheart + Vspleen)

## volume and surface area of vascular space in adipose
fadipose    <- 0.194                                             # Table S1 Weixiao Cheng ES&T 2017 A permeability-limited PBPK model for PFOA in male rats
V_Vadipose  <- Vadipose * fadipose                               # volume of vascular space 
V_Tadipose  <- Vadipose * (1 - fadipose)                         # volume of extravascular space 
S_adipose   <- 70                                                # cm2/g Table S1 Weixiao Cheng ES&T 2017 A permeability-limited PBPK model for PFOA in male rats


# L/kg
SDblood      <- 1.33  /100          # # per organ weight, Table 7, Scanes el al 2022 Quantitative morphometric 
SDheart      <- 0.04  /100          # # per organ weight, Table 7, Scanes el al 2022 Quantitative morphometric 
SDliver      <- 0.1   /100          # # per organ weight, Table 7, Scanes el al 2022 Quantitative morphometric 
SDkidney     <- 0.098 /100          # # per organ weight, Table 7, Scanes el al 2022 Quantitative morphometric 
SDspleen     <- 0.022 /100          # # per organ weight, Table 7, Scanes el al 2022 Quantitative morphometric 
SDgizzard    <- 1.10  /100          # # per organ weight, Table 7, Scanes el al 2022 Quantitative morphometric 
SDsmallint   <- 1.22  /100          # # per organ weight, Table 7, Scanes el al 2022 Quantitative morphometric 
SDovary      <- 0.21  /100          # # per organ weight, Table 7, Scanes el al 2022 Quantitative morphometric 
SDoviduct    <- 0.22  /100          # # per organ weight, Table 7, Scanes el al 2022 Quantitative morphometric 

SDbrain      <- 0.2  * 0.49 /100    # # Tissue weight fraction * CV (%/100) , Table 3 L.S. Lautz et al. 2020 (Environmental International)
SDlung       <- 0.8  * 0.27 /100    # # Tissue weight fraction * CV (%/100) , Table 3 L.S. Lautz et al. 2020 (Environmental International)
SDmuscle     <- 39.8 * 0.12 /100    # # Tissue weight fraction * CV (%/100) , Table 3 L.S. Lautz et al. 2020 (Environmental International)
SDadipose    <- 12.1 * 0.17 /100    # # Tissue weight fraction * CV (%/100) , Table 3 L.S. Lautz et al. 2020 (Environmental International)
SDint        <- 5.8  * 0.16 /100    # # Tissue weight fraction * CV (%/100) , Table 3 L.S. Lautz et al. 2020 (Environmental International)
SDbone       <- 20.3 * 0.36 /100    # Carcass # # Tissue weight fraction * CV (%/100) , Table 3 L.S. Lautz et al. 2020 (Environmental International)


###################################################################
###                   tissue flows (L/d/kg BW)                  ###  
###################################################################
# Table 3 L.S. Lautz et al. 2020 (Environmental International)
Qart        <- 0.34 * 60 * 24      # L/d/kg BW 20.4 * 24

# Cv can be 30% for all the organs according to Table 3 L.S. Lautz et al. 2020 (Environmental International); not used
Qovary      <- 1.9                 # 2.59%, derived from  Table 8, Scanes el al 2022 Quantitative morphometric 
Qoviduct    <- 5.88                # 7.94%, derived from  Table 8, Scanes el al 2022 Quantitative morphometric 

Qbrain      <- 0.4  / 100 * Qart   # % of cardiac output * cardiac output, Table 3 L.S. Lautz et al. 2020 (Environmental International)
Qliver      <- 10.01 / 100 * Qart  # % of cardiac output * cardiac output, average Table 16 Wang et al 2020 (13.43); Table 3 L.S. Lautz et al. 2020 (Environmental International); 6.6
Qkidney     <- 15.7 / 100 * Qart   # % of cardiac output * cardiac output, average of Table 25 Wang et al 2020 (20.12); Table 3 L.S. Lautz et al. 2020 (Environmental International) (11.4)
Qmuscle     <- 13.7 / 100 * Qart   # % of cardiac output * cardiac output, Average of Table 16 Wang et al 2020 (7.64); Table 3 L.S. Lautz et al. 2020 (Environmental International) 19.8
Qadipose    <- 1.5  / 100 * Qart   # % of cardiac output * cardiac output, Table 3 L.S. Lautz et al. 2020 (Environmental International)
Qheart      <- 5.5  / 100 * Qart   # % of cardiac output * cardiac output, Table 3 L.S. Lautz et al. 2020 (Environmental International)
Qrep        <- 14.3 / 100 * Qart   # % of cardiac output * cardiac output, Table 3 L.S. Lautz et al. 2020 (Environmental International) (14.3); Tabe 25 Wang 2021 (12.14)
Qspleen     <- 1.6  / 100 * Qart   # Table 8, Scanes el al 2022 Quantitative morphometric 
Qgizzard    <- 1.8  / 100 * Qart   # Table 8, Scanes el al 2022 Quantitative morphometric 
Qint        <- 17.7 / 100 * Qart   # intestine, % of cardiac output * cardiac output, Table 3 L.S. Lautz et al. 2020 (Environmental International)
#Qlung       <- 3
Qbone       <- 12.4 / 100 * Qart   # Skeletal Blood Flow in Bone Repair and Maintenance

Qgut        <- Qint + Qgizzard     # Wang 2021 13.57
Qrest_a     <- 100 - (Qbrain+Qliver+Qkidney+Qmuscle+Qadipose+Qint+Qheart+Qrep+Qspleen+Qgizzard+Qbone) * 100 / Qart  # 7%, critical equation to ensure flow balance
Qrest       <- Qart - (Qliver + Qkidney + Qmuscle + Qadipose + Qgut + Qrep)

CV_Qadipose   <- 30                # %, Table 3 L.S. Lautz et al. 2020 (Environmental International)
CV_Qbrain     <- 20                # %, Table 3 L.S. Lautz et al. 2020 (Environmental International)
CV_Qbone      <- 30                # %, Table 3 L.S. Lautz et al. 2020 (Environmental International)
CV_Qint       <- 39                # %, Table 3 L.S. Lautz et al. 2020 (Environmental International)
CV_Qheart     <- 33                # %, Table 3 L.S. Lautz et al. 2020 (Environmental International)
CV_Qkidney    <- 29                # %, Table 3 L.S. Lautz et al. 2020 (Environmental International)
CV_Qliver     <- 43                # %, Table 3 L.S. Lautz et al. 2020 (Environmental International)
CV_Qmuscle    <- 30                # %, Table 3 L.S. Lautz et al. 2020 (Environmental International)
CV_Qrep       <- 26                # %, Table 3 L.S. Lautz et al. 2020 (Environmental International)
#CV_Qlung      <- 30


###################################################################
###                       liver related                         ###  
###################################################################
# hepatic and intestinal clearance for parent; hepatic clearance for daughter (rat only)
liver.density           <-  1.05

#General hepatic parameter
HPGL_human              <-  117.5                           # million cells/g liver; Number of hepacytes per gram liver; Acibenzolar manuscript appendix.
MPPGL_human             <-  45                              # mg/g;microsomal protein content per gram liver; Acibenzolar manuscript appendix
MPPGGI_human            <-  3                               # mg/g tissue; microsomal protein content in the GI. # not used: table 2 of TK0648566-01 report

HPGL_rat                <-  108                             # million cells/g liver; Number of hepacytes per gram liver; Acibenzolar manuscript appendix.
MPPGL_rat               <-  45                              # mg/g tissue;microsomal protein content per gram liver; Acibenzolar manuscript appendix
MPPGGI_rat              <-  3                               # mg/g tissue; microsomal protein content in the GI. table 2 of TK0648566-01 report

MPPGL_mouse             <-  45                              # mg/g tissue;Miyoung Yoon 2019;Sakai C 2014
MPPGGI_mouse            <-  3                               # mg/g tissue;

MPPGL_chicken           <-  9.31                            # mg/g liver page 5 L.S. Lautz et al. 2020 (Environmental International)


###################################################################
###                           blood                             ###  
###################################################################               
pHIW                    <- 7.22                             # pH of intracellular of blood cells
pHEW                    <- 7.4                              # pH of extracellular water 
plasma.pH               <- 7.49                             # pH of plasma for chicken Hugo Chiodi, and James W. Terman; Arterial blood gases of the domestic hen


###################################################################
###                            egg                              ###  
###################################################################
##################### Yolk
## pH
pH_y          <- 6                        # pH is egg yolk
## fraction of total volume (sum equal to 1)
Fcell_y       <- 1                        # fraction of total volume (cells)
Fint_y        <- 0                        # fraction of total volume (interstitium) yolk
## fraction of cell volume (sum not necessary to be equal to 1)
FW_y          <- 0.52                     # fraction of water within cellular volume 
FL_y          <- 0.31                     # fraction of lipid within cellular volume
Fpr_y         <- 0.17                     # fraction of protein within cellular volume
## fraction of total lipid (sum close to 1)
Fn_L_y        <- 0.66 + 0.06              # fraction of neutral lipid for total lipid
Fn_PL_y       <- 0.28 * (0.75 + 0.15)     # fraction of neutral phospholipid 
Fa_PL_y       <- 0.28 * 0.1               # fraction of acidic phospholipid 

##################### White
## pH
pH_w          <- 7.6
## fraction of total volume (sum equal to 1)
#Fcell_w       <- 1                       # fraction of total volume (cells)
#Fint_w        <- 0                       # fraction of total volume (interstitium) yolk
## fraction of cell volume (sum not necessary to be equal to 1)
FW_w          <- 0.8897                   # fraction of water within egg white
FL_w          <- 0.03/100                 # fraction of lipid within egg white
Fpr_w         <- 0.11                     # fraction of protein within egg white
## fraction of total lipid (sum close to 1)
Fn_L_w        <- 0.85                     # fraction of neutral lipid for total lipid
Fn_PL_w       <- 0.15                     # fraction of neutral phospholipid 
Fa_PL_w       <- 0                        # fraction of acidic phospholipid 

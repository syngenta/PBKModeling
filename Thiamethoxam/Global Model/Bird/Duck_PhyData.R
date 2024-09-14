#######################################################################
###               tissue weight for mallard  (% bw)                 ### 
###               tissue volume for mallard   (L/kg)                ###  
#######################################################################
R_rat                   <- 0.11             # 0.2cm # radius of rat jejunum; https://github.com/wfsrqivive/rat_PBK/blob/main/model/generate_input_data.r
R_human                 <- 1.25  

# Table 7, Scanes el al 2022 Quantitative morphometric, no CV information provided (assumed 10%)
Hematocrit_fraction      <- 45.4     # Table 1, Scanes et al 2022 Quantitative morphometric
SD_Hematocrit_fraction   <- 1.82     # Table 1, Scanes et al 2022

# L/kg (organ weight as %bw)
# Table 7, Scanes el al 2022 Quantitative morphometric (L/kg)
Vblood      <- 6.3   /100   # Table 7, Scanes el al 2022 Quantitative morphometric %
#Vart       <- Vplasma / (1 - Hematocrit_fraction/100) /2             # L/kg
#Vven       <- Vplasma / (1 - Hematocrit_fraction/100) /2

Vart        <- Vblood * 0.16             # L/kg
Vven        <- Vblood * 0.595
Vplasma     <- Vblood * (1-Hematocrit_fraction/100)


Vheart      <- 0.45  /100   # Table 7, Scanes el al 2022 Quantitative morphometric
Vspleen     <- 0.16  /100   # Table 7, Scanes el al 2022 Quantitative morphometric

#Fsmallint   <- 3.39
Vliver      <- 3.7  /100    # Table 1, Baier et al 2022
Vkidney     <- 1.2   /100   # Table 1, Baier et al 2022
Vovary      <- 2.38  /100   # Table 7, Scanes el al 2022 Quantitative morphometric
Voviduct    <- 2.84  /100   # Table 7, Scanes el al 2022 Quantitative morphometric
Vgizzard    <- 3.21  /100   # Table 7, Scanes el al 2022 Quantitative morphometric

#Vskin       <- 12    /100   # Fereidoun et al. (2007) Mean Percentage of Skin and Visible Fat in 10 Chicken Carcass Weight 8-20%;  Table 1. Wang et al 2020 chicken and turkey
Vbrain      <- 0.2   /100   # # Table 3 L.S. Lautz et al. 2020 (Environmental International)
Vlung       <- 0.008 /0.74  # # Table 1, Baier et al 2022
Vmuscle     <- 43.1  /100   # # Table 1, Baier et al 2022
Vadipose    <- 3.92  /100   # # Table 1, Baier et al 2022
Vint        <- 5.8   /100   # # intestine, Table 3 L.S. Lautz et al. 2020 (Environmental International)
Vbone       <- 20.3  /100   # # Table 3 L.S. Lautz et al. 2020 (Environmental International) Carcass

Vrep        <- Vovary + Voviduct 
Vgut        <- Vint + Vgizzard # no crop for duck

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
Qart        <- 0.43 * 60 * 24      # L/d/kg BW 20.4 * 24 Table 1 Baier 2022 (642); Scanes et al 2022 (599)


# Cv can be 30% for all the organs according to Table 3 L.S. Lautz et al. 2020 (Environmental International); not used
Qovary      <- 1.9                 # 2.59%, derived from  Table 8, Scanes el al 2022 Quantitative morphometric 
Qoviduct    <- 5.88                # 7.94%, derived from  Table 8, Scanes el al 2022 Quantitative morphometric 

Qbrain      <- 0.4  / 100 * Qart   # % of cardiac output * cardiac output, Table 3 L.S. Lautz et al. 2020 (Environmental International)
Qliver      <- 0.58 * 24 * 60 * Vliver      # Table 1 Baier 2022; Scanes et al 2022
Qkidney     <- 1.08 * 24 * 60 * Vkidney     # Table 1 Baier 2022; Scanes et al 2022
Qmuscle     <- 0.55 * 24 * 60 * Vmuscle     # Table 1 Baier 2022; Scanes et al 2022
Qadipose    <- 0.158 * 24 * 60 * Vadipose   # Table 1 Baier 2022; Scanes et al 2022
Qint        <- 17.7 / 100 * Qart   # intestine, % of cardiac output * cardiac output, Table 3 L.S. Lautz et al. 2020 (Environmental International)
Qheart      <- 5.5  / 100 * Qart   # % of cardiac output * cardiac output, Table 3 L.S. Lautz et al. 2020 (Environmental International)
Qrep        <- 14.3 / 100 * Qart   # % of cardiac output * cardiac output, Table 3 L.S. Lautz et al. 2020 (Environmental International)
Qspleen     <- 1.6  / 100 * Qart   # Table 8, Scanes el al 2022 Quantitative morphometric 
Qgizzard    <- 1.8  / 100 * Qart   # Table 8, Scanes el al 2022 Quantitative morphometric 
#Qlung       <- 3
Qbone       <- 12.4 / 100 * Qart   # Skeletal Blood Flow in Bone Repair and Maintenance

Qgut        <- 13.57 / 100 * Qart  # Wang 2021 13.57
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



pHIW        <- 7.22             # pH of intracellular of blood cells
pHEW        <- 7.4              # pH of extracellular water 
plasma.pH   <- 7.49             # pH of plasma for chicken Hugo Chiodi, and James W. Terman; Arterial blood gases of the domestic hen


#########
# Yolk
## pH
pH_y          <- 6                        # pH is egg yolk
## fraction of total volume (sum equal to 1)
Fcell_y       <- 1                        # fraction of total volume (cells)
Fint_y        <- 0                        # fraction of total volume (interstitium) yolk
## fraction of cell volume (sum not necessary to be equal to 1)
FW_y          <- 0.45                     # fraction of water within cellular volume 
FL_y          <- 0.36                     # fraction of lipid within cellular volume
Fpr_y         <- 0.19                     # fraction of protein within cellular volume; CHI and TSENG 1998
## fraction of total lipid (sum close to 1)
Fn_L_y        <- 0.73                     # fraction of neutral lipid for total lipid
Fn_PL_y       <- 0.27 * 0.93              # fraction of neutral phospholipid 
Fa_PL_y       <- 0.27 * (1 - 0.93)        # fraction of acidic phospholipid 

# White
## pH
pH_w          <- 7.8
## fraction of total volume (sum equal to 1)
#Fcell_w       <- 1                       # fraction of total volume (cells)
#Fint_w        <- 0                       # fraction of total volume (interstitium) yolk
## fraction of cell volume (sum not necessary to be equal to 1)
FW_w          <- 0.875                    # fraction of water within duck egg white (Banerjee et al. 2011, Sun et al. 2019)
FL_w          <- 0.13/100                 # fraction of lipid within duck egg white
Fpr_w         <- 0.11                     # fraction of protein within egg white
## fraction of total lipid (sum close to 1)
Fn_L_w        <- 0.85                     # fraction of neutral lipid for total lipid
Fn_PL_w       <- 0.15                     # fraction of neutral phospholipid 
Fa_PL_w       <- 0                        # fraction of acidic phospholipid 

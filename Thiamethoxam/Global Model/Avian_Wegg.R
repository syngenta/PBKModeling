

# Reference
# PKNCA: https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk: https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk
# General hepatic parameter
HPGL_human              <-  117.5                           # million cells/g liver; Number of hepacytes per gram liver; Acibenzolar manuscript appendix.
MPPGL_human             <-  45                              # mg/g;microsomal protein content per gram liver; Acibenzolar manuscript appendix
MPPGGI_human            <-  3                               # mg/g tissue; microsomal protein content in the GI. # not used: table 2 of TK0648566-01 report

HPGL_rat                <-  108                             # million cells/g liver; Number of hepacytes per gram liver; Acibenzolar manuscript appendix.
MPPGL_rat               <-  45                              # mg/g tissue;microsomal protein content per gram liver; Acibenzolar manuscript appendix
MPPGGI_rat              <-  3                               # mg/g tissue; microsomal protein content in the GI. table 2 of TK0648566-01 report

MPPGL_mouse             <-  45                              # mg/g tissue;Miyoung Yoon 2019;Sakai C 2014
MPPGGI_mouse            <-  3                               # mg/g tissue;

kt_human                <-  0.02                             # h-1
kt_mouse                <-  0.11                             # h-1
kt_rat                  <-  0.054                            # h-1

BW_mouse                <- 0.02
BW_rat                  <- 0.25
BW_human                <- 70

#=========================================================================================
#                                 PBTK model equations                                   #
#=========================================================================================

pbtk5cpt_repeated <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    ## blood concentrations of the chemical in vein compartment
    Cgut_parent         =  Agut_parent / Vgut                    # tissue concentration uM (umol/kg / L/kg)
    Cliver_parent       =  Aliver_parent / Vliver                # tissue concentration uM
    Ckidney_parent      =  Akidney_parent / Vkidney              # tissue concentration uM
    Cadipose_parent     =  Aadipose_parent / Vadipose       # adipose concentration in the (extravascular) tissue uM
    Cmuscle_parent      =  Amuscle_parent / Vmuscle              # muscle concentration uM
    Crest_parent        =  Arest_parent / Vrest                  # tissue concentration uM
    Crep_parent         =  Arep_parent / Vrep                    # reproduction (ovary & oviduct) concentration uM
    Cven_parent         =  Aven_parent / Vven                    # blood concentration in vein uM
    Clung_parent        =  Alung_parent / Vlung                  # tissue concentration uM
    Cart_parent         =  Aart_parent / Vart                    # blood concentration in artery uM
    
    Cgut_daughter         =  Agut_daughter / Vgut                    # tissue concentration uM (ug/kg / L/kg)
    Cliver_daughter       =  Aliver_daughter / Vliver                # tissue concentration uM
    Ckidney_daughter      =  Akidney_daughter / Vkidney              # tissue concentration uM
    Cadipose_daughter     =  Aadipose_daughter / Vadipose       # adipose concentration in the (extravascular) tissue uM
    Cmuscle_daughter      =  Amuscle_daughter / Vmuscle              # muscle concentration uM
    Crest_daughter        =  Arest_daughter / Vrest                  # tissue concentration uM
    Crep_daughter         =  Arep_daughter / Vrep                    # reproduction (ovary & oviduct) concentration uM
    Cven_daughter         =  Aven_daughter / Vven                    # blood concentration in vein uM
    Clung_daughter        =  Alung_daughter / Vlung                  # tissue concentration uM
    Cart_daughter         =  Aart_daughter / Vart                    # blood concentration in artery uM
    
    ## Time for egg lay in nine egg compartments, assuming continues egg lays (everyday egg lay)
    ## assumed that egg was laid around 9 am daily; Choi JH et al 2004
    tlay1       <- ifelse(((t-0.38)/9)%%1 == 0, t, floor((t-0.38)/9)*9 + 9.38)          # egg lay start at day 9
    tlay2       <- ifelse(((t-1.38)/9)%%1 == 0, t, floor((t-1.38)/9)*9 + 10.38)         # egg lay start at day 10
    tlay3       <- ifelse(((t-2.38)/9)%%1 == 0, t, floor((t-2.38)/9)*9 + 11.38)         # egg lay start at day 11
    tlay4       <- ifelse(((t-3.38)/9)%%1 == 0, t, floor((t-3.38)/9)*9 + 12.38)         # egg lay start at day 12
    tlay5       <- ifelse(((t-4.38)/9)%%1 == 0, t, floor((t-4.38)/9)*9 + 13.38)         # egg lay start at day 13
    tlay6       <- ifelse(((t-5.38)/9)%%1 == 0, t, floor((t-5.38)/9)*9 + 14.38)         # egg lay start at day 14
    tlay7       <- ifelse(((t-6.38)/9)%%1 == 0, t, floor((t-6.38)/9)*9 + 15.38)         # egg lay start at day 15
    tlay8       <- ifelse(((t-7.38)/9)%%1 == 0, t, floor((t-7.38)/9)*9 + 16.38)         # egg lay start at day 16
    tlay9       <- ifelse(((t-8.38)/9)%%1 == 0, t, floor((t-8.38)/9)*9 + 17.38)         # egg lay start at day 17
    
    #import      <- input(t) 
    
    
    #############          equations for different compartment (parent)              ##########
    #####             equations for different compartment                #####
    # Gutlumen
    #dAgutlumen   =   import - ka * Agutlumen - kt * Agutlumen                           # umol/d kg bw
    kt                 <- 0
    dAgutlumen_parent   =   - ka * Agutlumen_parent - kt * Agutlumen_parent                     # umol/d kg bw, kt = 0
    
    # Gut 
    Cgutblood_parent    =   Rblood2plasma_parent / (Kgut2pu_parent * fub_parent) * Cgut_parent
    dAgut_parent        =   ka * fa * Agutlumen_parent + Qgut * (Cart_parent - Cgutblood_parent)
    
    # Liver
    Cliverblood_parent  =   Rblood2plasma_parent / (Kliver2pu_parent * fub_parent) * Cliver_parent
    Relim_parent        =   Clint_parent / Kliver2pu_parent * Cliver_parent                               # L/d/kg bw * uM -> umol/d/kg bw
    dAliver_parent      =   Qliver * Cart_parent + Qgut * Cgutblood_parent  - (Qliver + Qgut) * Cliverblood_parent - Relim_parent
    
    # Kidney
    Ckidneyblood_parent =   Rblood2plasma_parent / (Kkidney2pu_parent * fub_parent) * Ckidney_parent
    Rgfr_parent         =   Qgfr_parent / Kkidney2pu_parent * Ckidney_parent                     # L/d/kg BW  * uM  -->umol/d/kg BW
    dAkidney_parent     =   Qkidney * Cart_parent - Qkidney * Ckidneyblood_parent - Rgfr_parent
    
    # Adipose
    ##  vascular
    Cadiposeblood_parent  =   Rblood2plasma_parent / (Kadipose2pu_parent * fub_parent) * Cadipose_parent
    dAadipose_parent      =   Qadipose * (Cart_parent - Cadiposeblood_parent) 

    
    # Muscle
    Cmuscleblood_parent   =   Rblood2plasma_parent / (Kmuscle2pu_parent * fub_parent) * Cmuscle_parent
    dAmuscle_parent       =   Qmuscle * (Cart_parent - Cmuscleblood_parent)
    
    # Rest of body
    Crestblood_parent    =   Rblood2plasma_parent / (Krest2pu_parent * fub_parent) * Crest_parent
    dArest_parent        =   Qrest * (Cart_parent - Crestblood_parent)
    
    #############  reproduction (Overy & Oviduct)  ############
    ## yolk
    Crepplasma_parent    =   1 / (Krep2pu_parent * fub_parent) * Crep_parent
 
    ### Egg 1 
    Wy1                  =   egg_yolk * exp(-(t-(tlay1-tsig-tlag))) / ((1+exp(-(t-(tlay1-tsig-tlag))))^2) /1000        # // growth rate of yolk kg/d
    Ryolk1_parent        =   ky_parent * Crepplasma_parent * Wy1 / BW                                                          # uM * kg/d kg BW  --> umol/d kg bw
    ### Egg 2
    Wy2                  =   egg_yolk * exp(-(t-(tlay2-tsig-tlag))) / ((1+exp(-(t-(tlay2-tsig-tlag))))^2) /1000        # // growth rate kg/d
    Ryolk2_parent        =   ky_parent * Crepplasma_parent * Wy2 / BW  
    ### Egg 3
    Wy3                  =   egg_yolk * exp(-(t-(tlay3-tsig-tlag))) / ((1+exp(-(t-(tlay3-tsig-tlag))))^2) /1000        # // growth rate kg/d
    Ryolk3_parent        =   ky_parent * Crepplasma_parent * Wy3 / BW  
    ### Egg 4
    Wy4                  =   egg_yolk * exp(-(t-(tlay4-tsig-tlag))) / ((1+exp(-(t-(tlay4-tsig-tlag))))^2) /1000        # // growth rate kg/d
    Ryolk4_parent        =   ky_parent * Crepplasma_parent * Wy4 / BW  
    ### Egg 5
    Wy5                  =   egg_yolk * exp(-(t-(tlay5-tsig-tlag))) / ((1+exp(-(t-(tlay5-tsig-tlag))))^2) /1000        # // growth rate kg/d
    Ryolk5_parent        =   ky_parent * Crepplasma_parent * Wy5 / BW  
    ### Egg 6
    Wy6                  =   egg_yolk * exp(-(t-(tlay6-tsig-tlag))) / ((1+exp(-(t-(tlay6-tsig-tlag))))^2) /1000        # // growth rate kg/d
    Ryolk6_parent        =   ky_parent * Crepplasma_parent * Wy6 / BW  
    ### Egg 7
    Wy7                  =   egg_yolk * exp(-(t-(tlay7-tsig-tlag))) / ((1+exp(-(t-(tlay7-tsig-tlag))))^2) /1000        # // growth rate kg/d
    Ryolk7_parent        =   ky_parent * Crepplasma_parent * Wy7 / BW  
    ### Egg 8
    Wy8                  =   egg_yolk * exp(-(t-(tlay8-tsig-tlag))) / ((1+exp(-(t-(tlay8-tsig-tlag))))^2) /1000        # // growth rate kg/d
    Ryolk8_parent        =   ky_parent * Crepplasma_parent * Wy8 / BW  
    ### Egg 9
    Wy9                  =   egg_yolk * exp(-(t-(tlay9-tsig-tlag))) / ((1+exp(-(t-(tlay9-tsig-tlag))))^2) /1000        # // growth rate kg/d
    Ryolk9_parent        =   ky_parent * Crepplasma_parent * Wy9 / BW                                                          # uM * kg (L)/d kg bw -> umol/d/kg
    
    ### Egg white (daily)
    dWw             =   egg_white / (10/24) * ((t-floor(t)-0.38) >= 0) * (((t-floor(t)-0.38)) <= talbumen) /1000        # // weight of the egg white (kg) at time t (d) [tlay-tlag, tlay-tlag+talbumen]
    #Rw             =   kw * Crepplasma * 3.4 * (t >= tlay-tlag) * (t <= tlay - tlag + talbumen) * t  
    
    # assuming egg laid at time 0h at each day
    Rw_parent            =   kw_parent * Crepplasma_parent * egg_white/(10/24) * ((t-floor(t)-0.38) >= 0) * (((t-floor(t)-0.38)) <= talbumen) /1000 / BW          # (t/1-floor(t/1)) * 1  
  
    Crepblood_parent     =   Rblood2plasma_parent / (Krep2pu_parent * fub_parent) * Crep_parent
    dArep_parent         =   Qrep * (Cart_parent - Crepblood_parent) - Ryolk1_parent - Ryolk2_parent - Ryolk3_parent - Ryolk4_parent -
                                Ryolk5_parent - Ryolk6_parent - Ryolk7_parent - Ryolk8_parent - Ryolk9_parent - Rw_parent
    
    dAyolk1_parent       =   Ryolk1_parent
    dAyolk2_parent       =   Ryolk2_parent
    dAyolk3_parent       =   Ryolk3_parent
    dAyolk4_parent       =   Ryolk4_parent
    dAyolk5_parent       =   Ryolk5_parent
    dAyolk6_parent       =   Ryolk6_parent
    dAyolk7_parent       =   Ryolk7_parent
    dAyolk8_parent       =   Ryolk8_parent
    dAyolk9_parent       =   Ryolk9_parent
    dAwhite_parent       =   Rw_parent                     # ug/kg
    dWyolk1              =   Wy1
    dWyolk2              =   Wy2
    dWyolk3              =   Wy3
    dWyolk4              =   Wy4
    dWyolk5              =   Wy5
    dWyolk6              =   Wy6
    dWyolk7              =   Wy7
    dWyolk8              =   Wy8
    dWyolk9              =   Wy9
    
    # Venous blood
    dAven_parent         =   (Qliver + Qgut) * Cliverblood_parent + Qkidney * Ckidneyblood_parent + 
                                  Qadipose * Cadiposeblood_parent + Qmuscle * Cmuscleblood_parent + 
                                  Qrep * Crepblood_parent + Qrest * Crestblood_parent - Qart * Cven_parent
    
    # Lung
    Clungblood_parent    =   Rblood2plasma_parent / (Klung2pu_parent * fub_parent) * Clung_parent
    dAlung_parent        =   Qart * (Cven_parent - Clungblood_parent)
    
    # Artery
    dAart_parent         =   Qart * (Clungblood_parent - Cart_parent)
    
    dAurine_parent       =  Rgfr_parent  
    
    dAUC_Cplasma_parent  =  Cven_parent / Rblood2plasma_parent
    dAUC_Cblood_parent   =  Cven_parent
    
    
    
    #############          equations for different compartment (daughter)              ##########
    #####             equations for different compartment                #####
    # Gut 
    Cgutblood_daughter    =   Rblood2plasma_daughter / (Kgut2pu_daughter * fub_daughter) * Cgut_daughter
    dAgut_daughter        =   Qgut * (Cart_daughter - Cgutblood_daughter)
    
    # Liver
    Cliverblood_daughter  =   Rblood2plasma_daughter / (Kliver2pu_daughter * fub_daughter) * Cliver_daughter
    Relim_daughter        =   Clint_daughter / Kliver2pu_daughter * Cliver_daughter                               # L/d/kg bw * uM -> ug/d/kg bw
    dAliver_daughter      =   Qliver * Cart_daughter + Qgut * Cgutblood_daughter + fraction_daughter * Relim_parent - (Qliver + Qgut) * Cliverblood_daughter - Relim_daughter
    
    # Kidney
    Ckidneyblood_daughter =   Rblood2plasma_daughter / (Kkidney2pu_daughter * fub_daughter) * Ckidney_daughter
    Rgfr_daughter         =   Qgfr_daughter / Kkidney2pu_daughter * Ckidney_daughter
    dAkidney_daughter     =   Qkidney * Cart_daughter - Qkidney * Ckidneyblood_daughter - Rgfr_daughter
    
    # Adipose
    ##  vascular
    Cadiposeblood_daughter  =   Rblood2plasma_daughter / (Kadipose2pu_daughter * fub_daughter) * Cadipose_daughter
    dAadipose_daughter      =   Qadipose * (Cart_daughter - Cadiposeblood_daughter) 
    
    # Muscle
    Cmuscleblood_daughter   =   Rblood2plasma_daughter / (Kmuscle2pu_daughter * fub_daughter) * Cmuscle_daughter
    dAmuscle_daughter       =   Qmuscle * (Cart_daughter - Cmuscleblood_daughter)
    
    # Rest of body
    Crestblood_daughter    =   Rblood2plasma_daughter / (Krest2pu_daughter * fub_daughter) * Crest_daughter
    dArest_daughter        =   Qrest * (Cart_daughter - Crestblood_daughter)
    
    #############  reproduction (Overy & Oviduct)  ############
    ## yolk
    Crepplasma_daughter    =   1 / (Krep2pu_daughter * fub_daughter) * Crep_daughter
    
    ### Egg 1 
    Ryolk1_daughter        =   ky_daughter * Crepplasma_daughter * Wy1 / BW                                                          # umol/L * kg/d kg BW  --> umol/d kg bw
    ### Egg 2
    Ryolk2_daughter        =   ky_daughter * Crepplasma_daughter * Wy2 / BW  
    ### Egg 3
    Ryolk3_daughter        =   ky_daughter * Crepplasma_daughter * Wy3 / BW  
    ### Egg 4
    Ryolk4_daughter        =   ky_daughter * Crepplasma_daughter * Wy4 / BW  
    ### Egg 5
    Ryolk5_daughter        =   ky_daughter * Crepplasma_daughter * Wy5 / BW  
    ### Egg 6
    Ryolk6_daughter        =   ky_daughter * Crepplasma_daughter * Wy6 / BW  
    ### Egg 7
    Ryolk7_daughter        =   ky_daughter * Crepplasma_daughter * Wy7 / BW  
    ### Egg 8
    Ryolk8_daughter        =   ky_daughter * Crepplasma_daughter * Wy8 / BW  
    ### Egg 9
    Ryolk9_daughter        =   ky_daughter * Crepplasma_daughter * Wy9 / BW                                                          # uM * kg (L)/d kg bw -> umol/d/kg
    
    # assuming egg laid at time 0h at each day
    Rw_daughter            =   kw_daughter * Crepplasma_daughter * egg_white/10 * 24 * ((t-floor(t)-0.38) >= 0) * (((t-floor(t)-0.38)) <= talbumen) /1000 / BW          # (t/1-floor(t/1)) * 1  
    
    Crepblood_daughter     =   Rblood2plasma_daughter / (Krep2pu_daughter * fub_daughter) * Crep_daughter
    dArep_daughter         =   Qrep * (Cart_daughter - Crepblood_daughter) - Ryolk1_daughter - Ryolk2_daughter - Ryolk3_daughter - Ryolk4_daughter -
                                Ryolk5_daughter - Ryolk6_daughter - Ryolk7_daughter - Ryolk8_daughter - Ryolk9_daughter - Rw_daughter
    
    dAyolk1_daughter       =   Ryolk1_daughter
    dAyolk2_daughter       =   Ryolk2_daughter
    dAyolk3_daughter       =   Ryolk3_daughter
    dAyolk4_daughter       =   Ryolk4_daughter
    dAyolk5_daughter       =   Ryolk5_daughter
    dAyolk6_daughter       =   Ryolk6_daughter
    dAyolk7_daughter       =   Ryolk7_daughter
    dAyolk8_daughter       =   Ryolk8_daughter
    dAyolk9_daughter       =   Ryolk9_daughter
    dAwhite_daughter       =   Rw_daughter                     # ug/kg
 
    
    # Venous blood
    dAven_daughter         =   (Qliver + Qgut) * Cliverblood_daughter + Qkidney * Ckidneyblood_daughter + 
                                Qadipose * Cadiposeblood_daughter + Qmuscle * Cmuscleblood_daughter + 
                                Qrep * Crepblood_daughter + Qrest * Crestblood_daughter - Qart * Cven_daughter
    
    # Lung
    Clungblood_daughter    =   Rblood2plasma_daughter / (Klung2pu_daughter * fub_daughter) * Clung_daughter
    dAlung_daughter        =   Qart * (Cven_daughter - Clungblood_daughter)
    
    # Artery
    dAart_daughter         =   Qart * (Clungblood_daughter - Cart_daughter)
    
    dAurine_daughter       =  Rgfr_daughter 
    
    dAUC_Cplasma_daughter  =  Cven_daughter / Rblood2plasma_daughter
    dAUC_Cblood_daughter   =  Cven_daughter
    
    ######################         Mass balance check          #######################
    ## Mass balance of parent
    dAgut_parent_in       = ka * fa * Agutlumen_parent 
    dAgut_parent_out      = kt * Agutlumen_parent 
    dAliver_parent_out    = Clint_parent / Kliver2pu_parent * Cliver_parent
    dAven_parent_out      = 0 * Vven *  Cven_parent + Qgfr_parent / Kkidney2pu_parent * Ckidney_parent
    dAint_parent_out      = 0  * Cgut_parent / (Kgut2pu_parent * fub_parent)
    dAart_parent_out      = 0  * Vart * Cart_parent
    
    Mass_parent_in           = Agut_parent_in                                                    # Total Amount Acibenzolar from IV Dose and Absorbed in Gut; umol/kg bw
    Mass_parent_stored       = Agut_parent + Aliver_parent + Akidney_parent + Aadipose_parent + 
                                 Amuscle_parent + Arest_parent + Aven_parent + Alung_parent + Aart_parent + Arep_parent           # Total Amount of Parent Remaining in the Body (??mol)
    Mass_parent_out          = Aliver_parent_out + Aven_parent_out + Aint_parent_out + Aart_parent_out +
                                 Ayolk1_parent + Ayolk2_parent + Ayolk3_parent +
                                 Ayolk4_parent + Ayolk5_parent + Ayolk6_parent + Ayolk7_parent +
                                 Ayolk8_parent + Ayolk9_parent + Awhite_parent           # Total Amount of Parent Excreted from the Body (??mol)
    Mass_parent_bal          = Mass_parent_in - Mass_parent_stored - Mass_parent_out
    
    ## Mass balance of Parent-acid
    dAurine_daughter_out  = Qgfr_daughter / Kkidney2pu_daughter * Ckidney_daughter
    dAliver_daughter_out  = Clint_daughter / Kliver2pu_daughter * Cliver_daughter 
    
    Mass_daughter_in      = fraction_daughter * Aliver_parent_out 
    Mass_daughter_stored  = Agut_daughter + Aliver_daughter + Akidney_daughter + Aadipose_daughter + 
                             Amuscle_daughter + Arest_daughter +  Aven_daughter + Alung_daughter + Aart_daughter + Arep_daughter   
    Mass_daughter_out     = Aliver_daughter_out + Aurine_daughter_out + Ayolk1_daughter + Ayolk2_daughter + Ayolk3_daughter +
                             Ayolk4_daughter + Ayolk5_daughter +  Ayolk6_daughter + Ayolk7_daughter +
                             Ayolk8_daughter + Ayolk9_daughter + Awhite_daughter 
    
    Mass_daughter_bal     = Mass_daughter_in  - Mass_daughter_stored - Mass_daughter_out
    
    ## total mass balance
    Mass_in                  = Agut_parent_in 
    Mass_stored_out          = Mass_parent_stored + Mass_daughter_stored + Aven_parent_out + Mass_daughter_out    # Total Amount of Parent and Parent-Acid in the Body and Parent-Acid Excreted (??mol)
    Mass_bal                 = Mass_in - Mass_stored_out
    
    ##### output parameters
    # Return the rate of change in the same ordering as the specification of the state variables
    
    list(c(dAgutlumen_parent, 
           dAgut_parent, 
           dAliver_parent,
           dAkidney_parent,
           dAadipose_parent,
           dAmuscle_parent,
           dArest_parent,
           dWw,
           dArep_parent,
           dAyolk1_parent,
           dAyolk2_parent,
           dAyolk3_parent,
           dAyolk4_parent,
           dAyolk5_parent,
           dAyolk6_parent,
           dAyolk7_parent,
           dAyolk8_parent,
           dAyolk9_parent,
           dAwhite_parent, 
           dWyolk1,
           dWyolk2,
           dWyolk3,
           dWyolk4,
           dWyolk5,
           dWyolk6,
           dWyolk7,
           dWyolk8,
           dWyolk9,
           dAven_parent,
           dAlung_parent,
           dAart_parent,
           dAurine_parent,
           dAUC_Cplasma_parent,
           dAUC_Cblood_parent,
           
           dAgut_daughter, 
           dAliver_daughter,
           dAkidney_daughter,
           dAadipose_daughter,
           dAmuscle_daughter,
           dArest_daughter,
           dArep_daughter,
           dAyolk1_daughter,
           dAyolk2_daughter,
           dAyolk3_daughter,
           dAyolk4_daughter,
           dAyolk5_daughter,
           dAyolk6_daughter,
           dAyolk7_daughter,
           dAyolk8_daughter,
           dAyolk9_daughter,
           dAwhite_daughter, 
           dAven_daughter,
           dAlung_daughter,
           dAart_daughter,
           dAurine_daughter,
           dAUC_Cplasma_daughter,
           dAUC_Cblood_daughter,
           
           dAgut_parent_in, 
           dAliver_parent_out,
           dAven_parent_out,
           dAint_parent_out,
           dAart_parent_out,
           
           dAurine_daughter_out,
           dAliver_daughter_out),
         
         "C_blood_parent"  = Cven_parent,                              # uM
         "C_plasma_parent" = Cven_parent / Rblood2plasma_parent,       # uM
         "Cyolk1_parent"   = Ayolk1_parent / (Wyolk1 / BW),            # umol/kg
         "Cyolk2_parent"   = Ayolk2_parent / (Wyolk2 / BW),            # umol/kg
         "Cyolk3_parent"   = Ayolk3_parent / (Wyolk3 / BW),            # umol/kg
         "Cyolk4_parent"   = Ayolk4_parent / (Wyolk4 / BW),            # umol/kg
         "Cyolk5_parent"   = Ayolk5_parent / (Wyolk5 / BW),            # umol/kg
         "Cyolk6_parent"   = Ayolk6_parent / (Wyolk6 / BW),            # umol/kg
         "Cyolk7_parent"   = Ayolk7_parent / (Wyolk7 / BW),            # umol/kg
         "Cyolk8_parent"   = Ayolk8_parent / (Wyolk8 / BW),            # umol/kg
         "Cyolk9_parent"   = Ayolk9_parent / (Wyolk9 / BW),            # umol/kg
         "C_white_parent"  = Awhite_parent / (Ww / BW),                # umol/kg
         
         "Cegg1_parent"   = (Ayolk1_parent + Awhite_parent)/ ((Wyolk1+Ww) / BW),            # umol/kg
         "Cegg2_parent"   = (Ayolk2_parent + Awhite_parent)/ ((Wyolk2+Ww) / BW),            # umol/kg
         "Cegg3_parent"   = (Ayolk3_parent + Awhite_parent)/ ((Wyolk3+Ww) / BW),            # umol/kg
         "Cegg4_parent"   = (Ayolk4_parent + Awhite_parent)/ ((Wyolk4+Ww) / BW),            # umol/kg
         "Cegg5_parent"   = (Ayolk5_parent + Awhite_parent)/ ((Wyolk5+Ww) / BW),            # umol/kg
         "Cegg6_parent"   = (Ayolk6_parent + Awhite_parent)/ ((Wyolk6+Ww) / BW),            # umol/kg
         "Cegg7_parent"   = (Ayolk7_parent + Awhite_parent)/ ((Wyolk7+Ww) / BW),            # umol/kg
         "Cegg8_parent"   = (Ayolk8_parent + Awhite_parent)/ ((Wyolk8+Ww) / BW),            # umol/kg
         "Cegg9_parent"   = (Ayolk9_parent + Awhite_parent)/ ((Wyolk9+Ww) / BW),            # umol/kg
      
         "C_adipose_parent"  = Aadipose_parent / Vadipose,         # uM
         "C_liver_parent"    = Aliver_parent / Vliver,                   # uM
         "C_kidney_parent"   = Akidney_parent / Vkidney,                 # uM
         "C_muscle_parent"   = Amuscle_parent / Vmuscle,
         
         "C_blood_daughter"  = Cven_daughter,                              # uM
         "C_plasma_daughter" = Cven_daughter / Rblood2plasma_daughter,     # uM
         "Cyolk1_daughter"   = Ayolk1_daughter / (Wyolk1 / BW),            # umol/kg
         "Cyolk2_daughter"   = Ayolk2_daughter / (Wyolk2 / BW),            # umol/kg
         "Cyolk3_daughter"   = Ayolk3_daughter / (Wyolk3 / BW),            # umol/kg
         "Cyolk4_daughter"   = Ayolk4_daughter / (Wyolk4 / BW),            # umol/kg
         "Cyolk5_daughter"   = Ayolk5_daughter / (Wyolk5 / BW),            # umol/kg
         "Cyolk6_daughter"   = Ayolk6_daughter / (Wyolk6 / BW),            # umol/kg
         "Cyolk7_daughter"   = Ayolk7_daughter / (Wyolk7 / BW),            # umol/kg
         "Cyolk8_daughter"   = Ayolk8_daughter / (Wyolk8 / BW),            # umol/kg
         "Cyolk9_daughter"   = Ayolk9_daughter / (Wyolk9 / BW),            # umol/kg
         "C_white_daughter"  = Awhite_daughter / (Ww / BW),                # umol/kg
         
         "Cegg1_daughter"   = (Ayolk1_daughter + Awhite_daughter)/ ((Wyolk1+Ww) / BW),            # umol/kg
         "Cegg2_daughter"   = (Ayolk2_daughter + Awhite_daughter)/ ((Wyolk2+Ww) / BW),            # umol/kg
         "Cegg3_daughter"   = (Ayolk3_daughter + Awhite_daughter)/ ((Wyolk3+Ww) / BW),            # umol/kg
         "Cegg4_daughter"   = (Ayolk4_daughter + Awhite_daughter)/ ((Wyolk4+Ww) / BW),            # umol/kg
         "Cegg5_daughter"   = (Ayolk5_daughter + Awhite_daughter)/ ((Wyolk5+Ww) / BW),            # umol/kg
         "Cegg6_daughter"   = (Ayolk6_daughter + Awhite_daughter)/ ((Wyolk6+Ww) / BW),            # umol/kg
         "Cegg7_daughter"   = (Ayolk7_daughter + Awhite_daughter)/ ((Wyolk7+Ww) / BW),            # umol/kg
         "Cegg8_daughter"   = (Ayolk8_daughter + Awhite_daughter)/ ((Wyolk8+Ww) / BW),            # umol/kg
         "Cegg9_daughter"   = (Ayolk9_daughter + Awhite_daughter)/ ((Wyolk9+Ww) / BW),            # umol/kg
         
         "C_adipose_daughter" = Aadipose_daughter / Vadipose,         # uM
         "C_liver_daughter"   = Aliver_daughter   / Vliver,                   # uM
         "C_kidney_daughter"  = Akidney_daughter  / Vkidney,                 # uM
         "C_muscle_daughter"  = Amuscle_daughter  / Vmuscle,
         
         'C_blood_total'     = Cven_parent + Cven_daughter,
         
         'Mass_parent_in'        = Mass_parent_in ,   
         'Mass_parent_stored'    = Mass_parent_stored ,    
         'Mass_parent_out'       = Mass_parent_out ,         
         'Mass_parent_bal'       = Mass_parent_bal,   
         
         'Mass_daughter_in'      = Mass_daughter_in,         
         'Mass_daughter_stored'  = Mass_daughter_stored,      
         'Mass_daughter_out'     = Mass_daughter_out,       
         'Mass_daughter_bal'     = Mass_daughter_bal,        
         
         'Mass_in'               = Mass_in,                 
         'Mass_stored_out'       = Mass_stored_out,         
         'Mass_bal'              = Mass_bal)                 # uM
    
  })
} 
  

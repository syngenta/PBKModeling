

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

kt_human                <-  0.02                            # h-1
kt_mouse                <-  0.11                            # h-1
kt_rat                  <-  0.054                           # h-1

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
    Cadipose_parent     =  Aadipose_parent / Vadipose            # adipose concentration in the (extravascular) tissue uM
    Cmuscle_parent      =  Amuscle_parent / Vmuscle              # muscle concentration uM
    Crest_parent        =  Arest_parent / Vrest                  # tissue concentration uM
    Crep_parent         =  Arep_parent / Vrep                    # reproduction (ovary & oviduct) concentration uM
    Cven_parent         =  Aven_parent / Vven                    # blood concentration in vein uM
    Clung_parent        =  Alung_parent / Vlung                  # tissue concentration uM
    Cart_parent         =  Aart_parent / Vart                    # blood concentration in artery uM
    
    Cgut_daughter         =  Agut_daughter / Vgut                    # tissue concentration uM (ug/kg / L/kg)
    Cliver_daughter       =  Aliver_daughter / Vliver                # tissue concentration uM
    Ckidney_daughter      =  Akidney_daughter / Vkidney              # tissue concentration uM
    Cadipose_daughter     =  Aadipose_daughter / Vadipose            # adipose concentration in the (extravascular) tissue uM
    Cmuscle_daughter      =  Amuscle_daughter / Vmuscle              # muscle concentration uM
    Crest_daughter        =  Arest_daughter / Vrest                  # tissue concentration uM
    Crep_daughter         =  Arep_daughter / Vrep                    # reproduction (ovary & oviduct) concentration uM
    Cven_daughter         =  Aven_daughter / Vven                    # blood concentration in vein uM
    Clung_daughter        =  Alung_daughter / Vlung                  # tissue concentration uM
    Cart_daughter         =  Aart_daughter / Vart                    # blood concentration in artery uM
    

    
    
    #############          equations for different compartment (parent)              ##########
    #####             equations for different compartment                #####
    # Gutlumen
    kt                 <- 0
    dAgutlumen_parent   =   - ka * Agutlumen_parent - kt * Agutlumen_parent                               # umol/d kg bw, kt = 0
    
    # Gut 
    Cgutblood_parent    =   Rblood2plasma_parent / (Kgut2pu_parent * fub_parent) * Cgut_parent
    dAgut_parent        =   ka * fa * Agutlumen_parent + Qgut * (Cart_parent - Cgutblood_parent)
    
    # Liver
    Cliverblood_parent  =   Rblood2plasma_parent / (Kliver2pu_parent * fub_parent) * Cliver_parent
    Relim_parent        =   Clint_parent / Kliver2pu_parent * Cliver_parent                               # L/d/kg bw * uM -> umol/d/kg bw
    dAliver_parent      =   Qliver * Cart_parent + Qgut * Cgutblood_parent  - (Qliver + Qgut) * Cliverblood_parent - Relim_parent
    
    # Kidney
    Ckidneyblood_parent =   Rblood2plasma_parent / (Kkidney2pu_parent * fub_parent) * Ckidney_parent
    Rgfr_parent         =   Qgfr_parent / Kkidney2pu_parent * Ckidney_parent                              # L/d/kg BW  * uM  -->umol/d/kg BW
    dAkidney_parent     =   Qkidney * Cart_parent - Qkidney * Ckidneyblood_parent - Rgfr_parent
    
    # Adipose
    Cadiposeblood_parent  =   Rblood2plasma_parent / (Kadipose2pu_parent * fub_parent) * Cadipose_parent
    dAadipose_parent      =   Qadipose * (Cart_parent - Cadiposeblood_parent) 
    
    # Muscle
    Cmuscleblood_parent   =   Rblood2plasma_parent / (Kmuscle2pu_parent * fub_parent) * Cmuscle_parent
    dAmuscle_parent       =   Qmuscle * (Cart_parent - Cmuscleblood_parent)
    
    # Rest of body
    Crestblood_parent    =   Rblood2plasma_parent / (Krest2pu_parent * fub_parent) * Crest_parent
    dArest_parent        =   Qrest * (Cart_parent - Crestblood_parent)
    
    #############  reproduction (Overy & Oviduct)  ############
    Crepplasma_parent    =   1 / (Krep2pu_parent * fub_parent) * Crep_parent
    Crepblood_parent     =   Rblood2plasma_parent / (Krep2pu_parent * fub_parent) * Crep_parent
    dArep_parent         =   Qrep * (Cart_parent - Crepblood_parent) 
    
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
    Crepblood_daughter     =   Rblood2plasma_daughter / (Krep2pu_daughter * fub_daughter) * Crep_daughter
    dArep_daughter         =   Qrep * (Cart_daughter - Crepblood_daughter) 
 
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
    
    Mass_parent_in           = Agut_parent_in                                                                                    # Total Amount from IV Dose and Absorbed in Gut; umol/kg bw
    Mass_parent_stored       = Agut_parent + Aliver_parent + Akidney_parent + Aadipose_parent + 
                                Amuscle_parent + Arest_parent + Aven_parent + Alung_parent + Aart_parent + Arep_parent           # Total Amount of Parent Remaining in the Body 
    Mass_parent_out          = Aliver_parent_out + Aven_parent_out + Aint_parent_out + Aart_parent_out                           # Total Amount of Parent Excreted from the Body 
    Mass_parent_bal          = Mass_parent_in - Mass_parent_stored - Mass_parent_out
    
    ## Mass balance of daughter
    dAurine_daughter_out  = Qgfr_daughter / Kkidney2pu_daughter * Ckidney_daughter
    dAliver_daughter_out  = Clint_daughter / Kliver2pu_daughter * Cliver_daughter 
    
    Mass_daughter_in      = fraction_daughter * Aliver_parent_out 
    Mass_daughter_stored  = Agut_daughter + Aliver_daughter + Akidney_daughter + Aadipose_daughter + 
                             Amuscle_daughter + Arest_daughter +  Aven_daughter + Alung_daughter + Aart_daughter + Arep_daughter   
    Mass_daughter_out     = Aliver_daughter_out + Aurine_daughter_out
    
    Mass_daughter_bal     = Mass_daughter_in  - Mass_daughter_stored - Mass_daughter_out
    
    ## total mass balance
    Mass_in                  = Agut_parent_in 
    Mass_stored_out          = Mass_parent_stored + Mass_daughter_stored + 
                                Aven_parent_out + Mass_daughter_out                                                              # Total Amount of Parent and daughter in the Body and Excreted
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
           dArep_parent,
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
         "C_adipose_parent" = Aadipose_parent / Vadipose,              # uM
         "C_liver_parent"  = Aliver_parent / Vliver,                   # uM
         "C_kidney_parent" = Akidney_parent / Vkidney,                 # uM
         "C_muscle_parent" = Amuscle_parent / Vmuscle,
         
         "C_blood_daughter"  = Cven_daughter,                              # uM
         "C_plasma_daughter" = Cven_daughter / Rblood2plasma_daughter,     # uM
         "C_adipose_daughter" = Aadipose_daughter / Vadipose,              # uM
         "C_liver_daughter"  = Aliver_daughter / Vliver,                   # uM
         "C_kidney_daughter" = Akidney_daughter / Vkidney,                 # uM
         "C_muscle_daughter" = Amuscle_daughter / Vmuscle,
         
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
         'Mass_bal'              = Mass_bal)                               # uM
    
        

  })
} 
  

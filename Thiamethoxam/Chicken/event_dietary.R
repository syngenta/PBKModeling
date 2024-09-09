# egg laying event
## event triggered if time %%24 = 0
## https://stat.ethz.ch/pipermail/r-sig-dynamic-models/2013q1/000160.html
## https://stat.ethz.ch/pipermail/r-sig-dynamic-models/2013q1/000158.html
# https://stat.ethz.ch/pipermail/r-sig-dynamic-models/2018/000606.html          continuous infusion
# https://cran.r-project.org/web/packages/pksensi/vignettes/pbtk1cpt.html
# https://rstudio-pubs-static.s3.amazonaws.com/557239_2e4541dc478b4ef69a2ff2664da352bd.html  (inhalation)
# https://tpetzoldt.github.io/deSolve-forcing/deSolve-forcing.html    approxfun (useful for infusion and inhalation)
# https://stackoverflow.com/questions/71009453/combining-root-and-non-root-event-functions-in-r-desolve

'root <- function(t, y, parms) {return(t%%1)}
## set state variable Ww = 0
event <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  return(y)
}'

# DOSING WAS FOR 5 CONSECUTIVE DAYS FOLLOWED BY A WITHDRAWAL PERIOD OF 7 DAYS
eventdose <- data.frame(var    = rep("Agutlumen_parent", 4),
                        time   = c(0.3, 1.3, 2.3, 3.3),
                        value  = rep(Oral_input, 4),
                        method = "add")

# egg laying event
eventW <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white parent
  y[51]      <- 0                                             # chemical mass in egg white daughter
  return(y)
}

event1 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white parent
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[10]      <- 0                                             # chemical mass in egg yolk 1 parent
  y[42]      <- 0                                             # chemical mass in egg yolk 1 daughter
  y[20]      <- 0.1/1000                                      # mass in egg yolk 1
  return(y)
}

event2 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[11]      <- 0                                             # chemical mass in egg yolk 2
  y[43]      <- 0                                             # chemical mass in egg yolk 2 daughter
  y[21]      <- 0.1/1000                                      # mass in egg yolk 2
  return(y)
}

event3 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[12]      <- 0                                             # chemical mass in egg yolk 3
  y[44]      <- 0                                             # chemical mass in egg yolk 3 daughter
  y[22]      <- 0.1/1000                                      # mass in egg yolk 3
  return(y)
}

event4 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[13]      <- 0                                             # chemical mass in egg yolk 4
  y[45]      <- 0                                             # chemical mass in egg yolk 41 daughter
  y[23]      <- 0.1/1000                                      # mass in egg yolk 4
  return(y)
}

event5 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[14]      <- 0                                             # chemical mass in egg yolk 5
  y[46]      <- 0                                             # chemical mass in egg yolk 5 daughter
  y[24]      <- 0.1/1000                                      # mass in egg yolk 5
  return(y)
}

event6 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[15]      <- 0                                             # chemical mass in egg yolk 6
  y[47]      <- 0                                             # chemical mass in egg yolk 6 daughter
  y[25]      <- 0.1/1000                                      # mass in egg yolk 6
  return(y)
}

event7 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[16]      <- 0                                             # chemical mass in egg yolk 7
  y[48]      <- 0                                             # chemical mass in egg yolk 7 daughter
  y[26]      <- 0.1/1000                                      # mass in egg yolk 7
  return(y)
}

event8 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[17]      <- 0                                             # chemical mass in egg yolk 8
  y[49]      <- 0                                             # chemical mass in egg yolk 8 daughter
  y[27]      <- 0.1/1000                                      # mass in egg yolk 8
  return(y)
}

event9 <- function(t, y, parms) {
  y[8]       <- 0.1/1000                                      # 0.1 mL mass of egg white
  y[19]      <- 0                                             # chemical mass in egg white
  y[51]      <- 0                                             # chemical mass in egg white daughter
  y[18]      <- 0                                             # chemical mass in egg yolk 9
  y[50]      <- 0                                             # chemical mass in egg yolk 9 daughter
  y[28]      <- 0.1/1000                                      # mass in egg yolk 9
  return(y)
}



# dose event
event_dose0 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0                               
  return(y)
}

event_dose1 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.002022362                              
  return(y)
}

event_dose2 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.004044725                               
  return(y)
}

event_dose3 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.002022362                                
  return(y)
}

event_dose4 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.026560358                                 
  return(y)
}

event_dose5 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.047265253                               
  return(y)
}

event_dose6 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input *  0.048256832                                
  return(y)
}

event_dose7 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.048752621                                
  return(y)
}

event_dose8 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.053875778                                
  return(y)
}

event_dose9 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.059164198                                 
  return(y)
}

event_dose10 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.063461039                                
  return(y)
}

event_dose11 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.062799987                                
  return(y)
}

event_dose12 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.059494724                                
  return(y)
}

event_dose13 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.059494724                               
  return(y)
}

event_dose14 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.064948407                                
  return(y)
}

event_dose15 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.063791565                               
  return(y)
}

event_dose16 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.073376826                               
  return(y)
}

event_dose17 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.089076823                                
  return(y)
}

event_dose18 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.087919981                                
  return(y)
}

event_dose19 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.069245248                                
  return(y)
}

event_dose20 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.009302867                                
  return(y)
}

event_dose21 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0.005123318                              
  return(y)
}

event_dose22 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0                              
  return(y)
}

event_dose23 <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input * 0                               
  return(y)
}

# egg lay start at day 9
event_NoDose <- function(t, y, parms) {
  y[1]       <- 0                               
  return(y)
}

root_repeated <- function(t, y, parms) {
  yroot      <- c(((t-0.38)%%1),   ((t-0.38)/9)%%1, ((t-1.38)/9)%%1, ((t-2.38)/9)%%1, 
                  ((t-3.38)/9)%%1, ((t-4.38)/9)%%1, ((t-5.38)/9)%%1, ((t-6.38)/9)%%1, 
                  ((t-7.38)/9)%%1, ((t-8.38)/9)%%1, (t-17.38),
                  (t-64.38),        (t-65.38),      (t-66.38),        (t-67.38),
                  (t-68.38),        (t-69.38),      (t-70.38),        (t-71.38),
                  (t-126.38),       (t-127.38), 
                  (t-123.38),       (t-124.38),     (t-125.38),                       # egg laying events
                  
                  (t-dt)%%1,     (t-0.04)%%1,           (t-0.08)%%1,     (t-0.12)%%1, 
                  (t-0.16)%%1,   (t-0.21)%%1,           (t-0.25)%%1,     (t-0.29)%%1, 
                  (t-0.33)%%1,   (t-0.37)%%1,           (t-0.42)%%1,     (t-0.46)%%1, 
                  (t-0.5)%%1,    (t-0.54)%%1,           (t-0.58)%%1,     (t-0.62)%%1, 
                  (t-0.67)%%1,   (t-0.71)%%1,           (t-0.75)%%1,     (t-0.79)%%1, 
                  (t-0.83)%%1,   (t-0.87)%%1,           (t-0.92)%%1,     (t-0.96)%%1,  # feeding events
                  (t-1.67),      
                  (t-8.29),              (t-8.37),        (t-8.62),
                  (t-9.29),              (t-9.37),        (t-9.62),
                  (t-10.29),              (t-10.37),        (t-10.62),
                  (t-11.29),              (t-11.37),        (t-11.62),
                  (t-12.29),              (t-12.37),        (t-12.62),
                  (t-13.29),              (t-13.37),        (t-13.62),
                  (t-14.29),              (t-14.37),        (t-14.62),
                  (t-15.29),              (t-15.37),        (t-15.62),
                  (t-16.33),              (t-16.58),
                  (t-17.33),              (t-17.58),
                  (t-18.33),              (t-18.58),
                  (t-19.33),              (t-19.58),
                  (t-20.33),              (t-20.58),
                  (t-21.33),              (t-21.58),
                  (t-22.33),              (t-22.58),
                  (t-23.33),              (t-23.58),
                  (t-24.33),              (t-24.58),
                  (t-25.33),              (t-25.58),
                  (t-26.33),              (t-26.58),
                  (t-27.33),              (t-27.58),
                  (t-9)
                  
  )
  return(yroot)
}

event_repeated <- function(t, y, parms) {
  ret     <- y
  # SEQUENCE MATTERS
  # dose event
  if(abs((t-dt)%%1) < 1e-6)          ret <- event_dose0(t, y, parms)
  if(abs((t-0.04)%%1) < 1e-6)        ret <- event_dose1(t, y, parms)
  if(abs((t-0.08)%%1) < 1e-6)        ret <- event_dose2(t, y, parms)
  if(abs((t-0.12)%%1) < 1e-6)        ret <- event_dose3(t, y, parms)
  if(abs((t-0.16)%%1) < 1e-6)        ret <- event_dose4(t, y, parms)
  if(abs((t-0.21)%%1) < 1e-6)        ret <- event_dose5(t, y, parms)
  if(abs((t-0.25)%%1) < 1e-6)        ret <- event_dose6(t, y, parms)
  if(abs((t-0.29)%%1) < 1e-6)        ret <- event_dose7(t, y, parms)
  if(abs((t-0.33)%%1) < 1e-6)        ret <- event_dose8(t, y, parms)
  if(abs((t-0.37)%%1) < 1e-6)        ret <- event_dose9(t, y, parms)
  if(abs((t-0.42)%%1) < 1e-6)        ret <- event_dose10(t, y, parms)
  if(abs((t-0.46)%%1) < 1e-6)        ret <- event_dose11(t, y, parms)
  if(abs((t-0.50)%%1) < 1e-6)        ret <- event_dose12(t, y, parms)
  if(abs((t-0.54)%%1) < 1e-6)        ret <- event_dose13(t, y, parms)
  if(abs((t-0.58)%%1) < 1e-6)        ret <- event_dose14(t, y, parms)
  if(abs((t-0.62)%%1) < 1e-6)        ret <- event_dose15(t, y, parms)
  if(abs((t-0.67)%%1) < 1e-6)        ret <- event_dose16(t, y, parms)
  if(abs((t-0.71)%%1) < 1e-6)        ret <- event_dose17(t, y, parms)
  if(abs((t-0.75)%%1) < 1e-6)        ret <- event_dose18(t, y, parms)
  if(abs((t-0.79)%%1) < 1e-6)        ret <- event_dose19(t, y, parms)
  if(abs((t-0.83)%%1) < 1e-6)        ret <- event_dose20(t, y, parms)
  if(abs((t-0.87)%%1) < 1e-6)        ret <- event_dose21(t, y, parms)
  if(abs((t-0.92)%%1) < 1e-6)        ret <- event_dose22(t, y, parms)
  if(abs((t-0.96)%%1) < 1e-6)        ret <- event_dose23(t, y, parms)
  
  if(abs(t-1.67) < 1e-6)             ret <- event_dose16(t, y, parms)
  if(abs(t-8.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-8.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-8.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-9.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-9.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-9.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-10.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-10.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-10.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-11.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-11.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-11.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-12.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-12.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-12.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-13.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-13.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-13.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-14.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-14.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-14.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-15.29) < 1e-6)             ret <- event_dose7(t, y, parms)
  if(abs(t-15.37) < 1e-6)             ret <- event_dose9(t, y, parms)
  if(abs(t-15.62) < 1e-6)             ret <- event_dose15(t, y, parms)
  
  if(abs(t-16.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-16.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-17.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-17.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-18.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-18.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-19.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-19.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-20.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-20.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-21.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-21.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-22.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-22.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-23.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-23.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-24.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-24.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-25.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-25.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-26.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-26.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if(abs(t-27.33) < 1e-6)             ret <- event_dose8(t, y, parms)
  if(abs(t-27.58) < 1e-6)             ret <- event_dose14(t, y, parms)
  
  if((t-9) < 1e-6)                    ret <- event_NoDose(t, y, parms)  # actual dose start at day 9 because egg start lay at day 9
  
  # egg laying event
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
  if(abs(t-125.38) < 1e-6)         ret <- event9(t, y, parms)   # egg laying events above
  
  return((ret))
}

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



event_dose <- function(t, y, parms) {
  y[1]       <- y[1] + Oral_input                                  
  return(y)
}

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



root_repeated <- function(t, y, parms) {
  yroot      <- c(((t-0.38)%%1),   ((t-0.38)/9)%%1, ((t-1.38)/9)%%1, ((t-2.38)/9)%%1, 
                  ((t-3.38)/9)%%1, ((t-4.38)/9)%%1, ((t-5.38)/9)%%1, ((t-6.38)/9)%%1, 
                  ((t-7.38)/9)%%1, ((t-8.38)/9)%%1, (t-17.38), 
                  (t-64.38),        (t-65.38),      (t-66.38),        (t-67.38),
                  (t-68.38),        (t-69.38),      (t-70.38),        (t-71.38),
                  (t-126.38),       (t-127.38), 
                  (t-123.38),       (t-124.38),     (t-125.38), 
                  (t-0.3),          (t-1.3),        (t-2.3),          (t-3.3))
  return(yroot)
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
  if(abs(t-1.3) < 1e-6)              ret <- event_dose(t, y, parms)
  if(abs(t-2.3) < 1e-6)              ret <- event_dose(t, y, parms)
  if(abs(t-3.3) < 1e-6)              ret <- event_dose(t, y, parms)
  #if(abs(t-4) < 1e-6)              ret <- event_dose(t, y, parms)
  
  return((ret))
}


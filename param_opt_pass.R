########################################################
#Autoregulation Model parameter optimization
########################################################
rm(list=ls())
source("libs.R")
source("parms.R")
source("funcs_3_glomSS.R")
train_data <- readRDS("train_data.RDS")
#Passive parameter optimization - C_pass0 and C_pass1

P_aa_input <- c(100,125,150)      #Data from Takenaka, 1994, diltiazem
D_aa_input <- train_data$D_dilt 
D_try_EA <- DEA_0*(1 + (D_aa_input/DAA_0 - 1)*r_dilt)

Tp_aa_output <- 0*D_aa_input
P_aa_output <- 0*D_aa_input
  
for (i in seq(P_aa_input)){
  DAA_m <- D_aa_input[i]
  DEA_m <- D_try_EA[i]
  source("prep_anat.R")
  
  Rinf_init <- 1/(k*D*L*pi*0.1)
  
  OP_cont   <-      run_glom(src = src,
                       trg = trg,
                       D = D,
                       L = L,
                       k = k*src/src,
                       Rinf_init = Rinf_init,
                       Paa_i = P_aa_input[i],
                       Pea_o = Pea_o,
                       P_bs = P_bs,
                       Ca_i = Ca_i,
                       H_t_sys = H_t_sys,
                       mu_plas = mu_plas,
                       Rra=Rra,
                       t = t,
                       hpod=hpod,
                       wpod=wpod,
                       in.src = in.nodes,
                       out.trg = out.nodes,
                       num.iter = 200,
                       mu.tol = 1e-4,
                       Rinf.tol = 1e-4,
                       beta=2)

  Tp_aa_output[i] <- (OP_cont$G$Pm[src%in%in.nodes]-P_ext)*OP_cont$G$D[src%in%in.nodes]/2
  P_aa_output[i] <- OP_cont$G$Pm[src%in%in.nodes]
}

y <- log(Tp_aa_output)
x <- D_aa_input/DAA_0 - 1

A <- (sum(y)*sum(x^2)-sum(x)*sum(x*y))/(length(y)*sum(x^2)-(sum(x))^2)
C_pass0 <- exp(A)

B <- (length(y)*sum(x*y)-sum(x)*sum(y))/(length(y)*sum(x^2)-(sum(x))^2)
C_pass1 <- B

err <- ((Tp_aa_output - C_pass0*exp(C_pass1*(D_aa_input/DAA_0-1))))/Tp_aa_output*100
err  #% Relative error in tension prediction

autoreg_params <- data.frame(C_pass0=C_pass0,
                             C_pass1=C_pass1)

saveRDS(autoreg_params,"autoreg_params.RDS")


########################################################
#Active parameter optimization - C_act0: C_act1 assumed equal to 0.54
#Based on control state conditions
#source("parms.R")

DEA_m <- DEA_0
DAA_m <- DAA_0
source("prep_anat.R")

Rinf_init <- 1/(k*D*L*pi*0.1)

OP_cont   <-      run_glom(src = src,
                           trg = trg,
                           D = D,
                           L = L,
                           k = k*src/src,
                           Paa_i = Paa_i,
                           Pea_o = Pea_o,
                           P_bs = P_bs,
                           Ca_i = Ca_i,
                           H_t_sys = H_t_sys,
                           mu_plas = mu_plas,
                           Rra=Rra,
                           Rinf_init=Rinf_init,
                           t = t,
                           hpod=hpod,
                           wpod=wpod,
                           in.src = in.nodes,
                           out.trg = out.nodes,
                           num.iter = 150,
                           mu.tol = 1e-5,
                           Rinf.tol = 1e-5,
                           beta=beta)

Tp_cont <- (OP_cont$G$Pm[src%in%in.nodes]-P_ext)*OP_cont$G$D[src%in%in.nodes]/2
C_act0 <- 2*(Tp_cont - C_pass0)       #Control state condition in which A=1/2
C_act1 <- 0.54

autoreg_params <- data.frame(autoreg_params,
                             C_act0=C_act0,
                             C_act1=C_act1,
                             Tp_cont=Tp_cont)

saveRDS(autoreg_params,"autoreg_params.RDS")

#####################################################
#Myogenic parameter optimization - C_act2, C_myo and C_tone by Newton's Method
########################################################

rm(list=ls())
source("libs.R")
source("parms.R")
source("funcs_3_glomSS.R")
autoreg_params <- readRDS("autoreg_params.RDS")
train_data <- readRDS("train_data.RDS")

C_pass0 <- autoreg_params$C_pass0
C_pass1 <- autoreg_params$C_pass1
C_act0 <- autoreg_params$C_act0
C_act1 <- autoreg_params$C_act1

P_aa_input <- c(100,125,150)      #Data from Takenaka, 1994, furosemide
D_aa_input <- train_data$D_furo 

Tp_aa_output <- 0*D_aa_input
P_aa_output <- 0*D_aa_input
dd0 <- D_aa_input/DAA_0 - 1


for (i in seq(P_aa_input)){
  DAA_m <- D_aa_input[i]
  DEA_m <- DEA_0
  source("prep_anat.R")
  Rinf_init <- 1/(k*D*L*pi*0.1)
  OP_cont   <-      run_glom(src = src,
                             trg = trg,
                             D = D,
                             L = L,
                             k = k*src/src,
                             Paa_i = P_aa_input[i],
                             Pea_o = Pea_o,
                             P_bs = P_bs,
                             Ca_i = Ca_i,
                             H_t_sys = H_t_sys,
                             mu_plas = mu_plas,
                             Rra=Rra,
                             Rinf_init=Rinf_init,
                             t = t,
                             hpod=hpod,
                             wpod=wpod,
                             in.src = in.nodes,
                             out.trg = out.nodes,
                             num.iter = 150,
                             mu.tol = 1e-5,
                             Rinf.tol = 1e-5,
                             beta=beta)
  
  Tp_aa_output[i] <- (OP_cont$G$Pm[src%in%in.nodes]-P_ext)*OP_cont$G$D[src%in%in.nodes]/2
  P_aa_output[i] <- OP_cont$G$Pm[src%in%in.nodes]
}

start.data <- data.frame(P_aa_input = P_aa_input,
                       Tp_aa_output=Tp_aa_output,
                       S_tone=-log(C_act0*exp(-(dd0/C_act1)^2)/
                                      (Tp_aa_output - C_pass0*exp(C_pass1*dd0)) - 1))

plot(start.data$Tp_aa_output,start.data$S_tone)

Smyo_P_line <- lm(S_tone ~ P_aa_input, data=start.data)
summary(Smyo_P_line)

bS <- Smyo_P_line$coefficients[[1]]
mS <- Smyo_P_line$coefficients[[2]]

S80 <- bS + mS*80
S180 <- bS + mS*180

Smyo_Tp_line <- lm( Tp_aa_output ~ S_tone, data=start.data)
summary(Smyo_Tp_line)

bTp <- Smyo_Tp_line$coefficients[[1]]
mTp <- Smyo_Tp_line$coefficients[[2]]

Tp80 <- bTp + mTp*S80
Tp180 <- bTp + mTp*S180

tot.data <- data.frame(Tp = c(Tp80, Tp_aa_output, Tp180),
                       S_tone = c(S80, start.data$S_tone, S180))
save.tot.data <- data.frame(tot.data, P_aa_input=c(80,P_aa_input,180))

saveRDS(save.tot.data, "Myogenic_tone_fit_data.RDS")

plot(tot.data)
Tp_mid <- (Tp80 + Tp180)/2
S_mid <- (S80 + S180)/2
Tp_try <- 200:500

S_myo_try <- function(Tp_try,C){
  out <- (S180-S80)/(1+exp(-C*(Tp_try - Tp_mid))) + S80
  return(out)
}

Myo_mod <- nls(S_tone ~ S_myo_try(Tp_try = Tp,C = C), 
              data = tot.data, 
              start=list(C=0.01),
              control=list(maxiter=5000,tol=5e-5))

summary(Myo_mod)
C_myo <- coef(Myo_mod)[[1]]
lines(Tp_try,S_myo_try(Tp_try,C_myo))

Splus_myo <- S180-S80
Smin_myo <- S80
Tp_myo <- Tp_mid

autoreg_params <- data.frame(autoreg_params,
                             Splus_myo=Splus_myo, 
                             Smin_myo=Smin_myo,
                             Tp_myo=Tp_myo,
                             C_myo=C_myo)

saveRDS(autoreg_params,"autoreg_params.RDS")

###################################################
#TGF parameters from Darwin Bell 1982 "Relationship ..."
rm(list=ls())
source("libs.R")
source("parms.R")
source("funcs_3_glomSS.R")
autoreg_params <- readRDS("autoreg_params_20230711.RDS")
train_data <- readRDS("train_data.RDS")

C_pass0 <- autoreg_params$C_pass0
C_pass1 <- autoreg_params$C_pass1
C_act0 <- autoreg_params$C_act0
C_act1 <- autoreg_params$C_act1

Splus_myo <- autoreg_params$Splus_myo
Smin_myo <- autoreg_params$Smin_myo
C_myo <- autoreg_params$C_myo
Tp_myo <- autoreg_params$Tp_myo

Tp_cont <- autoreg_params$Tp_cont

conc <- c(50,100,125)      #Data from Bell, 1982
SFP_delt <- c(-3, -10, -13) 
D_furo <- train_data$D_furo[1] #Reference diameter

D_aa_input <- c(D_furo,seq(from = 7.5, to = 5, by = -0.5))        #Takenaka control: 19
#Pbs_input <- seq(from = 37, to = 22, by = -3)  
SNGFR <- 0*D_aa_input #matrix(0,ncol = length(Pbs_input),nrow= length(D_aa_input))
P_GC_output <- 0*D_aa_input
P_aa_output <- 0*D_aa_input

for (i in seq(D_aa_input)){
  #for(j in seq(Pbs_input)){
  DEA_m <- DEA_0
  DAA_m <- D_aa_input[i]
  source("prep_anat.R")
  Rinf_init <- 1/(k*D*L*pi*0.1)
  OP_cont   <-      run_glom(src = src,
                             trg = trg,
                             D = D,
                             L = L,
                             k = k*src/src,
                             Paa_i = Paa_i,
                             Pea_o = Pea_o,
                             P_bs = P_bs,
                             Ca_i = Ca_i,
                             H_t_sys = H_t_sys,
                             mu_plas = mu_plas,
                             Rra=Rra,
                             Rinf_init = Rinf_init,
                             t = t,
                             hpod=hpod,
                             wpod=wpod,
                             in.src = in.nodes,
                             out.trg = out.nodes,
                             num.iter = 150,
                             mu.tol = 1e-5,
                             Rinf.tol = 1e-5,
                             beta=beta)
  
  #Tp_aa_output[i] <- (OP_cont$G$Pm[src%in%in.nodes]-P_ext)*OP_cont$G$D[src%in%in.nodes]/2
  P_aa_output[i] <- OP_cont$G$Pm[src%in%in.nodes]
  P_GC_output[i] <- mean(OP_cont$G$Pm[!src%in%in.nodes & !trg%in%out.nodes])
  SNGFR[i] <- OP_cont$SNGFR
#}
}
plot(D_aa_input,P_GC_output-P_GC_output[1])
#plot(D_aa_input,log(-(P_GC_output-P_GC_output[1])))

D_P_dat <- data.frame(D=D_aa_input,
                      P=P_GC_output-P_GC_output[1])
mod <- lm(P ~ D, data=D_P_dat)
summary(mod)

bDP <- mod$coefficients[1]
mDP <- mod$coefficients[2]
D_try <- 4:9

plot(D_aa_input,P_GC_output-P_GC_output[1])
lines(D_try, bDP + mDP*D_try)
#points(log(SFP_delt/(-cDP))/kDP,SFP_delt,pch=4)

D_est <- (SFP_delt-bDP)/mDP
plot(conc,D_est)

D_est_p <- c(D_furo[1],D_est)#,D_est[length(D_est)])
conc_p <- c(25,conc,150)
Tp_aa_output <- 0*D_est_p
P_aa_output <- 0*D_est_p
SNGFR <- 0*D_est_p
  for (i in seq(D_est_p)){
    DEA_m <- DEA_0
  DAA_m <- D_est_p[i]
  source("prep_anat.R")
  Rinf_init <- 1/(k*D*L*pi*0.1)
  OP_cont   <-      run_glom(src = src,
                             trg = trg,
                             D = D,
                             L = L,
                             k = k*src/src,
                             Paa_i = Paa_i,
                             Pea_o = Pea_o,
                             P_bs = P_bs,
                             Ca_i = Ca_i,
                             H_t_sys = H_t_sys,
                             mu_plas = mu_plas,
                             Rra=Rra,
                             Rinf_init = Rinf_init,
                             t = t,
                             hpod=hpod,
                             wpod=wpod,
                             in.src = in.nodes,
                             out.trg = out.nodes,
                             num.iter = 150,
                             mu.tol = 1e-5,
                             Rinf.tol = 1e-5,
                             beta=beta)
  
  Tp_aa_output[i] <- (OP_cont$G$Pm[src%in%in.nodes]-P_ext)*OP_cont$G$D[src%in%in.nodes]/2
  P_aa_output[i] <- OP_cont$G$Pm[src%in%in.nodes]
  SNGFR[i] <- OP_cont$SNGFR
}

S_myo_try <- function(Tp_try){
  out <- Splus_myo/(1+exp(-C_myo*(Tp_try - Tp_myo))) + Smin_myo
  return(out)
}

D_est_p <- c(D_est_p,D_est_p[length(D_est_p)])
dd0 <- D_est_p/DAA_0-1
Tp_aa_output <- c(Tp_aa_output,Tp_aa_output[length(Tp_aa_output)])
#Tp_aa_output <- Tp_aa_output[1:5]
start.data <- data.frame(conc = conc_p,
                         Tp_aa_output=Tp_aa_output,
                         S_tone= - log(C_act0*exp(-(dd0/C_act1)^2)/
                                           (Tp_aa_output - C_pass0*exp(C_pass1*dd0))- 1),
                         S_tone_m= - log(C_act0*exp(-(dd0/C_act1)^2)/
                                       (Tp_aa_output - C_pass0*exp(C_pass1*dd0)) - 1)
                                        - S_myo_try(Tp_aa_output),
                         S_tone_d= - log(C_act0*exp(-(dd0/C_act1)^2)/
                                           (Tp_aa_output - C_pass0*exp(C_pass1*dd0)) - 1)/
                         S_myo_try(Tp_aa_output))
plot(0,0,xlim=c(0,150),ylim=c(-1,1))
plot(start.data$conc,start.data$S_tone_m)
#points(start.data$conc,start.data$S_tone_d)
x <- 0:150
#lines(x,(50-x)^2/x^2 + min(start.data$S_tone_d))# - 0.001*x - 0.4)
#plot(start.data$Tp_aa_output,start.data$S_tone_m)

#points(start.data$conc,start.data$S_tone,pch=5)

C_MD <- 0*D_est

S_max <- max(start.data$S_tone_m)
S_min <- -0.5
C_MD_max <- 150
C_MD_min <- 25
C_MD_mid <- 40#(C_MD_max+C_MD_min)/2

S_TGF_try <- function(C_try,C){
  out <- (S_max-S_min)/(1+exp(-C*(C_try - C_MD_mid))) + S_min
  #out <- (C1-C_try)^2/C_try^2 + C2*C_try + S_min
  return(out)
}


TGF_mod <- nls(S_tone_m ~ S_TGF_try(C_try = conc,C = C), 
               data = start.data, 
               start=list(C=0.1),
               control=list(maxiter=5000,tol=5e-5))

summary(TGF_mod)
C_TGF <- coef(TGF_mod)[[1]]
C_try <- 0:300
lines(C_try,S_TGF_try(C_try = C_try,C=C_TGF))

autoreg_params <- data.frame(autoreg_params[1:9], 
                             C_TGF=C_TGF, 
                             Splus_TGF = S_max-S_min,
                             Smin_TGF = S_min,
                             C_MD_TGF = C_MD_mid)

saveRDS(autoreg_params,"autoreg_params_20230711.RDS")
writeMat("autoreg_params_20230711.mat",autoreg_params=autoreg_params)
mtfd <- readRDS("Myogenic_tone_fit_data.RDS")
writeMat("Myogenic_tone_fit_data.mat",mtfd=mtfd)
writeMat("TGF_tone_fit_data.mat",ttfd=start.data)

###################################################
#Myogenic mechanism modulation by TGF from Takenaka control data
rm(list=ls())
surr_glom_df <- readRDS("surr_glom_20230705.RDS")
source("libs.R")
source("parms.R")
source("funcs_4.R")
autoreg_params <- readRDS("autoreg_params_20230711.RDS")
train_data <- readRDS("train_data.RDS")

C_pass0 <- autoreg_params$C_pass0
C_pass1 <- autoreg_params$C_pass1
C_act0 <- autoreg_params$C_act0
C_act1 <- autoreg_params$C_act1

Splus_myo <- autoreg_params$Splus_myo
Smin_myo <- autoreg_params$Smin_myo
C_myo <- autoreg_params$C_myo
Tp_myo <- autoreg_params$Tp_myo

Splus_TGF <- autoreg_params$Splus_TGF
Smin_TGF <- autoreg_params$Smin_TGF
C_TGF <- autoreg_params$C_TGF
C_MD_TGF <- autoreg_params$C_MD_TGF

Tp_cont <- autoreg_params$Tp_cont

D_aa_input <- train_data$D_cont        #Takenaka control: 19
P_aa_input <- c(100,125,150)
#Pbs_input <- seq(from = 37, to = 22, by = -3)  
SNGFR <- 0*D_aa_input #matrix(0,ncol = length(Pbs_input),nrow= length(D_aa_input))
Tp_aa_output <- 0*D_aa_input
C_MD_output <- 0*D_aa_input

for (i in seq(D_aa_input)){
  CT <- check_takenaka(D=D_aa_input[i], 
                 Pa=P_aa_input[i], 
                 Ca=Ca_i,
                 surr_glom_df=surr_glom_df,
                 PT_frac = PT_frac*P_aa_input/P_aa_input,
                 C_0_Tub = C_0_Tub,
                 furo_yn = 0,
                 dilt_yn = 0,
                 int_yn=0)
  
  Tp_aa_output[i] <- CT$Tp
  C_MD_output[i] <- CT$C_MD
  SNGFR[i] <- CT$SNGFR
}

S_TGF <- Splus_TGF/(1+exp(-C_TGF*(C_MD_output - C_MD_TGF))) + Smin_TGF
S_myo <- Splus_myo/(1+exp(-C_myo*(Tp_aa_output - Tp_myo))) + Smin_myo

dd0 <- D_aa_input/DAA_0-1

start.data <- data.frame(conc = C_MD_output,
                         Tp_aa_output=Tp_aa_output,
                         S_myo=S_myo,
                         S_TGF=S_TGF,
                         S_tone= - log(C_act0*exp(-(dd0/C_act1)^2)/
                                         (Tp_aa_output - C_pass0*exp(C_pass1*dd0))- 1)
                                          - S_TGF)
plot(P_aa_input,start.data$conc)

plot(P_aa_input,start.data$S_tone)

lines(P_aa_input,start.data$S_myo)

S_myo_try <- function(Tp_try,alpha,beta){
  out <- Splus_myo/(1+exp(-beta*Tp_try + alpha)) + Smin_myo
  return(out)
}

S_myo_try_a <- function(Tp_try,alpha){
  out <- Splus_myo/(1+exp(-C_myo*Tp_try + alpha)) + Smin_myo
  return(out)
}

S_myo_try_b <- function(Tp_try,beta){
  out <- Splus_myo/(1+exp(-beta*(Tp_try - Tp_myo))) + Smin_myo
  return(out)
}

MyoTGF_mod <- nls(S_tone ~ S_myo_try(Tp_aa_output,alpha = alpha,beta=beta), 
               data = start.data, 
               start=list(alpha=C_myo*Tp_cont,beta=C_myo),
               control=list(maxiter=10000,tol=5e-5))
summary(MyoTGF_mod)

MyoTGF_mod_a <- nls(S_tone ~ S_myo_try_a(Tp_aa_output,alpha = alpha), 
                  data = start.data, 
                  start=list(alpha=C_myo*Tp_cont),
                  control=list(maxiter=10000,tol=5e-5))
summary(MyoTGF_mod_a)

lines(P_aa_input, S_myo_try(Tp_try=Tp_aa_output,alpha=coef(MyoTGF_mod)[1],beta=coef(MyoTGF_mod)[2]))

MyoTGF_mod_b <- nls(S_tone ~ S_myo_try_b(Tp_aa_output,beta = beta), 
                    data = start.data, 
                    start=list(beta=C_myo),
                    control=list(maxiter=10000,tol=5e-5))
summary(MyoTGF_mod_b)


C_int <- coef(MyoTGF_mod)[[2]]
Tp_TGF <- coef(MyoTGF_mod)[[1]]/C_int

autoreg_params <- data.frame(autoreg_params[1:13],Tp_TGF = Tp_TGF, C_int=C_int)

saveRDS(autoreg_params,"autoreg_params_20230711.RDS")

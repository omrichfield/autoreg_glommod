#Adjust Model to Fit Takenaka Data
###########################################
#Step 1: Get model's control Q, translate Takenaka data to model terms
rm(list=ls())
source("libs.R")
source("parms.R")
source("funcs_3_glomSS.R")

P_try <- Paa_i#c(100,125,150)
D_try <- DAA_0#seq(from = 22.2, to = 26.2, by = 2)
P_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))
Q_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))
mu_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))
SNGFR_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))
j<-1
i<-1

#for (i in seq(P_try)){
#for (j in 1){#seq(D_try)){
DEA_m <- DEA_0
DAA_m <- D_try[j]
source("prep_anat.R")

Rinf_init <- 1/(k*D*L*pi*0.1)

OP_cont   <-      run_glom(src = src,
                           trg = trg,
                           D = D,
                           L = L,
                           k = k*src/src,
                           Rinf_init = Rinf_init,
                           Paa_i = P_try[i],
                           Pea_o = Pea_o,
                           P_bs = P_bs,
                           Ca_i = Ca_i,
                           H_t_sys = H_t_sys,
                           mu_plas = mu_plas,
                           Rra=1.09,
                           t = t,
                           hpod=hpod,
                           wpod=wpod,
                           in.src = in.nodes,
                           out.trg = out.nodes,
                           num.iter = 150,
                           mu.tol = 1e-4,
                           Rinf.tol = 1e-4,
                           beta=2)
Q_cont_mod <- OP_cont$G$Q[src%in%in.nodes]
mean(OP_cont$G$Pm[!(src%in%in.nodes | trg%in%out.nodes)])
OP_cont$G$mu[src%in%in.nodes]
OP_cont$SNGFR
#}
#}



#Takenaka Data

#Control flow
Q_cont_100 <- 203.333
Q_cont_rel <- c(203.333,200,200)/Q_cont_100
Q_cont <- Q_cont_mod*Q_cont_rel

#Diltiazem flow

Q_dilt_rel <- c(380,533.333,713.333)/Q_cont_100
Q_dilt <- Q_cont_mod*Q_dilt_rel


#Furosemide flow

Q_furo_rel <- c(250,280,326.666)/Q_cont_100
Q_furo <- Q_cont_mod*Q_furo_rel


##################################################
#Step 2: Use Hayashi data to evaluate EA changes as AA changes, to fit to diltiazem data
#Then perform simulations to get Q's for when both AA and EA change

#rm(list=ls())
source("libs.R")
source("parms.R")
source("funcs_3_glomSS.R")

#Hayashi Data

#Diltiazem Efferent/Efferent Arteriole

r_AA <- 110.8 #Relative increase in AA diameter at 10^(-5)M dilt
r_EA <- 21.6  #Relative increase in EA diameter at 10^(-5)M dilt
r_dilt <- r_EA/r_AA

P_try <- c(100,125,150)
D_try <- 7:15
D_try_EA <- DEA_0*(1 + (D_try/DAA_0 - 1)*r_dilt)
P_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))
Q_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))
mu_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))
SNGFR_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))

for (i in seq(P_try)){
  for (j in seq(D_try)){
    DAA_m <- D_try[j]
    DEA_m <- D_try_EA[j]
    source("prep_anat.R")
    
    Rinf_init <- 1/(k*D*L*pi*0.1)
    
    OP_cont   <-      run_glom(src = src,
                               trg = trg,
                               D = D,
                               L = L,
                               k = k*src/src,
                               Rinf_init = Rinf_init,
                               Paa_i = P_try[i],
                               Pea_o = Pea_o,
                               P_bs = P_bs,
                               Ca_i = Ca_i,
                               H_t_sys = H_t_sys,
                               mu_plas = mu_plas,
                               Rra=1.09,
                               t = t,
                               hpod=hpod,
                               wpod=wpod,
                               in.src = in.nodes,
                               out.trg = out.nodes,
                               num.iter = 200,
                               mu.tol = 1e-4,
                               Rinf.tol = 1e-4,
                               beta=2)
    Q_keep[i,j] <- OP_cont$G$Q[src%in%in.nodes]
    P_keep[i,j] <- mean(OP_cont$G$Pm[!(src%in%in.nodes | trg%in%out.nodes)])
    mu_keep[i,j] <- OP_cont$G$mu[src%in%in.nodes]
    SNGFR_keep[i,j] <- OP_cont$SNGFR
  }
}
plot(D_try,Q_keep[1,],type='l')
plot(D_try,SNGFR_keep[1,],type='l')
plot(D_try,P_keep[1,],type='l')

D_dilt <- 0*Q_dilt

for (i in seq(Q_dilt)){
  Q_curve <- Q_keep[i,]
  SQ <- sort(abs(Q_curve-Q_dilt[i]),index.return=TRUE)
  indx <- SQ$ix[1:2]
  Qval <- Q_keep[i,indx]
  m <- (Qval[1]-Qval[2])/(D_try[indx[1]]-D_try[indx[2]])
  b <- Qval[1] - m*D_try[indx[1]]
  D_dilt[i] <- (Q_dilt[i]-b)/m
}

##################################################
#Step 3: Run simulations to map Q's for control and furosemide cases to diameters

#rm(list=ls())
source("libs.R")
source("parms.R")
source("funcs_3_glomSS.R")

P_try <- c(100,125,150)
D_try <- 7:15
#D_try_EA <- DEA_0*(1 + (D_try/DAA_0 - 1)*r_dilt)
P_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))
Q_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))
mu_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))
SNGFR_keep <- matrix(0,nrow=length(P_try),ncol=length(D_try))

for (i in seq(P_try)){
  for (j in seq(D_try)){
    DAA_m <- D_try[j]
    DEA_m <- DEA_0
    source("prep_anat.R")
    
    Rinf_init <- 1/(k*D*L*pi*0.1)
    
    OP_cont   <-      run_glom(src = src,
                               trg = trg,
                               D = D,
                               L = L,
                               k = k*src/src,
                               Rinf_init = Rinf_init,
                               Paa_i = P_try[i],
                               Pea_o = Pea_o,
                               P_bs = P_bs,
                               Ca_i = Ca_i,
                               H_t_sys = H_t_sys,
                               mu_plas = mu_plas,
                               Rra=1.09,
                               t = t,
                               hpod=hpod,
                               wpod=wpod,
                               in.src = in.nodes,
                               out.trg = out.nodes,
                               num.iter = 200,
                               mu.tol = 1e-4,
                               Rinf.tol = 1e-4,
                               beta=2)
    Q_keep[i,j] <- OP_cont$G$Q[src%in%in.nodes]
    P_keep[i,j] <- mean(OP_cont$G$Pm[!(src%in%in.nodes | trg%in%out.nodes)])
    mu_keep[i,j] <- OP_cont$G$mu[src%in%in.nodes]
    SNGFR_keep[i,j] <- OP_cont$SNGFR
  }
}
plot(D_try,Q_keep[1,],type='l')
plot(D_try,SNGFR_keep[1,],type='l')
plot(D_try,P_keep[1,],type='l')

D_furo <- 0*Q_furo

for (i in seq(Q_furo)){
  Q_curve <- Q_keep[i,]
  SQ <- sort(abs(Q_curve-Q_furo[i]),index.return=TRUE)
  indx <- SQ$ix[1:2]
  Qval <- Q_keep[i,indx]
  m <- (Qval[1]-Qval[2])/(D_try[indx[1]]-D_try[indx[2]])
  b <- Qval[1] - m*D_try[indx[1]]
  D_furo[i] <- (Q_furo[i]-b)/m
}

D_cont <- 0*Q_cont

for (i in seq(Q_cont)){
  Q_curve <- Q_keep[i,]
  SQ <- sort(abs(Q_curve-Q_cont[i]),index.return=TRUE)
  indx <- SQ$ix[1:2]
  Qval <- Q_keep[i,indx]
  m <- (Qval[1]-Qval[2])/(D_try[indx[1]]-D_try[indx[2]])
  b <- Qval[1] - m*D_try[indx[1]]
  D_cont[i] <- (Q_cont[i]-b)/m
}

train_data <- data.frame(D_cont=D_cont,Q_cont=Q_cont,
                         D_furo=D_furo,Q_furo=Q_furo,
                         D_dilt=D_dilt,Q_dilt=Q_dilt,
                         DEA_dilt=DEA_0*(1 + (D_dilt/DAA_0 - 1)*r_dilt))

saveRDS(train_data,"train_data.RDS")
writeMat("train_data.mat",train_data=train_data)
#############################################################












# #DILTIAZEM DATA
# 
# P_aa_dilt <- c(100,125,150,100,125,150)           #Takenaka, 1994
# #P_aa_dilt <- c(P_aa_dilt,100,125,150,100,125,150) #Feng, 2003
# #P_aa_dilt <- c(P_aa_dilt,rep(c(40,80,120,180),1)) #Sanchez-Ferrer, 1989
# 
# D_aa_dilt <- c(DAA_0/19*c(24,24.5,25.1),DAA_0/19*25.1/31*c(29,30,31))
# #D_aa_dilt <- c(D_aa_dilt,25.1/22*c(21,21.5,22),25.1/23.5*c(21.5,22.5,23.5))
# #D_aa_dilt <- c(D_aa_dilt,25/1.15*c(1,1.03,1.05,1.15))
# 
# PT_frac_dilt <- rep(PT_frac,6)
# #PT_frac_dilt <- c(PT_frac_dilt,rep(PT_frac,6))
# #PT_frac_dilt <- c(PT_frac_dilt,rep(PT_frac,4))
# 
# ############################################
# 
# #FUROSEMIDE DATA
# 
# P_aa_furo <- c(100,125,150)#,100,125,150)         #Takenaka, 1994
# #P_aa_furo <- c(P_aa_furo,c(100,148))            #Walker Thesis, 2001
# #P_aa_furo <- c(P_aa_furo,rep(c(80,120,180),2))            #Walker Thesis, 2001
# 
# D_aa_furo <- c(DAA_0/19*c(21,19.5,18.5))#,c(22.5,22,21))
# #D_aa_furo <- c(D_aa_furo,DAA_0/17*c(19,18))
# #D_aa_furo <- c(D_aa_furo,DAA_0/0.925*c(0.96,1.05,1.04),DAA_0/0.96*c(0.95,1.075,1.1))
# 
# PT_frac_furo <- rep(PT_frac,3)
# #PT_frac_furo <- c(PT_frac_furo,rep(PT_frac,2))
# #PT_frac_furo <- c(PT_frac_furo,rep(PT_frac,6))
# 
# ############################################
# 
# #CONTROL DATA
# 
# P_aa_cont <- c(100,125,150,100,125,150)               #Takenaka, 1994  
# P_aa_cont <- c(P_aa_cont,c(100,148,100,148))          #Walker Thesis, 2001
# P_aa_cont <- c(P_aa_cont,c(100,125,150))              #Feng, 2003
# P_aa_cont <- c(P_aa_cont,rep(c(100,140),8))           #Kulthinee, 2020
# P_aa_cont <- c(P_aa_cont,rep(c(80,120,180),3))     #sanchez-ferrer, 1989
# 
# D_aa_cont <- c(DAA_0/19*c(19,16.7,14),DAA_0/21*c(21,18,16))
# D_aa_cont <- c(D_aa_cont,DAA_0/17*c(17,14,15.5,11))
# D_aa_cont <- c(D_aa_cont,DAA_0/19*c(19,16,14))
# D_aa_cont <- c(D_aa_cont,DAA_0/13*c(13,9.5), DAA_0/13*c(13,9.5),
#                DAA_0/13*c(13,10),DAA_0/15*c(15,12),
#                DAA_0/15*c(15,12),DAA_0/14.5*c(14.5,11),
#                DAA_0/15*c(15,12),DAA_0/15*c(15,12))
# D_aa_cont <- c(D_aa_cont,DAA_0/0.975*c(1,0.95,0.82),
#                DAA_0/0.925*c(1,0.85,0.78),DAA_0/0.96*c(1,0.92,0.73))
# 
# PT_frac_cont <- rep(PT_frac,6)
# PT_frac_cont <- c(PT_frac_cont,rep(PT_frac,2),rep(PT_frac_ACZ,2))
# PT_frac_cont <- c(PT_frac_cont,rep(PT_frac,3))
# PT_frac_cont <- c(PT_frac_cont,rep(PT_frac,16))
# PT_frac_cont <- c(PT_frac_cont,rep(PT_frac,9))

############################################

#ONCOTIC DATA

# Ca_onc <- c(Ca_i,4)
# D_onc <- c(DAA_0,8.45)


# a1 <- 2.1#1.629
# a2 <- 0.16#0.2935
# a3 <- 0.009
# 
# Ca_TGF_test <- c(0.1, 2.82, 4.66, 5.72)                 #Schurek 1981
# Ca_TGF_test <- c(Ca_TGF_test, c(3.4, 5.8, 3.6))                       #Baylis 1977
# Ca_TGF_test <- c(Ca_TGF_test, c(4.3, 6, 4.4, 5.9))                    #Tucker 1981
#                  #Thomas 1974
# Ca_TGF_test <- c(Ca_TGF_test, Ca_i, Ca_i)                             #Blantz 1974
# 
# Q_TGF_test <- c(15.6, 24.1, 30.5, 32)/32
# Q_TGF_test <- c(Q_TGF_test, c(92.3, 170.5, 79.2)/170.5)
# Q_TGF_test <- c(Q_TGF_test, c(169, 285, 145, 272)/272)
# Q_TGF_test <- c(Q_TGF_test, c(1,1))
# 
# D_TGF_test <- 0*Q_TGF_test
# 
# for (i in seq(Ca_TGF_test)){
#   indx <- which(surr_glom_df$Ca.in < Ca_TGF_test[i] + 0.25 &
#                   surr_glom_df$Ca.in > Ca_TGF_test[i] - 0.25 &
#                   surr_glom_df$P.in == 100)
#   d_try <- surr_glom_df[indx,]
#   indx_top <- order(abs(Q_TGF_test[i]*130 - d_try$Q))[1:5]
# 
#   d_spec <- d_try[indx_top,]
#   weights <- 1/abs(Q_TGF_test[i]*130 - d_spec$Q)
#   D_TGF_test[i] <- sum(d_spec$D.in*weights)/sum(weights)
# }
# 
# Ca_onc <- Ca_TGF_test[Ca_TGF_test > 4]
# D_onc <- D_TGF_test[Ca_TGF_test > 4]
# plot(Ca_onc,D_onc)

# S6 <- surr_glom_df[surr_glom_df$P.in==Paa_i &
#                surr_glom_df$Ca.in==6 &
#                surr_glom_df$D.in==DAA_0,]$SNGFR
#
# S7 <- surr_glom_df[surr_glom_df$P.in==Paa_i &
#                surr_glom_df$Ca.in==5 &
#                surr_glom_df$D.in==DAA_0,]$SNGFR
#
# df_7 <- surr_glom_df[surr_glom_df$P.in==Paa_i &
#                      surr_glom_df$Ca.in==5,]
# df_7p <- df_7[order(abs(df_7$SNGFR - S6))[1:2],]
#
# weights <- 1/abs(df_7p$SNGFR - S6)^4
# D_Ca_test <- sum(weights*df_7p$D.in)/sum(weights)
#
# Ca_onc <- c(6,5)
# D_onc <- c(DAA_0, D_Ca_test)
##########################################################

#SHEAR STRESS DATA

# D_SS <- c(24,26,28)/24*DAA_0
# P_SS <- c(100,125,150)

# D_denom <- 0*Ca_TGF_test
# D_numer <- 0*Ca_TGF_test
#
# D_denom[Ca_TGF_test < 4] <- Q_TGF_test[Ca_TGF_test < 4]*130.6622 -
#                             int_Ql -
#                             Ca_TGF_test[Ca_TGF_test < 4]*
#                             (coef_Ca.in_Ql + coef_Ca.in.P.in_Ql*Paa_i) -
#                             coef_P.in_Ql*Paa_i
# 
# D_denom[Ca_TGF_test >= 4] <- Q_TGF_test[Ca_TGF_test >= 4]*130.6622 -
#                             int_Q -
#                             Ca_TGF_test[Ca_TGF_test >= 4]*
#                             (coef_Ca.in_Q + coef_Ca.in.P.in_Q*Paa_i) -
#                             coef_P.in_Q*Paa_i
# 
# D_numer[Ca_TGF_test < 4] <- coef_invD4_Ql  +
#                             Ca_TGF_test[Ca_TGF_test < 4]*
#                             coef_invD4.Ca.in.P.in_Ql*Paa_i +
#                             coef_invD4.P.in_Ql*Paa_i
# 
# D_numer[Ca_TGF_test >= 4] <- coef_invD4_Q  +
#                               Ca_TGF_test[Ca_TGF_test >= 4]*
#                               (coef_invD4.Ca.in.P.in_Q*Paa_i+
#                                  coef_invD4.Ca.in_Q) +
#                               coef_invD4.P.in_Q*Paa_i
# 
# D_TGF_test <- (D_numer/D_denom)^(1/4)

# D_denom_e <- Q_TGF_test_e*160.4786 - int_Qe - coef_Ca.in_Qe*Ca_TGF_test_e - 
#                   coef_Ca.in.P.in_Qe*Ca_TGF_test_e*(Paa_i+10) - coef_P.in_Qe*(Paa_i+10)  
# D_numer_e <- coef_invD4_Qe + coef_invD4.P.in_Qe*(Paa_i+10) + coef_invD4.Ca.in_Qe*Ca_TGF_test_e
# D_TGF_test_e <- (D_numer_e/D_denom_e)^(1/4)
# D_TGF_test_e

# D_onc <- D_TGF_test[is.na(D_TGF_test)==FALSE]
# Ca_onc <- Ca_TGF_test[is.na(D_TGF_test)==FALSE]
# 
# plot(Ca_TGF_test, Q_TGF_test)
# data <- data.frame(y=log(Q_TGF_test), x=Ca_TGF_test)
# 
# fits <- lm(y~x, data=data)
# summary(fits)
# Q_Ca_coef1 <- exp(as.numeric(fits$coefficients[1]))
# Q_Ca_coef2 <- as.numeric(fits$coefficients[2])
# 
# lines(seq(from=0, to=8, by=0.1), Q_Ca_coef1*exp(Q_Ca_coef2*seq(from=0, to=8, by=0.1)))
# Ca_TGF_test_e <- c(5.85/5.06*Ca_i) 
# Ca_TGF_test_e <- c(Ca_TGF_test_e,Ca_i*max(Re(cubic(c(a3,a2,a1,-23.4))))/max(Re(cubic(c(a3,a2,a1,-15.1)))))                             #Blantz 1974
# #Ca_TGF_test_e <- c(Ca_TGF_test_e, 7.59)    



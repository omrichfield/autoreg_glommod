rm(list=ls())
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
#C_TGF2 <- autoreg_params$C_TGF2
C_MD_TGF <- autoreg_params$C_MD_TGF

Tp_cont <- autoreg_params$Tp_cont

Tp_TGF <- autoreg_params$Tp_TGF
C_int <- autoreg_params$C_int
#C_0 <- autoreg_params$C_0

D_aa_furo <- train_data$D_furo
D_aa_cont <- train_data$D_cont
D_aa_dilt <- train_data$D_dilt

Try_Pf <- seq(100,150,2.5)
int_low <- 8
int_high <- 10.5
num_mins <- 5
move_sure <- 5
mins_tol <- 0.1
surr_glom_df <- readRDS("surr_glom_dilt_20230412.RDS")

#jpeg('Autoreg_Pf.jpg')

tak_err_dilt <- tak_err(
  P_aa_input=Try_Pf,
  D_aa_input=c(max(D_aa_dilt),min(D_aa_dilt),mean(D_aa_dilt)*rep(1,length(Try_Pf)-2)),
  Ca_input = Ca_i*Try_Pf/Try_Pf,
  surr_glom_df = surr_glom_df,
  PT_frac = PT_frac*Try_Pf/Try_Pf,
  C_0_Tub = C_0_Tub,
  furo_yn = 0,
  dilt_yn=1,
  int_yn=0,
  int_low=int_low,
  int_high=int_high)
D_pred_dilt <- tak_err_dilt$D_zer
Q_pred_dilt <- tak_err_dilt$Q
plot(Try_Pf,D_pred_dilt,type='l',ylim=c(5,11),xlab="Perfusion Pressure (mmHg)",ylab="AA Diameter, um")
points(c(100,125,150),D_aa_dilt)

#plot(Try_Pf,tak_err_furo$err)
#plot(Try_Pf,tak_err_furo$Tp,type='l')
#plot(Try_Pf,D_pred_furo,type='l',ylim=c(5,8),xlab="Perfusion Pressure (mmHg)",ylab="AA Diameter, um")
#points(c(100,125,150),D_aa_furo)

Try_Pf <- seq(100,150,2.5)
int_low <- 5
int_high <- 8
num_mins <- 5
move_sure <- 5
mins_tol <- 0.1
surr_glom_df <- readRDS("surr_glom_df_20220724.RDS")

tak_err_furo <- tak_err(
  P_aa_input=Try_Pf,
  D_aa_input=c(max(D_aa_furo),min(D_aa_furo),mean(D_aa_furo)*rep(1,length(Try_Pf)-2)),
  Ca_input = Ca_i*Try_Pf/Try_Pf,
  surr_glom_df = surr_glom_df,
  PT_frac = PT_frac*Try_Pf/Try_Pf,
  C_0_Tub = C_0_Tub,
  furo_yn = 1,
  dilt_yn=0,
  int_yn=0,
  int_low=int_low,
  int_high=int_high)
D_pred_furo <- tak_err_furo$D_zer
S_TGF_pred <- tak_err_furo$S_TGF

#plot(Try_Pf,tak_err_furo$err)
#plot(Try_Pf,tak_err_furo$Tp,type='l')
#plot(c(100,150),c(5,8))
lines(Try_Pf,D_pred_furo,type='l')
points(c(100,125,150),D_aa_furo)


tak_err_cont <- tak_err(
  P_aa_input=Try_Pf,
  D_aa_input=c(max(D_aa_cont),min(D_aa_cont),mean(D_aa_cont)*rep(1,length(Try_Pf)-2)),
  Ca_input = Ca_i*Try_Pf/Try_Pf,
  surr_glom_df = surr_glom_df,
  PT_frac = PT_frac*Try_Pf/Try_Pf,
  C_0_Tub = C_0_Tub,
  furo_yn = 0,
  dilt_yn=0,
  int_yn=0,
  int_low=int_low,
  int_high=int_high)
D_pred_cont <- tak_err_cont$D_zer
S_TGF_pred <- tak_err_cont$S_TGF

#plot(Try_Pf,tak_err_cont$err)
#plot(Try_Pf,tak_err_furo$Tp,type='l')
lines(Try_Pf,D_pred_cont,type='l',ylim=c(5,8))
points(c(100,125,150),D_aa_cont)

tak_err_cont_int <- tak_err(
  P_aa_input=Try_Pf,
  D_aa_input=c(max(D_aa_cont),min(D_aa_cont),mean(D_aa_cont)*rep(1,length(Try_Pf)-2)),
  Ca_input = Ca_i*Try_Pf/Try_Pf,
  surr_glom_df = surr_glom_df,
  PT_frac = PT_frac*Try_Pf/Try_Pf,
  C_0_Tub = C_0_Tub,
  furo_yn = 0,
  dilt_yn=0,
  int_yn=1,
  int_low=int_low,
  int_high=int_high)
D_pred_cont_int <- tak_err_cont_int$D_zer
S_TGF_pred <- tak_err_cont$S_TGF

#Get Steady-state mechanics (SS, HS, CSGFR/SA)

DEA_m <- DEA_0#*(1 + (D_aa_input[i_D]/DAA_0 - 1)*r_dilt)

DAA_m <- DAA_0
source("prep_anat.R")

P_aa_input <- seq(100,150,10)
glom_list <- list()
Rinf0 <- 1/(k*D*L*pi*0.1)
for (i_P in seq(length(P_aa_input))){
  DEA_m <- DEA_0#*(1 + (D_aa_input[i_D]/DAA_0 - 1)*r_dilt)
  DAA_m <- DAA_0
  source("prep_anat.R")
  OP_cont   <-      run_glom(src = src,
                             trg = trg,
                             D = D,
                             L = L,
                             k = k*src/src,
                             Rinf_init=Rinf0,
                             Paa_i = P_aa_input[i_P],
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
                             num.iter = 150,
                             mu.tol = 1e-3,
                             Rinf.tol = 1e-3,
                             beta=beta)
glom_list[[i_P]] <- OP_cont
}

df <- data.frame(P.in=P_aa_input, 
                 SNGFR=0*P_aa_input,
                 P.glom_est=0*P_aa_input,
                 Q=0*P_aa_input,
                 hoop_mn=0*P_aa_input,
                 hoop_sd=0*P_aa_input,
                 csgfr_sa_mn=0*P_aa_input,
                 csgfr_sa_sd=0*P_aa_input,
                 shear_mn=0*P_aa_input,
                 shear_sd=0*P_aa_input)

for (i in seq(P_aa_input)){
  G <- glom_list[[i]]$G
  
  G_SA <- (G$D[!src%in%in.nodes & !trg%in%out.nodes]*pi*G$L[!src%in%in.nodes & !trg%in%out.nodes])
  G_L <- G$L[!src%in%in.nodes & !trg%in%out.nodes]
  df$SNGFR[i] <- glom_list[[i]]$SNGFR
  df$P.glom_est[i] <- G$Pe[src%in%in.nodes]
  df$Q[i] <- G$Q[src%in%in.nodes]
  
  df$hoop_mn[i] <- mean(G$hoop[!src%in%in.nodes & !trg%in%out.nodes])#*G_L/sum(G_L))
  df$hoop_sd[i] <- sd(G$hoop[!src%in%in.nodes & !trg%in%out.nodes])#*G_L/sum(G_L))
  
  df$csgfr_sa_mn[i] <- mean(G$CSGFR[!src%in%in.nodes & !trg%in%out.nodes])#*G_SA/sum(G_SA))
  df$csgfr_sa_sd[i] <- sd(G$CSGFR[!src%in%in.nodes & !trg%in%out.nodes])#*G_SA/sum(G_SA))
  
  df$shear_mn[i] <- mean(G$shear[!src%in%in.nodes & !trg%in%out.nodes])#*G_SA/sum(G_SA))
  df$shear_sd[i] <- sd(G$shear[!src%in%in.nodes & !trg%in%out.nodes])#*G_SA/sum(G_SA))
  }

Q_shear_mn_lm <- lm(shear_mn ~ Q, data=df)
summary(Q_shear_mn_lm)
Q_shear_sd_lm <- lm(shear_sd ~ Q, data=df)
summary(Q_shear_sd_lm)

P_hoop_mn_lm <- lm(hoop_mn ~ P.glom_est, data=df)
summary(P_hoop_mn_lm)
P_hoop_sd_lm <- lm(hoop_sd ~ P.glom_est, data=df)
summary(P_hoop_sd_lm)

SNGFR_CSGFR_mn_lm <- lm(csgfr_sa_mn ~ SNGFR, data=df)
summary(SNGFR_CSGFR_mn_lm)
SNGFR_CSGFR_sd_lm <- lm(csgfr_sa_sd ~ SNGFR, data=df)
summary(SNGFR_CSGFR_sd_lm)

helper_func <- function(lm, x){
  coef1 <- lm$coefficients[1]
  coef2 <- lm$coefficients[2]
  return(coef1 + coef2*x)
}



#plot(Try_Pf,tak_err_cont$err)
#plot(Try_Pf,tak_err_furo$Tp,type='l')
lines(Try_Pf,D_pred_cont_int,type='l',ylim=c(5,8))
#points(c(100,125,150),D_aa_cont)

#dev.off()

# plot(Try_Pf, tak_err_cont$S_Myo,type='l')
# lines(Try_Pf, tak_err_cont$S_TGF,col='blue')
# 
# plot(Try_Pf,tak_err_cont$SNGFR,type='l')
# plot(tak_err_cont$SNGFR,tak_err_cont$S_TGF)
# 
# lines(Try_Pf,tak_err_furo$Tp,type='l')

plot(c(100,150),c(50,350))
lines(Try_Pf, tak_err_dilt$Q,type='l')
lines(Try_Pf, tak_err_furo$Q,type='l')
lines(Try_Pf, tak_err_cont$Q,type='l')
lines(Try_Pf, tak_err_cont_int$Q,type='l')

plot(c(100,150),c(-1,2.5))
lines(Try_Pf, tak_err_furo$S_Myo,type='l')
lines(Try_Pf, tak_err_cont$S_Myo,type='l')
lines(Try_Pf, tak_err_cont_int$S_Myo,type='l')

plot(c(100,150),c(-1.5,2))
lines(Try_Pf, tak_err_furo$S_TGF,type='l')
lines(Try_Pf, tak_err_cont$S_TGF,type='l')
lines(Try_Pf, tak_err_cont_int$S_TGF,type='l')

plot(c(100,150),c(10,80))
lines(Try_Pf, tak_err_furo$SNGFR,type='l')
lines(Try_Pf, tak_err_cont$SNGFR,type='l')
lines(Try_Pf, tak_err_cont_int$SNGFR,type='l')

S_myo_try <- function(Tp_try,beta,alpha){
  out <- Splus_myo/(1+exp(-beta*(Tp_try - alpha))) + Smin_myo
  return(out)
}

plot(c(100,500),c(-1,4))
lines(0:600, S_myo_try(Tp_try=0:600,beta=C_myo,alpha=Tp_myo) ,type='l')
lines(0:600, S_myo_try(Tp_try=0:600,beta=C_int,alpha=Tp_TGF) ,type='l')

S_TGF_try <- function(C_try){
  out <- Splus_TGF/(1+exp(-C_TGF*(C_try - C_MD_TGF))) + Smin_TGF
  #out <- (C1-C_try)^2/C_try^2 + C2*C_try + S_min
  return(out)
}

plot(c(0,200),c(-1,4))
lines(0:200, S_TGF_try(C_try=0:200) ,type='l')


lines(Try_Pf, tak_err_cont$SNGFR,type='l')


writeMat("Myo_TGF_model_curves_20231001.mat",D_pred_dilt=D_pred_dilt,
                                    D_pred_furo=D_pred_furo,
                                    D_pred_cont=D_pred_cont,
                                    D_pred_cont_int=D_pred_cont_int,
                                    Q_pred_dilt=Q_pred_dilt,
                                    Q_pred_furo=tak_err_furo$Q_out,
                                    Q_pred_cont=tak_err_cont$Q_out,
                                    Q_pred_cont_int=tak_err_cont_int$Q_out,
                                    SNGFR_dilt=tak_err_dilt$SNGFR,
                                    SNGFR_furo=tak_err_furo$SNGFR,
                                    SNGFR_cont=tak_err_cont$SNGFR,
                                    SNGFR_int=tak_err_cont_int$SNGFR,
         C_MD_dilt=tak_err_dilt$C_MD,
         C_MD_furo=tak_err_furo$C_MD,
         C_MD_cont=tak_err_cont$C_MD,
         C_MD_int=tak_err_cont_int$C_MD,
         vloop_dilt=tak_err_dilt$Qloop,
         vloop_furo=tak_err_furo$Qloop,
         vloop_cont=tak_err_cont$Qloop,
         vloop_int=tak_err_cont_int$Qloop,
         PG_dilt=tak_err_dilt$PG,
         PG_furo=tak_err_furo$PG,
         PG_cont=tak_err_cont$PG,
         PG_int=tak_err_cont_int$PG,
                                    S_TGF_furo=tak_err_furo$S_TGF,
                                    S_TGF_cont=tak_err_cont$S_TGF,
                                    S_TGF_cont_int=tak_err_cont_int$S_TGF,
                                    S_Myo_furo=tak_err_furo$S_Myo,
                                    S_Myo_cont=tak_err_cont$S_Myo,
                                    S_Myo_cont_int=tak_err_cont_int$S_Myo,
                                    Try_Pf=Try_Pf,
                                    Tp_show=100:600,
                                    C_MD_show=0:200,
                                    S_TGF_show=S_TGF_try(C_try=0:200),
                                    S_Myo_show_furo=S_myo_try(Tp_try=100:600,beta=C_myo,alpha=Tp_myo),
                                    S_Myo_show_int=S_myo_try(Tp_try=100:600,beta=C_int,alpha=Tp_TGF),
         CSGFR_mn_furo = helper_func(lm=SNGFR_CSGFR_mn_lm, x=tak_err_furo$SNGFR),
         CSGFR_sd_furo = helper_func(lm=SNGFR_CSGFR_sd_lm, x=tak_err_furo$SNGFR),
         
         CSGFR_mn_cont = helper_func(lm=SNGFR_CSGFR_mn_lm, x=tak_err_cont$SNGFR),
         CSGFR_sd_cont = helper_func(lm=SNGFR_CSGFR_sd_lm, x=tak_err_cont$SNGFR),
         
         CSGFR_mn_int = helper_func(lm=SNGFR_CSGFR_mn_lm, x=tak_err_cont_int$SNGFR),
         CSGFR_sd_int = helper_func(lm=SNGFR_CSGFR_sd_lm, x=tak_err_cont_int$SNGFR),
         
         CSGFR_mn_dilt = helper_func(lm=SNGFR_CSGFR_mn_lm, x=tak_err_dilt$SNGFR),
         CSGFR_sd_dilt = helper_func(lm=SNGFR_CSGFR_sd_lm, x=tak_err_dilt$SNGFR),
   
         shear_mn_furo = helper_func(lm=Q_shear_mn_lm, x=tak_err_furo$Q_out/0.6),
         shear_sd_furo = helper_func(lm=Q_shear_sd_lm, x=tak_err_furo$Q_out/0.6),
         
         shear_mn_cont = helper_func(lm=Q_shear_mn_lm, x=tak_err_cont$Q_out/0.6),
         shear_sd_cont = helper_func(lm=Q_shear_sd_lm, x=tak_err_cont$Q_out/0.6),
         
         shear_mn_int = helper_func(lm=Q_shear_mn_lm, x=tak_err_cont_int$Q_out/0.6),
         shear_sd_int = helper_func(lm=Q_shear_sd_lm, x=tak_err_cont_int$Q_out/0.6),
         
         shear_mn_dilt = helper_func(lm=Q_shear_mn_lm, x=tak_err_dilt$Q_out/0.6),
         shear_sd_dilt = helper_func(lm=Q_shear_sd_lm, x=tak_err_dilt$Q_out/0.6),
         
         hoop_mn_furo = helper_func(lm=P_hoop_mn_lm, x=tak_err_furo$PG),
         hoop_sd_furo = helper_func(lm=P_hoop_sd_lm, x=tak_err_furo$PG),
         
         hoop_mn_cont = helper_func(lm=P_hoop_mn_lm, x=tak_err_cont$PG),
         hoop_sd_cont = helper_func(lm=P_hoop_sd_lm, x=tak_err_cont$PG),
         
         hoop_mn_int = helper_func(lm=P_hoop_mn_lm, x=tak_err_cont_int$PG),
         hoop_sd_int = helper_func(lm=P_hoop_sd_lm, x=tak_err_cont_int$PG),
         
         hoop_mn_dilt = helper_func(lm=P_hoop_mn_lm, x=tak_err_dilt$PG),
         hoop_sd_dilt = helper_func(lm=P_hoop_sd_lm, x=tak_err_dilt$PG)
         )

###########################################################################
rm(list=ls())
source("libs.R")
source("parms.R")
source("funcs_4.R")
autoreg_params <- readRDS("autoreg_params_20230711.RDS")
train_data <- readRDS("train_data.RDS")


Try_Pf <- c(100,130)#seq(100,150,2.5)
int_low <- 5
int_high <- 8
num_mins <- 5
move_sure <- 5
mins_tol <- 0.1
surr_glom_df <- readRDS("surr_glom_df_20220724.RDS")

D_aa_furo <- train_data$D_furo
D_aa_cont <- train_data$D_cont
D_aa_dilt <- train_data$D_dilt

autoreg_param_modify <- rep(1,12)

Q_keep <- matrix(0, nrow=length(autoreg_param_modify), ncol=length(Try_Pf))

for (i in seq(autoreg_param_modify)){
  
  autoreg_param_modify <- rep(1,12)
  autoreg_param_modify[i] <- 1.01
  
C_pass0 <- autoreg_params$C_pass0*autoreg_param_modify[1]
C_pass1 <- autoreg_params$C_pass1*autoreg_param_modify[2]
C_act0 <- autoreg_params$C_act0*autoreg_param_modify[3]
C_act1 <- autoreg_params$C_act1*autoreg_param_modify[4]

Splus_myo <- autoreg_params$Splus_myo*autoreg_param_modify[5]
Smin_myo <- autoreg_params$Smin_myo*autoreg_param_modify[6]

C_myo <- autoreg_params$C_myo
Tp_myo <- autoreg_params$Tp_myo

Splus_TGF <- autoreg_params$Splus_TGF*autoreg_param_modify[7]
Smin_TGF <- autoreg_params$Smin_TGF*autoreg_param_modify[8]
C_TGF <- autoreg_params$C_TGF*autoreg_param_modify[9]
C_MD_TGF <- autoreg_params$C_MD_TGF*autoreg_param_modify[10]

Tp_cont <- autoreg_params$Tp_cont

Tp_TGF <- autoreg_params$Tp_TGF*autoreg_param_modify[11]
C_int <- autoreg_params$C_int*autoreg_param_modify[12]



tak_err_int <- tak_err(
  P_aa_input=Try_Pf,
  D_aa_input=c(max(D_aa_cont),min(D_aa_cont),mean(D_aa_cont)*rep(1,length(Try_Pf)-2)),
  Ca_input = Ca_i*Try_Pf/Try_Pf,
  surr_glom_df = surr_glom_df,
  PT_frac = PT_frac*Try_Pf/Try_Pf,
  C_0_Tub = C_0_Tub,
  furo_yn = 0,
  dilt_yn=0,
  int_yn=1,
  int_low=int_low,
  int_high=int_high)
Q_keep[i,] <- tak_err_int$Q

}

autoreg_param_modify <- rep(1,12)

C_pass0 <- autoreg_params$C_pass0*autoreg_param_modify[1]
C_pass1 <- autoreg_params$C_pass1*autoreg_param_modify[2]
C_act0 <- autoreg_params$C_act0*autoreg_param_modify[3]
C_act1 <- autoreg_params$C_act1*autoreg_param_modify[4]

Splus_myo <- autoreg_params$Splus_myo*autoreg_param_modify[5]
Smin_myo <- autoreg_params$Smin_myo*autoreg_param_modify[6]

C_myo <- autoreg_params$C_myo
Tp_myo <- autoreg_params$Tp_myo

Splus_TGF <- autoreg_params$Splus_TGF*autoreg_param_modify[7]
Smin_TGF <- autoreg_params$Smin_TGF*autoreg_param_modify[8]
C_TGF <- autoreg_params$C_TGF*autoreg_param_modify[9]
C_MD_TGF <- autoreg_params$C_MD_TGF*autoreg_param_modify[10]

Tp_cont <- autoreg_params$Tp_cont

Tp_TGF <- autoreg_params$Tp_TGF*autoreg_param_modify[11]
C_int <- autoreg_params$C_int*autoreg_param_modify[12]

tak_err_int <- tak_err(
  P_aa_input=Try_Pf,
  D_aa_input=c(max(D_aa_cont),min(D_aa_cont),mean(D_aa_cont)*rep(1,length(Try_Pf)-2)),
  Ca_input = Ca_i*Try_Pf/Try_Pf,
  surr_glom_df = surr_glom_df,
  PT_frac = PT_frac*Try_Pf/Try_Pf,
  C_0_Tub = C_0_Tub,
  furo_yn = 0,
  dilt_yn=0,
  int_yn=1,
  int_low=int_low,
  int_high=int_high)
Q_base <- tak_err_int$Q

NSC <- 0*Q_keep

for (i in seq(autoreg_param_modify)){
  NSC[i,] <- (Q_keep[i,]-Q_base)/Q_base/0.01
}

NSC









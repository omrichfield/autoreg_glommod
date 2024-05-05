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

Try_Pf <- c(100,130)#seq(100,150,2.5)
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
  TGF_only_yn=0,
  int_low=int_low,
  int_high=int_high)
D_pred_dilt <- tak_err_dilt$D_zer

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
  TGF_only_yn=0,
  int_low=int_low,
  int_high=int_high)
D_pred_furo <- tak_err_furo$D_zer

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
  TGF_only_yn=0,
  int_low=int_low,
  int_high=int_high)
D_pred_cont <- tak_err_cont$D_zer

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
  TGF_only_yn=0,
  int_low=int_low,
  int_high=int_high)
D_pred_int <- tak_err_int$D_zer

Try_Pf <- c(100,130)#seq(100,150,2.5)
int_low <- 6
int_high <- 10.5
num_mins <- 5
move_sure <- 5
mins_tol <- 0.1
surr_glom_df <- readRDS("surr_glom_20230412.RDS")


tak_err_TGF_only <- tak_err(
  P_aa_input=Try_Pf,
  D_aa_input=c(max(D_aa_cont),min(D_aa_cont),mean(D_aa_cont)*rep(1,length(Try_Pf)-2)),
  Ca_input = Ca_i*Try_Pf/Try_Pf,
  surr_glom_df = surr_glom_df,
  PT_frac = PT_frac*Try_Pf/Try_Pf,
  C_0_Tub = C_0_Tub,
  furo_yn = 0,
  dilt_yn=0,
  int_yn=0,
  TGF_only_yn=1,
  int_low=int_low,
  int_high=int_high)
D_pred_TGF <- tak_err_TGF_only$D_zer

D_AA_try <- rbind(D_pred_cont,
                  D_pred_furo,
                  D_pred_dilt,
                  D_pred_int,
                  D_pred_TGF)
D_EA_try <- rbind(DEA_0*D_pred_cont/D_pred_cont,
                  DEA_0*D_pred_cont/D_pred_cont,
                  DEA_0*(1 + (D_pred_dilt/DAA_0 - 1)*r_dilt)*D_pred_cont/D_pred_cont,
                  DEA_0*D_pred_cont/D_pred_cont,
                  DEA_0*D_pred_cont/D_pred_cont)

DEA_m <- DAA_0
DAA_m <- DEA_0
source("prep_anat.R")

delt_shear <- matrix(rep(0,dim(D_AA_try)[1]*length(src)),
                     ncol = length(src), 
                     nrow =  dim(D_AA_try)[1])

strain <- delt_shear

delt_CSGFR <- delt_shear

for (i in seq(dim(D_AA_try)[1])){
#strain/delt shear/ delt csgfr

DEA_m <- D_EA_try[i,1]
DAA_m <- D_AA_try[i,1]
source("prep_anat.R")
Rinf_init <- 1/(k*D*L*pi*0.1)

OP_cont_deflate   <-      compliant_glom(src=src,
                                  trg=trg,
                                  D=D,
                                  L=L,
                                  Rinf_init = Rinf_init,
                                  k=k,
                                  E_cmp=YM,
                                  Paa_i=Try_Pf[1],
                                  Paa_try=Try_Pf[1]-20,
                                  Pea_o= Pea_o,
                                  P_bs = P_bs,
                                  Ca_i= Ca_i,
                                  H_t_sys= H_t_sys,
                                  mu_plas= mu_plas,
                                  Rra= Rra,
                                  t=t,
                                  hpod=hpod,
                                  wpod=wpod,
                                  in.nodes=in.nodes,
                                  out.nodes=out.nodes,
                                  num.iter.g=150,
                                  num.iter.d=10,
                                  d.tol=1e-5,
                                  L.tol=1e-5,
                                  t.tol=1e-5,
                                  SNGFR.tol=1e-5,
                                  mu.tol=1e-5,
                                  Rinf.tol=1e-5,
                                  beta=3)

OP_cont_inflate   <-      compliant_glom(src=src,
                                    trg=trg,
                                    D=D,
                                    L=L,
                                    Rinf_init = Rinf_init,
                                    k=k,
                                    E_cmp=YM,
                                    Paa_i=Try_Pf[1],
                                    Paa_try=Try_Pf[1]+20,
                                    Pea_o= Pea_o,
                                    P_bs = P_bs,
                                    Ca_i= Ca_i,
                                    H_t_sys= H_t_sys,
                                    mu_plas= mu_plas,
                                    Rra= Rra,
                                    t=t,
                                    hpod=hpod,
                                    wpod=wpod,
                                    in.nodes=in.nodes,
                                    out.nodes=out.nodes,
                                    num.iter.g=150,
                                    num.iter.d=10,
                                    d.tol=1e-5,
                                    L.tol=1e-5,
                                    t.tol=1e-5,
                                    SNGFR.tol=1e-5,
                                    mu.tol=1e-5,
                                    Rinf.tol=1e-5,
                                    beta=3)

strain[i,] <- (OP_cont_inflate$OP$G$D-OP_cont_deflate$OP$G$D)/OP_cont_deflate$OP$G$D*100

delt_shear[i,] <- OP_cont_inflate$OP$G$shear - OP_cont_deflate$OP$G$shear

delt_CSGFR[i,] <- OP_cont_inflate$OP$G$CSGFR - OP_cont_deflate$OP$G$CSGFR

}

writeMat('transient_mech_100_20231012.mat',strain_100=strain, delt_shear_100=delt_shear, delt_CSGFR_100=delt_CSGFR)

D_AA_try <- rbind(D_pred_cont,
                  D_pred_furo,
                  D_pred_dilt,
                  D_pred_int,
                  D_pred_TGF)
D_EA_try <- rbind(DEA_0*D_pred_cont/D_pred_cont,
                  DEA_0*D_pred_cont/D_pred_cont,
                  DEA_0*(1 + (D_pred_dilt/DAA_0 - 1)*r_dilt)*D_pred_cont/D_pred_cont,
                  DEA_0*D_pred_cont/D_pred_cont,
                  DEA_0*D_pred_cont/D_pred_cont)

DEA_m <- DAA_0
DAA_m <- DEA_0
source("prep_anat.R")

delt_shear <- matrix(rep(0,dim(D_AA_try)[1]*length(src)),
                     ncol = length(src), 
                     nrow =  dim(D_AA_try)[1])

strain <- delt_shear

delt_CSGFR <- delt_shear

for (i in seq(dim(D_AA_try)[1])){
  #strain/delt shear/ delt csgfr
  
  DEA_m <- D_EA_try[i,2]
  DAA_m <- D_AA_try[i,2]
  source("prep_anat.R")
  Rinf_init <- 1/(k*D*L*pi*0.1)
  
  OP_cont_deflate   <-      compliant_glom(src=src,
                                           trg=trg,
                                           D=D,
                                           L=L,
                                           Rinf_init = Rinf_init,
                                           k=k,
                                           E_cmp=YM,
                                           Paa_i=Paa_i,
                                           Paa_try=Try_Pf[2]-20,
                                           Pea_o= Pea_o,
                                           P_bs = P_bs,
                                           Ca_i= Ca_i,
                                           H_t_sys= H_t_sys,
                                           mu_plas= mu_plas,
                                           Rra= Rra,
                                           t=t,
                                           hpod=hpod,
                                           wpod=wpod,
                                           in.nodes=in.nodes,
                                           out.nodes=out.nodes,
                                           num.iter.g=150,
                                           num.iter.d=10,
                                           d.tol=1e-5,
                                           L.tol=1e-5,
                                           t.tol=1e-5,
                                           SNGFR.tol=1e-5,
                                           mu.tol=1e-5,
                                           Rinf.tol=1e-5,
                                           beta=3)
  
  OP_cont_inflate   <-      compliant_glom(src=src,
                                           trg=trg,
                                           D=D,
                                           L=L,
                                           Rinf_init = Rinf_init,
                                           k=k,
                                           E_cmp=YM,
                                           Paa_i=Paa_i,
                                           Paa_try=Try_Pf[2]+20,
                                           Pea_o= Pea_o,
                                           P_bs = P_bs,
                                           Ca_i= Ca_i,
                                           H_t_sys= H_t_sys,
                                           mu_plas= mu_plas,
                                           Rra= Rra,
                                           t=t,
                                           hpod=hpod,
                                           wpod=wpod,
                                           in.nodes=in.nodes,
                                           out.nodes=out.nodes,
                                           num.iter.g=150,
                                           num.iter.d=10,
                                           d.tol=1e-5,
                                           L.tol=1e-5,
                                           t.tol=1e-5,
                                           SNGFR.tol=1e-5,
                                           mu.tol=1e-5,
                                           Rinf.tol=1e-5,
                                           beta=3)
  
  strain[i,] <- (OP_cont_inflate$OP$G$D-OP_cont_deflate$OP$G$D)/OP_cont_deflate$OP$G$D*100
  
  delt_shear[i,] <- OP_cont_inflate$OP$G$shear - OP_cont_deflate$OP$G$shear
  
  delt_CSGFR[i,] <- OP_cont_inflate$OP$G$CSGFR - OP_cont_deflate$OP$G$CSGFR
  
}

writeMat('transient_mech_130_20231012.mat',strain_130=strain, 
                                            delt_shear_130=delt_shear, 
                                            delt_CSGFR_130=delt_CSGFR)







rm(list=ls())
surr_glom_df <- readRDS("surr_glom_20230705.RDS")
source("libs.R")
source("parms.R")
source("funcs_4.R")
autoreg_params <- readRDS("autoreg_params.RDS")
train_data <- readRDS("train_data.RDS")


DAA_m <- DAA_0
DEA_m <- DEA_0
source("prep_anat.R")

Rinf_init <- 1/(k*D*L*pi*0.1)

OP_cont   <-      run_glom(src = src,
                           trg = trg,
                           D = D,
                           L = L,
                           k = k*src/src,
                           Rinf_init = Rinf_init,
                           Paa_i = Paa_i,
                           Pea_o = Pea_o,
                           P_bs = P_bs,
                           Ca_i = Ca_i,
                           H_t_sys = 0.4,
                           mu_plas = mu_plas,
                           Rra=1.09,
                           t = t,
                           hpod=hpod,
                           wpod=wpod,
                           in.src = in.nodes,
                           out.trg = out.nodes,
                           num.iter = 150,
                           mu.tol = 1e-5,
                           Rinf.tol = 1e-5,
                           beta=3)

writeMat('G_base_20230925.mat',G_base=OP_cont$G)


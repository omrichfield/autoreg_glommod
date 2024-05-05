rm(list=ls())
surr_glom_df <- readRDS("surr_glom_20230705.RDS")
source("libs.R")
source("parms.R")
source("funcs_4.R")
autoreg_params <- readRDS("autoreg_params.RDS")
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
Try_Pf <- seq(100,150,2.5)

Tubule_schem <- rbind(Try_Pf, 0*Try_Pf, 0*Try_Pf, 0*Try_Pf, 0*Try_Pf)

for (i in seq(Try_Pf)){
BB1 <- check_takenaka(D=DAA_0,
                      Pa=Try_Pf[i],
                      Ca = Ca_i,
                      surr_glom_df = surr_glom_df,
                      PT_frac=1,
                      C_0_Tub=C_0_Tub,
                      furo_yn=0,
                      dilt_yn=0,
                      int_yn=1,
                      TGF_only_yn=0)
Tubule_schem[2,i] <- BB1$Qloop
Tubule_schem[3,i] <- BB1$C_MD
Tubule_schem[4,i] <- BB1$SNGFR
Tubule_schem[5,i] <- BB1$Qloop*A_M*1e6*60
}

plot(Try_Pf,1-Tubule_schem[5,]/Tubule_schem[4,])

writeMat("Tubule_schem_P_mat.mat",Tubule_schem=Tubule_schem)

  BB1 <- check_takenaka(D=DAA_0,
                        Pa=Paa_i,
                        Ca = Ca_i,
                        surr_glom_df = surr_glom_df,
                        PT_frac=1,
                        C_0_Tub=C_0_Tub,
                        furo_yn=0,
                        dilt_yn=0,
                        int_yn=1)

jj <- readRDS("Tubule_schem_rds.RDS")
writeMat("Tubule_schem_L_mat.mat",Tubule_schem_L=jj)

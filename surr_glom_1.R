########################################################
#Surrogate Glomerular Model
########################################################
rm(list=ls())
source("libs.R")
source("parms.R")
source("funcs_3_glomSS.R") 

m <- c(0.05,0.1,0.15)
mm <- seq(0,0.8,0.2)


P_aa_input <- seq(from=77.5,to=92.5,by=2.5)  #Takenaka control: 100
D_aa_input <- seq(from = 5, to = 10.5, by = 0.05)        #Takenaka control: 19
Ca_input <- Ca_i#seq(from = 4.6, to = 6.6, by = 0.2)  #Takenaka control: 100

P_aa_output <- matrix(nrow = length(P_aa_input), ncol = length(D_aa_input))
P_ae_output <- matrix(nrow = length(P_aa_input), ncol = length(D_aa_input))
SNGFR_output <- matrix(nrow = length(P_aa_input), ncol = length(D_aa_input))
Q_output <- matrix(nrow = length(P_aa_input), ncol = length(D_aa_input))
Pi_output <- matrix(nrow = length(P_aa_input), ncol = length(D_aa_input))
mu_output <- matrix(nrow = length(P_aa_input), ncol = length(D_aa_input))
shear_output <- matrix(nrow = length(P_aa_input), ncol = length(D_aa_input))

r_AA <- 110.8 #Relative increase in AA diameter at 10^(-5)M dilt
r_EA <- 21.6  #Relative increase in EA diameter at 10^(-5)M dilt
r_dilt <- r_EA/r_AA

DEA_m <- DEA_0
DAA_m <- D_aa_input[1]
source("prep_anat.R")
tot_out <- list()
count_num_old <- 0

start_time <- Sys.time()
  
#for (i_Ca in seq(length(Ca_input))){
  Rinf0 <- 1/(k*D*L*pi*0.1)
for (i_D in seq(length(D_aa_input))){
for (i_P in seq(length(P_aa_input))){
  DEA_m <- DEA_0#*(1 + (D_aa_input[i_D]/DAA_0 - 1)*r_dilt)
  DAA_m <- D_aa_input[i_D]
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
  
  P_aa_output[i_P,i_D] <- OP_cont$G$Pm[src%in%in.nodes]
  P_ae_output[i_P,i_D] <- OP_cont$G$Pe[src%in%in.nodes]
  SNGFR_output[i_P,i_D] <- OP_cont$SNGFR
  Q_output[i_P,i_D] <- OP_cont$G$Q[1]-OP_cont$G$E[1]
  Pi_output[i_P,i_D] <- (OP_cont$G$Pi[src%in%in.nodes] + OP_cont$G$Pi[trg%in%out.nodes])/2
  mu_output[i_P,i_D] <- OP_cont$G$mu[src%in%in.nodes]
  shear_output[i_P,i_D] <- OP_cont$G$shear[src%in%in.nodes]
  #Rinf0 <- OP_cont$G$Rinf
  
}
}
  dat_Ca_i <- data.frame(Pa=melt(P_aa_output), 
                         SNGFR=melt(SNGFR_output), 
                         Q=melt(Q_output),
                         Pi=melt(Pi_output),
                         Pae=melt(P_ae_output),
                         mu=melt(mu_output),
                         shear=melt(shear_output)) 
  tot_out <- dat_Ca_i

end_time <- Sys.time()
#CHECK if you have small diameter, small Ca


surr_glom_1_df <- data.frame(P.in = c(),
                             D.in = c(),
                             Ca.in = c(),
                             Pavg = c(),
                             SNGFR = c(),
                             Q = c(),
                             mu = c(),
                             shear = c())
tot_out_i <- tot_out
surr_glom_1_df <- data.frame(P.in = P_aa_input[tot_out_i$Pa.Var1],
                             D.in = D_aa_input[tot_out_i$Pa.Var2],
                             Ca.in = rep(Ca_i,length(tot_out_i$Pa.Var1)),
                             Pavg = tot_out_i$Pa.value,
                             SNGFR = tot_out_i$SNGFR.value,
                             Q = tot_out_i$Q.value)


saveRDS(surr_glom_1_df,"surr_glom_low_20230703.rds")
##################################################
surr_glom_df_1 <- readRDS("surr_glom_20230412.RDS")
surr_glom_df_2 <- readRDS("surr_glom_low_20230703.rds")
surr_glom_df <- rbind(surr_glom_df_1,surr_glom_df_2)
surr_glom_df <- surr_glom_df[order(surr_glom_df),]
ss <- surr_glom_df[!duplicated(surr_glom_df[,1:3]),]

#ss <- aggregate(surr_glom_df,list(D=surr_glom_df$D.in, P=surr_glom_df$P.in),mean)
#surr_glom_df <- ss[,3:ncol(ss)
saveRDS(ss,"surr_glom_20230705.RDS")


#################################################################
#Test glomSS lookup on table
rm(list=ls())
source("libs.R")
source("parms.R")
source("funcs_3_glomSS.R") 
surr_glom_df1 <- readRDS("surr_glom_df_20220724.RDS")
surr_glom_df2 <- readRDS("surr_glom_tak_8_20220724")

surr_glom_df <- unique(rbind(surr_glom_df1,surr_glom_df2))
surr_glom_df <- surr_glom_df[order(surr_glom_df$D.in),]# surr_glom_df$P.in),]
surr_glom_df <- surr_glom_df[order(surr_glom_df$P.in),]

saveRDS(surr_glom_df,"surr_glom_df_20220724.RDS")



P1_try <- 111
P2_try <- 134

P_aa_input <- c(P1_try,P2_try)
colors <- c('red','blue')#,'black')#,'green','magenta')
plot(0,0,xlim=c(5,8), ylim=c(0,130))
D_try <- seq(5,8,0.05)
MG <- 0*D_try
SNGFR <- matrix(0, length(P_aa_input), length(D_try))
for (i in seq(P_aa_input)){
  for (j in seq(D_try)){
    print(D_try[j])
    MG[j] <- glom_SS(Pa=P_aa_input[i], D=D_try[j], Ca=Ca_i, surr_glom_df = surr_glom_df)$SNGFR
    SNGFR[i,j] <- MG[j]
  }
  #indx <- which(surr_glom_df$P.in==P_aa_input[i])
  #points(surr_glom_df$D.in[indx],surr_glom_df$SNGFR[indx],col=colors[i])
  # MG <- multi_glom_SS(Pa_vec=P_aa_input[i]*D_try/D_try,
  #                     Ca_vec=Ca_i*D_try/D_try,
  #                     D_vec=D_try,
  #                     surr_glom_df = surr_glom_df)
  lines(D_try,MG,col=colors[i])
}



plot(0,0,xlim=c(5,8), ylim=c(0,150))
for (i in seq(P_aa_input)){
  indx <- which(surr_glom_df$P.in==P_aa_input[i])
  points(surr_glom_df$D.in[indx],surr_glom_df$Pavg[indx],col=colors[i])
  MG <- multi_glom_SS(Pa_vec=P_aa_input[i]*D_try/D_try,
                      Ca_vec=Ca_i*D_try/D_try,
                      D_vec=D_try,
                      surr_glom_df = surr_glom_df)
  lines(D_try,MG$Pavg,col=colors[i])
}

D_aa_input <- seq(5,8,0.4)
Pavg_P1_check <- 0*D_aa_input
SNGFR_P1_check <- 0*D_aa_input

for (i_D in seq(length(D_aa_input))){
    DEA_m <- DEA_0
    DAA_m <- D_aa_input[i_D]
    Rinf0 <- 1/(k*D*L*pi*0.1)
    source("prep_anat.R")
    OP_cont   <-      run_glom(src = src,
                               trg = trg,
                               D = D,
                               L = L,
                               k = k*src/src,
                               Rinf_init=Rinf0,
                               Paa_i = P1_try,
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
    
    Pavg_P1_check[i_D] <- OP_cont$G$Pm[src%in%in.nodes]
    SNGFR_P1_check[i_D] <- OP_cont$SNGFR
  }

points(D_aa_input,Pavg_P1_check)



# #####################################################
# # # Post-processing
# 
# rm(list=ls())
# source("libs.R")
# source("parms.R")
# source("funcs_2_temp_1.R")
# 
# surr_glom_df <- readRDS("surr_glom_1_Pbs13")
# 
# # 
# # surr_glom_1_ca4to7 <- readRDS("surr_glom_1_ca4to7")
# # surr_glom_1_ca6 <- readRDS("surr_glom_1_ca6")
# # surr_glom_1_ca6lh <- readRDS("surr_glom_1_ca6lh")
# # 
# # surr_glom_1_ca5.8to6.2 <- list()
# # surr_glom_1_ca5.8to6.2[[1]] <- surr_glom_1_ca6lh[[1]]
# # surr_glom_1_ca5.8to6.2[[2]] <- surr_glom_1_ca6[[1]]
# # surr_glom_1_ca5.8to6.2[[3]] <- surr_glom_1_ca6lh[[2]]
# # 
# surr_glom_1_df <- data.frame(P.in = c(),
#                            D.in = c(),
#                            Ca.in = c(),
#                            Pavg = c(),
#                            SNGFR = c(),
#                            Q = c(),
#                            mu = c(),
#                            shear = c())
# 
# P_aa_input <- seq(from=80,to=180,by=10)  #Takenaka control: 100
# D_aa_input <- seq(from = 2, to = 30, by = 1)        #Takenaka control: 19
# Ca_input <- Ca_i#seq(from = 4.6, to = 6.6, by = 0.2)  #Takenaka control: 100
# 
# 
# # 
# # #Inputs according to low D, normal Ca
# # P_aa_input <- c(seq(from=80,to=180,by=5))#c(80,100,120,140,160,180)  #Takenaka control: 100
# # D_aa_input <- seq(from = 5, to = 15, by = 1)#c(10,15,20,25,30)        #Takenaka control: 19
# # Ca_input <- seq(from=5.8,to=6.2,by=0.2)#c(80,100,120,140,160,180)  #Takenaka control: 100
# # 
# # 
# for (i in seq(surr_glom_df)){
#   tot_out_i <- surr_glom_df[[i]]
#   surr_glom_temp <- data.frame(P.in = P_aa_input[tot_out_i$Pa.Var1],
#                                D.in = D_aa_input[tot_out_i$Pa.Var2],
#                                Ca.in = rep(Ca_input[i],length(tot_out_i$Pa.Var1)),
#                                Pavg = tot_out_i$Pa.value,
#                                SNGFR = tot_out_i$SNGFR.value,
#                                Q = tot_out_i$Q.value)
# 
#   surr_glom_1_df <- rbind(surr_glom_1_df,surr_glom_temp)
# 
# }
# 
# saveRDS(surr_glom_1_df, "surr_glom_df_20211108.rds")
# 
# #####################################################################
# 
# rm(list=ls())
# source("libs.R")
# source("parms.R")
# source("funcs_3_glomSS.R")
# 
# surr_glom_1_df <- readRDS("surr_glom_df_now.rds")
# 
# G <- glom_SS(Pa = 100,
#             Ca = Ca_i,
#             D = DAA_0,
#             surr_glom_df = surr_glom_1_df)
# G
# 
# # 
# # #Inputs according to high Ca vals
# # P_aa_input <- 100#c(seq(from=80,to=180,by=5))#c(80,100,120,140,160,180)  #Takenaka control: 100
# # D_aa_input <- seq(from = 5, to = 15, by = 1)#c(10,15,20,25,30)        #Takenaka control: 19
# # Ca_input <- seq(from=3.8,to=7,by=0.2)#c(80,100,120,140,160,180)  #Takenaka control: 100
# # 
# # for (i in seq(surr_glom_1_ca4to7)){
# #   tot_out_i <- surr_glom_1_ca4to7[[i]]
# #   surr_glom_temp <- data.frame(P.in = P_aa_input[tot_out_i$Pa.Var1], 
# #                                     D.in = D_aa_input[tot_out_i$Pa.Var2],
# #                                     Ca.in = rep(Ca_input[i],length(tot_out_i$Pa.Var1)),
# #                                     Pavg = tot_out_i$Pa.value,
# #                                     SNGFR = tot_out_i$SNGFR.value,
# #                                     Q = tot_out_i$Q.value,
# #                                     mu = tot_out_i$mu.value,
# #                                     shear = tot_out_i$shear.value)
# #   
# #   surr_glom_1_df <- rbind(surr_glom_1_df,surr_glom_temp)
# #   
# # }
# # 
# # saveRDS(surr_glom_1_df,"surr_glom_1_df")
# # 
# # 
# # ############################################################
# # #Testing
# # 
# # DAA_m <- DAA_0 
# # source("prep_anat.R")
# # OP_cont   <-      run_glom(src = src,
# #                            trg = trg,
# #                            D = D,
# #                            L = L,
# #                            k = k*src/src,
# #                            Paa_i = Paa_i,
# #                            Pea_o = Pea_o,
# #                            P_bs = P_bs,
# #                            Ca_i = Ca_i,
# #                            H_t_sys = H_t_sys,
# #                            mu_plas = mu_plas,
# #                            Rra=6,
# #                            t = t,
# #                            hpod=hpod,
# #                            wpod=wpod,
# #                            in.src = in.nodes,
# #                            out.trg = out.nodes,
# #                            num.iter = 150,
#                            mu.tol = 1e-3,
#                            Rinf.tol = 1e-3,
#                            beta=beta)
# 
# OP_cont$SNGFR
# OP_cont$G$Q[1]
# 
# 
# DAA_m <- 15
# source("prep_anat.R")
# OP_cont   <-      run_glom(src = src,
#                            trg = trg,
#                            D = D,
#                            L = L,
#                            k = k*src/src,
#                            Paa_i = 80,
#                            Pea_o = Pea_o,
#                            P_bs = P_bs,
#                            Ca_i = Ca_i,
#                            H_t_sys = H_t_sys,
#                            mu_plas = mu_plas,
#                            Rra=15.5,
#                            t = t,
#                            hpod=hpod,
#                            wpod=wpod,
#                            in.src = in.nodes,
#                            out.trg = out.nodes,
#                            num.iter = 150,
#                            mu.tol = 1e-5,
#                            Rinf.tol = 1e-5,
#                            beta=beta)
# 
# OP_cont$SNGFR
# mean(OP_cont$G$Pm[!src%in%in.nodes & !trg%in%out.nodes])
# OP_cont$G$Q[src%in%in.nodes]
# OP_cont$G$mu[1]
# 
# 
# L_AA <- 106.6698
# DAA_m <- 11*13/19
# source("prep_anat.R")
# OP_cont   <-      run_glom(src = src,
#                            trg = trg,
#                            D = D,
#                            L = L,
#                            k = k*src/src,
#                            Paa_i = 160,
#                            Pea_o = Pea_o,
#                            P_bs = P_bs,
#                            Ca_i = Ca_i,
#                            H_t_sys = H_t_sys,
#                            mu_plas = mu_plas,
#                            Rra=15.5,
#                            t = t,
#                            hpod=hpod,
#                            wpod=wpod,
#                            in.src = in.nodes,
#                            out.trg = out.nodes,
#                            num.iter = 150,
#                            mu.tol = 1e-5,
#                            Rinf.tol = 1e-5,
#                            beta=beta)
# 
# OP_cont$SNGFR
# mean(OP_cont$G$Pm[!src%in%in.nodes & !trg%in%out.nodes])
# OP_cont$G$Q[src%in%in.nodes]
# OP_cont$G$mu[1]
# 
# # 128*OP_cont$G$mu[1]*L_AA/(pi*DAA_m^4)*1e3/133.32239/60
# # 5/60
# # 
# # L_AA <- 3369.056
# # DAA_m <- 15
# # source("prep_anat.R")
# # OP_hp   <-      run_glom(src = src,
# #                            trg = trg,
# #                            D = D,
# #                            L = L,
# #                            k = k*src/src,
# #                            Paa_i = 160,
# #                            Pea_o = Pea_o,
# #                            P_bs = P_bs,
# #                            Ca_i = Ca_i,
# #                            H_t_sys = H_t_sys,
# #                            mu_plas = mu_plas,
# #                            Rra=6,
# #                            t = t,
# #                            hpod=hpod,
# #                            wpod=wpod,
# #                            in.src = in.nodes,
# #                            out.trg = out.nodes,
# #                            num.iter = 150,
# #                            mu.tol = 1e-3,
# #                            Rinf.tol = 1e-3,
# #                            beta=beta)
# # OP_hp$G$Q[1]
# # 128*OP_hp$G$mu[1]*L_AA/(pi*15^4)*1e3/133.32239/60
# # 60/(128*140/pi)/(OP_hp$G$mu[1]/15^4-OP_cont$G$mu[1]/DAA_0^4)/(1e3/133.32239/60)
# # 
# # 
# # 
# 
# 
# 

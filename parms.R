#Parameters 
#Created March 12, 2020
#Last updated March 17, 2020

#Time parameters
Ttot <- 60*2 #s
Nt <- 1e3*Ttot #steps
dt <- Ttot/Nt
time <- seq(from=dt,to=Ttot,by=dt)
time_in <- seq(time)
samp <- 50
samp_vec <- seq(from=1,to=Nt,by=samp)
Pa_swing <- 20
myo_del <- dt #s
myo_del_ind <- myo_del/dt

#Pressure boundary conditions
Paa_i <- 100
Pea_o <- 15
P_bs <- 13# 12.5
P_ext <- 5

#Afferent and Efferent Arteriole parameters
DAA_0 <- 7
DAA_pass <- 22 #feng2003ttype
L_AA <- 106.6698 #3369.056#4575#300
DEA_0 <- 7.3
L_EA <- 106.6698#430
Dmax <- 20
Rra <-  1.09#15.5 #23 #nl*mmHg/s 
tau_c <- 925*2/2
r_dilt <- 21.6/110.8

# C_int <- 0.002028783
# C_tone2 <- 0.7014563
# C_TGF <- 0.3308623
# C_tone3 <- 0.5276675
# C_tone4 <- 0.2143185
# C_myo_un <- C_myo - 0.004102097
# C_TGF_un <- 0.980031

#Glomerular parameters
k <- 3e-5#1.8028e-05#2.052441e-05#1.699908e-05#1/3*4e-5 
Ca_i <- 5.93868#5.7
ym_glom <- 2.5 # - 4 kPa
t <- 138 + 40 #nm
hpod <- 300
wpod <- 0.17
cmp_0 <- 0.02081642 #6.15114e-4
YM <- 14.36735#14.61202#13.836 #10.46967 #5.693673
YMn <- 15.8405#15.16006#5.102345 #2.788068#2.587312
V_0 <- 512371.7
V_0n <- 559191.6 #525199.7#

#Tubule parameters
Vmax <- 16.6 * 1e-6 *1e3   # nmol cm^-2 s^-1 --> mmol cm^-2 s^-1 -->  mM cm s^-1    #Layton, 1991
Km <- 70        #mM
r_Tub <- 10e-4     #cm
L_Tub <- 0.5   #cm
P_Tub <- 1.5e-5 # cm/s

#PT-DHL params
Co = 0.275e3                      #mosmol
C_ref <- 0.3e3
r = 1.25e-3                   #cm
A_M <- pi*r^2                 #cm^2
L_PT <- 0.5
L_DHL <- 0.4
alpha_PT <- 0.793362
RTLp <- alpha_PT*3.7e-8                #ml/s/cm^2/mosmol
sigma <- 0.7
h <- 2.4e-4                   #cm/s
N <- 1.1e-7                   #mosmol/s/cm
Jo <- alpha_PT*3.8e-7                  #ml/s/cm
Jso <- 1.8e-5                 #mosmol/s/cm
gamma <- 0.8
vo <- 30/A_M/1e6/60           #cm/s
ko <- 1.48                    #s^-1
k1 <- 0.05                    #s^-2
k2 <- 1.56                    #s^-1
x <- seq(0,1.4,1.4/1000)
Pv_PT <- 1.5e-5
Pv_DHL <- 1.5e-5
Vmax_PT <-  28 * 1e-6 * 1e3  #nmol cm^-2 s^-1 --> mmol cm^-2 s^-1 --> mM cm s^-1 
Vmax_DHL <-  0 * 1e-6 * 1e3 
Ce_L_Tub <- 0.150e3    #mosmol
A3_Tub <- 2
Ce_L_PT <- 0.65e3
A3_PT <-2


# P_Tub <- 1.5e-5 # converted from cm/s
# PT_let <- 1/3
# PT_let_ACZ <- 1 - 0.2318841#0.0856567#0.66*0.23/0.4
# PT_frac <- 1/5
# PT_frac_ACZ <- PT_let_ACZ*3/5
# C_0_Tub <- 275*PT_let/PT_frac      #mM
# C_0_Tub_ACZ <- 275*PT_let_ACZ/PT_frac_ACZ      #mM


#Systemic parameters
mu_plas <- 1.24 #cP
H_t_sys <- 0.40
C_salt <- Co

#Calibration parameters
num.iter.g <- 60
SNGFR.tol <- 1e-9
Rinf.tol <- 1e-9
k.tol <- 1e-12
num.iter.k <- 15
beta <- 3 






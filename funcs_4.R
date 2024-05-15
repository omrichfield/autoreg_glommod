#Pressure Resistance Network Glomerulus Functions
#Created February 2, 2020
#Last updated March 13, 2020

#Function to calculate resistances of filtering vessels using Starling's Law of Filtration
#Assumes: Linear profile of plasma protein concentration on the length of the vessel
get_Rinf1 <- function(src,
                      trg,
                      k,    #hydraulic conductivity
                      D,    #vessel diameters
                      L,    #vessel lengths
                      Ca,   #concentration of plasma protein at the source
                      Ce,   #concentration of plasma protein at the target
                      Pa,   #pressure at the source
                      Pe,   #pressure at the target
                      Pbs,   #Bowman's Space pressure
                      Qa,
                      Rinf_prev,
                      E,
                      mu,
                      Rra,
                      in.src,
                      out.trg,
                      num_quad=1000){

  R <- 128*mu*L/(pi*D^4)*1e3/133.32239/60
  R[1] <- R[1] + Rra/60
  a <- sqrt(R/(L^2*Rinf_prev))
  
  #p_diff_int <- (cosh(a*L)-1)/(a*sinh(a*L))*((Pa+Pe)/2 - Pbs)/L 
  p_diff_int <- tanh(a*L/2)/(a*L/2)*((Pa+Pe)/2 - Pbs)
  x <- matrix(rep(c(0:num_quad), length(L)), ncol=(num_quad+1), nrow=length(L), byrow=T)
  ax <- sweep(x,MARGIN = 1,a*L/num_quad,'*')
  
  pi_term <- sinh(ax) - sweep(cosh(ax), MARGIN = 1, cosh(a*L)/sinh(a*L),'*') 
  pj_term <- sweep(cosh(ax), MARGIN=1, 1/sinh(a*L),'*')
  pbs_term <- sinh(ax) + sweep(cosh(ax), MARGIN=1, (1-cosh(a*L))/sinh(a*L),'*')
  
  dpdx_x <- sweep(pi_term,MARGIN = 1,a*Pa,'*') + 
            sweep(pj_term,MARGIN = 1,a*Pe,'*') - 
            sweep(pbs_term,MARGIN = 1,a*Pbs,'*')  

  Q_x <- sweep(dpdx_x,MARGIN = 1,-L/R,'*')
  
  Q0 <- -L/R*a/sinh(a*L)*(-cosh(a*L)*Pa + Pe - (1-cosh(a*L))*Pbs)
  Q_ratio <- sweep(1/(Q_x-E),MARGIN = 1,(Q0-E),'*')

  C_x <- sweep(Q_ratio,MARGIN = 1,Ca,'*')
  
  a1 <- 2.1#1.629
  a2 <- 0.16#0.2935
  a3 <- 0.009
  Pi_x <- a1*C_x + a2*C_x^2 + a3*C_x^3
  Pi_int <- 1/(num_quad+1)*(0.5*(Pi_x[,1] + Pi_x[,(num_quad+1)]) + rowSums(Pi_x[,2:(num_quad)]))
  Rinf <- 1/(k*L*pi*D*(1 - Pi_int/p_diff_int))

  return(list(Rinf=Rinf,
              Pi_x=Pi_x,
              Pi_int=Pi_int,
              p_int=p_diff_int,
              x=x,
              Ce=C_x[,num_quad+1],
              Q_x=Q_x,
              Q0=Q0,
              C_x=C_x,
              k=k))
}
###############################################################################
#Calculate pressures at each node WITHOUT needing an explicit middle node
#Assumes: One input and one output
#Created: April 1, 2020
get_P_filt <- function(src,             #source nodes
                       trg,             #target nodes
                       L,               #edge lengths
                       D,               #edge diameters
                       mu,              #apparent viscosities
                       Rra,             #additional preafferent resistance
                       Rinf,            #resistance of filtering vessel
                       Pa,              #inlet pressure
                       Pe,              #outlet pressure
                       P_bs,            #Bowman's space pressure
                       in.src=1,        #in node
                       out.trg=max(trg) #out node
){
  
  #Determine unknown nodes and length of that list
  u.nodes <- sort(setdiff(c(src,trg),c(in.src,out.trg)))
  lu <- length(u.nodes)
  node.hash <- hashmap(u.nodes,seq(u.nodes))
  
  #Resistance vectors
  R <- 128*mu*L/(pi*D^4)*1e3/133.32239/60 #will give flow in nl/min
  R[which(src%in%in.src)] <- R[which(src%in%in.src)] + Rra/60 
  
  #Matrix Setup
  A <- matrix(0,nrow=lu, ncol=lu)
  b <- 0*u.nodes
  bm <- 0*u.nodes
  
  #Matrix Build
  for (i in seq(u.nodes)){
    ni <- u.nodes[i]                                    #each unknown node, corresponding to matrix rows
    j <- which(src==ni|trg==ni)                         #index corresponding to edges of which this node is a part
    
    i.trg <- which(trg==ni & !src%in%c(in.src,out.trg))   #which nodes point to ni, NOT including boundary
    i.src <- which(src==ni & !trg%in%c(in.src,out.trg))   #which nodes ni points to, NOT including boundary
    i.trg.bd <- which(trg==ni)                            #which nodes point to ni, including the boundary
    i.src.bd <- which(src==ni)                            #which nodes ni points to, including the boundary
    
    i.bd <- c(i.src.bd,i.trg.bd)
    i.nbd <- c(i.src,i.trg)
    
    c.n <- c(trg[i.src],src[i.trg])                         #get the name of the node that ni points to
    
    a.bd <- sqrt(R[i.bd]/(L[i.bd]^2*Rinf[i.bd]))             #calculate "a" for when ni is the src
    a.nbd <- sqrt(R[i.nbd]/(L[i.nbd]^2*Rinf[i.nbd]))             #calculate "a" for when ni is the src

    aL.bd <- a.bd*L[i.bd]
    aL.nbd <- a.nbd*L[i.nbd]
    
    #coeff <- -tanh(a.nbd*L[i.nbd]/2)/(a.nbd*L[i.nbd]/2)/R[i.nbd]   #calculate the coefficient for the connecting trg
    coeff <- a.nbd*L[i.nbd]/2/R[i.nbd]*(tanh(a.nbd*L[i.nbd]/2)^2 - 1)/tanh(a.nbd*L[i.nbd]/2)
    #coeff <- (tanh(aL.nbd/2)^2 - R[i.nbd]/cosh(aL.nbd/2) + R[i.nbd])/(2*aL.nbd*Rinf[i.nbd]*tanh(aL.nbd/2))
    
    #coeff.diag <- (2 -tanh(a.bd*L[i.bd]/2)/(a.bd*L[i.bd]/2))/R[i.bd] #coefficient for ni, from src
    coeff.diag <- a.bd*L[i.bd]/2/R[i.bd]*(tanh(a.bd*L[i.bd]/2)^2 + 1)/tanh(a.bd*L[i.bd]/2)
    #coeff.diag <- (tanh(aL.bd/2)^2 + R[i.bd]/cosh(aL.bd/2) - R[i.bd])/(2*aL.bd*Rinf[i.bd]*tanh(aL.bd/2))
    
    #coeff.bs <- (1 - tanh(a.bd*L[i.bd]/2)/(a.bd*L[i.bd]/2))*2/R[i.bd] #Bowman's space pressure coefficients
    coeff.bs <- tanh(aL.bd/2)/(Rinf[i.bd]*aL.bd)
    
    agg.coeff <- aggregate(coeff, list(c.n), FUN=sum)
    jc.s <- node.hash$find(agg.coeff$Group.1)               #get indices in A corresponding to those nodes
    A[i,jc.s] <- agg.coeff$x
    
    A[i,i] <- sum(coeff.diag)          #Diagonal entries
    b[i] <- P_bs*sum(coeff.bs)           #update vector b
  }

  #Update to take into account boundary condition
  i.in.s <- which(src==in.src) #indices for nodes connected to the in-node
  i.in.t <- which(trg==in.src)
  i.out.s <- which(src==out.trg) #indices for nodes connected to the out-node
  i.out.t <- which(trg==out.trg)
  
  in.conns <- c(trg[i.in.s],src[i.in.t]) #get the nodes associated with these connections
  out.conns <- c(trg[i.out.s],src[i.out.t])
  
  in.ind <- node.hash[[in.conns]] #get indices for those nodes in terms of unknown nodes
  out.ind <- node.hash[[out.conns]]
  
  i.in <- c(i.in.s,i.in.t)
  i.out <- c(i.out.s,i.out.t)
  
  a.in <- sqrt(R[i.in]/(L[i.in]^2*Rinf[i.in]))
  a.out <- sqrt(R[i.out]/(L[i.out]^2*Rinf[i.out]))
  aL.in <- a.in*L[i.in]
  aL.out <- a.out*L[i.out]
  
  #Update temp vector bm for input (Pa) and output (Pe) pressures
  #bm[in.ind] <- bm[in.ind] + Pa*tanh(a.in*L[i.in]/2)/(a.in*L[i.in]/2)/R[i.in]
  #bm[out.ind] <- bm[out.ind] + Pe*tanh(a.out*L[i.out]/2)/(a.out*L[i.out]/2)/R[i.out]
  bm[in.ind] <- bm[in.ind] - Pa*a.in*L[i.in]/2/R[i.in]*(tanh(a.in*L[i.in]/2)^2 - 1)/tanh(a.in*L[i.in]/2)
  bm[out.ind] <- bm[out.ind] - Pe*a.out*L[i.out]/2/R[i.out]*(tanh(a.out*L[i.out]/2)^2 - 1)/tanh(a.out*L[i.out]/2)
  # bm[in.ind] <- bm[in.ind] - Pa*(tanh(aL.in/2)^2 - R[i.in]/cosh(aL.in/2) + R[i.in])/
  #                                         (2*aL.in*Rinf[i.in]*tanh(aL.in/2))
  # bm[out.ind] <- bm[out.ind] - Pe*(tanh(aL.out/2)^2 - R[i.out]/cosh(aL.out/2) + R[i.out])/
  #                                         (2*aL.out*Rinf[i.out]*tanh(aL.out/2))

  bb <- b + bm #update b to take into account boundary condition
  p <- solve(A)%*%bb
  return(list(p=p,b=bb,A=A,n=u.nodes))
}
#########################################################################################
#Function to calculate erythrocyte flow through each segment given blood flow
#Created April 4, 2020
#Assumes tube hematocrit BEFORE filtered volume, discharge hematocrit AFTER filtered volume
get_C_filt <- function(src,
                       trg,
                       Ca_i,
                       Qa,
                       Qe,
                       in.src,
                       out.trg){
  Ca <- NaN*src
  Ce <- NaN*src
  
  in.indx <- which(src%in%in.src)
  Ca[in.indx] <- Ca_i
  
  input <- data.frame(sn = src, tn = trg)
  g <- graph_from_data_frame(input, directed=TRUE)
  node.list <- as.numeric(names(topo_sort(graph=g,mode="out")))
  
  for (ni in node.list){
    i.trg <- which(trg==ni)
    i.src <- which(src==ni)
    
    if (length(i.trg) & length(i.src)){
      Q.parents <- Qe[i.trg]
      C.parents <- Ce[i.trg]
      Ca[i.src] <- sum(C.parents*Q.parents)/sum(Qa[i.src])
      Ce[i.src] <- Ca[i.src]*Qa[i.src]/Qe[i.src]
    }
    if (!length(i.trg)){
      Ce[i.src] <- Ca[i.src]*Qa[i.src]/Qe[i.src]
    }
  }
  return(list(Ca=Ca,Ce=Ce))
}

#########################################################################################
#Function to solve for pressures and get concentrations for a set of Rinfs
#Created April 4, 2020

run_net_one <- function(src,
                        trg,
                        D,
                        L,
                        k,
                        Ca_i,
                        Paa_i,
                        Pea_o,
                        E,
                        mu_app,
                        Rra,
                        Rinf,
                        P_bs,
                        t,
                        hpod,
                        wpod,
                        in.src,
                        out.trg){
  
  P_filt <- get_P_filt(src=src,             #source nodes
                       trg=trg,             #target nodes
                       L=L,               #edge lengths
                       D=D,               #edge diameters
                       mu=mu_app,              #apparent viscosities
                       Rra=Rra,
                       Rinf=Rinf,            #resistance of filtering vessel
                       Pa=Paa_i,              #inlet pressure
                       Pe=Pea_o,              #outlet pressure
                       P_bs=P_bs,            #Bowman's space pressure
                       in.src=in.src,        #in node
                       out.trg=out.trg #out node
  )

  p.hash <- hashmap(c(in.nodes,P_filt$n,out.nodes),c(Paa_i,P_filt$p,Pea_o))
  p.src <- p.hash[[src]]
  p.trg <- p.hash[[trg]]
  R <- 128*L*mu_app/(pi*D^4)*1e3/133.32239/60 #will give flow in nl/min
  R[which(src%in%in.src)] <- R[which(src%in%in.src)] + Rra/60 
  Rinf <- Rinf
  a <- sqrt(R/(L^2*Rinf))
  
  src1 <- src
  trg1 <- trg
  src1[which(p.src < p.trg)] <- trg[which(p.src < p.trg)]
  trg1[which(p.src < p.trg)] <- src[which(p.src < p.trg)]
  src <- src1
  trg <- trg1
  
  p.src <- p.hash[[src]]
  p.trg <- p.hash[[trg]]
  p.mid <- P_bs + ((p.src+p.trg)/2-P_bs)/cosh(a*L/2)#
  p.avg <- P_bs + tanh(a*L/2)/(a*L/2)*((p.src+p.trg)/2 - P_bs)
  Q <- (p.src - p.trg)/R
  
  Q.src.a <- -L/R*a/sinh(a*L)*(-cosh(a*L)*p.src + p.trg - P_bs*(1-cosh(a*L)))
  Q.trg.a <- -L/R*a/sinh(a*L)*(-p.src + cosh(a*L)*p.trg + P_bs*(1-cosh(a*L)))
  Q.src <- 2*(p.src-p.avg)/R
  Q.trg <- 2*(p.avg-p.trg)/R
  
  #Qfilt <- (p.avg-P_bs)/(Rinf)
  Qfilt <- (p.avg-P_bs)/Rinf
  
  C_filt <- get_C_filt(src=src,
                       trg=trg,
                       Ca_i=Ca_i,
                       Qa=Q.src.a-E,
                       Qe=Q.trg.a-E,
                       in.src=in.src,
                       out.trg=out.trg
  )
  
  #c.hash <- hashmap(c(in.src, C_filt$n, out.trg), c(Ca_i,C_filt$c,C_filt$Ce))
  c.src <- C_filt$Ca #c.hash[[src]]
  c.trg <- C_filt$Ce #c.hash[[trg]]
  
  t_new <- t + hpod/2 #+ hpod/2*wpod/(pi^2*D*2)*sin(2*pi/wpod*(pi*D-wpod/2))
  #t_new <- t;
  
  shear <- 32*mu_app*Q/(pi*D^3)*1e3/6 #dyne/cm^2
  hoop <- (p.avg-P_bs)*D/(2*t_new)*133.32239 #kPa
  
  output <- data.frame(src=src,
                       trg=trg,
                       D=D,
                       L=L,
                       t=t_new,
                       mu=mu_app,
                       Rinf=Rinf,
                       k=k,
                       Pbs=rep(P_bs,length(src)),
                       Q=Q,
                       Qsm=Q.src,
                       Qmt=Q.trg,                       
                       Qsm.a=Q.src.a,
                       Qmt.a=Q.trg.a,
                       CSGFR=Qfilt,
                       Pa=p.src,
                       Pe=p.trg,
                       Pm=p.avg,
                       Ca=c.src,
                       Ce=c.trg,
                       shear=shear,
                       hoop=hoop)
  return(output)
}


#########################################################################################
#Function to calculate hematocrit and discharge hematocrit from erythrocyte volume and blood flow
#Created April 4, 2020

hemato_helper <- function(E,
                          QB,
                          D){
  H_t <- E/QB
  H_d <- H_t/(0.5*(1+exp(-0.633*(D/10.43-1))))
  return(list(H_t=H_t,H_d=H_d))
}

#########################################################################################
#Function to calculate erythrocyte flow through each segment given blood flow
#Created April 4, 2020
#Assumes tube hematocrit BEFORE filtered volume, discharge hematocrit AFTER filtered volume
get_E_filt <- function(src,
                       trg,
                       D,
                       Ea,
                       Qa,
                       Qe,
                       in.src,
                       out.trg){
  E <- NaN*src                                         #Define initial arrays to populate - erythrocytes
  H_d <- 0*src                                         #discharge hematocrit
  H_t <- 0*src                                         #tube hematocrit
  mu_app <- 0*src                                      #apparent viscosity
  
  in.indx <- which(src%in%in.src)
  E[in.indx] <- Ea
  Hx <- hemato_helper(E = E, QB = (Qa + Qe)/2, D = D)
  H_t <- Hx$H_t
  H_d <- Hx$H_d
  
  input <- data.frame(sn = src, tn = trg)
  g <- graph_from_data_frame(input, directed=TRUE)
  node.list <- as.numeric(names(topo_sort(graph=g,mode="out")))
  
  for (ni in node.list){
    i.trg <- which(trg==ni)
    i.src <- which(src==ni)
    
    if (length(i.trg)==0){
      Esum <- Ea
    }else{
      Esum <- sum(E[i.trg])
    }
    
    H_d.parent <- sum(H_d[i.trg]*Qe[i.trg])/sum(Qe[i.trg])
    D.parent <- sqrt(sum(D[i.trg]^2))
    
    if (length(i.src) == 1){
      E[i.src] <- Esum
    }
    
    if (length(i.src) == 2){
      X0 <- 0.4/D.parent
      FQ <- Qa[i.src]/sum(Qa[i.src])
      xt <- (FQ - X0)/(1 - 2*X0)
      
      if (any(xt < 0)){
        xt1 <- 0*xt
        xt1[which(xt > 0)] <- 1
        xt <- xt1
      }
      
      A <- c(0,0)
      A[1] <- -6.96*log(D[i.src[1]]/D[i.src[2]])/D.parent
      A[2] <- -6.96*log(D[i.src[2]]/D[i.src[1]])/D.parent
      B <- 1 + 6.98*(1-H_d.parent)/D.parent
      E[i.src] <- Esum/(1 + exp(A - B*log(xt/(1-xt))))
    }
    
    
    if (length(i.src) > 2){
      li <- length(i.src)-1
      i.s.t <- i.src
      for (i in seq(li)){
        D1 <- D[i.s.t[1]]
        D2 <- sqrt(sum(D[i.s.t[2:length(i.s.t)]]^2))
      
        FQ<-c(0,0)
        FQ[1] <- Qa[i.s.t[1]]/sum(Qa[i.s.t])
        FQ[2] <- sum(Qa[i.s.t[2:length(i.s.t)]])/sum(Qa[i.s.t])
      
        X0 <- 0.4/D.parent
        xt <- (FQ - X0)/(1 - 2*X0)

        if (any(xt < 0)){
          xt1 <- 0*xt
          xt1[which(xt > 0)] <- 1
          xt <- xt1
        }
      
        A <- c(0,0)
        A[1] <- -6.96*log(D1/D2)/D.parent
        A[2] <- -6.96*log(D2/D1)/D.parent
        B <- 1 + 6.98*(1-H_d.parent)/D.parent
        E.temp <- Esum/(1 + exp(A - B*log(xt/(1-xt))))
      
        E[i.s.t[1]] <- E.temp[1]
        Esum <- E.temp[2]
        i.s.t <- i.s.t[2:length(i.s.t)]
      }
      E[i.src[li+1]] <- E.temp[2]
      # Q0 <- sum(Qe[i.trg])
      # #B <- 2.016 - 2.18*H_d.parent
      # #H_t[i.src] <- (Q0/Qa[i.src])*(Qa[i.src]/Q0)^B/((1-Qa[i.src]/Q0)^B + (Qa[i.src]/Q0)^B)
      # E[i.src] <- Esum*Qa[i.src]/sum(Qa[i.src]) #H_t[i.src]*Qa[i.src] #
    }
    
    Hx <- hemato_helper(E = E, QB = (Qa + Qe)/2, D = D)
    H_t <- Hx$H_t
    H_d <- Hx$H_d
  }
  
  return(list(E = E, H_t = H_t, H_d = H_d))
}
#########################################################################################
#Function to calculate apparent viscosity based on discharge hematocrit in each vessel
#Created April 8, 2020
get_mu_app <- function(H_d,
                       D, 
                       mu_plas)
{
  
  H_d[which(H_d > 1)]<-0.99
  eta045 <- 220*exp(-1.3*D) + 3.2 - 2.44*exp(-0.06*D^0.645)
  gamma <- (0.8+exp(-0.075*D))*(-1 + 1/(1+D^12/1e11)) + 1/(1+D^12/1e11)
  xi <- 1 + (eta045-1)*((1-H_d)^gamma-1)/(0.55^gamma-1)
  mu_app <- xi
  return(mu_app)
}

#########################################################################################
#Function to compute convergent glomerulus function
#Created April 8, 2020

run_glom <- function(src,
                     trg,
                     D,
                     L,
                     k,
                     Rinf_init,
                     Paa_i,
                     Pea_o,
                     P_bs,
                     Ca_i,
                     H_t_sys,
                     mu_plas,
                     Rra,
                     t,
                     hpod,
                     wpod,
                     in.src,
                     out.trg,
                     num.iter=30,
                     mu.tol=1e-5,
                     Rinf.tol=1e-5,
                     beta=4){
  count <- 0
  SNGFR0 <- 0
  Rinf0 <- Rinf_init
  mu0 <- mu_plas
  Pa0 <- 100*src/src
    
  mu.err <- mu.tol*5
  Rinf.err <- Rinf.tol*5
  
  Q.track <- 0*seq(num.iter)
  SNGFR.err.track <- 0*seq(num.iter)
  SNGFR.track <- 0*seq(num.iter)
  Rinf.track <- 0*seq(num.iter)
  mu.track <- 0*seq(num.iter)
  C.track.each <- matrix(0,length(src),num.iter)
  P.track.each <- matrix(0,length(src),num.iter)
  Rinf.track.each <- matrix(0,length(src),num.iter)
  mu.track.each <- matrix(0,length(src),num.iter)
  Rinf.track.ind <- 0*seq(num.iter)
  alpha.track <- 0*seq(num.iter)
  p.node.err.track <- 0*seq(num.iter)
  
  G <- run_net_one(src=src,
                    trg=trg,
                    D=D,
                    L=L,
                    k=k,
                    Ca_i=Ca_i,
                    Paa_i=Paa_i,
                    Pea_o=Pea_o,
                    E=0,
                    mu_app=mu_plas,
                    Rra=Rra,
                    Rinf=Rinf0,
                    P_bs=P_bs,
                    t=t,
                    hpod=hpod,
                    wpod=wpod,
                    in.src=in.src,
                    out.trg=out.trg)
  
  while ((mu.err > mu.tol | Rinf.err > Rinf.tol) &  count < num.iter ){
    Pa0 <- G$Pa
    
    H_filt <- get_E_filt( src=G$src,
                          trg=G$trg,
                          D=D,
                          Ea=H_t_sys*G$Q[which(src%in%in.nodes)],
                          Qa=G$Qsm.a,
                          Qe=G$Qmt.a,
                          in.src=in.src,
                          out.trg=out.trg)
    
    mu_plas_func <- 0.274 + 0.177*(G$Ca*G$Qsm.a/G$Q)
    mu_plas_func[which(mu_plas_func < 0)] <- mu_plas

    mu_app_inf <- get_mu_app(H_d=H_filt$H_d,
                         D=D,
                         mu_plas=mu_plas_func)

    Rinf <- get_Rinf1(src=G$src,
                      trg=G$trg,
                      k=k,    #hydraulic conductivity
                      D=D,    #vessel diameters
                      L=L,    #vessel lengths
                      Ca=G$Ca,   #concentration of plasma protein at the source
                      Ce=G$Ce,   #concentration of plasma protein at the target
                      Pa=G$Pa,   #pressure at the source
                      Pe=G$Pe,   #pressure at the target
                      Pbs=G$Pbs,   #Bowman's Space pressure
                      Qa=G$Qsm.a,
                      Rinf_prev=G$Rinf,
                      E=H_filt$E,
                      mu=mu0,#_app,
                      Rra=Rra,
                      in.src=in.src,
                      out.trg=out.trg,
                      num_quad=1000)
    
    Rinf1 <- Rinf$Rinf
    Rinf_prev <- G$Rinf
    
    i.b <- which((Rinf1 < 0 | Rinf$Ce < 0) & Rinf_prev < 1e10)
    i.c <- which((Rinf1 < 0 | Rinf$Ce < 0) & Rinf_prev >= 1e10)
    #print(i.b)
    #print(i.c)
    
    alpha_R <- abs(Rinf1[Rinf1>0]-Rinf_prev[Rinf1>0])/Rinf_prev[Rinf1>0]*beta
    alpha_mu <- abs(mu_app_inf - mu0)/mu0*beta
    alpha <- max(alpha_R, alpha_mu, 1)
    
    mu_app <- mu0 + (mu_app_inf - mu0)/alpha
    #print(mu_app[src%in%in.nodes])
    
    Rinf_new <- Rinf_prev + (Rinf1 - Rinf_prev)/alpha
    Rinf_new[i.b] <- (1+1/beta)*Rinf_prev[i.b] 
    Rinf_new[i.c] <- Rinf_prev[i.c]
    Rinf_new[which(src%in%in.src | trg%in%out.trg)] <- 1e20
    
    G <- run_net_one(src=src,
                     trg=trg,
                     D=D,
                     L=L,
                     k=k,
                     Ca_i=Ca_i,
                     Paa_i=Paa_i,
                     Pea_o=Pea_o,
                     E=H_filt$E,
                     mu_app=mu_app,
                     Rra=Rra,
                     Rinf=Rinf_new,
                     P_bs=P_bs,
                     t=t,
                     hpod=hpod,
                     wpod=wpod,
                     in.src=in.nodes,
                     out.trg=out.nodes)
    
    SNGFR <- sum(G$CSGFR)
    
    SNGFR.err <- abs(SNGFR - SNGFR0)/SNGFR0
    Rinf.err <- max(abs(Rinf_new-Rinf0)/Rinf0)
    mu.err <- max(abs(mu_app-mu0)/mu0)
    
    SNGFR0 <- SNGFR
    Rinf0 <- Rinf_new
    mu0 <- mu_app
    
    count <- count + 1
    if (!count%%25){
      print(paste("it: ", count))
      if (length(i.b)){
        print(paste("Pi_int > p_diff_int at ",which(Rinf$Pi_int > Rinf$p_int & Rinf_prev < 1e10)))
      }
    }
    
    C.track.each[,count] <- G$Ce
    P.track.each[,count] <- Rinf$p_int
    mu.track.each[,count] <- G$mu
    Rinf.track.each[,count] <- G$Rinf
    
    SNGFR.track[count]<-SNGFR
    SNGFR.err.track[count]<-SNGFR.err
    Rinf.track.ind[count] <- Rinf.track.each[121,count]
    Rinf.track[count]<-Rinf.err
    mu.track[count]<-mu.err
    Q.track[count] <- G$Q[1]
    alpha.track[count] <- alpha
    p.node.err.track[count] <- max(abs(G$Pa-Pa0)/Pa0)
  }
  
  H_filt1 <- get_E_filt( src=G$src,
                        trg=G$trg,
                        D=D,
                        Ea=H_t_sys*G$Q[which(src%in%in.nodes)],
                        Qa=G$Qsm.a,
                        Qe=G$Qmt.a,
                        in.src=in.src,
                        out.trg=out.trg)
  
  if (count >= num.iter){
    print(paste("Count exceeded max alloted. SNGFR.err = ", SNGFR.err))
  }else{
    print(paste("Converged in ", count, " iterations."))
    print(paste("Final SNGFR = ", SNGFR))
  }
  OP_prime <- data.frame(G, Pi=Rinf$Pi_int, E1 = H_filt1$E, E = H_filt$E, H_t = H_filt$H_t, H_d = H_filt$H_d)

  return(list(G=OP_prime,
              SNGFR.err=SNGFR.err.track[1:count],
              Rinf.err=Rinf.track[1:count],
              mu.err=mu.track[1:count],
              SNGFR=SNGFR,#.track[1:count],
              Q=Q.track[1:count],
              FF=SNGFR/(OP_prime$Q[1]-OP_prime$E[1]),
              num.iter=count,
              Rinf.tol=Rinf.tol,
              SNGFR.tol=SNGFR.tol,
              C.t=C.track.each[,1:count],
              P.t=P.track.each[,1:count],
              R.t=Rinf.track.each[,1:count],
              mu.t=mu.track.each[,1:count],
              alpha.t=alpha.track[1:count],
              p.err=p.node.err.track[1:count]))
}

#########################################################################
k_calibrate <- function(src,
                        trg,
                        D,
                        L,
                        k0,
                        k1,
                        SNGFR_ex,
                        Paa_i,
                        Pea_o,
                        P_bs,
                        Ca_i,
                        H_t_sys,
                        mu_plas,
                        Rra,
                        t,
                        hpod,
                        in.nodes,
                        out.nodes,
                        num.iter.g,
                        num.iter.k,
                        k.tol,
                        SNGFR.tol,
                        mu.tol,
                        Rinf.tol,
                        beta){

k.err <- k.tol*5
SNGFR.err <- SNGFR.tol*5
k_keep <- 0*seq(num.iter.k)
S_keep <- 0*seq(num.iter.k)
count <- 1

G0   <-      run_glom(src = src,
                      trg = trg,
                      D = D,
                      L = L,
                      k = k0*src/src,
                      Paa_i = Paa_i,
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
                      num.iter = num.iter.g,
                      mu.tol = mu.tol,
                      Rinf.tol = Rinf.tol,
                      beta=beta)

G1   <-      run_glom(src = src,
                      trg = trg,
                      D = D,
                      L = L,
                      k = k1*src/src,
                      Paa_i = Paa_i,
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
                      num.iter = num.iter.g,
                      mu.tol = mu.tol,
                      Rinf.tol = Rinf.tol,
                      beta=beta)

while((k.err > k.tol | SNGFR.err > SNGFR.tol) & count < num.iter.k){

  S0 <- G0$SNGFR
  S1 <- G1$SNGFR
  
  k_new <- (k0 + k1)/2
  
  G_new   <-      run_glom(src = src,
                          trg = trg,
                          D = D,
                          L = L,
                          k = k_new*src/src,
                          Paa_i = Paa_i,
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
                          num.iter = num.iter.g,
                          mu.tol = mu.tol,
                          Rinf.tol = Rinf.tol,
                          beta=beta)
    S_new <- G_new$SNGFR
    
    if (S_new < SNGFR_ex){
      k0 <- k_new
    }else{
      k1 <- k_new
    }
  
  k_keep[count] <- k_new
  S_keep[count] <- S_new
  
  count <- count+1
  k.err <- abs(k1 - k0)/k0
  SNGFR.err <- abs(S_new - SNGFR_ex)/SNGFR_ex
}

return(list(OP=G_new,k_keep=k_keep,S_keep=S_keep,count=count))
}

#########################################################################
compliant_glom <- function(src,
                          trg,
                          D,
                          L,
                          Rinf_init,
                          k,
                          E_cmp,
                          Paa_i,
                          Paa_try,
                          Pea_o,
                          P_bs,
                          Ca_i,
                          H_t_sys,
                          mu_plas,
                          Rra,
                          t,
                          hpod,
                          wpod,
                          in.nodes,
                          out.nodes,
                          num.iter.g,
                          num.iter.d,
                          d.tol,
                          t.tol,
                          L.tol,
                          SNGFR.tol,
                          mu.tol,
                          Rinf.tol,
                          beta){
  
  d.err <- d.tol*5
  L.err <- L.tol*5
  t.err <- t.tol*5
  hpod.err <- t.tol*5
  
  SNGFR.err <- SNGFR.tol*5
  d_keep <- matrix(0, ncol = num.iter.d, nrow = length(src))
  L_keep <- matrix(0, ncol = num.iter.d, nrow = length(src))
  t_keep <- matrix(0, ncol = num.iter.d, nrow = length(src))
  hpod_keep <- matrix(0, ncol = num.iter.d, nrow = length(src))
  wpod_keep <- matrix(0, ncol = num.iter.d, nrow = length(src))
  
  P_keep <- matrix(0, ncol = num.iter.d, nrow = length(src))
  S_keep <- 0*seq(num.iter.d)
  count <- 1

  D_temp <- D
  D_temp[src%in%in.nodes] <- DAA_0
  
  OP0   <-      run_glom(src = src,
                        trg = trg,
                        D = D_temp,
                        L = L,
                        Rinf_init = Rinf_init,
                        k = k*src/src,
                        Paa_i = Paa_i,
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
                        num.iter = num.iter.g,
                        mu.tol = mu.tol,
                        Rinf.tol = Rinf.tol,
                        beta=beta)
    P_1 <- OP0$G$Pm
    d_1 <- D
    t_1 <- t*src/src # + hpod/2 + hpod/2*wpod/(pi^2*D*2)*sin(2*pi/wpod*(pi*D-wpod/2))
    hpod_1 <- hpod*src/src
    wpod_1 <- wpod*src/src
    L_1 <- L
    
    L0 <- L 
    d0 <- D
    t0 <- t*src/src #OP0$G$t# + hpod/2 + hpod/2*wpod/(pi^2*D*2)*sin(2*pi/wpod*(pi*D-wpod/2))
    hpod0 <- hpod*src/src
    wpod0 <- wpod*src/src
    S_1 <- OP0$SNGFR
    
    d_keep[,1] <- d_1
    L_keep[,1] <- L_1
    t_keep[,1] <- t_1
    hpod_keep[,1] <- hpod_1
    wpod_keep[,1] <- wpod_1
    S_keep[1] <- S_1
    P_keep[,1] <- P_1   
    
while((d.err > d.tol | L.err > L.tol | t.err > t.tol | hpod.err > t.tol | SNGFR.err > SNGFR.tol) 
          & (count < num.iter.d)){
    
    OP_new   <-      run_glom(src = src,
                             trg = trg,
                             D = d0,
                             L = L0,
                             Rinf_init = Rinf_init,
                             k = k*src/src,
                             Paa_i = Paa_try,
                             Pea_o = Pea_o,
                             P_bs = P_bs,
                             Ca_i = Ca_i,
                             H_t_sys = H_t_sys,
                             mu_plas = mu_plas,
                             Rra=Rra,
                             t = t0,
                             hpod=hpod0,
                             wpod=wpod0,
                             in.src = in.nodes,
                             out.trg = out.nodes,
                             num.iter = num.iter.g,
                             mu.tol = mu.tol,
                             Rinf.tol = Rinf.tol,
                             beta=beta)
    P_new <- OP_new$G$Pm
    S_new <- OP_new$SNGFR
    d_new <- 0*d_1
    
    E_cmp_mmHg <- E_cmp/0.000133322
    t_um <- OP0$G$t/1000
    
    for (i in seq(d_1)){
      c_1 <- 1
      c_2 <- d_1[i]
      c_3 <- -4*E_cmp_mmHg*t_um[i]*d_1[i]/(3*(P_new[i]-P_bs))
      c_4 <- 4*E_cmp_mmHg*t_um[i]*d_1[i]^2/(3*(P_new[i]-P_bs)) - 2*d_1[i]^3*(P_1[i]-P_bs)/(P_new[i]-P_bs)
      
      d_cub <- cubic(c(c_1,c_2,c_3,c_4))
      d_new[i]     <- Re(d_cub[3])#Re(d_cub[which(Im(d_cub)==0)])
    }
    
    d_new[src%in%in.nodes | trg%in%out.nodes] <- d_1[src%in%in.nodes | trg%in%out.nodes]
    L_new     <- L_1*(3*d_1^2*(P_1-P_bs) - 4*E_cmp_mmHg*t_um*d_1)/(3*d_new^2*(P_new-P_bs) - 4*E_cmp_mmHg*t_um*d_1)
    t_new     <- t_1*L_1*d_1/(L_new*d_new)#(3*d_new^2*(P_new-P_bs) - 4*E_cmp_mmHg*t_um*d_1)/(3*d_new*d_1*(P_1-P_bs) - 4*E_cmp_mmHg*t_um*d_1)
    hpod_new  <- hpod_1*L_1*d_1/(L_new*d_new)#(3*d_new^2*(P_new-P_bs) - 4*E_cmp_mmHg*t_um*d_1)/(3*d_new*d_1*(P_1-P_bs) - 4*E_cmp_mmHg*t_um*d_1)
    wpod_new <- wpod_1*d_new/d_1
    
    d_new[src%in%in.nodes | trg%in%out.nodes] <- d_1[src%in%in.nodes | trg%in%out.nodes]
    L_new[src%in%in.nodes | trg%in%out.nodes] <- L_1[src%in%in.nodes | trg%in%out.nodes]
    t_new[src%in%in.nodes | trg%in%out.nodes] <- t_1[src%in%in.nodes | trg%in%out.nodes]
    hpod_new[src%in%in.nodes | trg%in%out.nodes] <- hpod_1[src%in%in.nodes | trg%in%out.nodes]
    wpod_new[src%in%in.nodes | trg%in%out.nodes] <- wpod_1[src%in%in.nodes | trg%in%out.nodes]
    
    d_keep[,count+1] <- d_new
    S_keep[count+1] <- S_new
    P_keep[,count+1] <- P_new
    L_keep[,count+1] <- L_new
    t_keep[,count+1] <- t_new
    hpod_keep[,count+1] <- hpod_new
    wpod_keep[,count+1] <- wpod_new
    
    count <- count+1
    d.err <- max(abs(d_new - d0)/d0)
    L.err <- max(abs(L_new - L0)/L0)
    t.err <- max(abs(t_new - t0)/t0)
    hpod.err <- max(abs(hpod_new - hpod0)/hpod0)
    SNGFR.err <- abs(S_new - S_1)/S_1
 
    P0 <- P_new
    d0 <- d_new
    t0 <- t_new
    hpod0 <- hpod_new
    wpod0 <- wpod_new
    L0 <- L_new
    S_1 <- S_new
  }
    OP_new   <-      run_glom(src = src,
                              trg = trg,
                              D = d0,
                              L = L0,
                              Rinf_init = Rinf_init,
                              k = k*src/src,
                              Paa_i = Paa_try,
                              Pea_o = Pea_o,
                              P_bs = P_bs,
                              Ca_i = Ca_i,
                              H_t_sys = H_t_sys,
                              mu_plas = mu_plas,
                              Rra=Rra,
                              t = t0,
                              hpod=hpod0,
                              wpod=wpod0,
                              in.src = in.nodes,
                              out.trg = out.nodes,
                              num.iter = num.iter.g,
                              mu.tol = mu.tol,
                              Rinf.tol = Rinf.tol,
                              beta=beta)
    
  return(list(OP=OP_new,
              d_keep=d_keep,
              P_keep=P_keep,
              S_keep=S_keep,
              count=count,
              L_keep=L_keep,
              t_keep=t_keep,
              wpod_keep=wpod_keep,
              hpod_keep=hpod_keep))
}

#########################################################################





check_takenaka <- function(D, 
                           Pa, 
                           Ca,
                           surr_glom_df,
                           PT_frac,
                           C_0_Tub,
                           furo_yn,
                           dilt_yn,
                           int_yn,
                           TGF_only_yn){

  GG <- glom_SS(Pa=Pa, Ca=Ca, D=D, surr_glom_df = surr_glom_df)
  SNGFR <- GG$SNGFR
  Pavg <- GG$Pavg
  Q <- GG$Q
  
  Tp <- D*(Pavg - P_ext)/2
  Tp0 <- Tp_cont#DAA_0*(glom_SS(Pa=Paa_i, Ca=Ca_i, D=DAA_0,surr_glom_df = surr_glom_df)$Pavg - P_ext)/2

  T_pass <- C_pass0*exp(C_pass1*(D/DAA_0-1))
  T_act <- C_act0*exp(-((D/DAA_0-1)/C_act1)^2)
  x <- seq(0, L_DHL+L_PT, by = (L_DHL+L_PT)/1000)
  
  parameters_PT <- c(Pv=Pv_PT,
                     A_M=A_M,Vmax=Vmax_PT,r=r)
  parameters_DHL <- c(Pv=Pv_DHL,
                      A_M=A_M,Vmax=Vmax_DHL,r=r)
  
  vo <- SNGFR/A_M/1e6/60
  x_PT <- seq(0,(L_PT),(L_PT)/1000)
  out_PT <- ode(y = c(Co,vo), 
                times = x_PT, 
                func = ss_Tubule_Richfield,
                parms = parameters_PT)
  x_DHL <- seq(L_PT,(L_PT+L_DHL),(L_DHL)/1000)
  out_DHL <- ode(y = out_PT[length(x_PT),2:3], 
                 times = x_DHL, 
                 func = ss_Tubule_Richfield,
                 parms = parameters_DHL)  
  
  x_AL <- seq(L_PT+L_DHL,(L_PT+L_DHL+L_Tub),(L_Tub)/1000)
  parameters_AL <- c(Loop_Flow = out_DHL[,3][length(out_DHL)/3]*A_M)
  out_AL <- ode(y = out_DHL[,2][length(out_DHL)/3], 
                times = x_AL, 
                func = ss_Tubule,
                parms = parameters_AL)
  C <- out_AL[,2]
  C_MD <- C[length(C)]

  out_all <- rbind(out_PT,out_DHL,cbind(out_AL,rep(out_DHL[,3][length(out_DHL)/3], length(out_AL[,2]))))
  #saveRDS(out_all,file='Tubule_schem_rds.rds')
  
    
  S_TGF <- Splus_TGF/(1+exp(-C_TGF*(C_MD - C_MD_TGF))) + Smin_TGF
  
  if (furo_yn==0){
    S_TGF_func <- S_TGF # + S_0
  }else{
    S_TGF_func <- 0
  }
  
  if (int_yn==0){
    Tp_temp <- Tp_myo
    C_temp <- C_myo
  }else{
    Tp_temp <- Tp_TGF
    C_temp <- C_int
  }
  
  if (dilt_yn){
    T_act<-0
  }

  
  S_myo <- Splus_myo/(1+exp(-C_temp*(Tp - Tp_temp))) + Smin_myo
  S_myo_func <- S_myo
  
  if (TGF_only_yn){
    S_myo_func<-0
  }

  S_tone <- S_TGF_func + S_myo_func #+ C_MyoTGF*S_TGF_func*S_myo_func + S_0
  
  A <- 1/(1+exp(-S_tone))
  
  return(data.frame(D=D,
                    Tp=Tp,
                    Pa=Pavg,
                    T_pass=T_pass,
                    T_act=T_act,
                    S_TGF=S_TGF,
                    S_myo=S_myo,
                    SNGFR=SNGFR,
                    C_MD=C_MD,
                    A=A,
                    S_tone=S_tone,
                    err=Tp-T_pass-A*T_act,
                    Q=Q,
                    Qloop=out_DHL[,3][length(out_DHL)/3]))
    
}

tak_err <- function(P_aa_input,
                    D_aa_input,
                    Ca_input,
                    surr_glom_df,
                    PT_frac,
                    C_0_Tub,
                    furo_yn,
                    dilt_yn,
                    int_yn,
                    TGF_only_yn,
                    int_low,
                    int_high){
  
  D_zer <- 0*P_aa_input
  C_MD_out <- 0*P_aa_input
  Tp_out <- 0*P_aa_input
  S_tone_out <- 0*P_aa_input
  S_myo_out <- 0*P_aa_input
  S_TGF_out <- 0*P_aa_input
  err <- 0*P_aa_input
  SNGFR_out <- 0*P_aa_input
  Q_out <- 0*P_aa_input
  PG_out <- 0*P_aa_input
  Qloop_out <- 0*P_aa_input
  
  D_try <- seq(int_low,int_high,(int_high-int_low)/20)

  for (i in seq(P_aa_input)){
    err_x <- 0*D_try
      for (j in 1:length(D_try)){
      BB1 <- check_takenaka(D=D_try[j],
                         Pa=P_aa_input[i],
                         Ca = Ca_input[i],
                         surr_glom_df = surr_glom_df,
                         PT_frac=PT_frac[i],
                         C_0_Tub=C_0_Tub,
                         furo_yn=furo_yn,
                         dilt_yn=dilt_yn,
                         int_yn=int_yn,
                         TGF_only_yn=TGF_only_yn)
      err_x[j] <- BB1$err
      }
    
    #plot(spline(D_try,err_x))
    ff <- splinefun(D_try,err_x)
    
    if (ff(int_low)*ff(int_high)>0){
      int_high_neg_temp <- D_try[which.min(err_x)]
          if (1==0){#(int_high_neg_temp < int_low){
            int_high_neg <- int_low
            int_low <- int_high_neg_temp
          }else{
            int_high_neg <- int_high_neg_temp
          }
    }else{
      int_high_neg<-int_high
      #print('neg product')
    }

    optim_dummy <- uniroot(ff,interval=c(int_low,int_high_neg),maxiter=1e5)
    #print(optim_dummy)
    D_zer[i] <- optim_dummy$root   #BFfzero(dummy_tak,10,30) #
    #print(D_zer[i])
    
BB <- check_takenaka(D=D_zer[i],
                     Ca = Ca_input[i],
                     Pa=P_aa_input[i],
                     surr_glom_df=surr_glom_df,
                     PT_frac=PT_frac[i],
                     C_0_Tub=C_0_Tub,
                     furo_yn=furo_yn,
                     dilt_yn=dilt_yn,
                     int_yn=int_yn,
                     TGF_only_yn=TGF_only_yn)


A <- BB$A
Tp <- BB$Tp
T_act <- BB$T_act
T_pass <- BB$T_pass
err[i] <- Tp - A*T_act - T_pass
C_MD_out[i] <- BB$C_MD
S_myo_out[i] <- BB$S_myo
S_TGF_out[i] <- BB$S_TGF
Tp_out[i] <- BB$Tp
S_tone_out[i] <- BB$S_tone
SNGFR_out[i] <- BB$SNGFR
PG_out[i] <- 2*BB$Pa - Try_Pf[i]
Q_out[i] <- BB$Q
Qloop_out[i] <- BB$Qloop
}
D_zer_out <- as.numeric(D_zer)

return(list(err=err, D_zer = D_zer_out, C_MD = C_MD_out,S_Myo = S_myo_out, S_TGF=S_TGF_out, PG=PG_out,
            Tp = Tp_out, S_tone=S_tone_out,SNGFR=SNGFR_out,Q_out=Q_out,Qloop=Qloop_out))
}

check_bell <- function(    D, 
                           Pa, 
                           Ca,
                           C_MD,
                           surr_glom_df,
                           PT_frac,
                           C_0_Tub,
                           furo_yn,
                           dilt_yn,
                           int_yn,
                           TGF_only_yn){
  
  GG <- glom_SS(Pa=Pa, Ca=Ca, D=D, surr_glom_df = surr_glom_df)
  SNGFR <- GG$SNGFR
  Pavg <- GG$Pavg
  Q <- GG$Q
  
  Tp <- D*(Pavg - P_ext)/2
  Tp0 <- Tp_cont#DAA_0*(glom_SS(Pa=Paa_i, Ca=Ca_i, D=DAA_0,surr_glom_df = surr_glom_df)$Pavg - P_ext)/2
  
  T_pass <- C_pass0*exp(C_pass1*(D/DAA_0-1))
  T_act <- C_act0*exp(-((D/DAA_0-1)/C_act1)^2)

  S_TGF <- Splus_TGF/(1+exp(-C_TGF*(C_MD - C_MD_TGF))) + Smin_TGF
  
  if (furo_yn==0){
    S_TGF_func <- S_TGF # + S_0
  }else{
    S_TGF_func <- 0
  }
  
  if (int_yn==0){
    Tp_temp <- Tp_myo
    C_temp <- C_myo
  }else{
    Tp_temp <- Tp_TGF
    C_temp <- C_int
  }
  
  if (dilt_yn){
    T_act<-0
  }
  
  
  S_myo <- Splus_myo/(1+exp(-C_temp*(Tp - Tp_temp))) + Smin_myo
  S_myo_func <- S_myo
  
  if (TGF_only_yn){
    S_myo_func<-0
  }
  
  S_tone <- S_TGF_func + S_myo_func #+ C_MyoTGF*S_TGF_func*S_myo_func + S_0
  
  A <- 1/(1+exp(-S_tone))
  
  return(data.frame(D=D,
                    Tp=Tp,
                    Pa=Pavg,
                    T_pass=T_pass,
                    T_act=T_act,
                    S_TGF=S_TGF,
                    S_myo=S_myo,
                    SNGFR=SNGFR,
                    C_MD=C_MD,
                    A=A,
                    S_tone=S_tone,
                    err=Tp-T_pass-A*T_act,
                    Q=Q))
  
}

bell_err <- function(P_aa_input,
                    D_aa_input,
                    Ca_input,
                    C_MD_input,
                    surr_glom_df,
                    PT_frac,
                    C_0_Tub,
                    furo_yn,
                    dilt_yn,
                    int_yn,
                    TGF_only_yn,
                    int_low,
                    int_high){
  
  D_zer <- 0*P_aa_input
  C_MD_out <- 0*P_aa_input
  Tp_out <- 0*P_aa_input
  S_tone_out <- 0*P_aa_input
  S_myo_out <- 0*P_aa_input
  S_TGF_out <- 0*P_aa_input
  err <- 0*P_aa_input
  SNGFR_out <- 0*P_aa_input
  Q_out <- 0*P_aa_input
  PG_out <- 0*P_aa_input
  Qloop_out <- 0*P_aa_input
  
  D_try <- seq(int_low,int_high,(int_high-int_low)/20)
  
  for (i in seq(C_MD_input)){
    err_x <- 0*D_try
    for (j in 1:length(D_try)){
      BB1 <- check_bell(D=D_try[j],
                            Pa=P_aa_input[i],
                            Ca = Ca_input[i],
                            C_MD = C_MD_input[i],
                            surr_glom_df = surr_glom_df,
                            PT_frac=PT_frac[i],
                            C_0_Tub=C_0_Tub,
                            furo_yn=furo_yn,
                            dilt_yn=dilt_yn,
                            int_yn=int_yn,
                            TGF_only_yn=TGF_only_yn)
      err_x[j] <- BB1$err
    }
    
    #plot(spline(D_try,err_x))
    ff <- splinefun(D_try,err_x)
    
    if (ff(int_low)*ff(int_high)>0){
      int_high_neg_temp <- D_try[which.min(err_x)]
      if (1==0){#(int_high_neg_temp < int_low){
        int_high_neg <- int_low
        int_low <- int_high_neg_temp
      }else{
        int_high_neg <- int_high_neg_temp
      }
    }else{
      int_high_neg<-int_high
      #print('neg product')
    }
    
    optim_dummy <- uniroot(ff,interval=c(int_low,int_high_neg),maxiter=1e5)
    #print(optim_dummy)
    D_zer[i] <- optim_dummy$root   #BFfzero(dummy_tak,10,30) #
    #print(D_zer[i])
    
    BB <- check_bell(D=D_zer[i],
                         Ca = Ca_input[i],
                         Pa=P_aa_input[i],
                         C_MD=C_MD_input[i],
                         surr_glom_df=surr_glom_df,
                         PT_frac=PT_frac[i],
                         C_0_Tub=C_0_Tub,
                         furo_yn=furo_yn,
                         dilt_yn=dilt_yn,
                         int_yn=int_yn,
                         TGF_only_yn=TGF_only_yn)
    
    
    A <- BB$A
    Tp <- BB$Tp
    T_act <- BB$T_act
    T_pass <- BB$T_pass
    err[i] <- Tp - A*T_act - T_pass
    C_MD_out[i] <- BB$C_MD
    S_myo_out[i] <- BB$S_myo
    S_TGF_out[i] <- BB$S_TGF
    Tp_out[i] <- BB$Tp
    S_tone_out[i] <- BB$S_tone
    SNGFR_out[i] <- BB$SNGFR
    PG_out[i] <- 2*BB$Pa - Try_Pf[i]
    Q_out[i] <- BB$Q
  }
  D_zer_out <- as.numeric(D_zer)
  
  return(list(err=err, D_zer = D_zer_out, C_MD = C_MD_out,S_Myo = S_myo_out, S_TGF=S_TGF_out, PG=PG_out,
              Tp = Tp_out, S_tone=S_tone_out,SNGFR=SNGFR_out,Q_out=Q_out))
}


ss_Tubule <- function(x, Css, parameters){
  Ce = Ce_calc(x)
  Q <- as.numeric(parameters["Loop_Flow"])
  
  # if (Q <= 0){
  #   Q <- 1/1e6/60
  # }

  dCssdx <- -2*pi*r_Tub*(Vmax*Css/(Km+Css) + P_Tub*(Css - Ce))/Q
  as.list(dCssdx)
}

Ce_calc <- function(x){
  Ce <- 0*x
  if (x < L_PT){
    Ce <- Co
  }
  if (x >= L_PT & x <= (L_PT+L_DHL)){
    A1_PT <- (1-Ce_L_PT/Co)/(1-exp(A3_PT))
    A2_PT <- 1 - A1_PT
    Ce <- Co*(A1_PT*exp(A3_PT*(x-L_PT)/L_DHL) + A2_PT)
  }
  if (x > (L_PT+L_DHL) & x <= (L_PT+L_DHL+L_Tub)){
    A1_Tub <- (1-Ce_L_Tub/Ce_calc(x=L_PT+L_DHL))/(1-exp(-A3_Tub))
    A2_Tub <- 1 - A1_Tub
    Ce <- Ce_calc(x=L_PT+L_DHL) * (A1_Tub*(exp(-A3_Tub*(x-L_PT-L_DHL)/L_Tub)) + A2_Tub)
  }
  return(Ce)
}

ss_Tubule_Richfield <- function(x,BB,parameters){
  C <- BB[1]
  v <- BB[2]
  
  r <- as.numeric(parameters["r"])
  A_M <- as.numeric(parameters["A_M"])
  Pv <- as.numeric(parameters["Pv"])
  #Pq <- as.numeric(parameters["Pq"])
  Vmax <- as.numeric(parameters["Vmax"])
  
  Pv_func <- Pv/(1 + exp(-(v-0.03)*350))
  #Pv_vec <- Pv*rep(1,length(v))

  J_v <-  Pv_func*(C - Ce_calc(x))/A_M
  J_s <- -J_v*C/v - Vmax*C/(C + Km)/v/A_M
  
  dvdx <- 2*pi*r*J_v
  dCdx <- 2*pi*r*J_s
  
  list(c(dCdx,dvdx))
}

Tubule_vx <- function(x,v,parameters){
  
  r <- as.numeric(parameters["r"])
  A_M <- as.numeric(parameters["A_M"])
  Pv <- as.numeric(parameters["Pv"])
  #Pq <- as.numeric(parameters["Pq"])
  Vmax <- as.numeric(parameters["Vmax"])
  
  Pv_func <- Pv/(1 + exp(-(v-0.03)*350))
  #Pv_vec <- Pv*rep(1,length(v))
  
  J_v <-  Pv_func*(C - Ce_calc(x))/A_M

  dvdx <- 2*pi*r*J_v

  list(c(dvdx))
}

#plot velocity profile
# ss_Tubule_Weinstein <- function(x,BB,parameters){
#   C <- BB[1]
#   v <- BB[2]
#   
#   Jo <- as.numeric(parameters["Jo"])
#   A_M <- as.numeric(parameters["A_M"])
#   sigma <- as.numeric(parameters["sigma"])
#   RTLp <- as.numeric(parameters["RTLp"])
#   C_ref <- as.numeric(parameters["C_ref"])
#   h <- as.numeric(parameters["h"])
# 
#   J_v <- Jo - RTLp*sigma*(C - Ce_calc(x))
#   J_s <- N + h*(C - Ce_calc(x)) + J_v*(1-sigma)*C_ref
#   
#   dvdx <- -J_v/A_M
#   dCdx <- -J_s/A_M
#   
#   list(c(dCdx,dvdx))
# }

tubule_param_choose_PT <- function(alpha,ss=5.55/32.67){
  
  Jo <- alpha*Jo
  RTLp <- alpha*RTLp
  
  parameters_PT <- c(Jo=Jo,
                     A_M=A_M,
                     sigma=sigma,
                     RTLp=RTLp,
                     C_ref=C_ref,
                     h=h)
  x_PT <- seq(0,(L_DHL + L_PT),(L_DHL + L_PT)/1000)
  out_PT <- ode(y = c(Co,vo), 
                times = x_PT, 
                func = ss_Tubule_Weinstein,
                parms = parameters_PT)
  
  return(out_PT[,3][length(out_PT[,3])]/vo - ss)
}

glom_SS <- function(Pa,Ca,D,surr_glom_df){
  ee <- 1e-6
  indx1 <- (surr_glom_df$P.in < (Pa + 2.5) & surr_glom_df$P.in > (Pa - 2.5))
  indx2 <- (round(surr_glom_df$D.in - D,2) < 0.05 & round(surr_glom_df$D.in - D,2) > -0.05)
  indx <- which(indx1 & indx2)
  
    if(!length(indx)){
    print("I can't find that entry!")
    print("Pa, Ca, D=")
    print(Pa)
    print(Ca)
    print(D)
    return(0)
  }
    
  df_try <- surr_glom_df[indx,]
  #df_try <- distinct(df_try)
  #print(df_try)
  
  Ca_findr <- df_try$Ca.in-Ca
  Ca_findr_sort <- sort(Ca_findr)
  #df_try <- df_try[which(df_try$Ca.in - Ca <= Ca_findr_sort[2]),]
  
  #print(length(indx))
  
  if(dim(df_try)[1] == 1){
    Q <- df_try$Q
    SNGFR <- df_try$SNGFR
    Pavg <- df_try$Pavg
  }
  
  if(dim(df_try)[1] == 2){
    
    input_vec <- c(Pa,D,Ca)
    
    u1 <- round(abs(df_try[,1][1] - df_try[,1][2]),2)
    u2 <- round(abs(df_try[,2][1] - df_try[,2][2]),2)
    u3 <- round(abs(df_try[,3][1] - df_try[,3][2]),2)
    indx_test <- which(c(u1,u2,u3) != 0)
    
    xi <- input_vec[indx_test]
    #print(xi)
    Q <- interp1(x=df_try[,indx_test], y=df_try$Q, xi=xi)
    SNGFR <- interp1(x=df_try[,indx_test],y=df_try$SNGFR,xi=xi)
    Pavg <- interp1(x=df_try[,indx_test],y=df_try$Pavg,xi=xi)
  }

  if(dim(df_try)[1] == 4){
    input_vec <- c(Pa,D,Ca)
    l_length <- 0*c(1,1,1)
    l_length[1] <- length(unique(df_try[,1]))
    l_length[2] <- length(unique(df_try[,2]))
    l_length[3] <- length(unique(df_try[,3]))
    indx_test <- which(l_length > 1)
    #print(indx_test)
    
    xi <- input_vec[indx_test]
    
    x <- linspace(min(df_try[,indx_test[1]]), max(df_try[,indx_test[1]]), 2)
    y <- linspace(min(df_try[,indx_test[2]]), max(df_try[,indx_test[2]]), 2)
    
    zQ <- array(df_try$Q, dim=c(length(x),length(y)))
    zS <- array(df_try$SNGFR, dim=c(length(x),length(y)))
    zP <- array(df_try$Pavg, dim=c(length(x),length(y)))
    
    Q <- interp2(x = x, y = y, Z = zQ, xp = xi[1], yp = xi[2])
    SNGFR <- interp2(x = x, y = y, Z = zS, xp = xi[1], yp = xi[2])
    Pavg <- interp2(x = x, y = y, Z = zP, xp = xi[1], yp = xi[2])
  }  

  if(dim(df_try)[1] == 8){

    x <- unique(df_try[,1])
    y <- unique(df_try[,2])
    z <- unique(df_try[,3])
    
    fQ <- array(df_try$Q, dim=c(length(x),length(y),length(z)))
    fS <- array(df_try$SNGFR, dim=c(length(x),length(y),length(z)))
    fP <- array(df_try$Pavg, dim=c(length(x),length(y),length(z)))
    
    Q <- approx3d(x=x, y=y, z=z, f=fQ, xout = Pa, yout = D, zout = Ca)
    SNGFR <- approx3d(x=x, y=y, z=z, f=fS, xout = Pa, yout = D, zout = Ca)
    Pavg <- approx3d(x=x, y=y, z=z, f=fP, xout = Pa, yout = D, zout = Ca)
  }  
  #err <- sqrt((df_try$P.in - Pa)^2 + (df_try$Ca.in - Ca)^2 + (df_try$D.in - D)^2)

  # if (any(abs(err) < 1e-11)){
  #   Q <- mean(df_try$Q[abs(err) < 1e-11])
  #   Pavg <- mean(df_try$Pavg[abs(err) < 1e-11])
  #   SNGFR <- mean(df_try$SNGFR[abs(err) < 1e-11])
  # }else{
  #   
  # weights <- (1/err)/sum(1/err)
  # #print(weights)
  # Q <- sum(weights*df_try$Q/sum(weights))
  # Pavg <- sum(weights*df_try$Pavg/sum(weights))
  # SNGFR <- sum(weights*df_try$SNGFR/sum(weights))
  #  }
  return(list(Q=Q,Pavg=Pavg,SNGFR=SNGFR))
}

autoreg <- function(D, 
                    Pa, 
                    Ca,
                    surr_glom_df,
                    furo_yn,
                    dilt_yn,
                    int_yn,
                    C_myo_t){
  
  GG_base <- glom_SS(Pa=Paa_i, Ca=Ca_i, D=DAA_0,surr_glom_df = surr_glom_df)
  Tp0 <- DAA_0*(GG_base$Pavg - P_ext)/2
  SNGFR <- GG_base$SNGFR
  
  Tp_keep <- rep(0,Nt)
  D_keep <- rep(0,Nt)
  Q_keep <- rep(0,Nt)
  SNGFR_keep <- rep(0,Nt)
  S_Myo_keep <- rep(0,Nt)
  S_TGF_keep <- rep(0,Nt)
  
  C_PT_keep <- rep(0,Nt)
  C_DHL_keep <- rep(0,Nt)
  C_MD_keep <- rep(0,Nt)
  
  vo_keep <- rep(0,Nt)
  v_DHL_keep <- rep(0,Nt)
  v_AL_keep <- rep(0,Nt)
  
  P_avg_keep <- rep(0,Nt)
  P_G_keep <- rep(0,Nt)
  
  x_PT <- seq(0,(L_PT),(L_PT)/1000)
  x_DHL <- seq(L_PT,(L_PT+L_DHL),(L_DHL)/1000)
  
  parameters_PT <- c(Pv=Pv_PT,
                     A_M=A_M,Vmax=Vmax_PT,r=r)
  parameters_DHL <- c(Pv=Pv_DHL,
                      A_M=A_M,Vmax=Vmax_DHL,r=r)
  
  vo <- SNGFR/A_M/1e6/60
  x_PT <- seq(0,(L_PT),(L_PT)/1000)
  out_PT <- ode(y = c(Co,vo), 
                times = x_PT, 
                func = ss_Tubule_Richfield,
                parms = parameters_PT)
  x_DHL <- seq(L_PT,(L_PT+L_DHL),(L_DHL)/1000)
  out_DHL <- ode(y = out_PT[length(x_PT),2:3], 
                 times = x_DHL, 
                 func = ss_Tubule_Richfield,
                 parms = parameters_DHL)  
  
  x_AL <- seq(L_PT+L_DHL,(L_PT+L_DHL+L_Tub),(L_Tub)/1000)
  parameters_AL <- c(Loop_Flow = out_DHL[,3][length(out_DHL)/3]*A_M)
  out_AL <- ode(y = out_DHL[,2][length(out_DHL)/3], 
                times = x_AL, 
                func = ss_Tubule,
                parms = parameters_AL)
  C_PT <- out_PT[,2]
  C_DHL <- out_DHL[,2]
  C_AL <- out_AL[,2]
  C_MD <- C_AL[length(C_AL)]
  
  v_PT <- out_PT[,3]
  v_DHL <- out_DHL[,3]
  v_AL <- rep(v_DHL[length(v_DHL)], length(x_AL))

  print('done with steady state ########')
  print('percent progress:')
  
  for (t in time_in){
  
    if (t%%100==0){
      print(t/Nt*100)
    }
    
    GG <- glom_SS(Pa=Pa[t],Ca=Ca,D=D,surr_glom_df = surr_glom_df)
    SNGFR <- GG$SNGFR
    Pavg <- GG$Pavg
    Q <- GG$Q
  
    #print(GG)
    
    Tp <- D*(Pavg - P_ext)/2
  
    T_pass <- C_pass0*exp(C_pass1*(D/DAA_0-1))
    T_act <- C_act0*exp(-((D/DAA_0-1)/C_act1)^2)
    
    x <- seq(0, L_DHL+L_PT, by = (L_DHL+L_PT)/1000)
  
    vo <- SNGFR/A_M/1e6/60
    out_PT_int <- ode(y = c(Co,vo), 
                  times = x_PT, 
                  func = ss_Tubule_Richfield,
                  parms = parameters_PT,maxsteps=1e5)
    v_PT<-out_PT_int[,3]
    
    Pv_func <- Pv_PT/(1 + exp(-(v_PT-0.03)*350))

    J_v <-  Pv_func*(C_PT - unlist(lapply(x_PT,Ce_calc)))/A_M
    J_s <- -J_v*C_PT/v_PT - Vmax_PT*C_PT/(C_PT + Km)/v_PT/A_M
    
    dC_PTdx <- c(0,(C_PT[2:length(C_PT)] - C_PT[1:(length(C_PT)-1)]))/(x_PT[2]-x_PT[1])
    
    Ct_PT <- 2*pi*r*J_s - dC_PTdx*v_PT
    
    #print(v_PT)
    C_PT <- C_PT + dt*Ct_PT

    out_DHL_int <- ode(y = c(C_PT[length(C_PT)],v_PT[length(v_PT)]), 
                      times = x_DHL, 
                      func = ss_Tubule_Richfield,
                      parms = parameters_DHL)
    v_DHL<-out_DHL_int[,3]
    
    Pv_func <- Pv_DHL/(1 + exp(-(v_DHL-0.03)*350))
    
    J_v <-  Pv_func*(C_DHL - unlist(lapply(x_DHL,Ce_calc)))/A_M
    J_s <- -J_v*C_DHL/v_DHL - Vmax_DHL*C_DHL/(C_DHL + Km)/v_DHL/A_M
    
    dC_DHLdx <- c(C_DHL[1]-C_PT[length(C_PT)],
                  (C_DHL[2:length(C_DHL)] - C_DHL[1:(length(C_DHL)-1)]))/(x_DHL[2]-x_DHL[1])
    
    Ct_DHL <- 2*pi*r*J_s - dC_DHLdx*v_DHL
    
    #print(v_PT)
    C_DHL <- C_DHL + dt*Ct_DHL
    
    v_AL <- rep(v_DHL[length(v_DHL)]*A_M/(pi*r_Tub^2),length(x_AL))
    
    J_s <- - (Vmax*C_AL/(C_AL + Km) + P_Tub*(C_AL - unlist(lapply(x_AL,Ce_calc))))
    
    dC_ALdx <- c(C_AL[1]-C_DHL[length(C_DHL)],
                  (C_AL[2:length(C_AL)] - C_AL[1:(length(C_AL)-1)]))/(x_AL[2]-x_AL[1])
    
    Ct_AL <- 2*J_s/r_Tub - dC_ALdx*v_AL
    
    #print(v_PT)
    C_AL <- C_AL + dt*Ct_AL   
    
    C_MD <- C_AL[length(C_AL)]
  
    S_TGF <- Splus_TGF/(1+exp(-C_TGF*(C_MD - C_MD_TGF))) + Smin_TGF
  
    if (furo_yn==0){
      S_TGF_func <- S_TGF # + S_0
    }else{
      S_TGF_func <- 0
    }
  
    if (int_yn==0){
      Tp_temp <- Tp_myo
      C_temp <- C_myo
    }else{
      Tp_temp <- Tp_TGF
      C_temp <- C_int
    }
  
    #MYOGENIC DELAY
    if (t == 1){
      #Tp_m <- Tp0
     S_myo <- Splus_myo/(1+exp(-C_temp*(Tp0 - Tp_temp))) + Smin_myo
     S_myo_ss <- S_myo
    }else{
      #Tp_m <- Tp_keep[t-1]
      S_myo <- S_Myo_keep[t-1]
      #S_myo_ss <- Splus_myo/(1+exp(-C_temp*(Tp_keep[t - myo_del_ind] - Tp_temp))) + Smin_myo
      S_myo_ss <- Splus_myo/(1+exp(-C_temp*(Tp - Tp_temp))) + Smin_myo
    }
    
  dS_myodt <- -C_myo_t*(S_myo-S_myo_ss)
  
  #print(S_myo)
  
  S_myo_func <- S_myo + dS_myodt*dt
  
  S_tone <- S_TGF_func + S_myo_func #+ C_MyoTGF*S_TGF_func*S_myo_func + S_0
  
  A <- 1/(1+exp(-S_tone))
  
  Dt <- (Tp - T_pass-A*T_act)/tau_c
  
  
  
  D <- D + Dt*dt
  
 # print(D)
  
  # if (t==1){
  #   Tp_m <- Tp0
  #   S_myo <- S_myo_ss
  # }else{
  #   Tp_m <- Tp_keep[t-1]
  #   S_myo <- S_Myo_keep[t-1]
  # }
  
  Tp_keep[t] <- Tp 
  D_keep[t] <- D
  Q_keep[t] <- Q
  SNGFR_keep[t] <- SNGFR
  S_Myo_keep[t] <- S_myo_func
  S_TGF_keep[t] <- S_TGF
  
  C_PT_keep[t] <- C_PT[length(C_PT)]
  C_DHL_keep[t] <- C_DHL[length(C_DHL)]
  C_MD_keep[t] <- C_MD
  
  vo_keep[t] <- vo
  v_DHL_keep[t] <- v_DHL[1]
  v_AL_keep[t] <- v_AL[1]
  
  P_avg_keep[t] <- Pavg
  P_G_keep[t] <- 2*Pavg - Pa[t]
  
  }
  
  #plot(c(x_PT,x_DHL,x_AL),c(C_PT,C_DHL,C_AL))
  
  return(data.frame(Tp=Tp_keep,
                    D=D_keep,
                    Q=Q_keep,
                    SNGFR=SNGFR_keep,
                    S_Myo=S_Myo_keep,
                    S_TGF=S_TGF_keep,
                    vo_keep,
                    v_DHL_keep,
                    v_AL_keep,
                    C_PT_keep,
                    C_DHL_keep,
                    C_MD_keep,
                    P_avg_keep,
                    P_G_keep))
}



# multi_glom_SS <- function(Pa_vec,Ca_vec,D_vec,surr_glom_df){
#   Q_keep <- 0*Pa_vec
#   Pavg_keep <- 0*Pa_vec
#   SNGFR_keep <- 0*Pa_vec
#   
#   for (i in seq(Pa_vec)){
#     GG <- glom_SS(Pa = Pa_vec[i], Ca = Ca_vec[i], D = D_vec[i],surr_glom_df=surr_glom_df)
#     Q_keep[i] <- GG$Q
#     Pavg_keep[i] <- GG$Pavg
#     SNGFR_keep[i] <- GG$SNGFR
#   }
#   return(data.frame(Q=Q_keep, Pavg=Pavg_keep, SNGFR=SNGFR_keep))
# }


BFfzero <- function (f, a, b, num = 100, eps = 1e-05)
{
  h = abs(b - a)/num
  i = 0
  j = 0
  a1 = b1 = 0
  while (i <= num) {
    a1 = a + i * h
    b1 = a1 + h
    if (f(a1) == 0) {
      print(a1)
      print(f(a1))
    }
    else if (f(b1) == 0) {
      print(b1)
      print(f(b1))
    }
    else if (f(a1) * f(b1) < 0) {
      repeat {
        if (abs(b1 - a1) < eps)
          break
        x <- (a1 + b1)/2
        if (f(a1) * f(x) < 0)
          b1 <- x
        else a1 <- x
      }
      print(j + 1)
      j = j + 1
      print((a1 + b1)/2)
      print(f((a1 + b1)/2))
    }
    i = i + 1
  }
  if (j == 0)
    print("finding root is fail")
  else print("finding root is successful")
  return(x)
}

#########################################################################
# k_calibrate <- function(src,
#                         trg,
#                         D,
#                         L,
#                         k0,
#                         k1,
#                         SNGFR_ex,
#                         Paa_i,
#                         Pea_o,
#                         P_bs,
#                         Ca_i,
#                         H_t_sys,
#                         mu_plas,
#                         Rra,
#                         t,
#                         in.nodes,
#                         out.nodes,
#                         num.iter.g,
#                         num.iter.k,
#                         SNGFR.k.tol,
#                         SNGFR.tol,
#                         Rinf.tol){
#   
#   SNGFR.k.err <- SNGFR.k.tol*5
#   k_keep <- 0*seq(num.iter.k)
#   S_keep <- 0*seq(num.iter.k)
#   count <- 1
#   
#   G0   <-      run_glom(src = src,
#                         trg = trg,
#                         D = D,
#                         L = L,
#                         k = k0*src/src,
#                         Paa_i = Paa_i,
#                         Pea_o = Pea_o,
#                         P_bs = P_bs,
#                         Ca_i = Ca_i,
#                         H_t_sys = H_t_sys,
#                         mu_plas = mu_plas,
#                         Rra=Rra,
#                         t = t,
#                         in.src = in.nodes,
#                         out.trg = out.nodes,
#                         num.iter = num.iter.g,
#                         SNGFR.tol = SNGFR.tol,
#                         Rinf.tol = Rinf.tol)
#   
#   while(SNGFR.k.err > SNGFR.k.tol & count < num.iter.k){
#     
#     G1   <-      run_glom(src = src,
#                           trg = trg,
#                           D = D,
#                           L = L,
#                           k = k1*src/src,
#                           Paa_i = Paa_i,
#                           Pea_o = Pea_o,
#                           P_bs = P_bs,
#                           Ca_i = Ca_i,
#                           H_t_sys = H_t_sys,
#                           mu_plas = mu_plas,
#                           Rra=Rra,
#                           t = t,
#                           in.src = in.nodes,
#                           out.trg = out.nodes,
#                           num.iter = num.iter.g,
#                           SNGFR.tol = SNGFR.tol,
#                           Rinf.tol = Rinf.tol)
#     
#     S0 <- G0$SNGFR
#     S1 <- G1$SNGFR
#     
#     m <- (S1-S0)/(k1-k0)
#     k_new <- (SNGFR_ex - S1 + m*k1)/m
#     
#     k0 <- k1
#     G0 <- G1
#     k1 <- k_new
#     
#     k_keep[count] <- k1
#     S_keep[count] <- S1
#     
#     count <- count+1
#     SNGFR.k.err <- abs(k1 - k0)
#   }
#   
#   plot(k_keep[1:count])
#   
#   if (count>=num.iter.k & num.iter.k > 2){
#     k0 <- k_keep[which.min(abs(S_keep-SNGFR_ex))-1]
#     G1   <-      run_glom(src = src,
#                           trg = trg,
#                           D = D,
#                           L = L,
#                           k = k0*src/src,
#                           Paa_i = Paa_i,
#                           Pea_o = Pea_o,
#                           P_bs = P_bs,
#                           Ca_i = Ca_i,
#                           H_t_sys = H_t_sys,
#                           mu_plas = mu_plas,
#                           Rra=Rra,
#                           t = t,
#                           in.src = in.nodes,
#                           out.trg = out.nodes,
#                           num.iter = num.iter.g,
#                           SNGFR.tol = SNGFR.tol,
#                           Rinf.tol = Rinf.tol)
#     print("Iterations Exceeded! Choosing optimal k...")
#   }
#   
#   return(list(OP=G1,k=k0,k_keep=k_keep,S_keep=S_keep))
# }







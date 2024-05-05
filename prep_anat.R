#Prep anatomy
#created March 18, 2020
#last updated March 28, 2020

in.nodes <- 1
out.nodes <- 195

source("shea_anatomy_pressure.R")
src <- c(source_nodes1,194)
trg <- c(target_nodes1,out.nodes)
edge_diam1[which(src%in%in.nodes)] <- DAA_m
edge_lengths1[which(src%in%in.nodes)] <- L_AA
D <- c(edge_diam1,DEA_m)
L <- c(edge_lengths1,L_EA)

# src <- c(1,2,2,3,4,5)                                                       #original source nodes
# trg <- c(2,3,4,5,5,6)                                                       #original target nodes
# D <- 2*src/src#c(2,2,2,2)
# L <- 6*src/src#c(6,6,6,6)                                                         #length of each segment




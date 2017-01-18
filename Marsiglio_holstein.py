# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 01:44:04 2016

@author: Ben
"""

from numpy import *
from init_functions import *
from Functions import *
import time

################
###############
##############
##############  WHAT ABOUT SUBTRACTING THE REAL PART AT ZERO FREQUENCY TO PREVENT CHEMICAL POTENTIAL SHIFT!!!!????
###GONNA RUN AT HALF FILLING TO AVOID THIS FOR NOW


###Nk must be odd or else the momentum points do not 
###form a group under addition!
###Choose Nw to be even (number of Fermionic frequencies)

tstart = time.time()

[Nk,Nw,beta,g,omega,blank,superconductivity] = load("data_holstein/params.npy")
Nk = int(Nk)
Nk = 81
Nw = int(Nw)

N_split_omega = 8

iter_selfconsistency = 5

kxs, kys  = init_momenta(Nk)
iw_bose   = init_boson_freq(Nw, beta)
iw_fermi  = init_fermion_freq(Nw, beta)
band      = init_band(kxs, kys, Nk)
ws,Nr     = init_ws(N_split_omega, omega, iter_selfconsistency, band, Nk)

print "Nr ",Nr

omeg = N_split_omega

save("data_Marsiglio_holstein/params",asarray([Nk,Nw,beta,g,omega,0.0,superconductivity,Nr]))
save("data_Marsiglio_holstein/ws",    ws)

Gloc = load("data_holstein/Gloc.npy")
Term1 = zeros([Nr,2,2],dtype=complex)

##### check sign of Term1. Different than formula in JJ paper (absorbed into D?) #####
for ir in range(Nr):
    for iw in range(Nw):
        Term1[ir,:,:] -= 1./(Nk**2 * beta) * computeD(ws[ir]-iw_fermi[iw], omega) * einsum('ij,jk,kl->il',tau3,Gloc[iw,:,:],tau3)

print " "
print "Done computing Term1"
print " "

Sigma     = zeros([Nr,2,2], dtype=complex)

nB = 1./(exp(beta*omega) - 1.0)
nF = zeros(Nr)
for ir in range(Nr):
    nF[ir] = 1./(exp(beta*ws[ir]) + 1.0)

for myiter in range(iter_selfconsistency):

    Sigma_old = Sigma.copy()
    Gloc      = zeros([Nr,2,2], dtype=complex)
    
    #compute new Gloc
    for ik1 in range(Nk):
        for ik2 in range(Nk):
            for ir in range(Nr):
                Gloc[ir,:,:] += linalg.inv((ws[ir]+0.0001*1j)*tau0 - band[ik1,ik2]*tau3 - Sigma[ir,:,:])

    Term2     = zeros([Nr,2,2], dtype=complex)
    for ir in range(Nr):
        if ir-omeg>=0:
            Term2[ir,:,:] += (nB+1.-nF[ir-omeg]) * 1./Nk**2 * g**2 * einsum('ij,jk,kl->il',tau3,Gloc[ir-omeg,:,:],tau3)
        if ir+omeg<Nr:
            Term2[ir,:,:] += (nB+nF[ir+omeg]) * 1./Nk**2 * g**2 * einsum('ij,jk,kl->il',tau3,Gloc[ir+omeg,:,:],tau3)
            
    Sigma = Term1 + Term2

    change = zeros([2,2], dtype=complex)    
    change += sum(abs(Sigma-Sigma_old), axis=0)

    print "change ",change
    print "time elapsed ",time.time() - tstart
    
    save("data_Marsiglio_holstein/Gloc.npy", Gloc)
    save("data_Marsiglio_holstein/Sigma.npy", Sigma)    
    
    
                        
print "total run time ", time.time() - tstart
        

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

[Nk,Nw,beta,g,omega,q0,superconductivity] = load("data/params.npy")
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
print "Nk ",Nk
print "Nw ",Nw
print "g  ",g
print "omega ",omega
print "beta ",beta

omeg = N_split_omega

gofq      = init_gofq(kxs, kys, Nk, g, q0)
fft_gofq2 = fft.fft2(gofq**2)

save("data_Marsiglio_forward/params",asarray([Nk,Nw,beta,g,omega,0.0,superconductivity,Nr]))
save("data_Marsiglio_forward/ws",    ws)

G = load("data/G.npy")
Term1 = zeros([Nk,Nk,Nr,2,2],dtype=complex)
Conv = zeros([Nk,Nk,2,2], dtype=complex)

##### check sign of Term1. Different than formula in JJ paper (absorbed into D?) #####
for m in range(Nw):
    fft_G = fft.fft2(einsum('ij,abjk,kl->abil',tau3,G[:,:,m,:,:],tau3), axes=(0,1))

    for ir in range(Nr):
        Conv = 1./(Nk**2 * beta) * computeD(ws[ir]-iw_fermi[m], omega) * fft.ifft2( einsum('ij,ijab->ijab', fft_gofq2 , fft_G) , axes=(0,1))
        Conv = roll(Conv, -(Nk-1)/2, axis=0)
        Conv = roll(Conv, -(Nk-1)/2, axis=1)
        
        Term1[:,:,ir,:,:] -= Conv

print " "
print "Done computing Term1"
print " "

G = zeros([Nk,Nk,Nr,2,2], dtype=complex)
fft_G = zeros([Nk,Nk,Nr,2,2], dtype=complex)
Term2 = zeros([Nk,Nk,Nr,2,2], dtype=complex)
Sigma = Term1 + Term2

nB = 1./(exp(beta*omega) - 1.0)
nF = zeros(Nr)
for ir in range(Nr):
    nF[ir] = 1./(exp(beta*ws[ir]) + 1.0)

for myiter in range(iter_selfconsistency):

    print " "
    print myiter
    
    #compute new G
    for ik1 in range(Nk):
        for ik2 in range(Nk):
            for n in range(Nr):
                
                G[ik1,ik2,n,:,:] = linalg.inv((ws[n]-0.0001*1j)*tau0 - band[ik1,ik2]*tau3 - Sigma[ik1,ik2,n,:,:])

    print "G computed"
                
    #compute fft_G
    for n in range(Nr):
        fft_G[:,:,n,:,:] = fft.fft2(einsum('ij,abjk,kl->abil',tau3,G[:,:,n,:,:],tau3), axes=(0,1))

    print "fft_G computed"
        
    Sigma_old = Sigma.copy()
    
    Term2 = zeros([Nr,2,2], dtype=complex)
    for ir in range(Nr):

        Conv = fft.ifft2( einsum('ij,ijab->ijab', fft_gofq2 , fft_G) , axes=(0,1))
        if ir+omeg<Nr:
            Term2[ir+omeg,:,:] += (nB+1.-nF[ir]) * 1./Nk**2 * Conv
        if ir-omeg>=0:
            Term2[ir-omeg,:,:] += (nB+nF[ir]) * 1./Nk**2 * Conv
            
                        
    Sigma = Term1 + Term2

    change = zeros([2,2], dtype=complex)    
    change += sum(abs(Sigma-Sigma_old), axis=0)

    print "change ",change
    print "time elapsed ",time.time() - tstart
    
    save("data_Marsiglio_forward/Gloc.npy", Gloc)
    save("data_Marsiglio_forward/Sigma.npy", Sigma)    
    
                        
print "total run time ", time.time() - tstart
        

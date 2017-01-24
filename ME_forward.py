# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 01:44:04 2016

@author: Ben
"""

from numpy import *
from init_functions import *
from Functions import *
import time
import subprocess


###Nk must be odd or else the momentum points do not 
###form a group under addition!
###Choose Nw to be even (number of Fermionic frequencies)

#correcting for chemical potential shift!

tstart = time.time()

Nk    = 41
Nw    = 200
beta  = 4.8
g     = 1.0
omega = 1.2
q0 = 12345678.9
#q0 = 0.2
superconductivity = True

save("data_forward/params",asarray([Nk,Nw,beta,g,omega,q0,superconductivity]))

#assert 1==0

iter_selfconsistency = 8

kxs, kys  = init_momenta(Nk)
gofq      = init_gofq(kxs, kys, Nk, g, q0)
iw_bose   = init_boson_freq(Nw, beta)
iw_fermi  = init_fermion_freq(Nw, beta)
band      = init_band(kxs, kys, Nk)
D         = init_D(Nw, beta, omega, iw_bose)

#now do the same calculation but with the FFT
#G         = init_G(Nk, Nw, beta, omega, band, kxs, kys, iw_fermi, superconductivity)
#G         = load("data/G.npy")

G         = zeros([Nk,Nk,Nw,2,2], dtype=complex)
Conv      = zeros([Nk,Nk,2,2], dtype=complex)
Sigma     = zeros([Nk,Nk,Nw,2,2], dtype=complex)

fft_gofq2 = fft.fft2(gofq**2)

#selfconsistency loop
for myiter in range(iter_selfconsistency):

    #compute new G
    for ik1 in range(Nk):
        for ik2 in range(Nk):
            
            for n in range(Nw):
                iwn = iw_fermi[n]
                
                G[ik1,ik2,n,:,:] = linalg.inv(iwn*tau0 - band[ik1,ik2]*tau3 - Sigma[ik1,ik2,n,:,:])
    
    Sigma_old = Sigma.copy()
    Sigma_proc = zeros([Nk,Nk,Nw,2,2], dtype=complex)
    Sigma = zeros([Nk,Nk,Nw,2,2], dtype=complex)
    
    #compute new Sigma    
    change = zeros([2,2], dtype=complex)

    for m in range(Nw): ### check this bound. it was n
        fft_G = fft.fft2(einsum('ij,abjk,kl->abil',tau3,G[:,:,m,:,:],tau3), axes=(0,1))

        for n in range(Nw): 
            n_m = subtract_freqs(n,m,Nw)
                            
            if(n_m>=0 and n_m<Nw-1):

                Conv = 1.0/(Nk**2)/beta * D[n_m] * fft.ifft2( einsum('ij,ijab->ijab', fft_gofq2 , fft_G) , axes=(0,1))
                Conv = roll(Conv, -(Nk-1)/2, axis=0)
                Conv = roll(Conv, -(Nk-1)/2, axis=1)

                Sigma_proc[:,:,n,:,:] -= Conv
                                    
    Sigma = Sigma_proc
    
    change += sum(abs(Sigma-Sigma_old), axis=(0,1,2))

    print " "
    print "iteration ",myiter
    print "change ", change
    print "iteration time ",time.time() - tstart
    print "filling : ", sum(G[:,:,:,0,0], axis=(0,1,2))/Nk**2/beta
    
    save("data_forward/G.npy", G)
    save("data_forward/Sigma.npy", Sigma)    

        
print "total run time ", time.time() - tstart
        

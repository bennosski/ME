# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 01:44:04 2016

@author: Ben
"""

from numpy import *
from init_functions import *
from Functions import *
import time
from mpi4py import MPI

comm = MPI.COMM_WORLD
nprocs = comm.size
myrank = comm.rank


###Nk must be odd or else the momentum points do not 
###form a group under addition!
###Choose Nw to be even (number of Fermionic frequencies)

tstart = time.time()

Nk = 11
Nw = 200
beta = 80.0
g = sqrt(0.12)
omega = 0.2
q0 = 123456789.
superconductivity = True

if myrank==0:
    save("data/params.txt",asarray([Nk,Nw,beta,g,omega,q0,superconductivity]))

    
iter_selfconsistency = 15

kxs, kys  = init_momenta(Nk)
gofq      = init_gofq(kxs, kys, Nk, g, q0)
iw_bose   = init_boson_freq(Nw, beta)
iw_fermi  = init_fermion_freq(Nw, beta)
band      = init_band(kxs, kys, Nk)
G         = init_G(Nk, Nw, beta, omega, band, kxs, kys, iw_fermi, superconductivity)
D         = init_D(Nw, beta, omega, iw_bose)
Sigma     = zeros([Nk,Nk,Nw,2,2],dtype=complex)




#now do the same calculation but with the FFT
G         = init_G(Nk, Nw, beta, omega, band, kxs, kys, iw_fermi, superconductivity)
Conv      = zeros([Nk,Nk,Nw,2,2], dtype=complex)
Sigma     = zeros([Nk,Nk,Nw,2,2],dtype=complex)

fft_gofq2 = fft.fft2(gofq**2)

#selfconsistency loop
for myiter in range(iter_selfconsistency):

    Sigma_old = Sigma.copy()
    Sigma_proc = zeros([Nk,Nk,Nw,2,2], dtype=complex)
    Sigma = zeros([Nk,Nk,Nw,2,2], dtype=complex)
    
    #compute new Sigma    
    change = zeros([2,2], dtype=complex)

    for m in range(myrank,Nw,nprocs): ### check this bound. it was n
        fft_G = fft.fft2(einsum('ij,abjk,kl->abil',tau3,G[:,:,m,:,:],tau3), axes=(0,1))

        for n in range(Nw): 
            n_m = subtract_freqs(n,m,Nw)
                            
            if(n_m>=0 and n_m<Nw-1):

                Conv[:,:,n,:,:] = 1.0/(Nk**2)/beta * D[n_m] * fft.ifft2( einsum('ij,ijab->ijab', fft_gofq2 , fft_G) , axes=(0,1))
                
                Conv[:,:,n,:,:] = roll(Conv[:,:,n,:,:], -(Nk-1)/2, axis=0)
                Conv[:,:,n,:,:] = roll(Conv[:,:,n,:,:], -(Nk-1)/2, axis=1)

                Sigma_proc[:,:,n,:,:] -= Conv[:,:,n,:,:]
                    
                
    comm.Allreduce(Sigma_proc, Sigma, op=MPI.SUM)
    
    change += sum(abs(Sigma-Sigma_old), axis=(0,1,2))
                  
    #compute new G
    for ik1 in range(Nk):
        kx = kxs[ik1]
        for ik2 in range(Nk):
            ky = kys[ik2]
            
            for n in range(Nw):
                iwn = iw_fermi[n]
                
                G[ik1,ik2,n,:,:] = linalg.inv((iwn*tau0 - band[ik1,ik2]*tau3 - Sigma[ik1,ik2,n,:,:]))

    if myrank==0:
        print "change ", change
        print "iteration time ",time.time() - tstart

        
if myrank==0:
    print "total run time ", time.time() - tstart
    equaltime = sum(G, axis=2)
    save("data/Geqfft.npy", equaltime)
    save("data/G.npy", G[:,:,:,:,:])
    save("data/Sigma.npy", Sigma[:,:,:,:,:])

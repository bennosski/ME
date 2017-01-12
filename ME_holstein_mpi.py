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

################
##################
#################
###############
##############
##############  WHAT ABOUT SUBTRACTING THE REAL PART AT ZERO FREQUENCY TO PREVENT CHEMICAL POTENTIAL SHIFT!!!!????
###GONNA RUN AT HALF FILLING TO AVOID THIS FOR NOW


###Nk must be odd or else the momentum points do not 
###form a group under addition!
###Choose Nw to be even (number of Fermionic frequencies)

tstart = time.time()

Nk    = 41
Nw    = 200
beta  = 2.4
g     = 0.4
omega = 0.8
superconductivity = True

if myrank==0:
    save("data_holstein/params",asarray([Nk,Nw,beta,g,omega,0.0,superconductivity]))
    
iter_selfconsistency = 8

kxs, kys  = init_momenta(Nk)
iw_bose   = init_boson_freq(Nw, beta)
iw_fermi  = init_fermion_freq(Nw, beta)
band      = init_band(kxs, kys, Nk)
D         = init_D(Nw, beta, omega, iw_bose)

#now do the same calculation but with the FFT
Gloc      = zeros([Nw,2,2], dtype=complex)
Sigma     = zeros([Nw,2,2], dtype=complex)

if superconductivity:
    Sigma[:,0,1] = 0.01*1j
    Sigma[:,1,0] = 0.01*1j
    
#selfconsistency loop
for myiter in range(iter_selfconsistency):

    Gloc = zeros([Nw,2,2], dtype=complex)
    
    #compute new Gloc
    for ik1 in range(Nk):
        for ik2 in range(Nk):
            
            for n in range(Nw):
                iwn = iw_fermi[n]
                
                Gloc[n,:,:] += linalg.inv(iwn*tau0 - band[ik1,ik2]*tau3 - Sigma[n,:,:])

                
    Sigma_old = Sigma.copy()
    Sigma_proc = zeros([Nw,2,2], dtype=complex)
    Sigma = zeros([Nw,2,2], dtype=complex)
    
    #compute new Sigma    
    change = zeros([2,2], dtype=complex)

    for m in range(myrank,Nw,nprocs): ### check this bound. it was n
        for n in range(Nw): 
            n_m = subtract_freqs(n,m,Nw)
                            
            if(n_m>=0 and n_m<Nw-1):

                Sigma_proc[n,:,:] -= 1.0/(Nk**2)/beta * g**2 * D[n_m] * einsum('ij,jk,kl->il',tau3,Gloc[m,:,:],tau3)

    
    comm.Allreduce(Sigma_proc, Sigma, op=MPI.SUM)
    
    change += sum(abs(Sigma-Sigma_old), axis=0)
                  
    if myrank==0:
        print "change ", change
        print "iteration time ",time.time() - tstart

        save("data_holstein/Gloc.npy", Gloc)
        save("data_holstein/Sigma.npy", Sigma)    

        
print "total run time ", time.time() - tstart
if myrank==0:
    save("data_holstein/Gloc.npy", Gloc)
    save("data_holstein/Sigma.npy", Sigma)
        

# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 01:44:04 2016

@author: Ben
"""

from numpy import *
from init_functions import *
from Functions import *
import time
import sys

#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#nprocs = comm.size
#myrank = comm.rank
myrank = 0

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

def parse_line(f):
    line = f.readline()
    index = line.index('#')+1
    if '.' in line:
        return float(line[index:])
    else:
        return int(float(line[index:]))
    
with open(sys.argv[1],'r') as f:
    g_dqmc = parse_line(f)
    Nk     = parse_line(f)
    Nw     = parse_line(f)
    beta   = parse_line(f)
    omega  = parse_line(f)
    mu     = parse_line(f)
    superconductivity = parse_line(f)
f.close()

#savedir = sys.argv[2]
import os
savedir = 'qInf_omega%1.1f'%omega+'_g%1.3f'%g_dqmc+'/'
if not os.path.exists(savedir):
    os.makedirs(savedir)

#g_dqmc = 0.100
#Nk    = 40
#Nw    = 200
#beta  = 2.4
#omega = 0.6

g     = g_dqmc * Nk * 1./ sqrt(2. * omega)
#superconductivity = False

print ' g_dqmc ',g_dqmc
print ' Nk     ',Nk
print ' Nw     ',Nw
print ' beta   ',beta
print ' omega  ',omega
print ' superconductivty ',superconductivity


if myrank==0:
    save(savedir+"params",asarray([Nk,Nw,beta,g,omega,0.0,superconductivity]))
    
iter_selfconsistency = 30

kxs, kys  = init_momenta(Nk)
iw_bose   = init_boson_freq(Nw, beta)
iw_fermi  = init_fermion_freq(Nw, beta)
band      = init_band(kxs, kys, Nk, mu)
D         = init_D(Nw, beta, omega, iw_bose)

#now do the same calculation but with the FFT
Gloc      = zeros([Nw,2,2], dtype=complex)
Sigma     = zeros([Nw,2,2], dtype=complex)

if superconductivity:
    Sigma[:,0,1] = 0.01*1j
    Sigma[:,1,0] = 0.01*1j

change = ones([2,2], dtype=complex)
    
#selfconsistency loop
for myiter in range(iter_selfconsistency):
    if abs(change[0,0]) < 1e-8:
        break

    print ' '
    print 'iteration ',myiter
    
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
    for m in range(Nw):
        for n in range(Nw): 
            n_m = subtract_freqs(n,m,Nw)
                            
            if n_m>=0 and n_m<Nw-1:

                Sigma_proc[n,:,:] -= 1.0/(Nk**2)/beta * g**2 * D[n_m] * einsum('ij,jk,kl->il',tau3,Gloc[m,:,:],tau3)

    Sigma = Sigma_proc
    #comm.Allreduce(Sigma_proc, Sigma, op=MPI.SUM)
    
    change = sum(abs(Sigma-Sigma_old), axis=0)
                  
    if myrank==0:
        print "change ", change
        print "iteration time ",time.time() - tstart

        save(savedir+"Gloc.npy", Gloc[:,0,0])
        save(savedir+"Sigma.npy", Sigma)    

        
print "total run time ", time.time() - tstart


import subprocess

def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

bash_command('cp '+sys.argv[1]+' '+savedir)



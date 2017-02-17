# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 01:44:04 2016

@author: Ben
"""

import subprocess
def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

bash_command('source mlpython2.7.5.sh')
    
from numpy import *
from init_functions import *
from Functions import *
import time
import subprocess
import sys
import os

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
    superconductivity = parse_line(f)
    mu = parse_line(f)
    q0     = parse_line(f)
f.close()

print ' g_dqmc ',g_dqmc
print ' Nk     ',Nk
print ' Nw     ',Nw
print ' beta   ',beta
print ' omega  ',omega
print ' superconductivty ',superconductivity
print ' q0     ',q0

params = [g_dqmc,Nk,Nw,beta,omega,superconductivity,q0]

import os
from datetime import date
today = date.today()
yr = today.timetuple()[0]
mn = today.timetuple()[1]
dy = today.timetuple()[2]
savedir = 'data_%d'%mn+'_%d'%dy+'_%d'%yr+'/fig2/q%1.1f'%q0+'_omega%1.1f'%omega+'_g%1.3f'%g_dqmc+'_mu%1.3f'%abs(mu)+'_Nw%d'%Nw+'_Nk%d'%Nk + '_beta%1.1f'%beta+'/'
print 'savedir = ',savedir
if not os.path.exists(savedir):
    os.makedirs(savedir)

def get_mu(mu, myparams):
    myrank=0

    [g_dqmc,Nk,Nw,beta,omega,superconductivity,q0] = myparams
    
    ###Nk must be odd or else the momentum points do not 
    ###form a group under addition!
    ###Choose Nw to be even (number of Fermionic frequencies)

    #correcting for chemical potential shift!


    # g_ME = g_dqmc * sqrt(N_total_DQMC) / sqrt(2 omega)

    #bash_command('source setup.sh')

    tstart = time.time()

    #savedir = sys.argv[2]

    if myrank==0:
        bash_command('cp '+sys.argv[1]+' '+savedir)

    g     = g_dqmc * Nk * 1./ sqrt(2. * omega)
    q0    = 2*pi*q0

    iter_selfconsistency = 300

    kxs, kys  = init_momenta(Nk)
    gofq      = init_gofq(kxs, kys, Nk, g, q0)
    iw_bose   = init_boson_freq(Nw, beta)
    iw_fermi  = init_fermion_freq(Nw, beta)
    band      = init_band(kxs, kys, Nk, mu)
    D         = init_D(Nw, beta, omega, iw_bose)

    #now do the same calculation but with the FFT
    #G         = init_G(Nk, Nw, beta, omega, band, kxs, kys, iw_fermi, superconductivity)
    #G         = load("data/G.npy")

    G         = zeros([Nk,Nk,Nw,2,2], dtype=complex)
    G_proc    = zeros([Nk,Nk,Nw,2,2], dtype=complex)
    Conv      = zeros([Nk,Nk,2,2], dtype=complex)
    Sigma     = zeros([Nk,Nk,Nw,2,2], dtype=complex)

    fft_gofq2 = fft.fft2(gofq**2)

    change = ones([2,2], dtype=complex)

    if superconductivity:
        Sigma[:,:,:,0,1] = 0.01*1j
        Sigma[:,:,:,1,0] = 0.01*1j

    #selfconsistency loop
    for myiter in range(iter_selfconsistency):
        if abs(change[0,0]) < 1e-4:
            break

        #compute new G
        for ik1 in range(Nk):
            for ik2 in range(Nk):

                for n in range(Nw):
                    iwn = iw_fermi[n]

                    G_proc[ik1,ik2,n,:,:] = linalg.inv(iwn*tau0 - band[ik1,ik2]*tau3 - Sigma[ik1,ik2,n,:,:])

        G = zeros([Nk,Nk,Nw,2,2], dtype=complex)
        G = G_proc
        #comm.Allreduce(G_proc, G, op=MPI.SUM)


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
                    Conv = roll(Conv, -Nk/2, axis=0)
                    Conv = roll(Conv, -Nk/2, axis=1)

                    Sigma_proc[:,:,n,:,:] -= Conv

        Sigma = Sigma_proc
        #comm.Allreduce(Sigma_proc, Sigma, op=MPI.SUM)

        change += sum(abs(Sigma-Sigma_old), axis=(0,1,2))/Nk**2
        filling =  abs(1.0 + 2.0*sum(G[:,:,:,0,0], axis=(0,1,2))/Nk**2/beta)
        
        if myrank==0:
            print " "
            print "iteration ",myiter
            print "change ", change
            print "iteration time ",time.time() - tstart
            print "filling : ",filling

    if myrank==0:
        save(savedir+"GM.npy", G)
        save(savedir+"Sigma.npy", Sigma)    
        print "total run time ", time.time() - tstart


    return filling
        

filling = 0
desired_filling = float(sys.argv[2])
mu = -2.0*(1.0-desired_filling)

print '_____________________'
print 'searching for bounds '
filling0 = get_mu(mu, params)
mu0 = mu
if filling0 > desired_filling:
    dmu = -0.2
else:
    dmu = 0.2

done = False
while not done:
    mu = mu + dmu
    print '_____________________'
    print 'searching for bounds '
    filling = get_mu(mu, params)

    if filling < desired_filling and dmu < 0:
        filling1 = filling
        mu1 = mu
        break
    if filling > desired_filling and dmu > 0:
        filling1 = filling
        mu1 = mu
        break


# choose linearly interpolated point between the two fillings

while abs(filling1-desired_filling)>0.0001 and abs(filling0-desired_filling)>0.0001:

    mu = mu0*(desired_filling-filling1)/(filling0-filling1) + mu1*(desired_filling-filling0)/(filling1-filling0)
    print '___________________'
    print 'trying mu = ',mu
    print '     mu0       mu1',mu0,mu1
    print 'filling0, filling1',filling0,filling1
    filling = get_mu(mu, params)
    
    if filling > desired_filling:
        if filling0 < filling1:
            filling1 = filling
            mu1      = mu
        else:
            filling0 = filling
            mu0      = mu   
    else:
        if filling0 < filling1:
            filling0 = filling
            mu0      = mu
        else:
            filling1 = filling
            mu1      = mu


if abs(filling1-desired_filling)<abs(filling0-desired_filling):
    print 'final filling ',filling1
    print 'final mu      ',mu1
    
    savetxt(savedir+'mu',[mu1,filling1]) 
else:
    print 'final filling ',filling0
    print 'final mu      ',mu0
    savetxt(savedir+'mu',[mu0,filling0]) 



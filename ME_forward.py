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


myrank=0

###Nk must be odd or else the momentum points do not 
###form a group under addition!
###Choose Nw to be even (number of Fermionic frequencies)

#correcting for chemical potential shift!


# g_ME = g_dqmc * sqrt(N_total_DQMC) / sqrt(2 omega)

#bash_command('source setup.sh')

tstart = time.time()

def parse_line(f):
    line = f.readline()
    index = line.index('#')
    if '.' in line:
        return float(line[:index])
    else:
        return int(float(line[:index]))
    
with open(sys.argv[1],'r') as f:
    g_dqmc = parse_line(f)
    Nk     = parse_line(f)
    Nw     = parse_line(f)
    beta   = parse_line(f)
    omega  = parse_line(f)
    superconductivity = parse_line(f)
    q0     = parse_line(f)
f.close()

#savedir = sys.argv[2]
import os
savedir = 'data_2_6_17_small_cluster/q%1.1f'%q0+'_omega%1.1f'%omega+'_g%1.3f'%g_dqmc+'/'
if not os.path.exists(savedir):
    os.makedirs(savedir)

if myrank==0:
    bash_command('cp '+sys.argv[1]+' '+savedir)
    
g     = g_dqmc * Nk * 1./ sqrt(2. * omega)

print ' g_dqmc ',g_dqmc
print ' Nk     ',Nk
print ' Nw     ',Nw
print ' beta   ',beta
print ' omega  ',omega
print ' superconductivty ',superconductivity
print ' q0     ',q0

q0    = 2*pi*q0

iter_selfconsistency = 30

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
    if abs(change[0,0]) < 1e-8:
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

    if myrank==0:
        print " "
        print "iteration ",myiter
        print "change ", change
        print "iteration time ",time.time() - tstart
        print "filling : ", 1.0 + 2.0*sum(G[:,:,:,0,0], axis=(0,1,2))/Nk**2/beta

if myrank==0:
    save(savedir+"GM.npy", G)
    save(savedir+"Sigma.npy", Sigma)    
    print "total run time ", time.time() - tstart


# copy input file into the savedir


# do this processing in post
'''
print 'now making and saving Gk'
    
def plotME_k(folder_ME, k_index):
    [Nk,Nw,beta,g,omega,q0,sc] = load(folder_ME+"params.npy")
    taus = linspace(0, beta)

    print Nk
    print Nw
    print beta
    print g
    print omega
    
    Nw = int(Nw)
    Nk = int(Nk)
    iw_fermi = zeros(Nw, dtype=complex)
    Nw2 = int(Nw/2.)
    for n in range(Nw):
        iw_fermi[n] = 1j*(2.*(n - Nw2) + 1.)*pi/beta
    ws = imag(iw_fermi)

    Ntau = 50
    taus = linspace(0.042,beta-0.042,Ntau)
    
    N_selected_k = len(k_index)
    Gtau = zeros([Ntau, N_selected_k], dtype=complex)
    G = load(folder_ME+"GM.npy")

    print 'G shape ',shape(G)
    
    iky = Nk/2
    ikxs = k_index

    for ik in range(N_selected_k):
        for itau in range(Ntau):
            for iw in range(len(iw_fermi)):
                Gtau[itau, ik] += 1./beta * exp(-1j*taus[itau]*imag(iw_fermi[iw])) * G[ikxs[ik],iky,iw,0,0]                
                #plot(taus,-Gtau)

    return taus,Gtau

if myrank==0:
    #k_index = [0,10,20,30,40]

    k_index = []
    for i in range(5):
        k_index.append(i*Nk/8)
    print "k_index = ", k_index
        
    taus,Gtau = plotME_k("data_forward/", k_index)

    print "saving new Gtau"
    save("data_forward/taus.npy", taus)
    save("data_forward/Gtau.npy", Gtau)
'''

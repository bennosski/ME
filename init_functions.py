# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 20:51:30 2016

@author: Ben
"""

from numpy import *

tau0 = array([[1.0,0.0],[0.0,1.0]])
tau1 = array([[0.0,1.0],[1.0,0.0]])
tau3 = array([[1.0,0.0],[0.0,-1.0]])

def Ek(kx, ky):
    return -2.0*0.25*(cos(kx) + cos(ky))

def init_Sigma(Nk,Nw,superconductivity):
    Sigma = zeros([Nk,Nk,Nw,2,2], dtype=complex)


    return Sigma

def init_G(Nk, Nw, beta, omega, band, kxs, kys, iw_fermi, superconductivity):
    G = zeros([Nk,Nk,Nw,2,2], dtype=complex)

    if superconductivity:
        S0 = asarray([[0., 0.01], [0.01, 0.]])
    else:
        S0 = zeros([2,2])
        
    for ik1 in range(Nk):
        for ik2 in range(Nk):
            
            for n in range(Nw): 
                iwn = iw_fermi[n]

                G[ik1,ik2,n,:,:] = linalg.inv(iwn*tau0 - band[ik1,ik2]*tau3 - S0)
                
                
    return G 
    
def init_D(Nw, beta, omega, iw_bose):
    D = zeros(Nw-1, dtype=complex)    
    for n in range(Nw-1):
        iwn = iw_bose[n]
        D[n] = 2*omega/(iwn**2 - omega**2)
    return D

def init_boson_freq(Nw, beta):
    iw_bose = zeros(Nw-1, dtype=complex)    
    Nw2 = int(Nw/2.) - 1
    for n in range(Nw-1):
        iw_bose[n] = 1j*2.*(n - Nw2)*pi/beta        
    return iw_bose

def init_fermion_freq(Nw, beta):
    iw_fermi = zeros(Nw, dtype=complex)
    Nw2 = int(Nw/2.)
    for n in range(Nw):
        iw_fermi[n] = 1j*(2.*(n - Nw2) + 1.)*pi/beta         
    return iw_fermi
    
def init_momenta(Nk):
    kxs = zeros(Nk)
    kys = zeros(Nk)
    dk = 2.*pi/Nk
    for ik in range(Nk):
        kxs[ik] = -pi + dk*ik + dk/2.
        kys[ik] = -pi + dk*ik + dk/2.
    return kxs, kys

def init_gofq(kxs, kys, Nk, g, q0):
    gofq = zeros([Nk,Nk],dtype=complex)
    for ik1 in range(Nk):
        for ik2 in range(Nk):
            kx, ky = kxs[ik1], kys[ik2] 
            gofq[ik1,ik2] = g * exp(-sqrt(kx**2 + ky**2)/q0) 
    return gofq

def init_band(kxs, kys, Nk):
    band = zeros([Nk,Nk])
    for ik1 in range(Nk):
        for ik2 in range(Nk):
            band[ik1,ik2] = Ek(kxs[ik1], kys[ik2])
    return band

def progress_bar(i, imax):
    bar = '\rprogress : |'
    L = 40
    end = int(i*1./(imax-1)*L)
    for j in range(L):
        if(j<end):
            bar += '-'
        else:
            bar += ' '
    bar += '|'
    print bar,
        

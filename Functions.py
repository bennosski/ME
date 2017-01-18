# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 22:38:24 2016

@author: Ben
"""
from numpy import *

def computeD(w, omega):
    return 2.*omega/(w**2 - omega**2)

#computes k - q
def subtract_momenta(ik1,ik2,iq1,iq2,Nk):
    Nk2 = int(Nk/2.)
    return (ik1 - iq1 + Nk2)%Nk, (ik2 - iq2 + Nk2)%Nk
    
#computes k - q
def add_momenta(ik1,ik2,iq1,iq2,Nk):
    Nk2 = int(Nk/2.)
    return (ik1 + iq1 - Nk2)%Nk, (ik2 + iq2 - Nk2)%Nk
    
#not a group
def subtract_freqs(n,m,Nw):
    Nw2 = int(Nw/2.) - 1
    return n - m + Nw2
    
def inv2b2(M):
    out = zeros([2,2],dtype=complex)
    
    
    

# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 22:38:24 2016

@author: Ben
"""
from numpy import *

def computeD(w, omega):
    return 2.*omega/(w**2 - omega**2)
    
#not a group
def subtract_freqs(n,m,Nw):
    Nw2 = int(Nw/2.) - 1
    return n - m + Nw2
    
def inv2b2(M):
    out = zeros([2,2],dtype=complex)
    
    
    

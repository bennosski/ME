from numpy import *

def iw_to_tau(F, iw_fermi):
    
    tau_start = 0.1
    Ntau = 100
    taus = linspace(tau_start,beta-tau_start,Ntau)
    
    Gtau = zeros(Ntau, dtype=complex)
    
    for itau in range(Ntau):
        for iw in range(len(iw_fermi)):
            Gtau[itau] += 1./beta * exp(-1j*taus[itau]*imag(iw_fermi[iw])) * F[iw]
    return taus, Gtau

def get_ME_Gtau_k(folder, ikx, iky):
    num_w = 10000

    GM = load(folder+'GM.npy')
    Giw0 = GM[ikx,iky,:,0,0]
	
    Giw = zeros(num_w, dtype=complex)    
    iw_fermi = zeros(num_w, dtype=complex)
    Nw2 = int(num_w/2.)
    for n in range(num_w):
        iw_fermi[n] = 1j*(2.*(n - Nw2) + 1.)*pi/beta
    Giw = 1./iw_fermi

    save(folder+"Giw", Giw)
    
    for i in range(Nw):
        Giw[num_w/2-Nw/2+i] = Giw0[i]
        
    taus, Gtau = iw_to_tau(Giw, iw_fermi)
    return taus, Gtau   

import sys
import os
folder = sys.argv[1]

files = os.listdir(folder)

for afile in files:
    if "input" in afile:
        myfile = afile
        break

def parse_line(f):
    line = f.readline()
    index = line.index('#')
    if '.' in line:
        return float(line[:index])
    else:
        return int(float(line[:index]))

with open(folder+myfile,'r') as f:
    g_dqmc = parse_line(f)
    Nk     = parse_line(f)
    Nw     = parse_line(f)
    beta   = parse_line(f)
    omega  = parse_line(f)
    superconductivity = parse_line(f)
    q0     = parse_line(f)
f.close()

print 'Nk ',Nk
print 'g_dqmc ',g_dqmc

taus, Gtau = get_ME_Gtau_k(folder, Nk/2, 0)
print folder
save(folder+'taus', taus)
save(folder+'Gtau', Gtau)

print 'done'

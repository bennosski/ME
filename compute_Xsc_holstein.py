
from numpy import *
import sys
import os
import time
from init_functions import *

t0 = time.time()

folder = sys.argv[1]
Sigma = load(folder+'Sigma.npy')
#the momentum independent Self-energy if this Sigma was calculated using the forward scattering code:
Sigma = Sigma[0,0,:,0,0]

files = os.listdir(folder)
for myfile in files:
    if "input" in myfile:
        inputfile = myfile
        break
        
def parse_line(f):
    line = f.readline()
    index = line.index('#')+1
    if '.' in line:
        return float(line[index:])
    else:
        return int(float(line[index:]))
    
with open(folder+inputfile,'r') as f:
    g_dqmc = parse_line(f)
    Nk     = parse_line(f)
    Nw     = parse_line(f)
    beta   = parse_line(f)
    omega  = parse_line(f)
    superconductivity = parse_line(f)
    mu     = parse_line(f)
    q0     = parse_line(f)
f.close()

g     = g_dqmc * Nk * 1./ sqrt(2. * omega)

print ' g_dqmc ',g_dqmc
print ' Nk     ',Nk
print ' Nw     ',Nw
print ' beta   ',beta
print ' omega  ',omega
print ' superconductivty ',superconductivity
print ' q0     ',q0
print ' mu     ',mu


q0    = 2*pi*q0

NwOriginal = Nw
Nws = [20,30,40,50]
#Nws = [24]

for Nw in Nws:

     kxs, kys = init_momenta(Nk)
     band     = init_band(kxs, kys, Nk, mu)
     iwn      = init_fermion_freq(Nw, beta)
     wn = imag(iwn)

     mid = NwOriginal/2
     Z = 1.0 - Sigma[mid-Nw/2:mid+Nw/2]/iwn

     x0 = 0.
     for ik1 in range(Nk):
         for ik2 in range(Nk):
             for n in range(Nw):
                 x0 += 1./(beta*Nk**2) * 1./(Z[n]**2*wn[n]**2 + band[ik1,ik2]**2)

     print 'done with x0 \n'


     g2D = zeros([Nw,Nw], dtype=complex)
     for n1 in range(Nw):
         for n2 in range(Nw):
             g2D[n1,n2] = -2.*g**2*omega / ((wn[n1]-wn[n2])**2+omega**2)

     t = g2D.copy()
     change = 1.0

     iter_selfconsistency = 20
     for myiter in range(iter_selfconsistency):
         if change < 1e-4:
             break

         print 'iter ',myiter

         tnew = zeros([Nw,Nw], dtype=complex)

         for ik1 in range(Nk):
             for ik2 in range(Nk):
                 for n in range(Nw):
                     for np in range(Nw):
                         for npp in range(Nw):
                             tnew[np,n] -= 1./(beta*Nk**2)*1./(Z[npp]**2*wn[npp]**2 + band[ik1,ik2]**2)*g2D[npp,np] * t[npp,n]

         tnew += g2D 
         change = sum(abs(tnew-t))/Nw**2
         print ' '
         print change
         t = tnew.copy()

         print '\ntime elapsed ',time.time()-t0

     save(folder,'t')

     x = 0.
     for ik1 in range(Nk):
         for ik2 in range(Nk):
             for ip1 in range(Nk):
                 for ip2 in range(Nk):
                     for n in range(Nw):
                         for np in range(Nw):
                             x -= 1./(beta*Nk**2)**2 * 1./(Z[np]**2*wn[np]**2 + band[ip1,ip2]**2) * t[np, n] * 1./(Z[n]**2*wn[n]**2 + band[ik1,ik2]**2)  

     x += x0

     print '---------------------'
     print 'Xsc = ', x
     print '---------------------'

     print 'total time elapsed ',time.time()-t0

     savetxt(folder+'xsc_Nw%d'%Nw, [real(x)])
    

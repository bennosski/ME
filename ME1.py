# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 01:44:04 2016

@author: Ben
"""

from numpy import *
from matplotlib.pyplot import *
from init_functions import *
from Functions import *

###Nk must be odd or else the momentum points do not 
###form a group under addition!
###Choose Nw to be even (number of Fermionic frequencies)

Nk = 5
Nw = 40
beta = 1.0
g = 0.3
omega = 0.2
q0 = 123456789.

iter_selfconsistency = 3 

kxs, kys  = init_momenta(Nk)
gofq      = init_gofq(kxs, kys, Nk, g, q0)
iw_bose   = init_boson_freq(Nw, beta)
iw_fermi  = init_fermion_freq(Nw, beta)
band      = init_band(kxs, kys, Nk)
G         = init_G(Nk, Nw, beta, omega, band, kxs, kys, iw_fermi)
D         = init_D(Nw, beta, omega, iw_bose)
Sigma     = zeros([Nk,Nk,Nw,2,2], dtype=complex)


#test momenta
Nk2 = int(Nk/2.)
print Nk2
print "should be zero"
print kxs[Nk2],kys[Nk2]
ik1,ik2 = add_momenta(Nk2+1,Nk2,Nk2+2,Nk2,Nk)
print kxs[ik1], kxs[ik2]
ik1,ik2 = subtract_momenta(Nk2,Nk2,Nk2+1,Nk2,Nk)
print kxs[ik1], kxs[ik2]

print subtract_freqs(2,2,Nw)

#assert 1==0

#sum D(q) * g(q) * G(k - q)

#selfconsistency loop
for myiter in range(iter_selfconsistency):

    Sigma_old = Sigma.copy()
    Sigma = zeros([Nk,Nk,Nw,2,2], dtype=complex)

    #compute new Sigma    
    change = zeros([2,2], dtype=complex)
    for ik1 in range(Nk):
        progress_bar(ik1,Nk)
        
        kx = kxs[ik1]
        for ik2 in range(Nk):
            ky = kys[ik2]

            for iq1 in range(Nk):
                qx = kxs[iq1]
                for iq2 in range(Nk):
                    qy = kys[iq2]
            
                    #p = k - q
                    #ip1, ip2 = subtract_momenta(ik1,ik2,iq1,iq2,Nk)
                    ip1, ip2 = add_momenta(ik1,ik2,iq1,iq2,Nk)
                        
                    for n in range(Nw): 
                        for m in range(n):
                            n_m = subtract_freqs(n,m,Nw)
                            
                            if(n_m>=0 and n_m<Nw):
                                Sigma[ik1,ik2,n,:,:] -= gofq[iq1,iq2]**2/(Nk**2)/beta * D[n-m] * \
                                                        einsum('ij,jk,kl->il',tau3,G[ip1,ip2,m,:,:],tau3) 
                    
    #change += 1./(Nk**2*Nw)*abs(Sigma[ik1,ik2,n,:,:]-Sigma_old[ik1,ik2,n,:,:])
    change += sum(abs(Sigma-Sigma_old), axis=(0,1,2))
                  
    #compute new G
    for ik1 in range(Nk):
        kx = kxs[ik1]
        for ik2 in range(Nk):
            ky = kys[ik2]
            
            for n in range(Nw):
                iwn = iw_fermi[n]
                
                G[ik1,ik2,n,:,:] = linalg.inv((iwn*tau0 - band[ik1,ik2]*tau3 - Sigma[ik1,ik2,n,:,:]))

    print "change ", change


figure()
plot(real(D))
plot(imag(D))


equaltime = sum(G,axis=2)

figure()
imshow(real(equaltime[:,:,0,0]))
colorbar()
figure()
imshow(imag(equaltime[:,:,0,0]))
colorbar()

figure()
imshow(real(equaltime[:,:,1,1]))
colorbar()
figure()
imshow(imag(equaltime[:,:,1,1]))
colorbar()







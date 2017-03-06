from numpy import *

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

k_index = [0,10,20,30,40]
taus,Gtau = plotME_k("data_forward/", k_index)

print "saving new Gtau"
save("data_forward/taus.npy", taus)
save("data_forward/Gtau.npy", Gtau)

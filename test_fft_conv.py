
from numpy import *



x = random.random((5,5))
y = random.random((5,5))

c = zeros([5,5])

for ipx in range(5):
    for ipy in range(5):

        for iqx in range(5):
            for iqy in range(5):

                kx = (ipx - (iqx - 2) )%5
                ky = (ipy - (iqy - 2) )%5

                #kx = (ipx - iqx )%5
                #ky = (ipy - iqy )%5
                
                
                c[ipx,ipy] += x[iqx,iqy]*y[kx, ky]

print c


#xr = roll(x,2,axis=1)
#xr = roll(x,2,axis=0)
#yr = roll(y,2,axis=1)
#yr = roll(y,2,axis=0)

fft1 = fft.fft2(x)
fft2 = fft.fft2(y)

out = real(fft.ifft2(fft1*fft2))

out = roll(out,-(5-1)/2,axis=0)
out = roll(out,-(5-1)/2,axis=1)

print " "
print out

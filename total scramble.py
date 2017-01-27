import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import ascii
from scipy import random as ra

bins=np.arange(0, 91, 1)

#isotropic data through randomization of RA
totcounts=0.

start=time.time()
for i in range(0,100):
    rRA=ra.random(100)*360
    rDec=ra.random(100)*180-90
    recoord=coord.SkyCoord(ra=rRA, dec=rDec, unit='deg', frame='icrs')
    rblat=abs(recoord.galactic.b.deg)
    (lcounts,ignored)=np.histogram(rblat,bins)
    totcounts=totcounts+lcounts

end=time.time()
print('Time =',end-start)

meancounts=totcounts/(i+1)
meancounts=19*meancounts/sum(meancounts)
#print(meancounts)
#print('b<10 =',sum(meancounts[:9]),'b>50 =',sum(meancounts[50:]))

data = ascii.read('events.txt')

#classification
energy=data['col2']
declin=data['col6']
righta=data['col7']

#signal data
cuts=energy >100 #higher than 100TeV
ecoord=coord.SkyCoord(ra=righta[cuts], dec=declin[cuts], unit='deg', frame='icrs')
blat=abs(ecoord.galactic.b.deg) #b latitude from galactic coord
nrevents=righta[cuts].size

bins=np.arange(0, 91, 1)

#isotropic data through geometric interpretation
teta=(90-bins)*u.deg
st=2*np.pi*(np.cos(teta[1:])-np.cos(teta[:-1])) #solidangle in str for every bin
frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
geocounts=frac*nrevents

#plotting the isotropic random and geometric dist with the hist
x1=bins[0],bins[1]
y1=meancounts[0],meancounts[0]
y2=geocounts[0],geocounts[0]
for k in range(1,bins.size-1):
    x1=np.hstack((x1,(bins[k],bins[k+1])))
    y1=np.hstack((y1,(meancounts[k],meancounts[k])))
    y2=np.hstack((y2,(geocounts[k],geocounts[k])))

plt.hist(blat,bins)
plt.plot(x1,y1,lw=3,label='iso ran dist')
plt.plot(x1,y2,label='iso geo dist',color='r')
plt.xlabel('b')
plt.ylabel('counts')
plt.legend()
#plt.savefig('total scramble.pdf')
plt.show()

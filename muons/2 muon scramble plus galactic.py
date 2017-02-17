import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import ascii
from scipy import random as ra

mdata = ascii.read('events2.txt')

#classification
msprob=mdata['col9']
menergy=mdata['col8']
mdeclin=mdata['col4']
mrighta=mdata['col3']

energy=menergy
declin=mdeclin
righta=mrighta

#signal data
cuts=msprob>0.0
ecoord=coord.SkyCoord(ra=righta[cuts], dec=declin[cuts], unit='deg', frame='icrs')
blat=abs(ecoord.galactic.b.deg) #b latitude from galactic coord
nrevents=righta[cuts].size

bins=np.arange(0, 100, 10)

#isotropic data through geometric interpretation
teta=(90-bins)*u.deg
st=2*np.pi*(np.cos(teta[1:])-np.cos(teta[:-1])) #solidangle in str for every bin
frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
geocounts=frac*nrevents

#galactic template lifted from histogram
galcounts=np.array([7.54, 2.12, 0.89, 0.64, 0.39, 0.24, 0.07, 0.02, 0.])*2
galcounts=galcounts/sum(galcounts)*nrevents

#bestcounts=geocounts+galcounts

#isotropic data through scramble of events, randomize RA and smur dec
totcounts=0.
Dec=declin[cuts] #not so random
'''
start=time.time()
for i in range(0,100):     #1k=30s, 10k=5min, 100k=50min
    rRA=ra.random(nrevents)*360
    recoord=coord.SkyCoord(ra=rRA, dec=Dec, unit='deg', frame='icrs')
    rblat=abs(recoord.galactic.b.deg)
    (lcounts,ignored)=np.histogram(rblat,bins)
    totcounts=totcounts+lcounts

end=time.time()
print('Time =',end-start)

meancounts=totcounts/(i+1)
'''
meancounts=np.array([3.0526, 2.8865, 2.7748, 3.0953, 3.1701, 2.9682, 2.1178, 0.7268, 0.2079])

#plotting the isotropic random and geometric dist with the hist
x1=bins[0],bins[1]
y1=meancounts[0],meancounts[0]
y2=geocounts[0],geocounts[0]
y3=galcounts[0],galcounts[0]
#y4=bestcounts[0],bestcounts[0]
for k in range(1,bins.size-1):
    x1=np.hstack((x1,(bins[k],bins[k+1])))
    y1=np.hstack((y1,(meancounts[k],meancounts[k])))
    y2=np.hstack((y2,(geocounts[k],geocounts[k])))
    y3=np.hstack((y3,(galcounts[k],galcounts[k])))
#    y4=np.hstack((y4,(bestcounts[k],bestcounts[k])))

plt.hist(blat,bins)
plt.plot(x1,y1,lw=3,label='Scrambled signal')
plt.plot(x1,y2,label='Isotropic model',color='r')
plt.plot(x1,y3,label='Galactic model',color='c')
#plt.plot(x1,y4,label='total',lw=3,color='y')
plt.xlabel('$|b|$ (deg)')
plt.ylabel('Counts')
plt.legend()
#plt.savefig('smurxdeg.pdf')
plt.show()

#10000 --> [3.0526, 2.8865, 2.7748, 3.0953, 3.1701, 2.9682, 2.1178, 0.7268, 0.2079]
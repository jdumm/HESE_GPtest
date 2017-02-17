import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import ascii
from scipy import random as ra

#data = ascii.read('events_4year.txt')
data = ascii.read('events_6year.txt')
#data = ascii.read('events_new2.txt')

# Generate b (gal lat) distribution from text file and compare to scrambled and geometric models.

#classification
energy=data['col2']
declin=data['col6']
righta=data['col7']

#signal data
#cuts=energy > 100 #higher than 100TeV
cuts=energy > 96 #higher than 100TeV but adjusted downwards for energy recalibration

ecoord=coord.SkyCoord(ra=righta[cuts], dec=declin[cuts], unit='deg', frame='icrs')
blat=abs(ecoord.galactic.b.deg) #b latitude from galactic coord
nrevents=righta[cuts].size

bins=np.arange(0, 100, 10)
#bins=np.arange(0, 91, 1) #for future choice in bining

#isotropic data through geometric interpretation
teta=(90-bins)*u.deg
st=2*np.pi*(np.cos(teta[1:])-np.cos(teta[:-1])) #solidangle in str for every bin
frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
geocounts=frac*nrevents
#isotropic data through scramble of events, randomize RA and smur dec
totcounts=0.
Dec=declin[cuts] #not so random

start=time.time()
for i in range(0,100):     #1k=30s, 10k=5min, 100k=50min
    rRA=ra.random(nrevents)*360
    rDec_delta=ra.random(nrevents)*30-15
    rDec=Dec+rDec_delta
    rDec[rDec<-90]=Dec[rDec<-90]-rDec_delta[rDec<-90]
    rDec[rDec>90]=Dec[rDec>90]-rDec_delta[rDec>90]
    recoord=coord.SkyCoord(ra=rRA, dec=rDec, unit='deg', frame='icrs')
    rblat=abs(recoord.galactic.b.deg)
    (lcounts,ignored)=np.histogram(rblat,bins)
    totcounts=totcounts+lcounts

end=time.time()
print('Time =',end-start)

meancounts=totcounts/(i+1)
'''
'''
#meancounts=[3.2575, 3.3278, 3.3453, 2.8473, 2.2266, 1.6354, 1.2094, 0.8414, 0.3093]

#plotting the isotropic random and geometric dist with the hist
x1=bins[0],bins[1]
y1=meancounts[0],meancounts[0]
y2=geocounts[0],geocounts[0]
for k in range(1,bins.size-1):
    x1=np.hstack((x1,(bins[k],bins[k+1])))
    y1=np.hstack((y1,(meancounts[k],meancounts[k])))
    y2=np.hstack((y2,(geocounts[k],geocounts[k])))

plt.hist(blat,bins)
plt.plot(x1,y1,lw=3,label='Scrambled signal')
plt.plot(x1,y2,label='Geometric model',color='r')
plt.xlabel('$|b|$ (deg)')
plt.ylabel('Counts')
plt.legend()
plt.ylim(0,11)
#plt.savefig('2_scramble_with_smur_4year.pdf')
#plt.savefig('2_scramble_with_smur_6year.pdf')

print('b<10 =',meancounts[0],'b>50 =',sum(meancounts[5:]))
#print('b<10 =',sum(meancounts[:9]),'b>50 =',sum(meancounts[50:]))

#10000 15deg smur--> [3.2575, 3.3278, 3.3453, 2.8473, 2.2266, 1.6354, 1.2094, 0.8414, 0.3093]

plt.show()

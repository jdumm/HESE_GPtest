import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import ascii
from scipy import random as ra

# Generate a distribution of maximum likelihoods for scrambled trials. 

#data = ascii.read('events_4year.txt')
data = ascii.read('events_6year.txt')


#classification
energy=data['col2']
declin=data['col6']
righta=data['col7']

#signal data
#cuts=energy >100 #higher than 100TeV
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
geocounts=frac*nrevents*0.7

#galactic template lifted from histogram
galcounts=np.array([7.54, 2.12, 0.89, 0.64, 0.39, 0.24, 0.07, 0.02, 0.])*2
galcounts=galcounts/sum(galcounts)*nrevents*0.3

bestcounts=geocounts+galcounts
'''
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
#plotting the isotropic random and geometric dist with the hist
x1=bins[0],bins[1]
#y1=meancounts[0],meancounts[0]
y2=geocounts[0],geocounts[0]
y3=galcounts[0],galcounts[0]
y4=bestcounts[0],bestcounts[0]
for k in range(1,bins.size-1):
    x1=np.hstack((x1,(bins[k],bins[k+1])))
#    y1=np.hstack((y1,(meancounts[k],meancounts[k])))
    y2=np.hstack((y2,(geocounts[k],geocounts[k])))
    y3=np.hstack((y3,(galcounts[k],galcounts[k])))
    y4=np.hstack((y4,(bestcounts[k],bestcounts[k])))

plt.hist(blat,bins)
#plt.plot(x1,y1,lw=3,label='iso ran dist')
plt.plot(x1,y4,label='Total',lw=3,color='y')
plt.plot(x1,y2,label='Isotropic 70%',color='r')
plt.plot(x1,y3,label='Galactic 30%',color='c')
plt.ylim(0,11)
plt.xlabel('$|b|$ (deg)')
plt.ylabel('Counts')
plt.legend()
#plt.savefig('7_scramble_with_smur_plus_galactic_4year.pdf')
#plt.savefig('7_scramble_with_smur_plus_galactic_6year.pdf')
plt.show()

#print('b<10 =',meancounts[0],'b>50 =',sum(meancounts[5:]))
#print('b<10 =',sum(meancounts[:9]),'b>50 =',sum(meancounts[50:]))

#10000 15deg smur--> [3.2575, 3.3278, 3.3453, 2.8473, 2.2266, 1.6354, 1.2094, 0.8414, 0.3093]

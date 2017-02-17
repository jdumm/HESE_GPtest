import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import ascii
from scipy.stats import poisson
from scipy import random as ra

# Calculates the Poisson likelihood of the data for an isotropic (geometric) flux hypothesis

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
(counts,ignored)=np.histogram(blat,bins)

#isotropic data through geometric interpretation
teta=(90-bins)*u.deg
st=2*np.pi*(np.cos(teta[1:])-np.cos(teta[:-1])) #solidangle in str for every bin
frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
geocounts=frac*nrevents

pbins=poisson.pmf(counts,geocounts) #P for exact nr counts in every bin

mylikelihood=sum(np.log(pbins))*-1


#isotropic data through scramble of events, randomize RA and smur dec
Dec=declin[cuts] #not so random
likelihood=0.

start=time.time()
ntrials = 1000
for i in range(0,ntrials):
    rRA=ra.random(nrevents)*360
    rDec_delta=ra.random(nrevents)*30-15
    rDec=Dec+rDec_delta
    rDec[rDec<-90]=Dec[rDec<-90]-rDec_delta[rDec<-90]
    rDec[rDec>90]=Dec[rDec>90]-rDec_delta[rDec>90]
    recoord=coord.SkyCoord(ra=rRA, dec=rDec, unit='deg', frame='icrs')
    rblat=abs(recoord.galactic.b.deg)
    (lcounts,ignored)=np.histogram(rblat,bins)
    lpbins=poisson.pmf(lcounts,geocounts) #P for exact nr counts in every bin
    llikelihood=sum(np.log(lpbins))*-1
    likelihood=np.hstack((likelihood,llikelihood))
end=time.time()
print('Time =',end-start)

plt.hist(likelihood[1:],np.arange(min(likelihood[1:]),max(likelihood),0.2))
plt.xlabel('-log L')
#plt.savefig('5_P_likelihood.pdf')
probability=sum(likelihood>mylikelihood)/float(ntrials)

print('-log likelihood =',mylikelihood,'probability =',probability)

plt.show()

#1000 --> -log likelihood = 16.7107540235 probability = 0.076
#10000 --> -log likelihood = 16.7107540235 probability = 0.0674

#1000 for 4year: -log likelihood = 16.71075402354171, probability = 0.0599
#1000 for 6year: -log likelihood = 16.83599093861659, probability = 0.307


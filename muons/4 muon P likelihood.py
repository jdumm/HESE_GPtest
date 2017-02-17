import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import ascii
from scipy.stats import poisson
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
(counts,ignored)=np.histogram(blat,bins)

#isotropic data through geometric interpretation
teta=(90-bins)*u.deg
st=2*np.pi*(np.cos(teta[1:])-np.cos(teta[:-1])) #solidangle in str for every bin
frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
geocounts=frac*nrevents

pbins=poisson.pmf(counts,geocounts) #P for exact nr counts in every bin

mylikelihood=sum(np.log(pbins))*-1

#isotropic data through scramble of events, randomize RA
Dec=declin[cuts] #not so random
likelihood=0.

start=time.time()
for i in range(0,100):
    rRA=ra.random(nrevents)*360
    recoord=coord.SkyCoord(ra=rRA, dec=Dec, unit='deg', frame='icrs')
    rblat=abs(recoord.galactic.b.deg)
    (lcounts,ignored)=np.histogram(rblat,bins)
    lpbins=poisson.pmf(lcounts,geocounts) #P for exact nr counts in every bin
    llikelihood=sum(np.log(lpbins))*-1
    likelihood=np.hstack((likelihood,llikelihood))
end=time.time()
print('Time =',end-start)

plt.hist(likelihood[1:],np.arange(min(likelihood[1:]),max(likelihood),0.2))
plt.xlabel('-log L')
#plt.savefig('muon p likelihood.pdf')
plt.show()
probability=sum(likelihood>mylikelihood)/(i+1)

print('-log likelihood =',mylikelihood,'probability =',probability)

#1000 --> -log likelihood = 13.3274731832 probability = 0.81
#10000 --> probability = 0.83

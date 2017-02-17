import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import ascii
from scipy.stats import poisson
from scipy import random as ra

# Scans for the max likelihood fit to an isotropic + galactic components and reports log-likelihood ratio.

data = ascii.read('events_4year.txt')
#data = ascii.read('events_6year.txt')


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

#galactic template lifted from histogram
galcounts=np.array([7.54, 2.12, 0.89, 0.64, 0.39, 0.24, 0.07, 0.02, 0.])*2
galcounts=galcounts/sum(galcounts)*nrevents

maxptot=0.
nsteps = 21
for w in np.linspace(0,1,nsteps):
    #w=w/float(10.0)
    combcounts=w*geocounts+(1-w)*galcounts
    pbins=poisson.pmf(counts,combcounts)
    ptot=np.prod(pbins)
    if ptot > maxptot:
        wmin=w
        maxptot=ptot

combcounts=wmin*geocounts+(1-wmin)*galcounts
pbins=poisson.pmf(counts,combcounts) #P for exact nr counts in every bin
likelihoodfit=sum(np.log(pbins))*-1

lratio=2*(mylikelihood-likelihoodfit)

print('best fit for ',wmin,'geometric and ',1-wmin,'galactic')
print('best fit -log likeihood =',likelihoodfit,'and only iso -log likelihood =',mylikelihood,'gives a likelihood ratio of ',lratio)

#best fit for  0.3 geometric and  0.7 galactic
#best fit -log likeihood = 11.5 and only iso -log likelihood = 16.7
#gives a likelihood ratio of  10.4

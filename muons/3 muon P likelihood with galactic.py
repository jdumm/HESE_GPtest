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

#galactic template lifted from histogram
galcounts=np.array([7.54, 2.12, 0.89, 0.64, 0.39, 0.24, 0.07, 0.02, 0.01])*2
galcounts=galcounts/sum(galcounts)*nrevents

maxptot=0.
for w in range(0,11):
    w=w/10
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

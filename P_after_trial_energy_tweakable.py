import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import ascii
from scipy.stats import poisson
from scipy import random as ra

#data = ascii.read('events_4year.txt')
data = ascii.read('events_6year.txt')

#classification
energy=data['col2']
declin=data['col6']
righta=data['col7']

#signal data
#cuts=energy >100 #higher than 100TeV
cuts=energy > 96 #higher than 100TeV but adjusted downwards for energy recalibration
nrevents=righta[cuts].size
bins=np.arange(0, 100, 10)

ecoord=coord.SkyCoord(ra=righta[cuts], dec=declin[cuts], unit='deg', frame='icrs')
blat=abs(ecoord.galactic.b.deg) #b latitude from galactic coord

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
for w in range(0,11):
    w=w/float(10.0)
    combcounts=w*geocounts+(1-w)*galcounts
    pbins=poisson.pmf(counts,combcounts)
    ptot=np.prod(pbins)
    if ptot > maxptot:
        wmin=w
        maxptot=ptot

combcounts=wmin*geocounts+(1-wmin)*galcounts
pbins=poisson.pmf(counts,combcounts) #P for exact nr counts in every bin
likelihoodfit=sum(np.log(pbins))*-1

mylratio=2*(mylikelihood-likelihoodfit)

print('best fit for ',wmin,'geometric and ',1-wmin,'galactic')

#isotropic data through scramble of events, randomize RA and smur dec
Dec=declin[cuts] #not so random
lratio=0.

start=time.time()
for i in range(0,10000): #1k=½min, 10k=5min, 100k=50min
    rRA=ra.random(nrevents)*360
    rDec_delta=ra.random(nrevents)*30-15
    rDec=Dec+rDec_delta
    rDec[rDec<-90]=Dec[rDec<-90]-rDec_delta[rDec<-90]
    rDec[rDec>90]=Dec[rDec>90]-rDec_delta[rDec>90]
    recoord=coord.SkyCoord(ra=rRA, dec=rDec, unit='deg', frame='icrs')
    rblat=abs(recoord.galactic.b.deg)
    (lcounts,ignored)=np.histogram(rblat,bins)

    maxptot=0.
    for w in range(0,11):
        w=w/float(10.0)
        combcounts=w*geocounts+(1-w)*galcounts
        pbins=poisson.pmf(lcounts,combcounts)
        ptot=np.prod(pbins)
        if ptot > maxptot:
            wmin=w
            maxptot=ptot

    combcounts=wmin*geocounts+(1-wmin)*galcounts
    pbinsfit=poisson.pmf(lcounts,combcounts) #P for exact nr counts in every bin
    likelihoodfit=sum(np.log(pbinsfit))*-1

    pbinsgeo=poisson.pmf(lcounts,geocounts) #P for exact nr counts in every bin
    likelihoodgeo=sum(np.log(pbinsgeo))*-1

    llratio=2*(likelihoodgeo-likelihoodfit)
    lratio=np.hstack((lratio,llratio))
end=time.time()
print('Time =',end-start)

probability=sum(lratio[1:] > mylratio)/float(i+1)
print('likelihood ratio = ',mylratio,' with a probability of =',probability)

#4year result, 10k trials:
#('likelihood ratio = ', 10.385067025151599, ' with a probability of =', 0.00050000000000000001)
#6year result, 10k trials:
#('likelihood ratio = ', 4.4314130215871259, ' with a probability of =', 0.0124)



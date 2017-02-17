import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import ascii
from scipy.stats import poisson
from scipy import random as ra

# Generates a the distribution of maximized likelihoods for scrambled psuedo-experiments.
# Requires user to specify 'mylratio' to calculate a probability based on the distribution.


data = ascii.read('events_4year.txt')
#data = ascii.read('events_6year.txt')

#classification
energy=data['col2']
declin=data['col6']
righta=data['col7']

#signal data
cuts=energy >100 #higher than 100TeV
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

#isotropic data through scramble of events, randomize RA and smur dec
Dec=declin[cuts] #not so random
lratio=0.

start=time.time()
for i in range(0,100):
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

plt.hist(lratio[1:])
plt.xlabel('Likelihood ratio')
#plt.savefig('ratio from scramble and best fit.pdf')

mylratio=5.10827823029 #The likelihoodratio from best fit model

probability=sum(lratio[1:] > mylratio)/float(i+1)
print('likelihood ratio = ',mylratio,' with a probability of =',probability)

plt.show()

#ascii.write(lratio[1:,np.newaxis],output='likelihood ratio.txt')
#test=np.loadtxt('likelihood ratio.txt',skiprows=1)

#10000 --> likelihood ratio = 10.3850670252 with a probability of = 0.0009

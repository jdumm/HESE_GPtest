import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import ascii
from scipy import random as ra

data = ascii.read('events.txt')

#classification
energy=data['col2']
d_energyplus=data['col4']
d_energyminus=data['col3']
declin=data['col6']
righta=data['col7']

#signal data
cuts=energy >100 #higher than 100TeV
coordinates=coord.SkyCoord(ra=righta[cuts], dec=declin[cuts], unit='deg', frame='icrs')
bvalues=coordinates.galactic.b.deg #b values in galactic coord
bbins=np.arange(0, 100, 10)
plt.hist(abs(bvalues),bbins)
plt.xlabel('b')
plt.ylabel('counts')

#isotropic data through randomization of RA
sumb=0
nrevents=righta[cuts].shape

for i in range(0,50):
    fakeRA=ra.random(nrevents)*360
    fakeCo=coord.SkyCoord(ra=fakeRA, dec=declin[cuts], unit='deg', frame='icrs')
    fakeb=fakeCo.galactic.b.deg
    (counts,bins)=np.histogram(abs(fakeb),bbins)
    sumb=sumb+counts

medelcounts=sumb/(i+1)
#plt.plot(0.5*(bins[:-1]+bins[1:]),medelcounts, label='iso ran dist')

#isotropic data through geometric interpretation
teta=(90-bbins)*u.deg
st=2*np.pi*(np.cos(teta[1:])-np.cos(teta[:-1])) #solidangle in str for every bin
frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
geocounts=frac*nrevents
#plt.plot(0.5*(bins[:-1]+bins[1:]),geocounts, label='iso geo dist')

#plotting the isotropic random and geometric dist
x1=bbins[0],bbins[1]
y1=medelcounts[0],medelcounts[0]
y2=geocounts[0],geocounts[0]
for i in range(1,9):
    x1=np.hstack((x1,(bbins[i],bbins[i+1])))
    y1=np.hstack((y1,(medelcounts[i],medelcounts[i])))
    y2=np.hstack((y2,(geocounts[i],geocounts[i])))
plt.plot(x1,y1,label='iso ran dist')
plt.plot(x1,y2,label='iso geo dist')

plt.legend()

#errorbars
(dcounts,ignored)=np.histogram(abs(bvalues),bbins)
poissondist=ra.poisson(lam=dcounts , size=(100,dcounts.size) )
errorbar=np.std(poissondist,axis=0)
plt.errorbar(0.5*(bins[:-1]+bins[1:]),dcounts,yerr=errorbar,ls='none',ecolor='black')


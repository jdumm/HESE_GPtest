import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import ascii
from scipy import random as ra

# Draw the b distribution with alternative binnings and try another energy cut.

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

bins=np.arange(0, 93, 3)
#bins=np.arange(0, 91, 1) #for future choice in bining

#isotropic data through geometric interpretation
teta=(90-bins)*u.deg
st=2*np.pi*(np.cos(teta[1:])-np.cos(teta[:-1])) #solidangle in str for every bin
frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
geocounts=frac*nrevents

#plotting the isotropic random and geometric dist with the hist
x1=bins[0],bins[1]
y2=geocounts[0],geocounts[0]
for k in range(1,bins.size-1):
    x1=np.hstack((x1,(bins[k],bins[k+1])))
    y2=np.hstack((y2,(geocounts[k],geocounts[k])))

plt.subplot(221)
plt.title('$3^{\circ}$bins')
plt.hist(blat,bins)
plt.ylim(0,5)
plt.plot(x1,y2,color='r')
plt.xlabel('$|b|$ (deg)')
plt.ylabel('Counts')
plt.xticks(bins[::3])
plt.yticks([0,1,2,3])

bins=np.arange(0, 105, 15)
#bins=np.arange(0, 91, 1) #for future choice in bining

#isotropic data through geometric interpretation
teta=(90-bins)*u.deg
st=2*np.pi*(np.cos(teta[1:])-np.cos(teta[:-1])) #solidangle in str for every bin
frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
geocounts=frac*nrevents

#plotting the isotropic random and geometric dist with the hist
x1=bins[0],bins[1]
y2=geocounts[0],geocounts[0]
for k in range(1,bins.size-1):
    x1=np.hstack((x1,(bins[k],bins[k+1])))
    y2=np.hstack((y2,(geocounts[k],geocounts[k])))

plt.subplot(222)
plt.title('$15^{\circ}$bins')
plt.hist(blat,bins)
plt.ylim(0,12)
plt.plot(x1,y2,color='r')
plt.xlabel('$|b|$ (deg)')
plt.ylabel('Counts')
plt.xticks(bins)

bins=np.arange(0, 96, 6)
#bins=np.arange(0, 91, 1) #for future choice in bining

#isotropic data through geometric interpretation
teta=(90-bins)*u.deg
st=2*np.pi*(np.cos(teta[1:])-np.cos(teta[:-1])) #solidangle in str for every bin
frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
geocounts=frac*nrevents

#plotting the isotropic random and geometric dist with the hist
x1=bins[0],bins[1]
y2=geocounts[0],geocounts[0]
for k in range(1,bins.size-1):
    x1=np.hstack((x1,(bins[k],bins[k+1])))
    y2=np.hstack((y2,(geocounts[k],geocounts[k])))

plt.subplot(223)
plt.title('$6^{\circ}$bins')
plt.hist(blat,bins)
plt.ylim(0,8)
plt.plot(x1,y2,color='r')
plt.xlabel('$|b|$ (deg)')
plt.ylabel('Counts')
plt.xticks(bins[::2])

#signal data
cuts=energy >90 #higher than 100TeV
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

#plotting the isotropic random and geometric dist with the hist
x1=bins[0],bins[1]
y2=geocounts[0],geocounts[0]
for k in range(1,bins.size-1):
    x1=np.hstack((x1,(bins[k],bins[k+1])))
    y2=np.hstack((y2,(geocounts[k],geocounts[k])))

plt.subplot(224)
plt.title('Energy > 90TeV')
plt.hist(blat,bins)
plt.ylim(0,12)
plt.plot(x1,y2,color='r')
plt.xlabel('$|b|$ (deg)')
plt.ylabel('Counts')
plt.xticks(bins)
plt.yticks([0,2,4,6,8])

plt.tight_layout()
plt.savefig('2_4_diff_hist_4year.pdf')
#plt.savefig('2_4_diff_hist_6year.pdf')
plt.show()

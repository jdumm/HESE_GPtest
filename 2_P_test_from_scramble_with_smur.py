import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import ascii
from scipy import random as ra

# Simple counting experiment and p-value estimate from scrambled trials.

data = ascii.read('events_4year.txt')
#data = ascii.read('events_6year.txt')

#classification
energy=data['col2']
declin=data['col6']
righta=data['col7']

#signal data
#cuts=energy >100 #higher than 100TeV
cuts=energy > 96 #higher than 100TeV but adjusted downwards for energy recalibration

#isotropic data through scramble of events, randomize RA and smur dec
Dec=declin[cuts] #not so random
events=0

start=time.time()
for i in range(0,100):     #1k=30s, 10k=5min, 100k=50min
    rRA=ra.random(19)*360
    rDec_delta=ra.random(19)*30-15
    rDec=Dec+rDec_delta
    rDec[rDec<-90]=Dec[rDec<-90]-rDec_delta[rDec<-90]
    rDec[rDec>90]=Dec[rDec>90]-rDec_delta[rDec>90]
    recoord=coord.SkyCoord(ra=rRA, dec=rDec, unit='deg', frame='icrs')
    rblat=abs(recoord.galactic.b.deg)
    if sum(rblat<10) >= 9:
        if sum(rblat>50) == 0:
            events=events+1
            end=time.time()
            print('event found after ',end-start)

end=time.time()
print('Time =',end-start)

print('counted events =',events)

#10000 --> low 21 high 104
#10000 --> combined events 1

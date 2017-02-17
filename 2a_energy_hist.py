import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import ascii
from scipy import random as ra

# Draw the energy distribution

#data = ascii.read('events_4year.txt')
data = ascii.read('events_6year.txt')

#classification
energy=data['col2']
declin=data['col6']
righta=data['col7']

#signal data
bins=np.arange(1, 3.6, 0.5/3.)
plt.hist(np.log10(energy), bins)
plt.xlabel(r'log$_{10} E$ / TeV')
plt.ylabel('Counts')
plt.xlim(1.0,3.5)
plt.ylim(0,21)
#plt.savefig('2a_energy_hist.pdf')

plt.show()


#print('b<10 =',meancounts[0],'b>50 =',sum(meancounts[5:]))
#print('b<10 =',sum(meancounts[:9]),'b>50 =',sum(meancounts[50:]))

#10000 15deg smur--> [3.2575, 3.3278, 3.3453, 2.8473, 2.2266, 1.6354, 1.2094, 0.8414, 0.3093]

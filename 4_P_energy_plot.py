import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import ascii
from scipy.stats import poisson
from scipy.stats import combine_pvalues

# Scan for optimal energy cut for a fixed b_low and b_high.

#data = ascii.read('events_4year.txt')
data = ascii.read('events_6year.txt')


#classification
energy=data['col2']
declin=data['col6']
righta=data['col7']

blow=6.2
bhigh=47.1

pcomb=0.

for e in range(30,300,5):
    cuts=energy > e #higher than scan point
    ecoord=coord.SkyCoord(ra=righta[cuts], dec=declin[cuts], unit='deg', frame='icrs')
    blat=abs(ecoord.galactic.b.deg) #b latitude from galactic coord
    nrevents=righta[cuts].size

    #counts of events for b<10 from geometric model
    st10=2*np.pi*(np.cos((90-blow)*u.deg)-np.cos(90*u.deg)) #solidangle in str for b<10
    frac10=2*st10/(4*np.pi) #fraction of 2*str over 4pi
    geocounts10=frac10*nrevents

    if sum(blat<blow) == 0:
        plowb=1 #P for 0+ counts in bin 10
    else:
        plowb=1-sum(poisson.pmf(np.arange(0,sum(blat<blow)),geocounts10)) #P for nr+ counts in bin 10

    #counts of events for b>50 from geometric model
    st50_90=2*np.pi*(np.cos(0*u.deg)-np.cos((90-bhigh)*u.deg)) #solidangle in str for 50<b<90
    frac50_90=2*st50_90/(4*np.pi) #fraction of 2*str over 4pi
    geocounts50_90=frac50_90*nrevents

    if sum(blat>bhigh) == 0:
        phighb=poisson.pmf(0,geocounts50_90) #P for 0 counts in bin 50-90 deg
    else:
        phighb=sum(poisson.pmf(np.arange(0,sum(blat>bhigh)+1),geocounts50_90)) #P for nr counts in bin 50-90 deg

    lpcomb=combine_pvalues([plowb,phighb])[1]
    pcomb=np.hstack((pcomb,lpcomb))

plt.plot(np.arange(30,300,5),pcomb[1:],color='b')
plt.hlines((0.0455,0.0027,0.00006),0,300,linestyle='dashed')
plt.yscale('log')
plt.text(280,0.06,r'$2\sigma$',fontsize=13)
plt.text(280,0.003,r'$3\sigma$',fontsize=13)
plt.text(280,0.00007,r'$4\sigma$',fontsize=13)
plt.xlabel('Energy cut (TeV)')
plt.ylabel('P')
#plt.savefig('4_P_energy_plot_4year.pdf')
plt.savefig('4_P_energy_plot_6year.pdf')

print('min combined p=',min(pcomb[1:]),'at energy cut 100TeV')

plt.show()

#min combined p= 0.000102 3.8sigma at energy cut 100TeV

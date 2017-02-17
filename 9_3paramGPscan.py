import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import ascii
from scipy.stats import poisson
from scipy.stats import combine_pvalues
from scipy import random as ra

# First, find the best pre-trial p-value from a 3-parameter scan: b_low, b_high, and Energy.
# Then, generate pre-trial p-values from scrambled data sets in the same scan.
# Print out the trial-corrected p-value as well as dump trial results to a text file.
# The output file allows many jobs to be run in parallel and the final p-value determined after,
# since scans take ~30 sec/trial on one core.  


#data = ascii.read('events_4year.txt')
data = ascii.read('events_6year.txt')


#classification
energy=data['col2']
declin=data['col6']
righta=data['col7']
nrevents=energy.size

#f = open('part1_4year.out', 'w')
f = open('part1_6year.out', 'w')

# Find the best value for Blow, return it and the p-value
def findBlow(blat):
    plowb=0.
    minplowb=1.
    for blow in range(20,300):
        blow=float(blow/10.0)
        #counts of events for b<blow from geometric model (constant value)
        st=2*np.pi*(np.cos((90-blow)*u.deg)-np.cos(90*u.deg)) #solidangle in str for b<10
        frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
        geocounts=float(frac*blat.size)
        if sum(blat<blow) == 0:
            lplowb=1 #P for 0+ counts in bin 10
        else:
            lplowb=1-sum(poisson.pmf(np.arange(0,sum(blat<blow)),geocounts)) #P for nr+ counts in bin 10
        if lplowb<minplowb:
            minplowb=lplowb
            minblow=blow
    return minblow,minplowb

def findBhigh(blat):
    phighb=0.
    minphighb=1.
    for bhigh in range(400,900):
        bhigh=float(bhigh/10.0)
        st=2*np.pi*(np.cos(0*u.deg)-np.cos((90-bhigh)*u.deg)) #solidangle in str for 50<b<90
        frac=2*st/(4*np.pi) #fraction of 2*str over 4pi
        geocounts=float(frac*blat.size)
        if sum(blat>bhigh) == 0:
            lphighb=poisson.pmf(0,geocounts) #P for 0 counts in bin 50-90 deg
        else:
            lphighb=sum(poisson.pmf(np.arange(0,sum(blat>bhigh)+1),geocounts)) #P for nr counts in bin 50-90 deg
        if lphighb<minphighb:
            minphighb=lphighb
            minbhigh=bhigh
    return minbhigh,minphighb

# Real data
myminpcomb=1.0
for e in range(30,300,5):
    cuts=energy > e # energy cut scan
    ecoord=coord.SkyCoord(ra=righta[cuts], dec=declin[cuts], unit='deg', frame='icrs')
    blat=abs(ecoord.galactic.b.deg) #b latitude from galactic coord

    blow,plow   = findBlow(blat)
    bhigh,phigh = findBhigh(blat)

    pcomb=combine_pvalues([plow,phigh])[1] # combine using Fisher's method
    if (pcomb<myminpcomb):
        myminpcomb = pcomb
        mydiags = [e, blow, plow, bhigh, phigh]

print(myminpcomb)
print(mydiags)

# Monte Carlo (scrambled trials)
ntrials = 1000
minpcombs = np.ones(ntrials)
for i in range(0,ntrials):
    for e in range(30,300,5):
        rRA=ra.random(nrevents)*360
        rDec_delta=ra.random(nrevents)*30-15
        rDec=declin+rDec_delta
        # This edge effect is not ideal...
        rDec[rDec<-90]=declin[rDec<-90]-rDec_delta[rDec<-90]
        rDec[rDec>90]=declin[rDec>90]-rDec_delta[rDec>90]
        recoord=coord.SkyCoord(ra=rRA, dec=rDec, unit='deg', frame='icrs')
        rblat=abs(recoord.galactic.b.deg)

        blow,plow   = findBlow(rblat)
        bhigh,phigh = findBhigh(rblat)
        pcomb=combine_pvalues([plow,phigh])[1] # combine using Fisher's method
        #print(e, blow, plow, bhigh, phigh, pcomb)
        if (pcomb<minpcombs[i]):
            minpcombs[i] = pcomb
            diags = [e, blow, plow, bhigh, phigh]

    #print('{0:0.4e}'.format(minpcombs[i]))
    f.write('{0:0.4e}\n'.format(minpcombs[i]))
    #print(diags)
    

pvalue = len(minpcombs[minpcombs<=myminpcomb])/float(ntrials)
print('pvalue = {}'.format(pvalue))

#plt.hist(minpcombs)
#plt.xlabel('Pre-trial P')
#plt.savefig('pretrial_pvalues.pdf')
#probability=sum(likelihood>mylikelihood)/float(ntrials)

#plt.show()


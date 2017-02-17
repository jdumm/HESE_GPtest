import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import ascii
from scipy.stats import poisson
from scipy.stats import combine_pvalues

# Scan for b_low and b_high (separately) and then combine them with Fisher's method (note they are not truly independent for a given max number of events).

#data = ascii.read('events_4year.txt')
data = ascii.read('events_6year.txt')


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

plowb=0.
minplowb=1.
for blow in range(20,300):
    blow=float(blow/10.0)
    #counts of events for b<10 from geometric model (constant value)
    st10=2*np.pi*(np.cos((90-blow)*u.deg)-np.cos(90*u.deg)) #solidangle in str for b<10
    frac10=2*st10/(4*np.pi) #fraction of 2*str over 4pi
    geocounts10=frac10*nrevents

    if sum(blat<blow) == 0:
        lplowb=1 #P for 0+ counts in bin 10
    else:
        lplowb=1-sum(poisson.pmf(np.arange(0,sum(blat<blow)),geocounts10)) #P for nr+ counts in bin 10

    if lplowb<minplowb:
        minplowb=lplowb
        minblow=blow

    if blow==10:
        p10=lplowb
        geo10=geocounts10

    plowb=np.hstack((plowb,lplowb))

phighb=0.
minphighb=1.
for bhigh in range(400,900):
    bhigh=float(bhigh/10.0)
    #counts of events for b>50 from geometric model (constant value)
    st50_90=2*np.pi*(np.cos(0*u.deg)-np.cos((90-bhigh)*u.deg)) #solidangle in str for 50<b<90
    frac50_90=2*st50_90/(4*np.pi) #fraction of 2*str over 4pi
    geocounts50_90=frac50_90*nrevents

    if sum(blat>bhigh) == 0:
        lphighb=poisson.pmf(0,geocounts50_90) #P for 0 counts in bin 50-90 deg
    else:
        lphighb=sum(poisson.pmf(np.arange(0,sum(blat>bhigh)+1),geocounts50_90)) #P for nr counts in bin 50-90 deg

    if lphighb<minphighb:
        minphighb=lphighb
        minbhigh=bhigh

    if bhigh==50:
        p50=lphighb
        geo50=geocounts50_90

    phighb=np.hstack((phighb,lphighb))

plt.plot(np.arange(20,300)/10.0,plowb[1:],color='b',label='$b_{\mathrm{low}}$')
plt.plot(np.arange(400,900)/10.0,phighb[1:],color='g',label='$b_{\mathrm{high}}$')
plt.hlines((0.0455,0.0027),0,90,linestyle='dashed')
plt.yscale('log')
plt.text(80,0.06,r'$2\sigma$',fontsize=13)
plt.text(80,0.003,r'$3\sigma$',fontsize=13)
plt.xlabel('$b_{\mathrm{low}}$ , $b_{\mathrm{high}}$ (deg)')
plt.ylabel('P')
plt.legend(loc=9)
#plt.savefig('3_P_singel_plot_4year.pdf')
#plt.savefig('3_P_singel_plot_6year.pdf')

#combined p value
pcomb10_50=combine_pvalues([p10,p50])[1]
pcombmin=combine_pvalues([minplowb,minphighb])[1]

print('b=10 geocounts',geo10,'b=50 geocounts',geo50)
print('b=10 data counts p=',p10,'b=50 data counts p=',p50,'combined p=',pcomb10_50)
print('min blow at',minblow,'min p=',minplowb,'min bhigh at',minbhigh,'min p=',minphighb,'min combined p=',pcombmin)

plt.show()

#b=10 geocounts 3.2993 b=50 geocounts 4.4452
#b=10 9 counts p= 0.006903 b=50 0 counts p= 0.011735 combined p= 0.000844 3.3sigma
#min blow at 6.2 min p= 0.001287 min bhigh at 47.1 min p= 0.006209 min combined p= 0.000102 3.8sigma

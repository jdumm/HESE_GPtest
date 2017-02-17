import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import ascii
from scipy import random as ra

# Read in events, generate equitorial and galactic sky maps.

#data = ascii.read('events_4year.txt')
data = ascii.read('events_6year.txt')

#classification
energy=data['col2']
declin=data['col6']
righta=data['col7']

#signal data
#cuts=energy > 100 #higher than 100TeV
cuts=energy > 96 #higher than 100TeV but adjusted downwards for energy recalibration

gplanel=np.arange(-180,180,0.5)
gplanecoord=coord.SkyCoord(l=gplanel, b=0, unit='deg', frame='galactic')
gplaneRa=gplanecoord.icrs.ra.deg
gplanedec=gplanecoord.icrs.dec.deg

gcencoord=coord.SkyCoord(l=0, b=0, unit='deg', frame='galactic')
gcenRa=gcencoord.icrs.ra.deg
gcendec=gcencoord.icrs.dec.deg

gp10planecoord=coord.SkyCoord(l=gplanel, b=10, unit='deg', frame='galactic')
gp10planeRa=gp10planecoord.icrs.ra.deg
gp10planedec=gp10planecoord.icrs.dec.deg

gm10planecoord=coord.SkyCoord(l=gplanel, b=-10, unit='deg', frame='galactic')
gm10planeRa=gm10planecoord.icrs.ra.deg
gm10planedec=gm10planecoord.icrs.dec.deg

gp50planecoord=coord.SkyCoord(l=gplanel, b=50, unit='deg', frame='galactic')
gp50planeRa=gp50planecoord.icrs.ra.deg
gp50planedec=gp50planecoord.icrs.dec.deg

gm50planecoord=coord.SkyCoord(l=gplanel, b=-50, unit='deg', frame='galactic')
gm50planeRa=gm50planecoord.icrs.ra.deg
gm50planedec=gm50planecoord.icrs.dec.deg

eRa=righta[cuts]
edec=declin[cuts]

eplaneRa=np.arange(0,360,0.5)
eplanecoord=coord.SkyCoord(ra=eplaneRa, dec=0, unit='deg', frame='icrs')
eplaneb=eplanecoord.galactic.b.deg
eplanel=-eplanecoord.galactic.l.deg
eplanel[eplanel<-180]=eplanel[(eplanel<-180)]+360 #fixar disposition av -180 180 och 0 360

'''
etestcoord=coord.SkyCoord(ra=eplaneRa, dec=-60, unit='deg', frame='icrs')
gtestb=etestcoord.galactic.b.deg
gtestl=-etestcoord.galactic.l.deg
gtestl[gtestl<-180]=gtestl[(gtestl<-180)]+360 #fixar disposition av -180 180 och 0 360
'''

ecoord=coord.SkyCoord(ra=eRa, dec=edec, unit='deg', frame='icrs')
gb=ecoord.galactic.b.deg #b latitude from galactic coord
gl=-ecoord.galactic.l.deg
gl[gl<-180]=gl[(gl<-180)]+360 #fixar disposition av -180 180 och 0 360

#print(np.concatenate((eRa[:,np.newaxis],edec[:,np.newaxis],gl[:,np.newaxis],gb[:,np.newaxis],data['col9'][cuts][:,np.newaxis]),axis=1))

mcut=data['col9'][cuts]==1
mRa=eRa[mcut]
mdec=edec[mcut]
mecoord=coord.SkyCoord(ra=mRa, dec=mdec, unit='deg', frame='icrs')
mb=mecoord.galactic.b.deg
ml=-mecoord.galactic.l.deg
ml[ml<-180]=ml[(ml<-180)]+360

#fixar disposition av -180 180 och 0 360
eRa=180-eRa
mRa=180-mRa
gplaneRa=180-gplaneRa
gcenRa=180-gcenRa
gp10planeRa=180-gp10planeRa
gm10planeRa=180-gm10planeRa
gp50planeRa=180-gp50planeRa
gm50planeRa=180-gm50planeRa

plt.subplot(211,projection="mollweide")
plt.title('                                  Equatorial')
plt.plot(gplaneRa*u.deg.to(u.rad),gplanedec*u.deg.to(u.rad),',',color='r',alpha=1.0)
plt.plot(gcenRa*u.deg.to(u.rad),gcendec*u.deg.to(u.rad),'o',color='r',alpha=0.5)
plt.plot(gp10planeRa*u.deg.to(u.rad),gp10planedec*u.deg.to(u.rad),',',color='g',alpha=1.0)
plt.plot(gm10planeRa*u.deg.to(u.rad),gm10planedec*u.deg.to(u.rad),',',color='g',alpha=1.0)
plt.plot(gp50planeRa*u.deg.to(u.rad),gp50planedec*u.deg.to(u.rad),',',color='y',alpha=1.0)
plt.plot(gm50planeRa*u.deg.to(u.rad),gm50planedec*u.deg.to(u.rad),',',color='y',alpha=1.0)
plt.plot(eRa*u.deg.to(u.rad),edec*u.deg.to(u.rad),'.',color='b',ms=8)
plt.plot(mRa*u.deg.to(u.rad),mdec*u.deg.to(u.rad),'.',color='c',ms=8)
plt.grid(True)
plt.xticks(np.array([-180,-120,-60,0,60,120,180])*np.pi/180,('','360$^{\circ}$           ','','180$^{\circ}$','','                0$^{\circ}$','      RA'))
plt.yticks(np.array([-90,-60,-30,0,30,60,90])*np.pi/180,('-90$^{\circ}$','-60$^{\circ}$','-30$^{\circ}$','0$^{\circ}$','30$^{\circ}$','60$^{\circ}$','90$^{\circ}$'))
#plt.xlabel('Ra')
plt.ylabel('Dec')

plt.subplot(212,projection="mollweide")
plt.title('                               Galactic')
#plt.plot(eplanel*u.deg.to(u.rad),eplaneb*u.deg.to(u.rad),',',color='m',alpha=0.1)
#plt.plot(gtestl*u.deg.to(u.rad),gtestb*u.deg.to(u.rad),',')
plt.hlines(0*u.deg.to(u.rad),-180*u.deg.to(u.rad),180*u.deg.to(u.rad),color='r')
plt.plot(0,0,'o',color='r',alpha=0.5)
plt.hlines(10*u.deg.to(u.rad),-180*u.deg.to(u.rad),180*u.deg.to(u.rad),color='g')
plt.hlines(-10*u.deg.to(u.rad),-180*u.deg.to(u.rad),180*u.deg.to(u.rad),color='g')
plt.hlines(50*u.deg.to(u.rad),-180*u.deg.to(u.rad),180*u.deg.to(u.rad),color='y')
plt.hlines(-50*u.deg.to(u.rad),-180*u.deg.to(u.rad),180*u.deg.to(u.rad),color='y')
plt.plot(gl*u.deg.to(u.rad),gb*u.deg.to(u.rad),'.',color='b',ms=8)
plt.plot(ml*u.deg.to(u.rad),mb*u.deg.to(u.rad),'.',color='c',ms=8)
plt.grid(True)
plt.xticks(np.array([-180,-120,-60,0,60,120,180])*np.pi/180,('','180$^{\circ}$           ','','0$^{\circ}$','','           -180$^{\circ}$','   $l$'))
plt.yticks(np.array([-90,-60,-30,0,30,60,90])*np.pi/180,('-90$^{\circ}$','-60$^{\circ}$','-30$^{\circ}$','0$^{\circ}$','30$^{\circ}$','60$^{\circ}$','90$^{\circ}$'))
#plt.xlabel('l')
plt.ylabel('$b$')

plt.tight_layout()
plt.savefig('1_celestial_map_mollweide.pdf')
plt.show()


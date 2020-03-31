#plots the VEX magnetic field and the bow shock shape

from scipy import stats
import scipy.io
from matplotlib import cm
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import sunpy.time
import time
import pickle
import seaborn as sns
import math
import sys

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

plt.close('all')


sns.set_context("talk")     
sns.set_style("white")  
sns.set_style("ticks")
sns.set_style({"xtick.direction": "in","ytick.direction": "in"})


######################## Venus




print ('read VEX')
#get insitu data
vex= pickle.load( open( "DATACAT/VEX_2007to2014_VSO.p", "rb" ) )
print( 'read VEX done.')


rvenus=6052 

vex.x=vex.x/rvenus
vex.y=vex.y/rvenus
vex.z=vex.z/rvenus
vex.r=vex.r/rvenus

print('how much data is nan in the non-removed dataset?')
a=np.isnan(vex.btot) 
b=np.where(a == True)
print(np.size(b)/np.size(vex.btot)*100)


#units for plot are Venus radii	

xx=np.linspace(min(vex.x),max(vex.x),71)
yy=np.linspace(min(vex.y),max(vex.y),71)
zz=np.linspace(min(vex.z),max(vex.z),71)

btotaverage=np.nanmean(vex.btot)
background_wind=4

#whether the values for the grid are computed 
compute_grid_xy=1




################################################################## X-Y PLANE

print('X-Y plane')

if compute_grid_xy > 0:
 #the grid for the contours
 XX,YY=np.meshgrid(xx,yy)
 #this contains the values for the contours
 ZZ = np.zeros((np.size(xx),np.size(yy)))
 ZZ.fill(np.nan)

 print('start loop')
 insidez=np.where(np.logical_and(vex.z < 1.5, vex.z > -1.5))
 #go through each grid point
 for i in range(np.size(xx)-1):
   print(i/np.size(yy)*100)
   for j in range(np.size(yy)-1): 
     #check for x and y points that are inside the cell 
     insidex=np.where(np.logical_and(vex.x[insidez] < xx[i+1],vex.x[insidez] > xx[i]))
     insidey=np.where(np.logical_and(vex.y[insidez] < yy[j+1],vex.y[insidez] > yy[j]))
     #make an array of only those points that are inside the cell
     active=np.intersect1d(insidex,insidey)
     #minimum 10 datapoints
     if np.size(np.where(np.isnan(vex.btot[insidez][active]) == False)) > 10: ZZ[i,j]=np.nanmean(vex.btot[insidez][active])
 print('end loop')
 pickle.dump([XX,YY,ZZ], open( "plots/planets/VEX_XY_plane.p", "wb" ) )
 
 
if compute_grid_xy == 0: XX, YY, ZZ = pickle.load( open( "plots/planets/VEX_XY_plane.p", "rb" ) )

#contourf
#Transpose is used to make the orientation of the plot correct   
#plt.contourf(XX, YY, ZZ.transpose(), 100, alpha=1, cmap='inferno')
#plt.colorbar()
#plt.ylabel('Y VSO [km]')
#plt.xlabel('X VSO [km]')

ax = plt.axes()
plt.figtext(0.7,0.88,'to Sun', fontsize=15, ha='center')
ax.arrow(2.7, 2.7, 0.15, 0, head_width=0.1, head_length=0.1, fc='k', ec='k')

plt.figtext(0.68,0.13,'orbital motion ', fontsize=15, ha='center')
ax.arrow(2.3, -2.35, 0, -0.15, head_width=0.1, head_length=0.1, fc='k', ec='k')


plt.figure(1, figsize=(10,8), dpi=100, facecolor='w', edgecolor='w')
plt.subplot(111)

#Btotal is normalized to average solar wind field including ICMEs and Venus
#plt.imshow(np.log(ZZ.transpose()/np.nanmean(vex.btot)), cmap='inferno', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])
#plt.imshow(ZZ.transpose()/btotaverage, cmap='inferno', interpolation='nearest',extent=[min(vex.x), max(vex.x), min(vex.y),max(vex.y)])
#plt.imshow(ZZ.transpose(), cmap='inferno', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])

#plt.contourf(XX, YY, ZZ.transpose()/btotaverage, 100, alpha=1, cmap='inferno')
#plt.imshow(np.rot90(ZZ)/background_wind, cmap='inferno', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])

plt.imshow(np.log10(np.rot90(ZZ)), cmap='jet', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])


plt.figtext(0.02,0.01,r'Data from MAG instrument: Zhang et al. (2006). Shapes from Zhang et al. (2008). $C. M\ddot{o}stl$', fontsize=8, ha='left')


plt.colorbar()
plt.figtext(0.8,0.05,'log10(|B|)', fontsize=15, ha='center')
plt.ylabel('Y VSO [Venus radii]')
plt.xlabel('X VSO [Venus radii]')
plt.title('Total magnetic field      VEX at Venus 2007-2014')


#plot Venus as circle, with rvenus as 1
circle1=plt.Circle((0,0),1,color='k', fill=False, linewidth=2)
plt.gcf().gca().add_artist(circle1)
plt.axis('equal')

#shape model in polar coordinates used for removing Venus field
SZA=np.linspace(-110,110,100)*np.pi/180
R= 3.5  / (1+0.621*np.cos(SZA))
xshape,zshape=pol2cart(R,SZA)
#plt.plot(xshape,zshape, color='steelblue',linewidth=2)

#Zhang 2008 model
#solar zenith angle in rad
SZAless=np.linspace(-117,117,100)*np.pi/180
SZAtop=np.linspace(117,160,100)*np.pi/180

Rless= 2.14  / (1+0.621*np.cos(SZAless))
Rtop=2.364/np.sin(SZAtop+10.5*np.pi/180)

xless,yless=pol2cart(Rless,SZAless)
xtop,ytop=pol2cart(Rtop,SZAtop)

plt.plot(xless,yless, color='black',linewidth=2, linestyle='--')
plt.plot(xtop,ytop, color='black',linewidth=2, linestyle='--')
plt.plot(xtop,-ytop, color='black',linewidth=2, linestyle='--')	

plt.axis([-3,3,-3,3])

sns.despine()

plt.show()
plt.savefig('plots/planets/venus_shape_XY.pdf', format='pdf', dpi=300)
plt.savefig('plots/planets/venus_shape_XY.png', format='png', dpi=300)






sys.exit()

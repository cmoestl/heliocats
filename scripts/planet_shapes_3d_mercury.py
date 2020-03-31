#animation ICMECAT rose scatter plot with HI CME circles and in situ data	

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


######################## Mars



#mav= pickle.load( open( "DATACAT/MAVEN_2014to2016_MSO_removed_orbit.p", "rb" ) )
#print('how much data is valid in the removed interpolated dataset?')
#a=np.isnan(mav.btot) 
#b=np.where(a == False)
#print(np.size(b)/np.size(mav.btot)*100)

mes	= pickle.load( open( "DATACAT/MES_2011to2014_MSO.p", "rb" ) )
print('how much data is not NAN in the removed dataset?')

a=np.isnan(mes.btot) 
b=np.where(a == False)
print(np.size(b)/np.size(mes.btot)*100)


#Mars in km
rmerc=2440

mes.x=mes.x/rmerc
mes.y=mes.y/rmerc
mes.z=mes.z/rmerc
mes.r=mes.r/rmerc


#units for plot are Mars radii

xx=np.linspace(min(mes.x),max(mes.x),71)
yy=np.linspace(min(mes.y),max(mes.y),71)
zz=np.linspace(min(mes.z),max(mes.z),71)

btotaverage=np.nanmean(mes.btot)
background_wind=4


#valid for all 3
compute_grid_xy=0
compute_grid_yz=0
compute_grid_xz=0






################################################################## X-Y PLANE

print('X-Y plane')

if compute_grid_xy > 0:
 #the grid for the contours
 XX,YY=np.meshgrid(xx,yy)
 #this contains the values for the contours
 ZZ = np.zeros((np.size(xx),np.size(yy)))
 ZZ.fill(np.nan)

 print('start loop')
 insidez=np.where(np.logical_and(mes.z < 1.5, mes.z > -1.5))
 #go through each grid point
 for i in range(np.size(xx)-1):
   print(i/np.size(yy)*100)
   for j in range(np.size(yy)-1): 
     #check for x and y points that are inside the cell 
     insidex=np.where(np.logical_and(mes.x[insidez] < xx[i+1],mes.x[insidez] > xx[i]))
     insidey=np.where(np.logical_and(mes.y[insidez] < yy[j+1],mes.y[insidez] > yy[j]))
     #make an array of only those points that are inside the cell
     active=np.intersect1d(insidex,insidey)
     #minimum 10 datapoints
     if np.size(np.where(np.isnan(mes.btot[insidez][active]) == False)) > 10: ZZ[i,j]=np.nanmedian(mes.btot[insidez][active]) 
 print('end loop')
 pickle.dump([XX,YY,ZZ], open( "plots/planets/mes_XY_plane.p", "wb" ) )
 
 
if compute_grid_xy == 0: XX, YY, ZZ = pickle.load( open( "plots/planets/mes_XY_plane.p", "rb" ) )

#contourf
#Transpose is used to make the orientation of the plot correct   
#plt.contourf(XX, YY, ZZ.transpose(), 100, alpha=1, cmap='inferno')
#plt.colorbar()
#plt.ylabel('Y MSO [km]')
#plt.xlabel('X MSO [km]')



plt.figure(1, figsize=(10,8), dpi=100, facecolor='w', edgecolor='w')
plt.subplot(111)

#Btotal is normalized to average solar wind field including ICMEs and Mars
#plt.imshow(np.log(ZZ.transpose()/np.nanmean(mes.btot)), cmap='inferno', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])
#plt.imshow(ZZ.transpose()/btotaverage, cmap='inferno', interpolation='nearest',extent=[min(mes.x), max(mes.x), min(mes.y),max(mes.y)])
#plt.imshow(np.log10(np.rot90(ZZ)/background_wind), cmap='jet', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])
plt.imshow(np.log10(np.rot90(ZZ)), cmap='jet', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])


#plt.imshow((np.rot90(ZZ)), cmap='jet', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])


#plt.contourf(XX, YY, ZZ.transpose()/btotaverage, 100, alpha=1, cmap='inferno')


plt.colorbar()
#plt.figtext(0.8,0.05,'log10(|B|/<|Bsw|>)', fontsize=15, ha='center')
plt.figtext(0.8,0.05,'log10(|B|)', fontsize=15, ha='center')

#plt.figtext(0.8,0.05,'|B| [nT]', fontsize=15, ha='center')
ax = plt.axes()
plt.figtext(0.68,0.80,'to Sun', fontsize=15, ha='center')
ax.arrow(2.2, 2.1, 0.15, 0, head_width=0.1, head_length=0.1, fc='k', ec='k')
plt.figtext(0.68,0.13,'orbital motion ', fontsize=15, ha='center')
ax.arrow(2.3, -2.35, 0, -0.15, head_width=0.1, head_length=0.1, fc='k', ec='k')

plt.ylabel('Y MSO [Mercury radii]')
plt.xlabel('X MSO [Mercury radii]')
plt.title('Total magnetic field         MESSENGER at Mercury 2011-2014')

plt.figtext(0.02,0.01,r'Data from MAG instrument: Anderson et al. (2007). Shapes from Winslow et al. (2013). $C. M\ddot{o}stl$', fontsize=8, ha='left')

#plot Mars as circle, with rmars as 1
circle1=plt.Circle((0,0),1,color='k', fill=False, linewidth=2)
plt.gcf().gca().add_artist(circle1)
plt.axis('equal')

#shape model in polar coordinates used for removing Mars field
SZA=np.linspace(-110,110,100)*np.pi/180

#Bowshock
#winslow et al. 2013 section 3.1.2
Rbs= 0.5  +(2.75*1.04)/ (1+(1.04)*np.cos(SZA))
xbs,ybs=pol2cart(Rbs,SZA)
plt.plot(xbs,ybs, color='black',linewidth=2, linestyle='--')




#Magnetopause
#Winslow et al. 2013 with Shue et al. model v(see Shue et al. 1998 paper)
#X and Y should be abberrated, but here we use XMSO, YMSO; ZMSO is equal
#x=np.linspace(-3,3,100)
#y=np.linspace(-3,3,100)
#z=np.linspace(-3,3,100)
#rho=np.sqrt(y**2+(z-0.196)**2)
#theta=np.arctan(rho/x)

SZA=np.linspace(-140,140,100)*np.pi/180
Rmp=1.45*((2.0/(1+np.cos(SZA)))**0.5)

xmp,ymp=pol2cart(Rmp,SZA)
plt.plot(xmp,ymp, color='black',linewidth=2, linestyle='--')

#plt.axis([min(xx),max(xx),min(yy),max(yy)])
plt.axis([-3,3,-3,3])

sns.despine()

plt.show()
plt.savefig('plots/planets/Mercury_shape_XY.pdf', format='pdf', dpi=300)
plt.savefig('plots/planets/Mercury_shape_XY.png', format='png', dpi=300)


sys.exit()








################################################################## Y-Z PLANE

print('Y-Z plane')


insidex=0
insidey=0
insidez=0
XX=0
YY=0
ZZ=0


if compute_grid_yz > 0:
 #the grid for the contours
 YY,ZZ=np.meshgrid(yy,zz)
 #this contains the values for the contours
 XX = np.zeros((np.size(yy),np.size(zz)))
 XX.fill(np.nan)


 print('start loop')
 insidex=np.where(np.logical_and(mav.x < 1.5, mav.x > -1.5))
 #go through each grid point
 for i in range(np.size(yy)-1):
  print(i/np.size(yy)*100)
  for j in range(np.size(zz)-1): 
    #check for x and y points that are inside the cell 
    insidey=np.where(np.logical_and(mav.y[insidex] < yy[i+1],mav.y[insidex] > yy[i]))
    insidez=np.where(np.logical_and(mav.z[insidex] < zz[j+1],mav.z[insidex] > zz[j]))
    #make an array of only those points that are inside the cell
    active=np.intersect1d(insidey,insidez)
    #minimum 10 datapoints that are not NaN
    if np.size(np.where(np.isnan(mav.btot[insidex][active]) == False)) > 10: XX[i,j]=np.nanmean(mav.btot[insidex][active])
 print('end loop')
 
 pickle.dump([XX,YY,ZZ], open( "plots/planets/mav_YZ_plane.p", "wb" ) )
 
 
if compute_grid_yz == 0: XX, YY, ZZ = pickle.load( open( "plots/planets/mav_YZ_plane.p", "rb" ) )

#contourf
#Transpose is used to make the orientation of the plot correct   
#plt.contourf(XX, YY, ZZ.transpose(), 100, alpha=1, cmap='inferno')
#plt.colorbar()
#plt.ylabel('Y MSO [km]')
#plt.xlabel('X MSO [km]')


plt.figure(2, figsize=(10,8), dpi=100, facecolor='w', edgecolor='w')
plt.subplot(111)

#Btotal is normalized to average solar wind field including ICMEs
#plt.imshow(np.log(YY.transpose()/np.nanmean(mav.btot)), cmap='inferno', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])
#plt.imshow(XX.transpose()/np.nanmean(mav.btot), cmap='inferno', interpolation='nearest',extent=[min(mav.y), max(mav.y), min(mav.z),max(mav.z)])
#plt.imshow(ZZ.transpose(), cmap='inferno', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])
plt.imshow(np.log10(np.rot90(XX)/background_wind), cmap='jet', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])

plt.colorbar()
plt.figtext(0.8,0.05,'B/<B>', fontsize=15, ha='center')
plt.xlabel('Y MSO [Mars radii]')
plt.ylabel('Z MSO [Mars radii]')
plt.title('Total magnetic field      MAVEN    Y/Z plane')


#plot Mars as circle, with rmars as 1
circle1=plt.Circle((0,0),1,color='k', fill=False, linewidth=2)
plt.gcf().gca().add_artist(circle1)


Rbs= (0.55+0.12)  +(2.1 +0.09)

#plot Rbs at terminator with SZA=90
circle2=plt.Circle((0,0),Rbs,color='steelblue', fill=False, linewidth=2)
plt.gcf().gca().add_artist(circle2)


Rmpb= 0.86 + 0.90
#plot Rmbp at terminator
circle3=plt.Circle((0,0),Rmpb,color='orange', fill=False, linewidth=2)
plt.gcf().gca().add_artist(circle3)

plt.axis('equal')

plt.axis([min(yy),max(yy),min(zz),max(zz)])

plt.show()
plt.savefig('plots/planets/Mars_shape_YZ.pdf', format='pdf', dpi=300)
plt.savefig('plots/planets/Mars_shape_YZ.png', format='png', dpi=300)









#sys.exit()























################################################################## X-Z PLANE

print('X-Z plane')


insidex=0
insidey=0
insidez=0
XX=0
YY=0
ZZ=0


if compute_grid_xz > 0:
 #the grid for the contours
 XX,ZZ=np.meshgrid(xx,zz)
 #this contains the values for the contours
 YY = np.zeros((np.size(xx),np.size(zz)))
 YY.fill(np.nan)


 print('start loop')
 insidey=np.where(np.logical_and(mav.y < 1.5, mav.y > -1.5))
 #go through each grid point
 for i in range(np.size(xx)-1):
  print(i/np.size(yy)*100)
  for j in range(np.size(zz)-1): 
    #check for x and y points that are inside the cell 
    insidex=np.where(np.logical_and(mav.x[insidey] < xx[i+1],mav.x[insidey] > xx[i]))
    insidez=np.where(np.logical_and(mav.z[insidey] < zz[j+1],mav.z[insidey] > zz[j]))
    #make an array of only those points that are inside the cell
    active=np.intersect1d(insidex,insidez)
    #minimum 10 datapoints that are not NaN
    if np.size(np.where(np.isnan(mav.btot[insidey][active]) == False)) > 10:  YY[i,j]=np.nanmean(mav.btot[insidey][active])

 print('end loop')
 pickle.dump([XX,YY,ZZ], open( "plots/planets/mav_XZ_plane.p", "wb" ) )
 
if compute_grid_xz == 0: XX, YY, ZZ = pickle.load( open( "plots/planets/mav_XZ_plane.p", "rb" ) )

#contourf
#Transpose is used to make the orientation of the plot correct   
#plt.contourf(XX, YY, ZZ.transpose(), 100, alpha=1, cmap='inferno')
#plt.colorbar()
#plt.ylabel('Y MSO [km]')
#plt.xlabel('X MSO [km]')



plt.figure(3, figsize=(10,8), dpi=100, facecolor='w', edgecolor='w')
plt.subplot(111)

#Btotal is normalized to average solar wind field including ICMEs
#plt.imshow(np.log(YY.transpose()/np.nanmean(mav.btot)), cmap='inferno', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])
#plt.imshow(YY.transpose()/np.nanmean(mav.btot), cmap='inferno', interpolation='nearest',extent=[min(mav.x), max(mav.x), min(mav.z),max(mav.z)])

#plt.imshow(ZZ.transpose(), cmap='inferno', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])
plt.imshow(np.log10(np.rot90(YY)/background_wind), cmap='jet', interpolation='nearest',extent=[min(xx), max(xx), min(yy),max(yy)])

#plt.contourf(XX, ZZ, YY.transpose(), 100, alpha=1, cmap='inferno')

plt.colorbar()
plt.figtext(0.8,0.05,'B/<B>', fontsize=15, ha='center')
plt.xlabel('X MSO [Mars radii]')
plt.ylabel('Z MSO [Mars radii]')
plt.title('Total magnetic field      MAVEN    X/Z plane')


#plot Mars as circle, with rmars as 1
circle1=plt.Circle((0,0),1,color='k', fill=False, linewidth=2)
plt.gcf().gca().add_artist(circle1)
plt.axis('equal')


#shape model in polar coordinates used for removing Mars field
SZA=np.linspace(-110,110,100)*np.pi/180

#Bowshock
#Edberg epsilon 1.05 ± 0.04 L= 2.1 ± 0.09 X0=0.55 ± 0.12
#make the bow shock shape at maximum of error bars
Rbs= (0.55+0.12)  +(2.1 +0.09)/ (1+(1.05-0.04)*np.cos(SZA))
xbs,ybs=pol2cart(Rbs,SZA)
plt.plot(xbs,ybs, color='black',linewidth=2, linestyle='--')

#Magnetic pileup boundary
#Edberg 2008
#epsilon 0.92 ± 0.03 L= 0.90 ± 0.06 X0=0.86 ± 0.11
Rmpb= 0.86 + 0.90  / (1+0.92*np.cos(SZA))
xmpb,ympb=pol2cart(Rmpb,SZA)
plt.plot(xmpb,ympb, color='black',linewidth=2, linestyle='--')

plt.axis([min(xx),max(xx),min(zz),max(zz)])

plt.show()
plt.savefig('plots/planets/Mars_shape_XZ.pdf', format='pdf', dpi=300)
plt.savefig('plots/planets/Mars_shape_XZ.png', format='png', dpi=300)























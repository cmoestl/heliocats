'''
save spacecraft and planet trajectories in numpy incl. Bepi Colombo, PSP, Solar Orbiter

Author: C. Moestl, IWF Graz, Austria
twitter @chrisoutofspace, https://github.com/cmoestl
last update: March 2020

needs python 3.7 with sunpy, heliopy, numba (conda environment "helio"")

MIT LICENSE
Copyright 2019, Christian Moestl 
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, 
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
permit persons to whom the Software is furnished to do so, subject to the following 
conditions:
The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF mercuryHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EvenusT SHALL THE AUTHORS OR COPYRIGHT 
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''



#import scipy.io
import os
import datetime
from datetime import datetime, timedelta
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pdb
import sunpy.time
import pickle
import seaborn as sns
import sys
import heliopy.data.spice as spicedata
import heliopy.spice as spice
import astropy
import time
import numba
from numba import jit
import multiprocessing
from sunpy.time import parse_time
#ignore warnings
import warnings
warnings.filterwarnings('ignore')

#for server
#matplotlib.use('Agg')


############################################ SETTINGS


#Coordinate System
#frame='HCI'
#frame='HEEQ'
frame='HEEQ'

print(frame)


#Time resolution
res_hours=1


################################## FUNCTIONS #############################################

@jit(nopython=True)
def sphere2cart(r, phi, theta):
    x = r*np.cos(theta)*np.cos(phi)
    y = r*np.cos(theta)*np.sin(phi)
    z = r*np.sin(theta)
    return (x, y, z) 

@jit(nopython=True)
def cart2sphere(x,y,z):
    r = np.sqrt(x**2+ y**2 + z**2)            # r
    theta = np.arctan2(z,np.sqrt(x**2+ y**2))     # theta
    phi = np.arctan2(y,x)                        # phi
    return (r, theta, phi)


##########################################################################################
######################################### CODE START #####################################
##########################################################################################




start=time.time()
##########################################  PSP

starttime =datetime(2018, 8,13)
endtime = datetime(2025, 8, 31)
psp_time = []
while starttime < endtime:
    psp_time.append(starttime)
    starttime += timedelta(hours=res_hours)
psp_time_num=parse_time(psp_time).iso        

spice.furnish(spicedata.get_kernel('psp_pred'))
psp=spice.Trajectory('SPP')
psp.generate_positions(psp_time,'Sun',frame)
print('PSP pos')

psp.change_units(astropy.units.AU)  
[psp_r, psp_lat, psp_lon]=cart2sphere(psp.x,psp.y,psp.z)
print('PSP conv')

############################################## BepiColombo

starttime =datetime(2018, 10, 21)
endtime = datetime(2025, 11, 2)
bepi_time = []
while starttime < endtime:
    bepi_time.append(starttime)
    starttime += timedelta(hours=res_hours)
bepi_time_num=parse_time(bepi_time).iso 

spice.furnish(spicedata.get_kernel('bepi_pred'))
bepi=spice.Trajectory('BEPICOLOMBO MPO') # or BEPICOLOMBO MMO
bepi.generate_positions(bepi_time,'Sun',frame)
bepi.change_units(astropy.units.AU)  
[bepi_r, bepi_lat, bepi_lon]=cart2sphere(bepi.x,bepi.y,bepi.z)
print('Bepi')



#################################################### Solar Orbiter

starttime = datetime(2020, 3, 1)
endtime = datetime(2029, 12, 31)
solo_time = []
while starttime < endtime:
    solo_time.append(starttime)
    starttime += timedelta(hours=res_hours)
solo_time_num=parse_time(solo_time).iso  

spice.furnish(spicedata.get_kernel('solo_2020'))
solo=spice.Trajectory('Solar Orbiter')
solo.generate_positions(solo_time, 'Sun',frame)
solo.change_units(astropy.units.AU)
[solo_r, solo_lat, solo_lon]=cart2sphere(solo.x,solo.y,solo.z)
print('Solo')

'''
plt.figure(1, figsize=(12,9))
plt.plot_date(psp_time,psp_r,'-', label='R')
plt.plot_date(psp_time,psp_lat,'-',label='lat')
plt.plot_date(psp_time,psp_lon,'-',label='lon')
plt.ylabel('AU / RAD')
plt.legend()




plt.figure(2, figsize=(12,9))
plt.plot_date(bepi_time,bepi_r,'-', label='R')
plt.plot_date(bepi_time,bepi_lat,'-',label='lat')
plt.plot_date(bepi_time,bepi_lon,'-',label='lon')
plt.title('Bepi Colombo position '+frame)
plt.ylabel('AU / RAD')
plt.legend()


plt.figure(3, figsize=(12,9))
plt.plot_date(solo_time,solo_r,'-', label='R')
plt.plot_date(solo_time,solo_lat,'-',label='lat')
plt.plot_date(solo_time,solo_lon,'-',label='lon')
plt.title('Solar Orbiter position '+frame)
plt.ylabel('AU / RAD')
plt.legend()


######## R with all three
plt.figure(4, figsize=(16,10))
plt.plot_date(psp_time,psp.r,'-',label='PSP')
plt.plot_date(bepi_time,bepi.r,'-',label='Bepi Colombo')
plt.plot_date(solo_time,solo.r,'-',label='Solar Orbiter')
plt.legend()
plt.title('Heliocentric distance of heliospheric observatories')
plt.ylabel('AU')
plt.savefig('results/positions_plots/bepi_psp_solo_R.png')

##### Longitude all three
plt.figure(5, figsize=(16,10))
plt.plot_date(psp_time,psp_lon*180/np.pi,'-',label='PSP')
plt.plot_date(bepi_time,bepi_lon*180/np.pi,'-',label='Bepi Colombo')
plt.plot_date(solo_time,solo_lon*180/np.pi,'-',label='Solar Orbiter')
plt.legend()
plt.title(frame+' longitude')
plt.ylabel('DEG')
plt.savefig('results/positions_plots/bepi_psp_solo_longitude_'+frame+'.png')
'''

############# Earth for mercury, venusus, STA
#https://docs.heliopy.org/en/stable/data/spice.html


planet_kernel=spicedata.get_kernel('planet_trajectories')

starttime =datetime(2018, 1, 1)
endtime = datetime(2029, 12, 31)
earth_time = []
while starttime < endtime:
    earth_time.append(starttime)
    starttime += timedelta(hours=res_hours)
earth_time_num=parse_time(earth_time).iso  


earth=spice.Trajectory('399')  #399 for Earth, not barycenter (because of moon)
earth.generate_positions(earth_time,'Sun',frame)
earth.change_units(astropy.units.AU)  
[earth_r, earth_lat, earth_lon]=cart2sphere(earth.x,earth.y,earth.z)
print('Earth')

################ mercury
mercury_time_num=earth_time_num

mercury=spice.Trajectory('1')  #barycenter
mercury.generate_positions(earth_time,'Sun',frame)  
mercury.change_units(astropy.units.AU)  
[mercury_r, mercury_lat, mercury_lon]=cart2sphere(mercury.x,mercury.y,mercury.z)
print('mercury') 

################# venus
venus_time_num=earth_time_num
venus=spice.Trajectory('2')  
venus.generate_positions(earth_time,'Sun',frame)  
venus.change_units(astropy.units.AU)  
[venus_r, venus_lat, venus_lon]=cart2sphere(venus.x,venus.y,venus.z)
print('venus') 


############### Mars
mars_time_num=earth_time_num
mars=spice.Trajectory('4')  
mars.generate_positions(earth_time,'Sun',frame)  
mars.change_units(astropy.units.AU)  
[mars_r, mars_lat, mars_lon]=cart2sphere(mars.x,mars.y,mars.z)
print('mars') 

############### Jupiter
jupiter_time_num=earth_time_num
jupiter=spice.Trajectory('5')  
jupiter.generate_positions(earth_time,'Sun',frame)  
jupiter.change_units(astropy.units.AU)  
[jupiter_r, jupiter_lat, jupiter_lon]=cart2sphere(jupiter.x,jupiter.y,jupiter.z)
print('jupiter') 


############### saturn
saturn_time_num=earth_time_num
saturn=spice.Trajectory('6')  
saturn.generate_positions(earth_time,'Sun',frame)  
saturn.change_units(astropy.units.AU)  
[saturn_r, saturn_lat, saturn_lon]=cart2sphere(saturn.x,saturn.y,saturn.z)
print('saturn') 

############### uranus
uranus_time_num=earth_time_num
uranus=spice.Trajectory('7')  
uranus.generate_positions(earth_time,'Sun',frame)  
uranus.change_units(astropy.units.AU)  
[uranus_r, uranus_lat, uranus_lon]=cart2sphere(uranus.x,uranus.y,uranus.z)
print('uranus') 

############### neptune
neptune_time_num=earth_time_num
neptune=spice.Trajectory('8')  
neptune.generate_positions(earth_time,'Sun',frame)  
neptune.change_units(astropy.units.AU)  
[neptune_r, neptune_lat, neptune_lon]=cart2sphere(neptune.x,neptune.y,neptune.z)
print('neptune') 




#############stereo-A
sta_time_num=earth_time_num
spice.furnish(spicedata.get_kernel('stereo_a_pred'))
sta=spice.Trajectory('-234')  
sta.generate_positions(earth_time,'Sun',frame)  
sta.change_units(astropy.units.AU)  
[sta_r, sta_lat, sta_lon]=cart2sphere(sta.x,sta.y,sta.z)
print('STEREO-A') 









#save positions with date as datetime


psp_time= mdates.date2num(psp_time)   
bepi_time= mdates.date2num(bepi_time)   
solo_time= mdates.date2num(solo_time)   
earth_time= mdates.date2num(earth_time)   


psp=np.rec.array([psp_time,psp_r,psp_lon,psp_lat, psp.x, psp.y,psp.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
bepi=np.rec.array([bepi_time,bepi_r,bepi_lon,bepi_lat,bepi.x, bepi.y,bepi.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
solo=np.rec.array([solo_time,solo_r,solo_lon,solo_lat,solo.x, solo.y,solo.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
sta=np.rec.array([earth_time,sta_r,sta_lon,sta_lat,sta.x, sta.y,sta.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
earth=np.rec.array([earth_time,earth_r,earth_lon,earth_lat, earth.x, earth.y,earth.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
venus=np.rec.array([earth_time,venus_r,venus_lon,venus_lat, venus.x, venus.y,venus.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
mars=np.rec.array([earth_time,mars_r,mars_lon,mars_lat, mars.x, mars.y,mars.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
mercury=np.rec.array([earth_time,mercury_r,mercury_lon,mercury_lat,mercury.x, mercury.y,mercury.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])

jupiter=np.rec.array([earth_time,jupiter_r,jupiter_lon,jupiter_lat,jupiter.x, jupiter.y,jupiter.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
saturn=np.rec.array([earth_time,saturn_r,saturn_lon,saturn_lat,saturn.x, saturn.y,saturn.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
uranus=np.rec.array([earth_time,uranus_r,uranus_lon,uranus_lat,uranus.x, uranus.y,uranus.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
neptune=np.rec.array([earth_time,neptune_r,neptune_lon,neptune_lat,neptune.x, neptune.y,neptune.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])


pickle.dump([psp, bepi, solo, sta, earth, venus, mars, mercury,jupiter, saturn, uranus, neptune,frame], open( 'results/positions_'+frame+'_1hr.p', "wb" ) )
sys.exit()


#save positions with date as char
psp=np.rec.array([psp_time_num,psp_r,psp_lon,psp_lat, psp.x, psp.y,psp.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
bepi=np.rec.array([bepi_time_num,bepi_r,bepi_lon,bepi_lat,bepi.x, bepi.y,bepi.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
solo=np.rec.array([solo_time_num,solo_r,solo_lon,solo_lat,solo.x, solo.y,solo.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
sta=np.rec.array([sta_time_num,sta_r,sta_lon,sta_lat,sta.x, sta.y,sta.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
earth=np.rec.array([earth_time_num,earth_r,earth_lon,earth_lat, earth.x, earth.y,earth.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
venus=np.rec.array([venus_time_num,venus_r,venus_lon,venus_lat, venus.x, venus.y,venus.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
mars=np.rec.array([mars_time_num,mars_r,mars_lon,mars_lat, mars.x, mars.y,mars.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
mercury=np.rec.array([mercury_time_num,mercury_r,mercury_lon,mercury_lat,mercury.x, mercury.y,mercury.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])

jupiter=np.rec.array([jupiter_time_num,jupiter_r,jupiter_lon,jupiter_lat,jupiter.x, jupiter.y,jupiter.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
saturn=np.rec.array([saturn_time_num,saturn_r,saturn_lon,saturn_lat,saturn.x, saturn.y,saturn.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
uranus=np.rec.array([uranus_time_num,uranus_r,uranus_lon,uranus_lat,uranus.x, uranus.y,uranus.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
neptune=np.rec.array([neptune_time_num,neptune_r,neptune_lon,neptune_lat,neptune.x, neptune.y,neptune.z],dtype=[('time','U16'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])



pickle.dump([psp, bepi, solo, sta, earth, venus, mars, mercury,jupiter, saturn, uranus, neptune,frame], open( 'results/positions_vr_'+frame+'.p', "wb" ) )



#a[1].strftime('%Y %m %d  %H')        
#psp_time

#make array
#psptxt=np.zeros(np.size(psp),dtype=[('time',object),('r', float),('lat', float), ('lon', float)])   


np.savetxt('results/positions_ascii/psp_'+frame+'.txt',psp,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/bepi_'+frame+'.txt',bepi,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/solo_'+frame+'.txt',solo,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/sta_'+frame+'.txt',sta,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/earth_'+frame+'.txt',earth,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/venus_'+frame+'.txt',venus,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/mars_'+frame+'.txt',mars,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/mercury_'+frame+'.txt',mercury,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/jupiter_'+frame+'.txt',jupiter,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/saturn_'+frame+'.txt',saturn,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/uranus_'+frame+'.txt',uranus,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')
np.savetxt('results/positions_ascii/neptune_'+frame+'.txt',neptune,header='time, r (AU), lon (rad), lat (rad), x,y,z (AU), frame:'+frame,fmt='%16s %.18e %.18e %.18e %.18e %.18e %.18e')






end=time.time()
print( 'generate position took time in seconds:', round((end-start),1) )


#################################################### MAIN ###############################


#make_positions()

#[psp, bepi, solo, sta, earth, venus, mars, mercury,jupiter, saturn, uranus, neptune,frame]=pickle.load( open( 'results/positions_vr_HEEQ.p', "rb" ) )















#for updating data every day for Wind and STEREO-A
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy
import copy
import matplotlib.dates as mdates
import datetime
import urllib
import json
import os
import pdb
from sunpy.time import parse_time
import scipy.io
import pickle
import sys
import cdflib
import matplotlib.pyplot as plt
import heliosat
from numba import njit
from astropy.time import Time
import heliopy.data.cassini as cassinidata
import heliopy.data.helios as heliosdata
import heliopy.data.spice as spicedata
import heliopy.spice as spice
import astropy




@njit
def cart2sphere(x,y,z):
    r = np.sqrt(x**2+ y**2 + z**2)           
    theta = np.arctan2(z,np.sqrt(x**2+ y**2))
    phi = np.arctan2(y,x)                    
    return (r, theta, phi)
    


@njit
def sphere2cart(r,theta,phi):
    x = r * np.sin( theta ) * np.cos( phi )
    y = r * np.sin( theta ) * np.sin( phi )
    z = r * np.cos( theta )
    return (x, y,z)
    





def save_wind_data_update(path,file):
    
    print('start wind update')
    wind_sat = heliosat.WIND()
    t_start = datetime.datetime(2018, 1, 1)
    #t_end = datetime.datetime(2018, 2, 1)
    t_end = datetime.datetime.utcnow()   

    
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    tm, mag = wind_sat.get_data_raw(t_start, t_end, "wind_mfi_h0")
    tp, pro = wind_sat.get_data_raw(t_start, t_end, "wind_swe_h1")
    print('download complete')
    
    tm=parse_time(tm,format='unix').datetime 
    tp=parse_time(tp,format='unix').datetime 
    
    
    #convert to matplotlib time for linear interpolation
    tm_mat=mdates.date2num(tm) 
    tp_mat=mdates.date2num(tp) 
    
    print('time convert done')
        
    print('position start')
    frame='HEEQ'
    planet_kernel=spicedata.get_kernel('planet_trajectories')
    earth=spice.Trajectory('399')  #399 for Earth, not barycenter (because of moon)
    earth.generate_positions(time,'Sun',frame)
    earth.change_units(astropy.units.AU)
    #*****with astropy lagrange points exact value? L1 position with 0.01 AU 
    [r, lat, lon]=cart2sphere(earth.x-0.01*astropy.units.AU,earth.y,earth.z)
    print('position end ')
        
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
        
    den = np.interp(time_mat, tp_mat, pro[:,0])
    vt = np.interp(time_mat, tp_mat, pro[:,1])
    tp = np.interp(time_mat, tp_mat, pro[:,2])
    #p3 = np.interp(time_mat, tp_mat, pro[:,3])
    #p4 = np.interp(time_mat, tp_mat, pro[:,4])

    
    #make array
    win=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('np', float),('vt', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    win = win.view(np.recarray)  

    #fill with data
    win.time=time
    win.bx=bx
    win.by=by
    win.bz=bz 
    win.bt=bt

    win.x=earth.x
    win.y=earth.y
    win.z=earth.z
    
    win.r=r
    win.lat=np.rad2deg(lat)
    win.lon=np.rad2deg(lon)

    
    win.np=den
    win.vt=vt
    #https://en.wikipedia.org/wiki/Thermal_velocity convert from km/s to K
    from astropy.constants import m_p,k_B
    win.tp=np.pi*m_p*((tp*1e3)**2)/(8*k_B) 
        
    #win.p3=p3
    #win.p4=p4
    
    header='Wind magnetic field (MAG instrument) and plasma data (SWE), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/wind/  '+ \
    'Timerange: '+win.time[0].strftime("%Y-%b-%d %H:%M")+' to '+win.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(win.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by win.time, win.bx, win.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(win.size)+'. '+\
    'Units are btxyz [nT, GSE], vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    pickle.dump([win,header], open(path+file, "wb"))
    
    
    print('wind update done')
    print()
    




def save_stereoa_beacon_data_update(path,file):

    print('start STA')
    sta_sat = heliosat.STA()
    t_start = datetime.datetime(2018, 1, 1)
    t_end = datetime.datetime(2018, 2, 1)

    #t_end = datetime.datetime.utcnow()  

    #t_end = datetime.datetime(2019, 12,31)
     
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    #tm, mag = sta_sat.get_data_raw(t_start, t_end, "sta_impact_l1")
    #tp, pro = sta_sat.get_data_raw(t_start, t_end, "sta_plastic_l2")

    tm, mag = sta_sat.get_data_raw(t_start, t_end, "mag_beacon")
    tp, pro = sta_sat.get_data_raw(t_start, t_end, "proton_beacon")

    print('download complete')
   
    tm=parse_time(tm,format='unix').datetime 
    tp=parse_time(tp,format='unix').datetime 

    #convert to matplotlib time for linear interpolation
    tm_mat=mdates.date2num(tm) 
    tp_mat=mdates.date2num(tp) 
    
    print('time convert done')
    
    
    print('position start')
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('stereo_a_pred'))
    statra=spice.Trajectory('-234') #STEREO-A SPICE NAIF code
    statra.generate_positions(time,'Sun',frame)
    statra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(statra.x,statra.y,statra.z)
    print('position end ')
    
    
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
      
      
      
    #add speed!!!!!!!!!!!!!!!!    
    den = np.interp(time_mat, tp_mat, pro[:,0])
    tp = np.interp(time_mat, tp_mat, pro[:,1])
    #p2 = np.interp(time_mat, tp_mat, pro[:,2])
    #p3 = np.interp(time_mat, tp_mat, pro[:,3])
    #p4 = np.interp(time_mat, tp_mat, pro[:,4])

    
    #make array
    sta=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    sta = sta.view(np.recarray)  

    #fill with data
    sta.time=time
    sta.bx=bx
    sta.by=by
    sta.bz=bz 
    sta.bt=bt



    sta.x=statra.x
    sta.y=statra.y
    sta.z=statra.z
    
    sta.r=r
    sta.lat=np.rad2deg(lat)
    sta.lon=np.rad2deg(lon)
    

    
    sta.np=den
    sta.tp=tp    
    #sta.p2=p2
    #sta.p3=p3
    #sta.p4=p4
    
       

    
    
    header='BEACON STEREO-A magnetic field (IMPACT instrument) and plasma data (PLASTIC), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/  '+ \
    'Timerange: '+sta.time[0].strftime("%Y-%b-%d %H:%M")+' to '+sta.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(sta.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by sta.time, sta.bx, sta.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(sta.size)+'. '+\
    'Units are btxyz [nT, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    #'Units are btxyz [nT, RTN], vtxyz [km/s, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\

    pickle.dump([sta,header], open(path+file, "wb"))
    
    print('done sta')
    print()



############### MAIN

data_path='/nas/helio/data/insitu_python/'


#filesta="sta_2018_now.p" 

#save_stereoa_beacon_data_update(data_path,filesta)

filewin="wind_2018_now.p" 

save_wind_data_update(data_path,filewin)






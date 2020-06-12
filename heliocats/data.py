#data.py
#load and save data for heliocats
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
import scipy.signal
import pickle
import time
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
import requests
import math


from config import data_path
#data_path='/nas/helio/data/insitu_python/'

heliosat_data_path='/nas/helio/data/heliosat/data/'

data_path_sun='/nas/helio/data/SDO_realtime/'

'''
MIT LICENSE
Copyright 2020, Christian Moestl, Rachel L. Bailey 
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, 
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
permit persons to whom the Software is furnished to do so, subject to the following 
conditions:
The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''


####################################### get new data ####################################


def save_wsa_hux(filein):
    #load wsa hux 

    windraw = np.loadtxt('data/wsa_hux_mars_aug2014_jan2018.txt', dtype=[('time','<U30'),('time2','<U30'),('time_mat', float),('vt', float)] )
    windraw = windraw.view(np.recarray)  

    wind=np.zeros(len(windraw),dtype=[('time',object),('vt', float)])   
    wind=wind.view(np.recarray) 


    for i in np.arange(len(windraw)):
        wind_time_str=windraw.time[i][8:12]+'-'+windraw.time[i][4:7]+'-'+windraw.time[i][1:3]+' '+windraw.time2[i][0:8]
        wind.time[i]=(parse_time(wind_time_str).datetime)

    wind.vt=windraw.vt
    
    fileout='wsa_hux_mars_aug2014_jan2018.p'      
    pickle.dump(wind, open(data_path+fileout, "wb"))
    
    return 0



def load_mars_wsa_hux():

    file='wsa_hux_mars_aug2014_jan2018.p'  
    rad=pickle.load(open(data_path+file, "rb"))
    
    return rad





def load_maven_sir_huang():
    
    
    #Huang et al. 2019 APJ convert PDF to excel with https://pdftoxls.com
    
    
   
    mavensir='sircat/sources/Huang_2019_SIR_MAVEN_table_1.xlsx'
    
    print('load MAVEN Huang SIR catalog from ', mavensir)
    ms=pd.read_excel(mavensir)
    ms=ms.drop(index=[0,1,2])
    ms_num=np.array(ms['No.'])
    ms_start=np.array(ms['Start'])
    ms_end=np.array(ms['End'])
    ms_si=np.array(ms['SI'])
    
    
    
    ms=np.zeros(len(ms_num),dtype=[('start',object),('end',object),('si',object)])   
    ms=ms.view(np.recarray) 
        
        
    #make correct years for start time 
    ms_num[np.where(ms_num< 7)[0]]=2014
    ms_num[np.where(ms_num< 27)[0]]=2015
    ms_num[np.where(ms_num< 64)[0]]=2016
    ms_num[np.where(ms_num< 83)[0]]=2017
    ms_num[np.where(ms_num< 127)[0]]=2018
    
    #make correct years for end and si time
    ms_num2=copy.deepcopy(ms_num)
    ms_num2[3]=2015
    ms_num2[62]=2017

    #transform date of start time
    for t in np.arange(0,len(ms_start)):
        #check for nans in between time strings
        if pd.isna(ms_start[t])==False:                    
            
            
            ####################### start time
            #year
            year=str(ms_num[t])
            
            #month
            datetimestr=ms_start[t]    
            datestr=datetimestr[0:2]
            
            monthfloat=float(datestr)
            month=str(int(np.floor(monthfloat)))
            
            #day            
            if int(month) < 10: day=datetimestr[2:4]
            if int(month) > 9: day=datetimestr[3:5]           
            
            #time
            timestr=datetimestr[-5:]                        
            
            #construct year month day
            datetimestrfin=str(ms_num[t])+'-'+month+'-'+day                    
            #remove white spaces at the end and add time
            finaldatetime=datetimestrfin.strip()+' '+timestr
            #print(ms_start[t])
            #print(finaldatetime)            
            ms.start[t]=parse_time(finaldatetime).datetime
            
            
            
            
            ################### end time            

            #year
            year=str(ms_num2[t])

            #month
            datetimestr=ms_end[t]    
            datestr=datetimestr[0:2]
            
            monthfloat=float(datestr)
            month=str(int(np.floor(monthfloat)))
            
            #day            
            if int(month) < 10: day=datetimestr[2:4]
            if int(month) > 9: day=datetimestr[3:5]           
            
            #time
            timestr=datetimestr[-5:]                        
            
            #construct year month day
            datetimestrfin=str(ms_num2[t])+'-'+month+'-'+day                    
            #remove white spaces at the end and add time
            finaldatetime=datetimestrfin.strip()+' '+timestr
            #print(ms_end[t])
            #print(finaldatetime)            
            ms.end[t]=parse_time(finaldatetime).datetime
            
            
                        
            ############# stream interface time            

            #year
            year=str(ms_num2[t])

            #month
            datetimestr=ms_si[t]    
            datestr=datetimestr[0:2]
            
            monthfloat=float(datestr)
            month=str(int(np.floor(monthfloat)))
            
            #day            
            if int(month) < 10: day=datetimestr[2:4]
            if int(month) > 9: day=datetimestr[3:5]           
            
            #time
            timestr=datetimestr[-5:]                        
            
            #construct year month day
            datetimestrfin=str(ms_num2[t])+'-'+month+'-'+day                    
            #remove white spaces at the end and add time
            finaldatetime=datetimestrfin.strip()+' '+timestr
            #print(ms_si[t])
            #print(finaldatetime)            
            ms.si[t]=parse_time(finaldatetime).datetime
            #print()


    #get rid of zeros where the years where stated in the original data        
    ms2 = ms[np.argwhere(ms)]

    return ms2





def save_msl_rad():

    #convert file 
    # year, doy, sol, doseE hourly [uGy/day], doseE sol-filtered [uGy/day]
    raw=np.loadtxt('data/doseE_sol_filter_2019.dat')
    
    rad=np.zeros(len(raw),dtype=[('time',object),('sol', float),('dose_hour', float),('dose_sol', float)])   
    rad = rad.view(np.recarray) 
    
    rad.sol=raw[:,2]
    rad.dose_hour=raw[:,3]
    rad.dose_sol=raw[:,4]


    #make datetime array from year and doy
    for i in np.arange(len(rad)):
        rad[i].time=parse_time(str(int(raw[i,0]))+'-01-01 00:00').datetime+datetime.timedelta(days=raw[i,1]-1)
        print(rad[i].time)
        
    file='msl_2012_2019_rad.p'  
    pickle.dump(rad, open(data_path+file, "wb"))
    
    return 0



def load_msl_rad():

    file='msl_2012_2019_rad.p'  
    rad=pickle.load(open(data_path+file, "rb"))
    
    return rad








def save_psp_data(path, file, sceq):

    t_start = datetime.datetime(2018, 10, 6)

    #t_start = datetime.datetime(2019, 4, 22)


    #t_start = datetime.datetime(2019, 3, 1) # 
    t_end = datetime.datetime(2019, 4, 24) #  UNTIL ERROR on Apr 25
    psp1=get_psp_data(t_start,t_end)

    t_start = datetime.datetime(2019, 4, 26)    
    #t_end = datetime.datetime(2019, 4, 30)    
    #t_end = datetime.datetime(2019, 10, 15)
    
    t_end = datetime.datetime(2019, 12, 31)    
    psp2=get_psp_data(t_start,t_end)

    #add both
    
    psp=np.zeros(np.size(psp1.time)+np.size(psp2.time),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('vx', float),('vy', float),('vz', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    psp = psp.view(np.recarray)  
    psp.time=np.hstack((psp1.time,psp2.time))
    psp.bx=np.hstack((psp1.bx,psp2.bx))
    psp.by=np.hstack((psp1.by,psp2.by))
    psp.bz=np.hstack((psp1.bz,psp2.bz))
    psp.bt=np.hstack((psp1.bt,psp2.bt))
    psp.vt=np.hstack((psp1.vt,psp2.vt))
    psp.vx=np.hstack((psp1.vx,psp2.vx))
    psp.vy=np.hstack((psp1.vy,psp2.vy))
    psp.vz=np.hstack((psp1.vz,psp2.vz))
    psp.np=np.hstack((psp1.np,psp2.np))
    psp.tp=np.hstack((psp1.tp,psp2.tp))
    psp.x=np.hstack((psp1.x,psp2.x))
    psp.y=np.hstack((psp1.y,psp2.y))
    psp.z=np.hstack((psp1.z,psp2.z))
    psp.r=np.hstack((psp1.r,psp2.r))
    psp.lon=np.hstack((psp1.lon,psp2.lon))
    psp.lat=np.hstack((psp1.lat,psp2.lat))

    print('Merging done')
    
    
    #convert magnetic field to SCEQ
    coord='RTN'
    if sceq==True:
        coord='SCEQ'
        psp=convert_RTN_to_SCEQ(psp,'PSP')
    
    
    header='PSP magnetic field (FIELDS instrument) and plasma data (SWEAP), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/psp/  '+ \
    'Timerange: '+psp.time[0].strftime("%Y-%b-%d %H:%M")+' to '+psp.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(psp.time)).seconds)+' seconds. '+\
    'The data are put in a numpy recarray, fields can be accessed by psp.time, psp.bx, psp.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(psp.size)+'. '+\
    'Units are btxyz [nT,'+coord+'], vtxyz [km/s, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([psp,header], open(path+file, "wb"))




def get_psp_data(t_start,t_end):
     
    print('start PSP')
     
    psp_sat = heliosat.PSP()
    
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    tm, mag = psp_sat.get_data_raw(t_start, t_end, "psp_fields_l2")#,return_datetimes=True)
    tp, pro = psp_sat.get_data_raw(t_start, t_end, "psp_spc_l3")#,return_datetimes=True)
    
    
    tm=parse_time(tm,format='unix').datetime 
    tp=parse_time(tp,format='unix').datetime 
    
    print('download complete')
    
    print('start nan or interpolate')
    
    print('field')
    #round first each original time to full minutes   original data at 30sec
    tround=copy.deepcopy(tm)
    format_str = '%Y-%m-%d %H:%M'  
    for k in np.arange(np.size(tm)):
         tround[k] = datetime.datetime.strptime(datetime.datetime.strftime(tm[k], format_str), format_str) 
    tm_mat=parse_time(tround).plot_date
    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    
   
    #which values are not in original data compared to full time range
    isin=np.isin(time_mat,tm_mat)      
    setnan=np.where(isin==False)
    #set to to nan that is not in original data
    bx[setnan]=np.nan
    by[setnan]=np.nan
    bz[setnan]=np.nan

    bt = np.sqrt(bx**2+by**2+bz**2)
    

    
    print('plasma')
        
    #for plasma round first each original time to full minutes
    tround=copy.deepcopy(tp)
    format_str = '%Y-%m-%d %H:%M'  
    for k in np.arange(np.size(tp)):
         tround[k] = datetime.datetime.strptime(datetime.datetime.strftime(tp[k], format_str), format_str) 
    tp_mat=mdates.date2num(tround) 
    
    isin=np.isin(time_mat,tp_mat)      
    setnan=np.where(isin==False)
        
    den = np.interp(time_mat, tp_mat, pro[:,0])
    vx = np.interp(time_mat, tp_mat, pro[:,1])
    vy = np.interp(time_mat, tp_mat, pro[:,2])
    vz = np.interp(time_mat, tp_mat, pro[:,3])
    temp = np.interp(time_mat, tp_mat, pro[:,4])
    
    den[setnan]=np.nan
    temp[setnan]=np.nan
    vx[setnan]=np.nan
    vy[setnan]=np.nan
    vz[setnan]=np.nan
  
    vt=np.sqrt(vx**2+vy**2+vz**2)

    print('end nan or interpolate')

        
    print('position start')
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('psp_pred'))
    psptra=spice.Trajectory('SPP')
    psptra.generate_positions(time,'Sun',frame)
    psptra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(psptra.x,psptra.y,psptra.z)
    print('PSP pos')    
    print('position end')

        
    #make array
    psp=np.zeros(np.size(bx),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('np', float),('tp', float),('x', float),('y', float),\
                ('z', float),('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    psp = psp.view(np.recarray)  

    #fill with data
    psp.time=time
    psp.bx=bx
    psp.by=by
    psp.bz=bz 
    psp.bt=bt
    
    psp.x=psptra.x
    psp.y=psptra.y
    psp.z=psptra.z
    
    psp.r=r
    psp.lat=np.degrees(lat)
    psp.lon=np.degrees(lon)
    
    psp.vt=vt
    psp.vx=vx    
    psp.vy=vy  
    psp.vz=vz
    psp.np=den
    #https://en.wikipedia.org/wiki/Thermal_velocity convert from km/s to K
    from astropy.constants import m_p,k_B
    psp.tp=np.pi*m_p*((temp*1e3)**2)/(8*k_B) 
    
    #remove spikes
    psp.tp[np.where(psp.tp > 1e10)]=np.nan
    
    print('done get psp')
    print()
    
    return psp

    
    
    
      
    
    
    

    
def save_psp_data_non_merged(path, file):
    '''
    *** TO DO
    save PSP data as pickle file with 3 separate arrays for orbit, magnetic field and plasma data    
    '''
     
    print('start PSP')
     
    psp_sat = heliosat.PSP()
    t_start = datetime.datetime(2018, 10, 14,14,14, 30)
    #t_end = datetime.datetime(2018, 12, 12,23,59,30)
    t_end = datetime.datetime(2019, 4, 23,23,59,30)

    #t_end = datetime.datetime(2019, 5, 31,23,59,30)
    #t_end = datetime.datetime(2019, 5, 1,23,59,30)

    
    timeb, mag = psp_sat.get_data_raw(t_start, t_end, "psp_fields_l2")
    timep, pro = psp_sat.get_data_raw(t_start, t_end, "psp_spc_l3")
    print('download complete')
    
    #create an array with 1 minute resolution between t start and end for position
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    
    
    #make separate arrays for orbit plasma and mag
    psp_orbit=np.zeros(np.size(time),dtype=[('time',object),('x', float),('y', float),\
                ('z', float),('r', float),('lat', float),('lon', float)])   
    #convert to recarray
    psp_orbit = psp_orbit.view(np.recarray)  

    psp_mag=np.zeros(np.size(timeb),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float)])   
    psp_mag = psp_mag.view(np.recarray)  

    psp_plasma=np.zeros(np.size(timep),dtype=[('time',object),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('np', float),('tp', float)])
    psp_plasma = psp_plasma.view(np.recarray)  


    psp_orbit.time=time
    psp_mag.time=parse_time(timeb,format='unix').datetime 
    psp_plasma.time=parse_time(timep,format='unix').datetime 
    print('time convert done')


    print('position start')
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('psp_pred'))
    psptra=spice.Trajectory('SPP')
    psptra.generate_positions(time,'Sun',frame)
    psptra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(psptra.x,psptra.y,psptra.z)
    psp_orbit.x=psptra.x
    psp_orbit.y=psptra.y
    psp_orbit.z=psptra.z
    psp_orbit.r=r
    psp_orbit.lat=np.rad2deg(lat)
    psp_orbit.lon=np.rad2deg(lon)
    print('position end')
   
    #fields
    psp_mag.bx = mag[:,0] 
    psp_mag.by = mag[:,1] 
    psp_mag.bz = mag[:,2] 
    psp_mag.bt = np.sqrt(psp_mag.bx**2+psp_mag.by**2+psp_mag.bz**2)
     
    #sweap
    from astropy.constants import m_p,k_B
    psp_plasma.np = pro[:,0]
    #https://en.wikipedia.org/wiki/Thermal_velocity convert from km/s to K

    psp_plasma.tp = np.pi*m_p*((pro[:,4]*1e3)**2)/(8*k_B) 

    psp_plasma.vx = pro[:,1]
    psp_plasma.vy = pro[:,2]
    psp_plasma.vz = pro[:,3]
    psp_plasma.vt=np.sqrt(psp_plasma.vx**2+psp_plasma.vy**2+psp_plasma.vz**2)

    
    header='PSP magnetic field (FIELDS instrument) and plasma data (SWEAP), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/psp/  '+ \
    'Timerange: '+psp_orbit.time[0].strftime("%Y-%b-%d %H:%M")+' to '+psp_orbit.time[-1].strftime("%Y-%b-%d %H:%M")+\
    '. The data are put in 3 numpy recarrays, fields can be accessed by psp_plasma.timep (for plasma), psp.vt etc.; psp_mag.timeb (for magnetic field),psp.bt, etc.; psp_orbit.time (for position) psp.r, psp.lon, ... '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(psp_plasma.size)+' (plasma), '+str(psp_mag.size)+' (mag). '+\
    'Units are btxyz [nT, RTN], vtxyz [km/s, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_psp_data_non_merged (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([psp_orbit,psp_mag,psp_plasma,header], open(path+file, "wb"))

    print('done psp')
    print()









def save_wind_data(path,file,start_date,end_date,heeq):
    
    '''
    description of data sources used in heliosat:
    https://cdaweb.sci.gsfc.nasa.gov/misc/NotesW.html
    '''
  
    print('start wind update')
    wind_sat = heliosat.WIND()
    t_start = start_date
    t_end = end_date
    
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=2*n) for n in range(int ((t_end - t_start).days*30*24))]  
    time_mat=mdates.date2num(time) 
    
    #print(parse_time(time[0:10]).iso)
    
    tm, mag = wind_sat.get_data_raw(t_start, t_end, "wind_mfi_k0", extra_columns=["PGSE"])#, return_datetimes=True)
    #tm, mag = wind_sat.get_data_raw(t_start, t_end, "wind_mfi_h0") k0
    tp, pro = wind_sat.get_data_raw(t_start, t_end, "wind_swe_h1")#,return_datetimes=True)
    
    tm=parse_time(tm,format='unix').datetime 
    tp=parse_time(tp,format='unix').datetime 
   
    
    #print(parse_time(tm[0:10]).iso)
    #print(parse_time(tp[0:10]).iso)
    
    tm=parse_time(tm).plot_date
    tp=parse_time(tp).plot_date
    
    print('download complete')
     
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, tm, mag[:,0] )
    by = np.interp(time_mat, tm, mag[:,1] )
    bz = np.interp(time_mat, tm, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
        
    den = np.interp(time_mat, tp, pro[:,0])
    vt = np.interp(time_mat, tp, pro[:,1])
    tp = np.interp(time_mat, tp, pro[:,2])
        
    #interpolate the GSE position over full data range
    x_gse = np.interp(time_mat, tm, mag[:,3])*6378.1/149597870.7*astropy.units.AU #earth radii to km to AU
    y_gse = np.interp(time_mat, tm, mag[:,4])*6378.1/149597870.7*astropy.units.AU
    z_gse = np.interp(time_mat, tm, mag[:,5])*6378.1/149597870.7*astropy.units.AU
    
    print('position start')    
    frame='HEEQ'
    planet_kernel=spicedata.get_kernel('planet_trajectories')
    earth=spice.Trajectory('399')  #399 for Earth
    earth.generate_positions(time,'Sun',frame)
    earth.change_units(astropy.units.AU)       #from km to AU    
    #add gse position to Earth position
    x=earth.x-x_gse  #earth radii to km
    y=earth.y-y_gse
    z=earth.z+z_gse
    [r, lat, lon]=cart2sphere(x,y,z)    
    print('position end ')
    
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

    win.x=x
    win.y=y
    win.z=z
    
    win.r=r
    win.lat=np.degrees(lat)
    win.lon=np.degrees(lon)
    
    win.np=den
    win.vt=vt
    #https://en.wikipedia.org/wiki/Thermal_velocity convert from km/s to K
    from astropy.constants import m_p,k_B
    win.tp=np.pi*m_p*((tp*1e3)**2)/(8*k_B)
    
    
    ############ spike removal
            
    #plasma    
    win.np[np.where(win.np> 500)]=1000000
    #get rid of all single spikes with scipy signal find peaks
    peaks, properties = scipy.signal.find_peaks(win.np, height=500,width=(1, 250))
    #go through all of them and set to nan according to widths
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        win.np[peaks[i]-width-2:peaks[i]+width+2]=np.nan

    win.tp[np.where(win.tp> 1e8)]=1e11
    #get rid of all single spikes with scipy signal find peaks
    peaks, properties = scipy.signal.find_peaks(win.tp, height=1e8,width=(1, 250))
    #go through all of them and set to nan according to widths
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        win.tp[peaks[i]-width-2:peaks[i]+width+2]=np.nan

    win.vt[np.where(win.vt> 3000)]=1e11
    #get rid of all single spikes with scipy signal find peaks
    peaks, properties = scipy.signal.find_peaks(win.vt, height=1e8,width=(1, 250))
    #go through all of them and set to nan according to widths
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        win.vt[peaks[i]-width-2:peaks[i]+width+2]=np.nan

        
        
        
    #magnetic field    
    peaks, properties = scipy.signal.find_peaks(win.bt, prominence=30,width=(1, 10))
    #go through all of them and set to nan according to widths
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths'])[i])
        #remove data
        win.bt[peaks[i]-width-5:peaks[i]+width+5]=np.nan    

    peaks, properties = scipy.signal.find_peaks(abs(win.bx), prominence=30,width=(1, 10))
    for i in np.arange(len(peaks)):
        width=int(np.ceil(properties['widths'])[i])
        win.bx[peaks[i]-width-5:peaks[i]+width+5]=np.nan    

    peaks, properties = scipy.signal.find_peaks(abs(win.by), prominence=30,width=(1, 10))
    for i in np.arange(len(peaks)):
        width=int(np.ceil(properties['widths'])[i])
        win.by[peaks[i]-width-5:peaks[i]+width+5]=np.nan    

    peaks, properties = scipy.signal.find_peaks(abs(win.bz), prominence=30,width=(1, 10))
    for i in np.arange(len(peaks)):
        width=int(np.ceil(properties['widths'])[i])
        win.bz[peaks[i]-width-5:peaks[i]+width+5]=np.nan    



    #manual spike removal for magnetic field
    if t_start < datetime.datetime(2018, 7, 19, 16, 25):    
        if t_end > datetime.datetime(2018, 7, 19, 16, 25):         

            remove_start=datetime.datetime(2018, 7, 19, 16, 25)
            remove_end=datetime.datetime(2018, 7, 19, 17, 35)
            remove_start_ind=np.where(remove_start<win.time)[0][0]
            remove_end_ind=np.where(remove_end<win.time)[0][0] 

            win.bt[remove_start_ind:remove_end_ind]=np.nan
            win.bx[remove_start_ind:remove_end_ind]=np.nan
            win.by[remove_start_ind:remove_end_ind]=np.nan
            win.bz[remove_start_ind:remove_end_ind]=np.nan

    if t_start < datetime.datetime(2018, 8, 29, 19, 00):    
        if t_end > datetime.datetime(2018, 8, 29, 19, 00):         

            remove_start=datetime.datetime(2018, 8, 29, 19, 00)
            remove_end=datetime.datetime(2018,8, 30, 5, 00)
            remove_start_ind=np.where(remove_start<win.time)[0][0]
            remove_end_ind=np.where(remove_end<win.time)[0][0] 

            win.bt[remove_start_ind:remove_end_ind]=np.nan
            win.bx[remove_start_ind:remove_end_ind]=np.nan
            win.by[remove_start_ind:remove_end_ind]=np.nan
            win.bz[remove_start_ind:remove_end_ind]=np.nan

            
    if t_start < datetime.datetime(2019, 8, 8, 22, 45):    
        if t_end > datetime.datetime(2019, 8, 8, 22, 45):         

            remove_start=datetime.datetime(2019, 8, 8, 22, 45)
            remove_end=datetime.datetime(2019,   8, 9, 17, 00)
            remove_start_ind=np.where(remove_start<win.time)[0][0]
            remove_end_ind=np.where(remove_end<win.time)[0][0] 

            win.bt[remove_start_ind:remove_end_ind]=np.nan
            win.bx[remove_start_ind:remove_end_ind]=np.nan
            win.by[remove_start_ind:remove_end_ind]=np.nan
            win.bz[remove_start_ind:remove_end_ind]=np.nan

    if t_start < datetime.datetime(2019, 8, 21, 22, 45):    
        if t_end > datetime.datetime(2019, 8, 21, 22, 45):         

            remove_start=datetime.datetime(2019, 8, 20, 18, 0)
            remove_end=datetime.datetime(2019,   8, 21, 12, 0)
            remove_start_ind=np.where(remove_start<win.time)[0][0]
            remove_end_ind=np.where(remove_end<win.time)[0][0] 

            win.bt[remove_start_ind:remove_end_ind]=np.nan
            win.bx[remove_start_ind:remove_end_ind]=np.nan
            win.by[remove_start_ind:remove_end_ind]=np.nan
            win.bz[remove_start_ind:remove_end_ind]=np.nan            

    if t_start < datetime.datetime(2019, 8, 21, 22, 45):    
        if t_end > datetime.datetime(2019, 8, 21, 22, 45):         

            remove_start=datetime.datetime(2019, 8, 22, 1, 0)
            remove_end=datetime.datetime(2019,   8, 22, 9, 0)
            remove_start_ind=np.where(remove_start<win.time)[0][0]
            remove_end_ind=np.where(remove_end<win.time)[0][0] 

            win.bt[remove_start_ind:remove_end_ind]=np.nan
            win.bx[remove_start_ind:remove_end_ind]=np.nan
            win.by[remove_start_ind:remove_end_ind]=np.nan
            win.bz[remove_start_ind:remove_end_ind]=np.nan            


    
     
    coord='GSE'
    #convert magnetic field to SCEQ
    if heeq==True:
        win=convert_GSE_to_HEEQ(win)
        coord='HEEQ'
   
    
        
    ###################### 

    
    header='Wind magnetic field (MAG instrument) and plasma data (SWE), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/wind/  '+ \
    'Timerange: '+win.time[0].strftime("%Y-%b-%d %H:%M")+' to '+win.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(win.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by win.time, win.bx, win.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(win.size)+'. '+\
    'Units are btxyz [nT, '+coord+'], vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_wind_data (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    pickle.dump([win,header], open(path+file, "wb"))
    
    
    print('wind update done')
    print()
    



    
    
    
    
    
    


def load_stereoa_science_1min():

    varnames = ['Epoch', 'Vp', 'Vr_Over_V_RTN', 'Np', 'Tp', 'BFIELDRTN']
    alldata = {k: [] for k in varnames}
    if not os.path.exists(heliosat_data_path+'sta_magplasma_outside_heliosat'):
        os.mkdir(heliosat_data_path+'sta_magplasma_outside_heliosat')
    for year in ['2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014','2015','2016','2017','2018','2019']:
        print('get STEREO-A yearly 1min data file for ',year)
        cdf_write = heliosat_data_path+'sta_magplasma_outside_heliosat/STA_L2_MAGPLASMA_1m_{}_V01.cdf'.format(year)
        if not os.path.exists(cdf_write):
            cdf_url = ("https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/ahead/magplasma/STA_L2_MAGPLASMA_1m_{}_V01.cdf".format(year))
            cdf_file = requests.get(cdf_url)
            open(cdf_write, 'wb').write(cdf_file.content)
        cdf = cdflib.CDF(cdf_write)
        #cdf.cdf_info() shows all variable names and attributes
        for var in varnames:
            data = cdf[var][...]
            data[np.where(data < cdf.varattsget(var)['VALIDMIN'][0])] = np.NaN
            data[np.where(data > cdf.varattsget(var)['VALIDMAX'][0])] = np.NaN
            alldata[var].append(data)
    arrays = {}
    for var in varnames:
        arrays[var] = np.concatenate(alldata[var])
        
    return arrays



def save_all_stereoa_science_data(path,file,sceq):
    ''' 
    saves all STEREO-Ahead science data btxyz 
    vt np tp x y z r lat lon 1 min resolution as pickle
    sceq=True -> convert RTN to SCEQ coordinates for magnetic field components
    
    filesta_all='stereoa_2007_2020_sceq.p'
    hd.save_all_stereoa_science_data(data_path, filesta_all,sceq=True)

    filesta_all='stereoa_2007_2020.p'
    hd.save_all_stereoa_science_data(data_path, filesta_all,sceq=False)    
    
    
    [sta_t,hsta_t]=pickle.load(open(data_path+filesta_all, "rb" ) )
    '''
    
    #load all data with function
    s1=load_stereoa_science_1min()    
    print('download complete')

    #make array
    sta=np.zeros(len(s1['Epoch']),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
    
    #convert to recarray
    sta = sta.view(np.recarray)  
    
    #parse Epoch time to datetime objects
    sta.time=parse_time(s1['Epoch'],format='cdf_epoch').datetime  
    print('time conversion complete')

    sta.bx=s1['BFIELDRTN'][:,0]
    sta.by=s1['BFIELDRTN'][:,1]
    sta.bz=s1['BFIELDRTN'][:,2]
    sta.bt=np.sqrt(sta.bx**2+sta.by**2+sta.bz**2)
    
    sta.vt=s1['Vp']    
    sta.np=s1['Np']
    sta.tp=s1['Tp']    
    
    print('parameters into array complete')
    
    print('position start')
    frame='HEEQ'
    kernels = spicedata.get_kernel('stereo_a')
    kernels += spicedata.get_kernel('stereo_a_pred')
    spice.furnish(kernels)
    statra=spice.Trajectory('-234') #STEREO-A SPICE NAIF code
    statra.generate_positions(sta.time,'Sun',frame)
    statra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(statra.x,statra.y,statra.z)
    
    sta.x=statra.x
    sta.y=statra.y
    sta.z=statra.z
    
    sta.r=r
    sta.lat=np.degrees(lat)
    sta.lon=np.degrees(lon)

    print('position end ')
    
    
    
    #remove spike in magnetic field in 2015
    spike_ind=np.where(sta.bt >300)[0]
    if len(spike_ind) > 0:
        sta.bt[spike_ind[0]-77:spike_ind[-1]+5]=np.nan
        sta.bx[spike_ind[0]-77:spike_ind[-1]+5]=np.nan
        sta.by[spike_ind[0]-77:spike_ind[-1]+5]=np.nan
        sta.bz[spike_ind[0]-77:spike_ind[-1]+5]=np.nan



    
    coord='RTN'
    #convert magnetic field to SCEQ
    if sceq==True:
        print('convert RTN to SCEQ ')
        coord='SCEQ'
        sta=convert_RTN_to_SCEQ(sta,'STEREO-A')
    
    header='STEREO-A magnetic field (IMPACT instrument) and plasma data (PLASTIC, science), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/ahead/magplasma   '+ \
    'Timerange: '+sta.time[0].strftime("%Y-%b-%d %H:%M")+' to '+sta.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', with an average time resolution of '+str(np.mean(np.diff(sta.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by sta.time, sta.bx, sta.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(sta.size)+'. '+\
    'Units are btxyz [nT, '+coord+', vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats '+\
    'and https://github.com/heliopython/heliopy. '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, R. L. Bailey and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    print('save pickle file')
    pickle.dump([sta,header], open(path+file, "wb"))
    
    print('done sta')
    print()

   








def load_stereob_science_1min():

    varnames = ['Epoch', 'Vp', 'Vr_Over_V_RTN', 'Np', 'Tp', 'BFIELDRTN']
    alldata = {k: [] for k in varnames}
    if not os.path.exists(heliosat_data_path+'stb_magplasma_outside_heliosat'):
        os.mkdir(heliosat_data_path+'stb_magplasma_outside_heliosat')
    for year in ['2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014']:
    #for year in ['2007']:
        print('get STEREO-B yearly 1min data file for ',year)
        cdf_write = heliosat_data_path+'stb_magplasma_outside_heliosat/STB_L2_MAGPLASMA_1m_{}_V01.cdf'.format(year)
        if not os.path.exists(cdf_write):
            cdf_url = ("https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/behind/magplasma/STB_L2_MAGPLASMA_1m_{}_V01.cdf".format(year))
            cdf_file = requests.get(cdf_url)
            open(cdf_write, 'wb').write(cdf_file.content)
        cdf = cdflib.CDF(cdf_write)
        #cdf.cdf_info() shows all variable names and attributes
        for var in varnames:
            data = cdf[var][...]
            #fillval = cdf[var].attrs['FILLVAL']
            #fillval=cdf.varattsget(var)['FILLVAL'][0]
            data[np.where(data < cdf.varattsget(var)['VALIDMIN'][0])] = np.NaN
            data[np.where(data > cdf.varattsget(var)['VALIDMAX'][0])] = np.NaN
            alldata[var].append(data)
    arrays = {}
    for var in varnames:
        arrays[var] = np.concatenate(alldata[var])
        
    return arrays


def save_all_stereob_science_data(path,file,sceq):
    ''' 
    saves all STEREO-Behind science data btxyz 
    vt np tp x y z r lat lon 1 min resolution as pickle
    sceq=True -> convert RTN to SCEQ coordinates for magnetic field components
    
    use as:
    filestb_all='stereob_2007_2014_sceq.p'
    hd.save_all_stereob_science_data(data_path, filestb_all,sceq=True)
    
    filestb_all='stereob_2007_2014.p'
    hd.save_all_stereob_science_data(data_path, filestb_all,sceq=False)  
    
    [stb_t,hstb_t]=pickle.load(open(data_path+filestb_all, "rb" ) )    
    '''
    
    #load all data with function
    s1=load_stereob_science_1min()    
    print('download complete')

    #make array
    stb=np.zeros(len(s1['Epoch']),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
    
    #convert to recarray
    stb = stb.view(np.recarray)  
   
    
    #parse Epoch time to datetime objects
    stb.time=parse_time(s1['Epoch'],format='cdf_epoch').datetime  
    
    print('time conversion complete')

    stb.bx=s1['BFIELDRTN'][:,0]
    stb.by=s1['BFIELDRTN'][:,1]
    stb.bz=s1['BFIELDRTN'][:,2]
    stb.bt=np.sqrt(stb.bx**2+stb.by**2+stb.bz**2)

    
    stb.vt=s1['Vp']    
    stb.np=s1['Np']
    stb.tp=s1['Tp']    
    
    print('parameters into array complete')
    
    print('position start')
    frame='HEEQ'
    kernels = spicedata.get_kernel('stereo_b')
    spice.furnish(kernels)
    stbtra=spice.Trajectory('-235') #STEREO-A SPICE NAIF code
    stbtra.generate_positions(stb.time,'Sun',frame)
    stbtra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(stbtra.x,stbtra.y,stbtra.z)
    
    stb.x=stbtra.x
    stb.y=stbtra.y
    stb.z=stbtra.z
    
    stb.r=r
    stb.lat=np.rad2deg(lat)
    stb.lon=np.rad2deg(lon)

    print('position end ')

    
    coord='RTN'
    #convert magnetic field to SCEQ
    if sceq==True:
        print('convert RTN to SCEQ ')
        coord='SCEQ'
        stb=convert_RTN_to_SCEQ(stb,'STEREO-B')
    
    
    header='STEREO-B magnetic field (IMPACT instrument) and plasma data (PLASTIC, science), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/behind/magplasma   '+ \
    'Timerange: '+stb.time[0].strftime("%Y-%b-%d %H:%M")+' to '+stb.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', with a an average time resolution of '+str(np.mean(np.diff(stb.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by stb.time, stb.bx, stb.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(stb.size)+'. '+\
    'Units are btxyz [nT, '+coord+'], vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats '+\
    'and https://github.com/heliopython/heliopy. '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, R. L. Bailey and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
  
    print('save pickle file')
    pickle.dump([stb,header], open(path+file, "wb"))
    
    print('done stb')
    print()


 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 

def get_sdo_realtime_image():
    """Downloads latest SDO image."""


    sdo_latest='https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_0193.jpg'
    #PFSS
    #sdo_latest='https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_0193pfss.jpg'
    try: urllib.request.urlretrieve(sdo_latest,data_path_sun+'latest_1024_0193.jpg')
    except urllib.error.URLError as e:
        print('Failed downloading ', sdo_latest,' ',e)

    print('saved ',data_path_sun+'latest_1024_0193.jpg')    
    
    
    sdo_latest='https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_HMIB.jpg'        
    try: urllib.request.urlretrieve(sdo_latest,data_path_sun+'latest_1024_HMIB.jpg')
    except urllib.error.URLError as e:
        print('Failed downloading ', sdo_latest,' ',e)
        
    print('saved ',data_path_sun+'latest_1024_HMIB.jpg')    
        
    '''    
    #convert to png
    #check if ffmpeg is available locally in the folder or systemwide
    if os.path.isfile('ffmpeg'):
        os.system('./ffmpeg -i latest_1024_0193.jpg latest_1024_0193.png -loglevel quiet -y')
        ffmpeg_avail=True
        logger.info('downloaded SDO latest_1024_0193.jpg converted to png')
        os.system('rm latest_1024_0193.jpg')
    else:
        os.system('ffmpeg -i latest_1024_0193.jpg latest_1024_0193.png -loglevel quiet -y')
        os.system('rm latest_1024_0193.jpg')
    '''    

 
def save_noaa_rtsw_data(data_path,noaa_path,filenoaa):


    print(' ')
    print('convert NOAA real time solar wind archive to pickle file')
    items=os.listdir(noaa_path)  
    newlist = [] 
    for names in items: 
       if names.endswith(".json"):
            newlist.append(names)
    #print(newlist)

    a=sorted(newlist) #sort so that mag and plasma and dates are separated
    #print(a)
    nr_of_files=int(np.size(a)/2)#******************
    mag=a[0:nr_of_files]  
    pla=a[nr_of_files:-1]  

    #make array for 10 years
    noaa=np.zeros(5000000,dtype=[('time',object),('bx', float),('by', float),\
                    ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                    ('x', float),('y', float),('z', float),\
                    ('r', float),('lat', float),('lon', float)])   

    k=0
    
    

    for i in np.arange(nr_of_files)-1:

        #read in data of corresponding files
        #print(noaa_path+mag[i])
        m1=open(noaa_path+mag[i],'r')
        p1=open(noaa_path+pla[i],'r')
        try: 
            d1=get_noaa_realtime_data(m1, p1)
            #save in large array
            noaa[k:k+np.size(d1)]=d1
            k=k+np.size(d1) 
        except: 
            print(mag[i], ' ', pla[i], ' json not working')

    

    #cut zeros, sort, convert to recarray, and find unique times and data

    noaa_cut=noaa[0:k]
    noaa_cut.sort()
         
    nu=noaa_cut.view(np.recarray)
    [dum,ind]=np.unique(nu.time,return_index=True)  
    nf=nu[ind]

    header='Real time solar wind magnetic field and plasma data from NOAA, ' + \
        'obtained daily from https://services.swpc.noaa.gov/products/solar-wind/  '+ \
        'Timerange: '+nf.time[0].strftime("%Y-%b-%d %H:%M")+' to '+nf.time[-1].strftime("%Y-%b-%d %H:%M")+\
        ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(nf.time)).seconds)+' seconds. '+\
        'The data are available in a numpy recarray, fields can be accessed by nf.time, nf.bx, nf.vt etc. '+\
        'Total number of data points: '+str(nf.size)+'. '+\
        'Units are btxyz [nT, RTN], vt  [km s^-1], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
        'Made with https://github.com/cmoestl/heliocats save_noaa_rtsw_data  '+\
        'By C. Moestl (twitter @chrisoutofspace) and R. Bailey. File creation date: '+\
        datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'


    pickle.dump([nf,header], open(data_path+filenoaa, "wb"))
    
    #to read file
    #import pickle
    #filenoaa='noaa_rtsw_2020.p'
    #data_path='/nas/helio/data/insitu_python/'
    #[n,hn]=pickle.load(open(data_path+filenoaa, "rb" ) ) 
    
    print('NOAA done')        



def get_noaa_realtime_data(magfile, plasmafile):
    """
    Downloads and returns noaa real time solar wind data 
    data from http://services.swpc.noaa.gov/products/solar-wind/
    if needed replace with ACE
    http://legacy-www.swpc.noaa.gov/ftpdir/lists/ace/
    get 3 or 7 day data
    url_plasma='http://services.swpc.noaa.gov/products/solar-wind/plasma-3-day.json'
    url_mag='http://services.swpc.noaa.gov/products/solar-wind/mag-3-day.json'
    
    Author: R. Bailey, modified for heliocats by C. Moestl
    
    Parameters
    ==========
    None
    Returns: recarray with interpolated data
    =======
    """
    
    # Read plasma data:
    dp = json.loads(plasmafile.read())
    dpn = [[np.nan if x == None else x for x in d] for d in dp]     # Replace None w NaN
    dtype=[(x, 'float') for x in dp[0]]
    datesp = [datetime.datetime.strptime(x[0], "%Y-%m-%d %H:%M:%S.%f")  for x in dpn[1:]]
    #convert datetime to matplotlib times
    mdatesp=mdates.date2num(datesp)
    dp_ = [tuple([d]+[float(y) for y in x[1:]]) for d, x in zip(mdatesp, dpn[1:])] 
    DSCOVR_P = np.array(dp_, dtype=dtype)

    
    # Read magnetic field data:
    dm = json.loads(magfile.read())
    dmn = [[np.nan if x == None else x for x in d] for d in dm]     # Replace None w NaN
    dtype=[(x, 'float') for x in dmn[0]]
    datesm = [datetime.datetime.strptime(x[0], "%Y-%m-%d %H:%M:%S.%f")  for x in dmn[1:]]
    mdatesm=mdates.date2num(datesm)
    dm_ = [tuple([d]+[float(y) for y in x[1:]]) for d, x in zip(mdatesm, dm[1:])] 
    DSCOVR_M = np.array(dm_, dtype=dtype)
    
    

    #first_timestep = np.max([mdatesp[-1], mdatesm[-1]])
    #last_timestep = np.min([mdatesp[-1], mdatesm[-1]])
    #nminutes = int((num2date(last_timestep)-num2date(first_timestep)).total_seconds()/60.)
    #itime = np.asarray([date2num(num2date(first_timestep) + timedelta(minutes=i)) for i in range(nminutes)], dtype=np.float64)
    
    #use mag for times
    t_start=datesm[0]
    t_end=datesm[-1]

    #1 minute res
    itime = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int (((t_end - t_start).days+1)*60*24))]   #*******BUG everywhere with this line for last day

        
    itimeint=mdates.date2num(itime)
    
    rbtot_m = np.interp(itimeint, DSCOVR_M['time_tag'], DSCOVR_M['bt'])
    rbxgsm_m = np.interp(itimeint, DSCOVR_M['time_tag'], DSCOVR_M['bx_gsm'])
    rbygsm_m = np.interp(itimeint, DSCOVR_M['time_tag'], DSCOVR_M['by_gsm'])
    rbzgsm_m = np.interp(itimeint, DSCOVR_M['time_tag'], DSCOVR_M['bz_gsm'])
    rpv_m = np.interp(itimeint, DSCOVR_P['time_tag'], DSCOVR_P['speed'])
    rpn_m = np.interp(itimeint, DSCOVR_P['time_tag'], DSCOVR_P['density'])
    rpt_m = np.interp(itimeint, DSCOVR_P['time_tag'], DSCOVR_P['temperature'])

    #make array
    dscovr_data=np.zeros(np.size(rbtot_m),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
                        
    #convert to recarray
    dscovr_data = dscovr_data.view(np.recarray)                             

    dscovr_data.time=itime
    dscovr_data.bt=rbtot_m
    dscovr_data.bx=rbxgsm_m
    dscovr_data.by=rbygsm_m
    dscovr_data.bz=rbzgsm_m

    dscovr_data.vt=rpv_m
    dscovr_data.np=rpn_m
    dscovr_data.tp=rpt_m
    
    
    
    #print('position start')
    frame='HEEQ'
    planet_kernel=spicedata.get_kernel('planet_trajectories')
    earth=spice.Trajectory('399')  #399 for Earth, not barycenter (because of moon)
    earth.generate_positions(itime,'Sun',frame)
    #from km to AU
    earth.change_units(astropy.units.AU)
    #add gse position to Earth position
    x=earth.x-1.5*1e6*astropy.units.km
    y=earth.y
    z=earth.z
    [r, lat, lon]=cart2sphere(x,y,z)
    #*****with astropy lagrange points exact value? L1 position with 0.01 AU 
    #[r, lat, lon]=cart2sphere(earth.x-0.01*astropy.units.AU,earth.y,earth.z)
    #print('position end ')
       
    
    
    dscovr_data.x=x
    dscovr_data.y=y
    dscovr_data.z=z
    
    dscovr_data.r=r
    dscovr_data.lat=np.rad2deg(lat)
    dscovr_data.lon=np.rad2deg(lon)
       
    
    
    print('NOAA data read completed for file with end time: ',itime[-1])
    
    return dscovr_data





def save_stereob_beacon_data(path,file,start_time,end_time):

    print('start STB')
    stb_sat = heliosat.STB()
    t_start = start_time
    t_end = end_time
   
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
  
    tm, mag = stb_sat.get_data_raw(t_start, t_end, "mag_beacon")
    tp, pro = stb_sat.get_data_raw(t_start, t_end, "proton_beacon")

    print('download complete')
   
    tm=parse_time(tm,format='unix').datetime 
    tp=parse_time(tp,format='unix').datetime 

    #convert to matplotlib time for linear interpolation
    tm_mat=mdates.date2num(tm) 
    tp_mat=mdates.date2num(tp) 
    
    print('time convert done')
    
    print('position start')
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('stereo_b'))
    stbtra=spice.Trajectory('-235') #STEREO-A SPICE NAIF code
    stbtra.generate_positions(time,'Sun',frame)
    stbtra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(stbtra.x,stbtra.y,stbtra.z)
    print('position end ')
    
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
      
    den = np.interp(time_mat, tp_mat, pro[:,0])
    vt = np.interp(time_mat, tp_mat, pro[:,1])
    tp = np.interp(time_mat, tp_mat, pro[:,2])
    
    #make array
    stb=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    stb = stb.view(np.recarray)  

    #fill with data
    stb.time=time
    stb.bx=bx
    stb.by=by
    stb.bz=bz 
    stb.bt=bt

    stb.x=stbtra.x
    stb.y=stbtra.y
    stb.z=stbtra.z
    
    stb.r=r
    stb.lat=np.rad2deg(lat)
    stb.lon=np.rad2deg(lon)
    
    stb.np=den
    stb.tp=tp    
    stb.vt=vt  
        
    #remove spikes from plasma data
    #median filter
    stb.vt=scipy.signal.medfilt(stb.vt,9)
    #set nans to a high number
    stb.vt[np.where(np.isfinite(stb.vt) == False)]=1e5
    #get rid of all single spikes with scipy signal find peaks (cannot use nan)
    peaks,properties = scipy.signal.find_peaks(stb.vt, prominence=200,width=(1,200))
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        stb.vt[peaks[i]-width-2:peaks[i]+width+2]=np.nan
    #set nan again
    stb.vt[np.where(stb.vt == 1e5)]=np.nan     
    stb.tp[np.where(np.isfinite(stb.vt) == False)]=np.nan
    stb.np[np.where(np.isfinite(stb.vt) == False)]=np.nan   
    
    #remove spikes from magnetic field data
    #median filter
    #set nans to a high number
    stb.bt[np.where(np.isfinite(stb.bt) == False)]=1e5
    #get rid of all single spikes with scipy signal find peaks (cannot use nan)
    peaks,properties = scipy.signal.find_peaks(stb.bt, height=40,width=(1,20))
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths'])[i])
        #remove data
        stb.bt[peaks[i]-width-2:peaks[i]+width+2]=np.nan
        stb.bx[peaks[i]-width-2:peaks[i]+width+2]=np.nan
        stb.by[peaks[i]-width-2:peaks[i]+width+2]=np.nan
        stb.bz[peaks[i]-width-2:peaks[i]+width+2]=np.nan

    #set nan again
    stb.bt[np.where(stb.bt == 1e5)]=np.nan     

    #manual spike removal for speed
    remove_start=datetime.datetime(2007, 7, 18,22, 00)
    remove_end=datetime.datetime(2007, 7, 19, 16, 00)
    remove_start_ind=np.where(remove_start==stb.time)[0][0]
    remove_end_ind=np.where(remove_end==stb.time)[0][0] 

    stb.vt[remove_start_ind:remove_end_ind]=np.nan
    
    
    header='BEACON STEREO-B magnetic field (IMPACT instrument) and plasma data (PLASTIC), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/beacon/behind/  '+ \
    'Timerange: '+stb.time[0].strftime("%Y-%b-%d %H:%M")+' to '+stb.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(stb.time)).seconds)+' seconds. '+\
    'A median filter has been applied (plasma data only) and then spikes were removed with scipy.signal.find_peaks (plasma and field). '+\
    'The data are available in a numpy recarray, fields can be accessed by stb.time, stb.bx, stb.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(stb.size)+'. '+\
    'Units are btxyz [nT, RTN], vt  [km s^-1], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_stereob_beacon_data (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([stb,header], open(path+file, "wb"))
    
    print('done stb')
    print()





 
def save_stereoa_beacon_data(path,file,start_time,end_time,sceq):

    print('start STA')
    sta_sat = heliosat.STA()
    t_start = start_time
    t_end = end_time
 
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 

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

    den = np.interp(time_mat, tp_mat, pro[:,0])
    vt = np.interp(time_mat, tp_mat, pro[:,1])
    tp = np.interp(time_mat, tp_mat, pro[:,2])
    
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
    sta.vt=vt  
    
    
    #remove spikes from plasma data
    #median filter
    sta.vt=scipy.signal.medfilt(sta.vt,9)
    #set nans to a high number
    sta.vt[np.where(np.isfinite(sta.vt) == False)]=1e5
    #get rid of all single spikes with scipy signal find peaks (cannot use nan)
    peaks,properties = scipy.signal.find_peaks(sta.vt, prominence=200,width=(1,200))
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        sta.vt[peaks[i]-width-2:peaks[i]+width+2]=np.nan
    #set nan again
    sta.vt[np.where(sta.vt == 1e5)]=np.nan     
    sta.tp[np.where(np.isfinite(sta.vt) == False)]=np.nan
    sta.np[np.where(np.isfinite(sta.vt) == False)]=np.nan   
    
    
    #manual spike removal for magnetic field
    #remove_start=datetime.datetime(2018, 9, 23, 11, 00)
    #remove_end=datetime.datetime(2018, 9, 25, 00, 00)
    #remove_start_ind=np.where(remove_start==sta.time)[0][0]
    #remove_end_ind=np.where(remove_end==sta.time)[0][0] 

    #sta.bt[remove_start_ind:remove_end_ind]=np.nan
    #sta.bx[remove_start_ind:remove_end_ind]=np.nan
    #sta.by[remove_start_ind:remove_end_ind]=np.nan
    #sta.bz[remove_start_ind:remove_end_ind]=np.nan
    
        
    coord='RTN'
    #convert magnetic field to SCEQ
    if sceq==True:
        sta=convert_RTN_to_SCEQ(sta,'STEREO-A')
        coord='SCEQ'
   
    
    
    header='BEACON STEREO-A magnetic field (IMPACT instrument) and plasma data (PLASTIC), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/  '+ \
    'Timerange: '+sta.time[0].strftime("%Y-%b-%d %H:%M")+' to '+sta.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(sta.time)).seconds)+' seconds. '+\
    'A median filter has been applied as scipy.signal.medfilt(sta.vt,9) and then spikes were removed with '+\
    'scipy.signal.find_peaks(sta.vt, prominence=200,width=(1,200)). '+\
    'The data are available in a numpy recarray, fields can be accessed by sta.time, sta.bx, sta.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(sta.size)+'. '+\
    'Units are btxyz [nT, '+coord+'], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_stereoa_beacon_data (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    

    pickle.dump([sta,header], open(path+file, "wb"))
    
    print('done sta')
    print()






def save_stereoa_science_data(path,file,t_start, t_end,sceq):
    
    '''** TO DO 
    '''

    print('start STA')
    sta_sat = heliosat.STA()
     
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    tp, pro = sta_sat.get_data_raw(t_start, t_end, "sta_plastic_l2")
    #tm, mag = sta_sat.get_data_raw(t_start, t_end, "sta_impact_beacon")
    #tm, mag = sta_sat.get_data_raw(t_start, t_end, "sta_impact_l1")

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
      
      
    den = np.interp(time_mat, tp_mat, pro[:,0])
    vt = np.interp(time_mat, tp_mat, pro[:,1])
    tp = np.interp(time_mat, tp_mat, pro[:,2])
    
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
    
    sta.vt=vt    
    sta.np=den
    sta.tp=tp    
    
     
     
    #remove spikes from plasma data
    #median filter
    sta.vt=scipy.signal.medfilt(sta.vt,9)
    #set nans to a high number
    sta.vt[np.where(np.isfinite(sta.vt) == False)]=1e5
    #get rid of all single spikes with scipy signal find peaks (cannot use nan)
    peaks,properties = scipy.signal.find_peaks(sta.vt, prominence=200,width=(1,200))
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        sta.vt[peaks[i]-width-2:peaks[i]+width+2]=np.nan
    #set nan again
    sta.vt[np.where(sta.vt == 1e5)]=np.nan     
    sta.tp[np.where(np.isfinite(sta.vt) == False)]=np.nan
    sta.np[np.where(np.isfinite(sta.vt) == False)]=np.nan   
    
    #manual spike removal for magnetic field
    remove_start=datetime.datetime(2018, 9, 23, 11, 00)
    remove_end=datetime.datetime(2018, 9, 25, 00, 00)
    remove_start_ind=np.where(remove_start==sta.time)[0][0]
    remove_end_ind=np.where(remove_end==sta.time)[0][0] 

    sta.bt[remove_start_ind:remove_end_ind]=np.nan
    sta.bx[remove_start_ind:remove_end_ind]=np.nan
    sta.by[remove_start_ind:remove_end_ind]=np.nan
    sta.bz[remove_start_ind:remove_end_ind]=np.nan

    
    
    #convert magnetic field to SCEQ
    if sceq==True:
        sta=convert_RTN_to_SCEQ(sta,'STEREO-A')
    
    
    
    
    
    header='STEREO-A magnetic field (IMPACT instrument, beacon) and plasma data (PLASTIC, science), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/ins_data/   '+ \
    'Timerange: '+sta.time[0].strftime("%Y-%b-%d %H:%M")+' to '+sta.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(sta.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by sta.time, sta.bx, sta.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(sta.size)+'. '+\
    'Units are btxyz [nT, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
  
    pickle.dump([sta,header], open(path+file, "wb"))
    
    print('done sta')
    print()







def omni_loader(overwrite):
   '''
   downloads all omni2 data into the "data" folder
   '''

   if overwrite>0: 
         print('download OMNI2 again')
         os.remove('data/omni2_all_years.dat')
  
   if not os.path.exists('data/omni2_all_years.dat'):
      #see http://omniweb.gsfc.nasa.gov/html/ow_data.html
      print('OMNI2 .dat file not in "data" directory, so download OMNI2 data from')
      omni2_url='https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
      print(omni2_url)
      try: urllib.request.urlretrieve(omni2_url, 'data/omni2_all_years.dat')
      except urllib.error.URLError as e:
          print(' ', omni2_url,' ',e.reason)
          sys.exit()


def save_omni_data(path,file,overwrite):
    '''
    save variables from OMNI2 dataset as pickle

    documentation https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2.text
    omni2_url='https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
    '''        
  
    print('start omni')
    
    omni_loader(overwrite)
    #check how many rows exist in this file
    f=open('data/omni2_all_years.dat')
    dataset= len(f.readlines())
    
    #make array
    o=np.zeros(dataset,dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bygsm', float),('bzgsm', float),('bt', float),\
                ('vt', float),('np', float),('tp', float),('alpha', float),\
                ('dst', float),('kp', float),('spot', float),\
                ('ae', float),('ap', float),('f107', float),\
                ('pcn', float),('al', float),('au', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])
                
    o=o.view(np.recarray)  
    print(dataset, ' datapoints')   #for reading data from OMNI file
    
    j=0
    with open('data/omni2_all_years.dat') as f:
        for line in f:
            line = line.split() # to deal with blank 
            
            #time - need to convert from year doy hour to datetime object
            o.time[j]=datetime.datetime(int(line[0]), 1, 1) + datetime.timedelta(int(line[1]) - 1) \
                              + datetime.timedelta(hours=int(line[2])) 

            #25 is bulkspeed F6.0, in km/s
            o.vt[j]=line[24]
            if o.vt[j] == 9999: o.vt[j]=np.NaN
            
            #24 in file, index 23 proton density /ccm
            o.np[j]=line[23]
            if o.np[j] == 999.9: o.np[j]=np.NaN
            
            #23 in file, index 22 Proton temperature  /ccm
            o.tp[j]=line[22]
            if o.tp[j] == 9999999.: o.tp[j]=np.NaN

            #28 in file, index 27 alpha to proton ratio
            o.alpha[j]=line[27]
            if o.alpha[j] == 9.999: o.alpha[j]=np.NaN
           
            #9 is total B  F6.1 also fill ist 999.9, in nT
            o.bt[j]=line[9]
            if o.bt[j] == 999.9: o.bt[j]=np.NaN

            #GSE components from 13 to 15, so 12 to 14 index, in nT
            o.bx[j]=line[12]
            if o.bx[j] == 999.9: o.bx[j]=np.NaN
            o.by[j]=line[13]
            if o.by[j] == 999.9: o.by[j]=np.NaN
            o.bz[j]=line[14]
            if o.bz[j] == 999.9: o.bz[j]=np.NaN
          
            #GSM
            o.bygsm[j]=line[15]
            if o.bygsm[j] == 999.9: o.bygsm[j]=np.NaN
          
            o.bzgsm[j]=line[16]
            if o.bzgsm[j] == 999.9: o.bzgsm[j]=np.NaN    
          

            o.kp[j]=line[38]
            if o.kp[j] == 99: o.kp[j]=np.nan
           
            o.spot[j]=line[39]
            if o.spot[j] ==  999: o.spot[j]=np.nan

            o.dst[j]=line[40]
            if o.dst[j] == 99999: o.dst[j]=np.nan

            o.ae[j]=line[41]
            if o.ae[j] ==  9999: o.ae[j]=np.nan
            
            o.ap[j]=line[49]
            if o.ap[j] ==  999: o.ap[j]=np.nan
            
            o.f107[j]=line[50]
            if o.f107[j] ==  999.9  : o.f107[j]=np.nan
            
            o.pcn[j]=line[51]
            if o.pcn[j] ==  999.9  : o.pcn[j]=np.nan

            o.al[j]=line[52]
            if o.al[j] ==   99999  : o.al[j]=np.nan

            o.au[j]=line[53]
            if o.au[j] ==  99999 : o.au[j]=np.nan
          
          
            j=j+1     
            
    print('position start')
    frame='HEEQ'
    planet_kernel=spicedata.get_kernel('planet_trajectories')
    earth=spice.Trajectory('399')  #399 for Earth, not barycenter (because of moon)
    earth.generate_positions(o.time,'Sun',frame)
    earth.change_units(astropy.units.AU)
    [r, lat, lon]=cart2sphere(earth.x,earth.y,earth.z)
    print('position end ')   
    
    
    o.x=earth.x
    o.y=earth.y
    o.z=earth.z
    
    o.r=r
    o.lat=np.rad2deg(lat)
    o.lon=np.rad2deg(lon)
    
    header='Near Earth OMNI2 1 hour solar wind and geomagnetic indices data since 1963. ' + \
    'Obtained from https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/  '+ \
    'Timerange: '+o.time[0].strftime("%Y-%b-%d %H:%M")+' to '+o.time[-1].strftime("%Y-%b-%d %H:%M")+'. '+\
    'The data are available in a numpy recarray, fields can be accessed by o.time, o.bx, o.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(o.size)+'. '+\
    'For units and documentation see: https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2.text, the '+\
    'heliospheric position of Earth was added and is given in x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_omni_data (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    
    pickle.dump([o,header], open(path+file, 'wb') )
    
    print('done omni')
    print()

    




def convert_MAVEN_mat_original(data_path,filename):

    print('load MAVEN from MAT')
    
    file=data_path+'input/MAVEN_2014to2018_cyril.mat'
    mavraw = scipy.io.loadmat(file)
    
    #make array 
    mav=np.zeros(np.size(mavraw['BT']),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('tp', float),('np', float),('r', float),('lat', float),\
                ('lon', float),('x', float),('y', float),('z', float),\
                ('ro', float), ('lato', float), ('lono', float),\
                ('xo', float), ('yo', float), ('zo', float)])   
    #convert to recarray
    mav = mav.view(np.recarray)  
    #convert time from matlab to python
    t=mavraw['timeD'][:,0]    
    
    for p in np.arange(np.size(t)):
        mav.time[p]= datetime.datetime.fromordinal(t[p].astype(int) ) + \
        datetime.timedelta(days=t[p]%1) - datetime.timedelta(days = 366) 
        

    mav.bx=mavraw['Bx'][:,0]       
    mav.by=mavraw['By'][:,0] 
    mav.bz=mavraw['Bz'][:,0]      
    mav.bt=mavraw['BT'][:,0]      

    mav.vx=mavraw['Vx'][:,0]      
    mav.vy=mavraw['Vy'][:,0]      
    mav.vz=mavraw['Vz'][:,0]      
    mav.vt=mavraw['VT'][:,0] 

    mav.tp=mavraw['Tp'][:,0]*(1.602176634*1e-19)/(1.38064852*1e-23)        #from ev to K     
    mav.np=mavraw['np'][:,0]      
    
    #add position with respect to Mars center in km in MSO
    print('orbit position start')
    
    insertion=datetime.datetime(2014,9,22,2,24,0)
    #these are the indices of the times for the cruise phase     
    tc=np.where(mdates.date2num(mav.time) < mdates.date2num(insertion))       

    mars_radius=3389.5
    mav.xo=mavraw['Xsc'][:,0]*mars_radius
    mav.yo=mavraw['Ysc'][:,0]*mars_radius
    mav.zo=mavraw['Zsc'][:,0]*mars_radius  
    
    #set to nan for cruise phase
    mav.xo[tc]=np.nan
    mav.yo[tc]=np.nan
    mav.zo[tc]=np.nan
    
    [mav.ro,mav.lato,mav.lono]=cart2sphere(mav.xo,mav.yo,mav.zo)
    mav.lono=np.rad2deg(mav.lono)   
    mav.lato=np.rad2deg(mav.lato)

    print('HEEQ position start')
    frame='HEEQ'
    
    #add position in HEEQ for cruise phase and orbit
    
    #cruise phase
    #use heliopy to load own bsp spice file from MAVEN 
    #obtained through https://naif.jpl.nasa.gov/pub/naif/pds/pds4/maven/maven_spice/spice_kernels/spk/
    spice.furnish(data_path+'input/maven_cru_rec_131118_140923_v1.bsp') 
    cruise=spice.Trajectory('MAVEN') #or NAIF CODE -202 
    cruise.generate_positions(mav.time[tc],'Sun',frame)     
    cruise.change_units(astropy.units.AU)  
    mav.x[tc]=cruise.x
    mav.y[tc]=cruise.y
    mav.z[tc]=cruise.z
    [mav.r[tc], mav.lat[tc], mav.lon[tc]]=cart2sphere(mav.x[tc],mav.y[tc],mav.z[tc])

    #times in orbit
    to=np.where(mdates.date2num(mav.time) > mdates.date2num(insertion))       

    planet_kernel=spicedata.get_kernel('planet_trajectories')
    mars=spice.Trajectory('MARS BARYCENTER')
    mars.generate_positions(mav.time[to],'Sun',frame)
    mars.change_units(astropy.units.AU)  
    mav.x[to]=mars.x
    mav.y[to]=mars.y
    mav.z[to]=mars.z
    [mav.r[to], mav.lat[to], mav.lon[to]]=cart2sphere(mav.x[to],mav.y[to],mav.z[to])

    #convert to degree
    mav.lon=np.rad2deg(mav.lon)   
    mav.lat=np.rad2deg(mav.lat)
    print('position end ')
  
        
    print('save MAVEN as pickle')

    header='MAVEN merged magnetic field and plasma data, obtained from Toulouse and C. Simon Wedlund. '+\
    'Timerange: '+mav.time[0].strftime("%Y-%b-%d %H:%M")+' to '+mav.time[-1].strftime("%Y-%b-%d %H:%M")+'.'+\
    'Mean time resolution: '+str(np.mean(np.diff(mav.time)).seconds)+' seconds. '+\
    'The data are put in a numpy recarray, fields can be accessed by mav.time, mav.bx, mav.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(mav.size)+'. '+\
    'Units are btxyz [nT, MSO], vtxyz [km/s, MSO], np[cm^-3], tp [K], orbital position: '+ \
    'xo/yo/zo/ro/lono/lato [km, degree, MSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    'Made with https://github.com/cmoestl/heliocats heliocats.data.convert_MAVEN_mat_original (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    

    pickle.dump([mav,header], open(data_path+filename, "wb"))

   
    



def convert_MAVEN_mat_removed(data_path,filename):

    print('load MAVEN from MAT')
    
    file=data_path+'input/MAVEN_2014to2018_removed_cyril.mat'
    mavraw = scipy.io.loadmat(file)
    
    #make array 
    mav=np.zeros(np.size(mavraw['BT']),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('tp', float),('np', float),('r', float),('lat', float),\
                ('lon', float),('x', float),('y', float),('z', float),\
                ('ro', float), ('lato', float), ('lono', float),\
                ('xo', float), ('yo', float), ('zo', float)])   
    #convert to recarray
    mav = mav.view(np.recarray)  
    #convert time from matlab to python
    t=mavraw['timeD'][:,0]
    for p in np.arange(np.size(t)):
        mav.time[p]= datetime.datetime.fromordinal(t[p].astype(int) ) + \
        datetime.timedelta(days=t[p]%1) - datetime.timedelta(days = 366) 
        

    mav.bx=mavraw['Bx'][:,0]       
    mav.by=mavraw['By'][:,0] 
    mav.bz=mavraw['Bz'][:,0]      
    mav.bt=mavraw['BT'][:,0]      

    mav.vx=mavraw['Vx'][:,0]      
    mav.vy=mavraw['Vy'][:,0]      
    mav.vz=mavraw['Vz'][:,0]      
    mav.vt=mavraw['VT'][:,0] 

    mav.tp=mavraw['Tp'][:,0]*(1.602176634*1e-19)/(1.38064852*1e-23)        #from ev to K     
    mav.np=mavraw['np'][:,0]      
  
    
    #add position with respect to Mars center in km in MSO
 
    print('orbit position start')
    
    insertion=datetime.datetime(2014,9,22,2,24,0)
    #these are the indices of the times for the cruise phase     
    tc=np.where(mdates.date2num(mav.time) < mdates.date2num(insertion))       

    mars_radius=3389.5
    mav.xo=mavraw['Xsc'][:,0]*mars_radius
    mav.yo=mavraw['Ysc'][:,0]*mars_radius
    mav.zo=mavraw['Zsc'][:,0]*mars_radius  
    
    #set to nan for cruise phase
    mav.xo[tc]=np.nan
    mav.yo[tc]=np.nan
    mav.zo[tc]=np.nan
    
    [mav.ro,mav.lato,mav.lono]=cart2sphere(mav.xo,mav.yo,mav.zo)
    mav.lono=np.rad2deg(mav.lono)   
    mav.lato=np.rad2deg(mav.lato)


    print('HEEQ position start')
    frame='HEEQ'
    
    #add position in HEEQ for cruise phase and orbit
    
    #cruise phase
    #use heliopy to load own bsp spice file from MAVEN 
    #obtained through https://naif.jpl.nasa.gov/pub/naif/pds/pds4/maven/maven_spice/spice_kernels/spk/
    spice.furnish(data_path+'input/maven_cru_rec_131118_140923_v1.bsp') 
    cruise=spice.Trajectory('MAVEN') #or NAIF CODE -202 
    cruise.generate_positions(mav.time[tc],'Sun',frame)     
    cruise.change_units(astropy.units.AU)  
    mav.x[tc]=cruise.x
    mav.y[tc]=cruise.y
    mav.z[tc]=cruise.z
    [mav.r[tc], mav.lat[tc], mav.lon[tc]]=cart2sphere(mav.x[tc],mav.y[tc],mav.z[tc])

    #times in orbit
    to=np.where(mdates.date2num(mav.time) > mdates.date2num(insertion))       

    planet_kernel=spicedata.get_kernel('planet_trajectories')
    mars=spice.Trajectory('MARS BARYCENTER')
    mars.generate_positions(mav.time[to],'Sun',frame)
    mars.change_units(astropy.units.AU)  
    mav.x[to]=mars.x
    mav.y[to]=mars.y
    mav.z[to]=mars.z
    [mav.r[to], mav.lat[to], mav.lon[to]]=cart2sphere(mav.x[to],mav.y[to],mav.z[to])

    #convert to degree
    mav.lon=np.rad2deg(mav.lon)   
    mav.lat=np.rad2deg(mav.lat)
    print('position end ')
    
    
    header='MAVEN merged magnetic field and plasma data, obtained from Toulouse. ' + \
    'The magnetosphere is removed with the Gruesbeck et al. 3D model (by C. Simon Wedlund). '+ \
    'Timerange: '+mav.time[0].strftime("%Y-%b-%d %H:%M")+' to '+mav.time[-1].strftime("%Y-%b-%d %H:%M")+'. '+\
    'Mean time resolution: '+str(np.mean(np.diff(mav.time)).seconds)+' seconds. '+\
    'The data are put in a numpy recarray, fields can be accessed by mav.time, mav.bx, mav.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(mav.size)+'. '+\
    'Units are btxyz [nT, MSO], vtxyz [km/s, MSO], np[cm^-3], tp [K], orbital position: '+ \
    'xo/yo/zo/ro/lono/lato [km, degree, MSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    'Made with https://github.com/cmoestl/heliocats heliocats.data.convert_MAVEN_mat_removed (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
     
    pickle.dump([mav,header], open(data_path+filename, "wb"))







def MAVEN_smooth_orbit(data_path,filename):


    filemav=data_path+'maven_2014_2018_removed.p'
    [mav,hmav]=pickle.load(open(filemav, 'rb' ) )
    print('loaded ',filemav)
    
    
    ############# smooth over each orbit to extract solar wind signal
    #determine apogees
    #get rid of nans
    mav.ro[np.where(mav.ro)==np.nan]=-1e3 
    peaks,properties = scipy.signal.find_peaks(mav.ro,height=5000,width=(100))
    #bring nans back
    mav.ro[np.where(mav.ro)==-1e3]=np.nan
    
    print('Nr. of orbits in dataset: ',len(peaks))


    #make array 
    mavs=np.zeros(np.size(peaks),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('tp', float),('np', float),('r', float),('lat', float),\
                ('lon', float),('x', float),('y', float),('z', float),\
                ('ro', float), ('lato', float), ('lono', float),\
                ('xo', float), ('yo', float), ('zo', float)])   
     
    #convert to recarray
    mavs = mavs.view(np.recarray)  

                
    #2h on each side
    window=121
    for i in np.arange(len(peaks)):
          mavs.bt[i]=np.nanmedian(mav.bt[peaks[i]-window:peaks[i]+window])
          mavs.bx[i]=np.nanmedian(mav.bx[peaks[i]-window:peaks[i]+window])
          mavs.by[i]=np.nanmedian(mav.by[peaks[i]-window:peaks[i]+window])
          mavs.bz[i]=np.nanmedian(mav.bz[peaks[i]-window:peaks[i]+window])

          mavs.vt[i]=np.nanmedian(mav.vt[peaks[i]-window:peaks[i]+window])
          mavs.vx[i]=np.nanmedian(mav.vx[peaks[i]-window:peaks[i]+window])
          mavs.vy[i]=np.nanmedian(mav.vy[peaks[i]-window:peaks[i]+window])
          mavs.vz[i]=np.nanmedian(mav.vz[peaks[i]-window:peaks[i]+window])

          mavs.np[i]=np.nanmedian(mav.np[peaks[i]-window:peaks[i]+window])
          mavs.tp[i]=np.nanmedian(mav.tp[peaks[i]-window:peaks[i]+window])
         
          mavs.time[i]=mav.time[peaks[i]] 
          
          mavs.r[i]=mav.r[peaks[i]]
          mavs.lat[i]=mav.lat[peaks[i]]
          mavs.lon[i]=mav.lon[peaks[i]]

          mavs.x[i]=mav.x[peaks[i]]
          mavs.y[i]=mav.y[peaks[i]]
          mavs.z[i]=mav.z[peaks[i]]

          mavs.ro[i]=mav.ro[peaks[i]]
          mavs.lato[i]=mav.lato[peaks[i]]
          mavs.lono[i]=mav.lono[peaks[i]]

          mavs.xo[i]=mav.xo[peaks[i]]
          mavs.yo[i]=mav.yo[peaks[i]]
          mavs.zo[i]=mav.zo[peaks[i]]

    
    '''
    for testing:
    plt.figure(1)
    ax1 = plt.subplot(121)
    #ax1.plot_date(mav.time,mav.bt,'bo') 

    #ax1.plot_date(mav.time,mav.bt,'-r') 
    ax1.plot_date(mav.time[peaks],bt1,'-r') 
    ax1.plot_date(mav.time[peaks],vt1,'-k') 

    #ax1.plot_date(mav.time[peaks],bt2,'bo') 


    #ax1.plot_date(mav.time,g,'-b') 

    #ax1.set_xlim(timeset-days_window*10,timeset+days_window*10)

    ax2 = plt.subplot(122)
    ax2.plot_date(mav.time,mav.ro,'-k') 
    ax2.plot_date(mav.time[peaks],mav.ro[peaks],'bo') 
    #ax2.set_xlim(timeset-days_window,timeset+days_window)

    ax2.set_ylim(7000,10000)
    plt.show()


    plt.figure(2)
    plt.plot_date(mav.time,mav.vt,'-r') 

    plt.plot_date(mav.time[peaks],vt1,'ob') 
    '''
    #pickle.dump(mavs, open(data_path+filename, "wb"))

    
    
    header='MAVEN solar wind merged magnetic field and plasma data. ' + \
    'The magnetosphere was removed with the Gruesbeck et al. 3D model (by C. Simon Wedlund), '+\
    'and a +/-2h median filter around the apogee is used for 1 data point per orbit. '+ \
    'Timerange: '+mavs.time[0].strftime("%Y-%b-%d %H:%M")+' to '+mavs.time[-1].strftime("%Y-%b-%d %H:%M")+'. '+\
    'Mean time resolution: '+str(np.mean(np.diff(mavs.time)).seconds)+' seconds. '+\
    'The data are put in a numpy recarray, fields can be accessed by mav.time, mav.bx, mav.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(mavs.size)+'. '+\
    'Units are btxyz [nT, MSO], vtxyz [km/s, MSO], np[cm^-3], tp [K], orbital position: '+ \
    'xo/yo/zo/ro/lono/lato [km, degree, MSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    'Made with https://github.com/cmoestl/heliocats heliocats.data.convert_MAVEN_mat_removed_orbit (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([mavs,header], open(data_path+filename, "wb"))

    


    
    
    
    



########################################## load HISTORIC DATA ############################


def save_helios_data(file):
    
    '''
    **TO DO
    '''

    print('start Helios')
    t_start = datetime.datetime(1975, 1, 1)
    t_end = datetime.datetime(1980, 12, 31)
     
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  


    #h1=heliosdata.corefit(1,t_start,t_end)
    #h2=heliosdata.corefit(2,t_start,t_end)
    h1=heliosdata.merged(1,t_start,t_end)  

    print('end Helios')




def save_cassini_data(file):

    '''
    **TO DO
    '''

    print('start Cassini')
    t_start = datetime.datetime(1999, 8, 16)
    t_end = datetime.datetime(2016, 12, 31)
     
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  


    coords='RTN'

    #Cassini Orbiter Magnetometer Calibrated MAG data in 1 minute averages available 
    #covering the period 1999-08-16 (DOY 228) to 2016-12-31 (DOY 366). 
    #The data are provided in RTN coordinates throughout the mission, with Earth, Jupiter, 
    #and Saturn centered coordinates for the respective flybys of those planets.
    cas=cassinidata.mag_hires(t_start,t_end, coords)




def save_ulysses_data(data_path):

   
    print('read Ulysses data from cdf and convert to pickle')   

    datacat_path='/nas/helio/data/DATACAT/'
    #load cdf
    ulycdf = cdflib.CDF(datacat_path+'ulysses_merged_1990_2009_CDAWEB.cdf') 
    #check variables
    #ulycdf.cdf_info()

    #time conversion to datetime      
    time=ulycdf.varget('Epoch')
    t=parse_time(time,format='cdf_epoch').datetime  
    
    
    #cutoff time and later data so that it starts with available position on Oct 6 1990
    t=t[6696:-1]

    
    print('Ulysses position start')
    #position starts on Oct 6 1990
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('ulysses'))
    upos=spice.Trajectory('-55')
    upos.generate_positions(t,'Sun',frame)
    upos.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(upos.x,upos.y,upos.z)
    print('position end ')
    
    #make custom array
    uly=np.zeros(len(t),dtype=[('time',object),('bx', float),('by', float), \
     ('bz', float), ('bt', float),('vt', float),('np', float),('tp', float), \
      ('x', float),('y', float), ('z', float),('r', float),('lat', float), \
      ('lon', float)])   
    #convert to recarray
    uly = uly.view(np.recarray)  
    
    uly.time=t
    
    uly.bx=ulycdf.varget('BR')[6696:-1]
    uly.by=ulycdf.varget('BT')[6696:-1]
    uly.bz=ulycdf.varget('BN')[6696:-1]
    uly.bt=ulycdf.varget('ABS_B')[6696:-1]
    uly.vt=ulycdf.varget('plasmaFlowSpeed')[6696:-1]
    uly.np=ulycdf.varget('protonDensity')[6696:-1]
    uly.tp=ulycdf.varget('protonTempLarge')[6696:-1]
    
    
    uly.x=upos.x
    uly.y=upos.y
    uly.z=upos.z
    
    uly.r=r
    uly.lat=np.rad2deg(lat)
    uly.lon=np.rad2deg(lon)
        
    badmag=np.where(uly.bt < -10000)
    uly.bt[badmag]=np.nan  
    uly.bx[badmag]=np.nan  
    uly.by[badmag]=np.nan  
    uly.bz[badmag]=np.nan  
    
    badv=np.where(uly.vt < -100000)
    uly.vt[badv]=np.nan  
    
    badn=np.where(uly.np < -100000)
    uly.np[badn]=np.nan  
    
    badt=np.where(uly.tp < -100000)
    uly.tp[badt]=np.nan  
    
 
    header='Ulysses merged magnetic field and plasma data, obtained from CDAWEB. '+ \
    'Timerange: '+uly.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+uly.time[-1].strftime("%d-%b-%Y %H:%M:%S") +\
    '. Units are btxyz [nT, RTN], vt [km/s], np [cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ].\
Made with https://github.com/cmoestl/heliocats and https://github.com/heliopython/heliopy. By C. Moestl (twitter @chrisoutofspace) and D. Stansby.'
            

    file=data_path+'ulysses_1990_2009_rtn.p'
    pickle.dump([uly,header], open(file, "wb"))
    
    
    print('Ulysses done')




      






############################# HELCATS DATA into single file ###############################



#for reading .sav files
def getcat(filename):

    print('reading CAT '+filename)
    cat=scipy.io.readsav(filename, verbose='true')  
    print('done reading CAT')
    return cat  
  
def convert_sav_to_p():
    start=time.time()
    wind=getcat('/nas/helio/data/DATACAT/WIND_2007to2016_HEEQ.sav')
    end=time.time()
    print('load wind end. time in minutes:', (end-start)/60)
    #save time and components as pickle
    pickle.dump(wind.wind, open( "/nas/helio/data/insitu_python/WIND_2007to2016_HEEQ.p", "wb" ) )


def despike_helcats_speed_wind(vt,vx,vy,vz):


    #set all nan to 0 in the v gradient array
    v1=copy.deepcopy( abs( np.gradient(vt)   ))
    vnan_ind=np.where(np.isnan(v1)==True)[0]
    v1[vnan_ind]=0
 
    peaks, properties = scipy.signal.find_peaks(v1, prominence=40, width=(1, 10))
    vt[np.where(v1>100)]=np.nan
    vx[np.where(v1>100)]=np.nan
    vy[np.where(v1>100)]=np.nan
    vz[np.where(v1>100)]=np.nan


    for i in np.arange(len(peaks)):
         width=int(np.ceil(properties['widths'])[i])
         #print(width)   
         vt[peaks[i]-width:peaks[i]+width]=np.nan
         vx[peaks[i]-width:peaks[i]+width]=np.nan
         vy[peaks[i]-width:peaks[i]+width]=np.nan
         vz[peaks[i]-width:peaks[i]+width]=np.nan


    #print(properties['widths'])
    #plt.plot(win.vt[0:80000],'-b')
    #plt.plot(win.vt[0:80000],'-k',linewidth=1)
    #plt.plot(v,'-r',linewidth=1)
    #plt.plot(v1,'-b',linewidth=1)
    return vt,vx,vy,vz
 
def despike_helcats_density_wind(den):


    #set all nan to 0 in the v gradient array
    den1=copy.deepcopy( abs( np.gradient(den)   ))
    den1nan_ind=np.where(np.isnan(den1)==True)[0]
    den1[den1nan_ind]=0
 
    peaks, properties = scipy.signal.find_peaks(den1, prominence=10, width=(1, 10))
    den[np.where(den1>10)]=np.nan
    den1[np.where(den1>10)]=np.nan

    for i in np.arange(len(peaks)):
         width=int(np.ceil(properties['widths'])[i])
         #print(width)   
         den[peaks[i]-width:peaks[i]+width]=np.nan
            
            
    #print(properties['widths'])
    #plt.plot(win.vt[0:80000],'-b')
    #plt.plot(win.np[1200000:1500000]+50,'-g',linewidth=5)
    #plt.plot(den+1,'-k',linewidth=1)
    #plt.plot(den1,'-b',linewidth=1)

    return den
 

def despike_helcats_temperature_wind(den):

    den=den/1e6    
    #set all nan to 0 in the v gradient array
    den1=copy.deepcopy( abs( np.gradient(den)   ))
    den1nan_ind=np.where(np.isnan(den1)==True)[0]
    den1[den1nan_ind]=0

    peaks, properties = scipy.signal.find_peaks(den1, prominence=0.2, width=(1, 10))
    #den[np.where(den>100)[0]]=np.nan
    den[np.where(den1>0.2)]=np.nan
    den1[np.where(den1>0.2)]=np.nan

    for i in np.arange(len(peaks)):
         width=int(np.ceil(properties['widths'])[i])
         #print(width)   
         den[peaks[i]-width:peaks[i]+width]=np.nan

    return den*1e6



    
    
    
def save_helcats_datacat(data_path,removed):  
    ''' to save all of helcats DATACAT into a single file
    use: hd.save_helcats_datacat(data_path,removed=True)
    '''
      
    print('save all helcats DATACAT into single file')
    datacat_path='/nas/helio/data/DATACAT/'
    print('all data in', datacat_path)
      
      
      
      
      
    print( 'read Wind')
    winin= pickle.load( open(datacat_path+ "WIND_2007to2018_HEEQ.p", "rb" ) )
    winin_time=parse_time(winin.time,format='utime').datetime
    winin=winin.astype([('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('vt', 'float64'), ('vx', 'float64'), ('vy', 'float64'), \
                    ('vz', 'float64'), ('tp', 'float64'), ('np', 'float64'), \
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64')])  
    #make new array with xyz                
    win=np.zeros(np.size(winin),[('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('vt', 'float64'), ('vx', 'float64'), ('vy', 'float64'), \
                    ('vz', 'float64'), ('tp', 'float64'), ('np', 'float64'),\
                    ('x', 'float64'),('y', 'float64'), ('z', 'float64'),\
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64')])      
                    
    win = win.view(np.recarray)                                                       
    
    win.time=winin_time
    win.bx=winin.bx
    win.by=winin.by
    win.bz=winin.bz
    win.bt=winin.bt

    win.vt=winin.vt
    win.vx=winin.vx
    win.vy=winin.vy
    win.vz=winin.vz

    win.np=winin.np
    win.tp=winin.tp
 
    [win.vt,win.vx,win.vy,win.vz]=despike_helcats_speed_wind(win.vt,win.vx,win.vy,win.vz)
    win.np=despike_helcats_density_wind(win.np)
    win.tp=despike_helcats_temperature_wind(win.tp)
    






    win.r=winin.r/(astropy.constants.au.value/1e3)
    win.lat=winin.lat
    win.lon=winin.lon
     
        
    [win.x, win.y, win.z]=sphere2cart(win.r,np.abs(win.lat-np.radians(90)),win.lon)
    win.lon=np.rad2deg(win.lon)   
    win.lat=np.rad2deg(win.lat)

    del(winin)
    
    hwin='Wind merged magnetic field and plasma data, obtained from HELCATS (A. Isavnin). '+ \
    'Timerange: '+win.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+win.time[-1].strftime("%d-%b-%Y %H:%M:%S") +\
    'Units are btxyz [nT, SCEQ], vtxyz [km/s, SCEQ], np [#/cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    
    pickle.dump([win,hwin], open(data_path+ "helcats/wind_2007_2018_helcats.p", "wb" ) )
    print( 'convert Wind done.')   
     
    
    
    
    

    print( 'read STEREO-A')
    stain= pickle.load( open(datacat_path+ "STA_2007to2015_SCEQ.p", "rb" ) )
    stain_time=parse_time(stain.time,format='utime').datetime
    stain=stain.astype([('time', 'object'), ('bt', 'float'),\
                    ('bx', 'float'), ('by', 'float'), ('bz', 'float'), \
                    ('vt', 'float'), ('vx', 'float'), ('vy', 'float'), \
                    ('vz', 'float'), ('tp', 'float'), ('np', 'float'), \
                    ('r', 'float'),('lat', 'float'), ('lon', 'float')])     
    sta=np.zeros(np.size(stain),[('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('vt', 'float64'), ('vx', 'float64'), ('vy', 'float64'), \
                    ('vz', 'float64'), ('tp', 'float64'), ('np', 'float64'),\
                    ('x', 'float64'),('y', 'float64'), ('z', 'float64'),\
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64')])      
    sta = sta.view(np.recarray)                                                       
    
    sta.time=stain_time
    sta.bx=stain.bx
    sta.by=stain.by
    sta.bz=stain.bz
    sta.bt=stain.bt

    sta.vt=stain.vt
    sta.vx=stain.vx
    sta.vy=stain.vy
    sta.vz=stain.vz

    sta.np=stain.np
    sta.tp=stain.tp
 
    sta.r=stain.r/(astropy.constants.au.value/1e3)
    sta.lat=stain.lat
    sta.lon=stain.lon
     
        

    [sta.x, sta.y, sta.z]=sphere2cart(sta.r,np.abs(sta.lat-np.radians(90)),sta.lon)
    sta.lon=np.rad2deg(sta.lon)   
    sta.lat=np.rad2deg(sta.lat)

    del(stain)
  
    hsta='STEREO-A merged magnetic field and plasma data, obtained from HELCATS (A. Isavnin). '+ \
    'Timerange: '+sta.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+sta.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
    'Units are btxyz [nT, SCEQ], vtxyz [km/s, SCEQ], np [#/cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    pickle.dump([sta,hsta], open(data_path+ "helcats/stereoa_2007_2015_helcats.p", "wb" ) )
    print( 'read STA done.')
    
  
  
    
    
    print( 'read STEREO-B')
    stbin= pickle.load( open(datacat_path+ "STB_2007to2014_SCEQ.p", "rb" ) )
    stbin_time=parse_time(stbin.time,format='utime').datetime
    stbin=stbin.astype([('time', 'object'), ('bt', 'float'),\
                    ('bx', 'float'), ('by', 'float'), ('bz', 'float'), \
                    ('vt', 'float'), ('vx', 'float'), ('vy', 'float'), \
                    ('vz', 'float'), ('tp', 'float'), ('np', 'float'), \
                    ('r', 'float'),('lat', 'float'), ('lon', 'float')])     
    stb=np.zeros(np.size(stbin),[('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('vt', 'float64'), ('vx', 'float64'), ('vy', 'float64'), \
                    ('vz', 'float64'), ('tp', 'float64'), ('np', 'float64'),\
                    ('x', 'float64'),('y', 'float64'), ('z', 'float64'),\
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64')])      
    stb = stb.view(np.recarray)                                                       
    
    stb.time=stbin_time
    stb.bx=stbin.bx
    stb.by=stbin.by
    stb.bz=stbin.bz
    stb.bt=stbin.bt

    stb.vt=stbin.vt
    stb.vx=stbin.vx
    stb.vy=stbin.vy
    stb.vz=stbin.vz

    stb.np=stbin.np
    stb.tp=stbin.tp
 
    stb.r=stbin.r/(astropy.constants.au.value/1e3)
    stb.lat=stbin.lat
    stb.lon=stbin.lon
    
    
    
    
     
    [stb.x, stb.y, stb.z]=sphere2cart(stb.r,np.abs(stb.lat-np.radians(90)),stb.lon)
    stb.lon=np.rad2deg(stb.lon)   
    stb.lat=np.rad2deg(stb.lat)
    
    
    #replace missing 2014 plasma data for STEREO-B    
    filestb2='stereob_2013_2014.p'
    [stb2,hstb2]=pickle.load(open(data_path+filestb2, "rb" ) )   
    
    stb_time_mat=parse_time(stb.time).plot_date
    stb2_time_mat=parse_time(stb2.time).plot_date
    
    #interpolate times onto stb.time
    dumvt=np.interp(stb_time_mat, stb2_time_mat,stb2.vt)
    dumtp=np.interp(stb_time_mat, stb2_time_mat,stb2.tp)
    dumnp=np.interp(stb_time_mat, stb2_time_mat,stb2.np)
    #get indices of 1-1-2014 to end
    begin=np.where(stb_time_mat > parse_time('2014-1-1').plot_date )[0][0]
    end=np.size(stb.vt)
    stb.vt[begin:end]=dumvt[begin:end]   
    stb.tp[begin:end]=dumtp[begin:end]   
    stb.np[begin:end]=dumnp[begin:end]   


    

    del(stbin)
    del(dumvt)
    del(dumtp)
    del(dumnp)



    hstb='STEREO-B merged magnetic field and plasma data, obtained from HELCATS (A. Isavnin). '+ \
    'Timerange: '+stb.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+stb.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
    'Units are btxyz [nT, SCEQ], vtxyz [km/s, SCEQ], np [#/cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    pickle.dump([stb,hstb], open(data_path+ "helcats/stereob_2007_2014_helcats.p", "wb" ) )
    print( 'read STB done.')
    
     
      
      
      

    print( 'read MESSENGER')
    #get insitu data from helcats, converted from IDL .sav to pickle
    if removed == True:   
       mesin= pickle.load( open( datacat_path+"MES_2007to2015_SCEQ_removed.p", "rb" ) )
    #non removed dataset
    if removed == False:  
       mesin= pickle.load( open( datacat_path+"MES_2007to2015_SCEQ_non_removed.p", "rb" ) )
    
    
    #time conversion
    mesin_time=parse_time(mesin.time,format='utime').datetime
    #replace mes.time with datetime object
    #new variable names                  
  
    mes=np.zeros(np.size(mesin),[('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('x', 'float64'),('y', 'float64'), ('z', 'float64'),\
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64')])
                    
    #convert to recarray
    mes = mes.view(np.recarray)  
    #set time new  
    mes.time=mesin_time
    mes.bx=mesin.bx
    mes.by=mesin.by
    mes.bz=mesin.bz
    mes.bt=mesin.btot
    
    
    #convert distance from Sun from km to AU, astropy constant is given in m
    mes.r=mesin.mes_radius_in_km_heeq/(astropy.constants.au.value/1e3)
    mes.lon=np.degrees(mesin.mes_longitude_in_radians_heeq.astype('float64'))   
    mes.lat=np.degrees(mesin.mes_latitude_in_radians_heeq.astype('float64'))
    [mes.x, mes.y, mes.z]=sphere2cart(mes.r,np.radians(np.abs(mes.lat-90)),np.radians(mes.lon))
    
    
    if removed == True:   
      hmes='MESSENGER magnetic field data, obtained from NASA PDS. '+ \
      'Timerange: '+mes.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+mes.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
      '. The magnetosphere is removed with a manual magnetopause crossings list (Lydia Philpott, Reka Winslow, Brian Anderson). '+ \
      'Units are btxyz [nT, SCEQ], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
      pickle.dump([mes,hmes], open(data_path+ "helcats/messenger_2007_2015_helcats_removed.p", "wb" ) )
     
    
    if removed == False:  
       hmes='MESSENGER magnetic field data, obtained from NASA PDS. '+ \
       'Timerange: '+mes.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+mes.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
       '. The magnetosphere is removed with a manual magnetopause crossings list (Lydia Philpott, Reka Winslow, Brian Anderson). '+ \
       'Units are btxyz [nT, SCEQ], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
       pickle.dump([mes,hmes], open(data_path+ "helcats/messenger_2007_2015_helcats.p", "wb" ) )
    
    print('convert MESSENGER done.')    
    
        


    print ('read VEX')
    if removed == True:   
         vexin= pickle.load( open(datacat_path+ "VEX_2007to2014_SCEQ_removed.p", "rb" ) )
    if removed == False:  
         vexin= pickle.load( open(datacat_path+ "VEX_2007to2014_SCEQ.p", "rb" ) )
    #time conversion
    vexin_time=parse_time(vexin.time,format='utime').datetime
    vex=np.zeros(np.size(vexin),[('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('x', 'float64'),('y', 'float64'), ('z', 'float64'),\
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64'),\
                    ('xo', 'float64'),('yo', 'float64'), ('zo', 'float64'),\
                    ('ro', 'float64'),('lato', 'float64'), ('lono', 'float64')])      
    vex = vex.view(np.recarray)                                  
 
 
    vex.r=vexin.vex_radius_in_km_heeq/(astropy.constants.au.value/1e3)
    vex.lon=np.rad2deg(vexin.vex_longitude_in_radians_heeq)   
    vex.lat=np.rad2deg(vexin.vex_latitude_in_radians_heeq)

    
    [vex.x, vex.y, vex.z]=sphere2cart(vex.r,np.radians(np.abs(vex.lat-90)),np.radians(vex.lon))
    #convert to degree

    vex.time=vexin_time
    vex.bx=vexin.bx
    vex.by=vexin.by
    vex.bz=vexin.bz
    vex.bt=vexin.btot

       
    #add orbit position 
    #https://www.cosmos.esa.int/web/spice/spice-for-vex
    #or from VEX_2007to2014_VSO.p
    vex2in=pickle.load( open( datacat_path+"VEX_2007to2014_VSO.p", "rb" ) )
  
    vex.xo=vex2in.x
    vex.yo=vex2in.y
    vex.zo=vex2in.z
    
    [vex.ro, vex.lato, vex.lono]=cart2sphere(vex.xo,vex.yo,vex.zo)
    vex.lono=np.rad2deg(vex.lono)   
    vex.lato=np.rad2deg(vex.lato)
    
    del(vexin)
    del(vex2in)

  
    
    if removed == True:   
     hvex='VEX magnetic field data, obtained from the VEX magnetometer PI T. Zhang IWF Graz, Austria. '+ \
     'Timerange: '+vex.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+vex.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
     '. The magnetosphere was removed with the model from Zhang et al. (2008), see Moestl et al. (2017, doi: 10.1002/2017SW001614) for details. '+ \
     'Units are btxyz [nT, SCEQ], orbital position: '+ \
     'xo/yo/zo/ro/lono/lato [km, degree, VSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
     pickle.dump([vex,hvex], open(data_path+ "helcats/vex_2007_2014_helcats_removed.p", "wb" ) )
     
    
    if removed == False:  
     hvex='VEX magnetic field data, obtained from the VEX magnetometer PI T. Zhang IWF Graz, Austria. '+ \
     'Timerange: '+vex.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+vex.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
     '. Units are btxyz [nT, SCEQ], orbital position: '+ \
     'xo/yo/zo/ro/lono/lato [km, degree, VSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
     pickle.dump([vex,hvex], open(data_path+ "helcats/vex_2007_2014_helcats.p", "wb" ) )
    
    print( 'convert VEX done.')
    
    
    
    

    
    #the Ulysses file has been generated by selecting the merged Ulysses data in CDAWEB
    #and then saved as one cdf 2.7 file
    print('read Ulysses from CDAWEB cdf')    
    save_ulysses_data(data_path)
    #fileuly=data_path+'ulysses_1990_2009_helcats.p'
    #[uly,huly]=pickle.load(open(fileuly, 'rb' ) )
  

    #if removed==True: 
    #    pickle.dump([vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly], open(data_path+ "helcats_all_data_removed.p", "wb" ) )
    #    print('saved as ' +data_path+ 'helcats_all_data_removed.p')

    #if removed==False: 
    #    pickle.dump([vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly], open(data_path+ "helcats_all_data_non_removed.p", "wb" ) )
    #    print('saved as ' +data_path+ 'helcats_all_data_non_removed.p')




  
def load_helcats_datacat(file):  
    ''' to load all of helcats DATACAT from a single file'''
    
    print('load all helcats DATACAT from single file: ', file)
    [vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly]=pickle.load( open(file, "rb" ) )
    print('Use vex,win,sta,stb,mes,uly to access data and position, hvex,hwin, hmes, hsta, hstb, huly for headers.')
    return [vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly]
    
    
    
    
def recarray_to_numpy_array(rec):    
    
    '''convert data recarray to numpy structured array with matplotlib time '''

    #recarr.time=parse_time(recarr.time).plot_date
    #numarr = pd.DataFrame(recarr).to_numpy()
    
    num=copy.deepcopy(np.array(rec))
    num['time']=parse_time(num['time']).plot_date
    
    return num
    
    
  
    
    
    
    
    
#################################### MATH ################################################



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
    

   
def convert_GSE_to_HEEQ(sc):
    '''
    for Wind magnetic field components: convert GSE to HEE to HAE to HEEQ
    '''

    print('conversion GSE to HEEQ start')                                

    jd=np.zeros(len(sc))
    mjd=np.zeros(len(sc))
        

    for i in np.arange(0,len(sc)):

        jd[i]=parse_time(sc.time[i]).jd
        mjd[i]=float(int(jd[i]-2400000.5)) #use modified julian date    
        
        #GSE to HEE
        #Hapgood 1992 rotation by 180 degrees, or simply change sign in bx by    
        #rotangle=np.radians(180)
        #c, s = np.cos(rotangle), np.sin(rotangle)
        #T1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
        #[bx_hee,by_hee,bz_hee]=T1[sc.bx[i],sc.by[i],sc.bz[i]]        
        b_hee=[-sc.bx[i],-sc.by[i],sc.bz[i]]
        
        #HEE to HAE        
        
        #define T00 and UT
        T00=(mjd[i]-51544.5)/36525.0          
        dobj=sc.time[i]
        UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   

        #lambda_sun in Hapgood, equation 5, here in rad
        M=np.radians(357.528+35999.050*T00+0.04107*UT)
        LAMBDA=280.460+36000.772*T00+0.04107*UT        
        lambda_sun=np.radians( (LAMBDA+(1.915-0.0048*T00)*np.sin(M)+0.020*np.sin(2*M)) )
        
        #S-1 Matrix equation 12 hapgood 1992, change sign in lambda angle
        c, s = np.cos(-(lambda_sun+np.radians(180))), np.sin(-(lambda_sun+np.radians(180)))
        Sm1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
        b_hae=np.dot(Sm1,b_hee)

        #HAE to HEEQ
        
        iota=np.radians(7.25)
        omega=np.radians((73.6667+0.013958*((mjd[i]+3242)/365.25)))                      
        theta=np.arctan(np.cos(iota)*np.tan(lambda_sun-omega))                       
                      
        #quadrant of theta must be opposite lambda_sun minus omega; Hapgood 1992 end of section 5   
        #get lambda-omega angle in degree mod 360   
        
        #************** CHECK
                      
        lambda_omega_deg=np.degrees(np.mod(np.degrees(lambda_sun)-np.degrees(omega),360))
        ##get theta_node in deg
        theta_node_deg=np.degrees(theta)
        ##if in same quadrant, then theta_node = theta_node +pi   
        if abs(lambda_omega_deg-theta_node_deg) < 180: theta=theta+np.pi                                                            
        
        #rotation around Z by theta
        c, s = np.cos(theta), np.sin(theta)
        S2_1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))

        #rotation around X by iota  
        iota=np.radians(7.25)
        c, s = np.cos(iota), np.sin(iota)
        S2_2 = np.array(( (1,0,0), (0,c, s), (0, -s, c)) )
                
        #rotation around Z by Omega  
        c, s = np.cos(omega), np.sin(omega)
        S2_3 = np.array( ((c,s, 0), (-s, c, 0), (0, 0, 1)) )
        
        #matrix multiplication to go from HAE to HEEQ components                
        [bx_heeq,by_heeq,bz_heeq]=np.dot(  np.dot(   np.dot(S2_1,S2_2),S2_3), b_hae) 
        
        sc.bx[i]=bx_heeq
        sc.by[i]=by_heeq
        sc.bz[i]=bz_heeq

    
    print('conversion GSE to HEEQ done')                                
    return sc

    
   

def convert_RTN_to_SCEQ(sc,name):
    '''
    for STEREO-A, B and Parker Solar Probe   
    sc is the input recarray for the data, e.g. psp, sta    
    '''
    
    sc_len=len(sc)
   
    #HEEQ - make to SCEQ for each timestep,   solar rotation axis at 0, 0, 1 in HEEQ
    X_heeq=[1,0,0]
    Y_heeq=[0,1,0]
    Z_heeq=[0,0,1]
        
    # go through all data points    
    for i in np.arange(0,sc_len):        
        
        #rotation for HEEQ to SCEQ vectors, depending on current longitude (correct c, s)
        rotangle=np.radians(sc.lon[i])
        c, s = np.cos(rotangle), np.sin(rotangle)
        #rotation matrix around Z
        R = np.array(((c,-s, 0), (s, c, 0), (0, 0, 1)))
        
        #rotate X and Y  
        X_sceq=np.dot(R,X_heeq)
        Y_sceq=np.dot(R,Y_heeq)
         
        #make normalized RTN vectors
        Xrtn=[sc.x[i], sc.y[i],sc.z[i]]/np.linalg.norm([sc.x[i], sc.y[i],sc.z[i]])
        Yrtn=np.cross(Z_heeq,Xrtn)/np.linalg.norm(np.cross(Z_heeq,Xrtn))
        Zrtn=np.cross(Xrtn, Yrtn)/np.linalg.norm(np.cross(Xrtn, Yrtn))
        
        #project into new system
        sc.bx[i]=np.dot(np.dot(sc.bx[i],Xrtn)+np.dot(sc.by[i],Yrtn)+np.dot(sc.bz[i],Zrtn),X_sceq)
        sc.by[i]=np.dot(np.dot(sc.bx[i],Xrtn)+np.dot(sc.by[i],Yrtn)+np.dot(sc.bz[i],Zrtn),Y_sceq)
        sc.bz[i]=np.dot(np.dot(sc.bx[i],Xrtn)+np.dot(sc.by[i],Yrtn)+np.dot(sc.bz[i],Zrtn),Z_heeq)
        
        #project into new system
        if name=='PSP':
            sc.vx[i]=np.dot(np.dot(sc.vx[i],Xrtn)+np.dot(sc.vy[i],Yrtn)+np.dot(sc.vz[i],Zrtn),X_sceq)
            sc.vy[i]=np.dot(np.dot(sc.vx[i],Xrtn)+np.dot(sc.vy[i],Yrtn)+np.dot(sc.vz[i],Zrtn),Y_sceq)
            sc.vz[i]=np.dot(np.dot(sc.vx[i],Xrtn)+np.dot(sc.vy[i],Yrtn)+np.dot(sc.vz[i],Zrtn),Z_heeq)


    
    return sc



    

##########################################################################################################
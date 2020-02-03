#data.py
#load data for heliocats
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

data_path='/nas/helio/data/insitu_python/'

'''
MIT LICENSE
Copyright 2020, Christian Moestl 
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
        d1=get_noaa_realtime_data(m1, p1)
    
        #save in large array
        noaa[k:k+np.size(d1)]=d1
        k=k+np.size(d1) 


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
    dp = json.loads (plasmafile.read())
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





 
def save_stereoa_beacon_data(path,file,start_time,end_time):

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

    
    header='BEACON STEREO-A magnetic field (IMPACT instrument) and plasma data (PLASTIC), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/  '+ \
    'Timerange: '+sta.time[0].strftime("%Y-%b-%d %H:%M")+' to '+sta.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(sta.time)).seconds)+' seconds. '+\
    'A median filter has been applied as scipy.signal.medfilt(sta.vt,9) and then spikes were removed with '+\
    'scipy.signal.find_peaks(sta.vt, prominence=200,width=(1,200)). '+\
    'The data are available in a numpy recarray, fields can be accessed by sta.time, sta.bx, sta.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(sta.size)+'. '+\
    'Units are btxyz [nT, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_stereoa_beacon_data (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    #'Units are btxyz [nT, RTN], vtxyz [km/s, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\

    pickle.dump([sta,header], open(path+file, "wb"))
    
    print('done sta')
    print()







def save_wind_data(path,file,start_date,end_date):
    
    '''
    description of data sources used in heliosat:
    https://cdaweb.sci.gsfc.nasa.gov/misc/NotesW.html
    '''
    
    print('start wind update')
    wind_sat = heliosat.WIND()
    t_start = start_date
    t_end = end_date
    
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    tm, mag = wind_sat.get_data_raw(t_start, t_end, "wind_mfi_k0", extra_columns=["PGSE"], return_datetimes=True)
    #tm, mag = wind_sat.get_data_raw(t_start, t_end, "wind_mfi_h0")
    tp, pro = wind_sat.get_data_raw(t_start, t_end, "wind_swe_h1",return_datetimes=True)

    print('download complete')
    
    #tm=parse_time(tm,format='unix').datetime 
    #tp=parse_time(tp,format='unix').datetime 
    
    #convert to matplotlib time for linear interpolation
    tm_mat=mdates.date2num(tm) 
    tp_mat=mdates.date2num(tp) 
    print('time convert done')
    
     
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
        
    den = np.interp(time_mat, tp_mat, pro[:,0])
    vt = np.interp(time_mat, tp_mat, pro[:,1])
    tp = np.interp(time_mat, tp_mat, pro[:,2])
    
    
    #interpolate the GSE position over full data range
    x_gse = np.interp(time_mat, tm_mat, mag[:,3])*6378.1/149597870.7*astropy.units.AU #earth radii to km to AU
    y_gse = np.interp(time_mat, tm_mat, mag[:,4])*6378.1/149597870.7*astropy.units.AU
    z_gse = np.interp(time_mat, tm_mat, mag[:,5])*6378.1/149597870.7*astropy.units.AU
    
    
    
    #set nan over linear interpolated
    #..........    
        
    print('position start')
    
    frame='HEEQ'
    planet_kernel=spicedata.get_kernel('planet_trajectories')
    earth=spice.Trajectory('399')  #399 for Earth, not barycenter (because of moon)
    earth.generate_positions(time,'Sun',frame)

    #from km to AU
    earth.change_units(astropy.units.AU)
    
    #add gse position to Earth position
    x=earth.x-x_gse  #earth radii to km
    y=earth.y-y_gse
    z=earth.z+z_gse
    [r, lat, lon]=cart2sphere(x,y,z)
    
    #*****with astropy lagrange points exact value? L1 position with 0.01 AU 
    #[r, lat, lon]=cart2sphere(earth.x-0.01*astropy.units.AU,earth.y,earth.z)
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
    win.lat=np.rad2deg(lat)
    win.lon=np.rad2deg(lon)

    
    win.np=den
    win.vt=vt
    #https://en.wikipedia.org/wiki/Thermal_velocity convert from km/s to K
    from astropy.constants import m_p,k_B
    win.tp=np.pi*m_p*((tp*1e3)**2)/(8*k_B) 
        
    #win.p3=p3
    #win.p4=p4
    
    
    
    ###################### get rid of spikes
    win.vt[np.where(win.vt> 3000)]=100000
    #get rid of all single spikes with scipy signal find peaks
    peaks, properties = scipy.signal.find_peaks(win.vt, height=3000,width=(1, 250))
    #go through all of them and set to nan according to widths
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        win.vt[peaks[i]-width-2:peaks[i]+width+2]=np.nan
        
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
        
    #magnetic field    
    peaks, properties = scipy.signal.find_peaks(win.bt, height=50,width=(1, 10))
    #go through all of them and set to nan according to widths
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        win.bt[peaks[i]-width-2:peaks[i]+width+2]=np.nan    
        win.bx[peaks[i]-width-2:peaks[i]+width+2]=np.nan    
        win.by[peaks[i]-width-2:peaks[i]+width+2]=np.nan    
        win.bz[peaks[i]-width-2:peaks[i]+width+2]=np.nan    
        
        
        
    #manual spike removal for magnetic field
    remove_start=datetime.datetime(2018, 7, 19, 18, 25)
    remove_end=datetime.datetime(2018, 7, 19, 19, 30)
    remove_start_ind=np.where(remove_start==win.time)[0][0]
    remove_end_ind=np.where(remove_end==win.time)[0][0] 

    win.bt[remove_start_ind:remove_end_ind]=np.nan
    win.bx[remove_start_ind:remove_end_ind]=np.nan
    win.by[remove_start_ind:remove_end_ind]=np.nan
    win.bz[remove_start_ind:remove_end_ind]=np.nan
    

        
    ###################### 

    
    header='Wind magnetic field (MAG instrument) and plasma data (SWE), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/wind/  '+ \
    'Timerange: '+win.time[0].strftime("%Y-%b-%d %H:%M")+' to '+win.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(win.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by win.time, win.bx, win.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(win.size)+'. '+\
    'Units are btxyz [nT, GSE], vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_wind_data (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    pickle.dump([win,header], open(path+file, "wb"))
    
    
    print('wind update done')
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

    
   
  


def save_stereoa_science_data(path,file):

    print('start STA')
    sta_sat = heliosat.STA()
    t_start = datetime.datetime(2018, 1, 1)
    t_end = datetime.datetime(2018, 1,30)
     
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    tp, pro = sta_sat.get_data_raw(t_start, t_end, "sta_plastic_l2")
    tm, mag = sta_sat.get_data_raw(t_start, t_end, "sta_impact_beacon")

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
      
      
      
    #add speed!!!!!!!!!!!!!!!!  check parameters
    #STA_L1_MAG_RTN_20180101_V06.cdf
    #STA_L2_PLA_1DMax_1min_20181011_V11.cdf  
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
    
    #'Units are btxyz [nT, RTN], vtxyz [km/s, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\

    pickle.dump([sta,header], open(path+file, "wb"))
    
    print('done sta')
    print()



    
def save_psp_data_non_merged(path, file):
    '''
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




def save_psp_data(path, file):
     
    print('start PSP')
     
    psp_sat = heliosat.PSP()
    #!!!if you change these times, see below for another time set!
    t_start = datetime.datetime(2018, 10, 14,14,14, 30)
    t_end = datetime.datetime(2019, 4, 23,23,59,30)
    #t_end = datetime.datetime(2018, 12, 12,23,59,30)
    #t_end = datetime.datetime(2018, 11, 23,23,59,30)
    #t_end = datetime.datetime(2019, 5, 31,23,59,30)    
    #!!!if you change these times, see below for another time set!
    
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    tm, mag = psp_sat.get_data_raw(t_start, t_end, "psp_fields_l2")
    tp, pro = psp_sat.get_data_raw(t_start, t_end, "psp_spc_l3")

    print('download complete')

    tm=parse_time(tm,format='unix').datetime 
    tp=parse_time(tp,format='unix').datetime 
    
    #convert to matplotlib time for linear interpolation
    tm_mat=mdates.date2num(tm) 
    tp_mat=mdates.date2num(tp) 
    
    print('time convert done')

      
    print('position start')
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('psp_pred'))
    psptra=spice.Trajectory('SPP')
    psptra.generate_positions(time,'Sun',frame)
    psptra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(psptra.x,psptra.y,psptra.z)
    print('PSP pos')    
    print('position end')
   
    

    #linear interpolation to time_mat times 
    
    #which values are not in original data?
    isin=np.isin(time_mat,tm_mat)      
    setnan=np.where(isin==False)
  
    #linear interpolation to time_mat times now with new array at full minutes
    t_start = datetime.datetime(2018, 10, 14,14,15, 0)
    #t_end = datetime.datetime(2018, 11, 23,23,59,0)
    t_end = datetime.datetime(2019, 4, 23,23,59,0)
    #t_end = datetime.datetime(2019, 5, 31,23,59,0)
    
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
   
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )

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
    tround_mat=mdates.date2num(tround) 

    #same as above for magnetic field
    isin=np.isin(time_mat,tround_mat)      
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
    psp.lat=np.rad2deg(lat)
    psp.lon=np.rad2deg(lon)
    
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
    
    
    header='PSP magnetic field (FIELDS instrument) and plasma data (SWEAP), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/psp/  '+ \
    'Timerange: '+psp.time[0].strftime("%Y-%b-%d %H:%M")+' to '+psp.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(psp.time)).seconds)+' seconds. '+\
    'The data are put in a numpy recarray, fields can be accessed by psp.time, psp.bx, psp.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(psp.size)+'. '+\
    'Units are btxyz [nT, RTN], vtxyz [km/s, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([psp,header], open(path+file, "wb"))

    print('done psp')
    print()

    ########## for debugging
    '''
    plt.figure(1)
    plt.plot_date(tm_mat,mag[:,0],'r-')
    plt.plot_date(time_mat,bx,'g-')
    plt.figure(2)
    plt.plot_date(tp_mat,pro[:,2],'r-')
    plt.plot_date(time_mat,vx,'g-')
    plt.figure(3)
    plt.plot_date(time_mat,bx,'r-')
    plt.plot_date(time_mat,vx,'g-')
    plt.show()
    sys.exit()
    '''










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



    print('start Cassini')
    t_start = datetime.datetime(1999, 8, 16)
    t_end = datetime.datetime(2016, 12, 31)
     
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  


    coords='RTN'

    #Cassini Orbiter Magnetometer Calibrated MAG data in 1 minute averages available covering the period 1999-08-16 (DOY 228) to 2016-12-31 (DOY 366). The data are provided in RTN coordinates throughout the mission, with Earth, Jupiter, and Saturn centered coordinates for the respective flybys of those planets.
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
    uly.vp=ulycdf.varget('plasmaFlowSpeed')[6696:-1]
    uly.np=ulycdf.varget('protonDensity')[6696:-1]
    uly.tp=ulycdf.varget('protonTempLarge')[6696:-1]
    
    
    uly.x=upos.x
    uly.y=upos.y
    uly.z=upos.z
    
    uly.r=r
    uly.lat=lat
    uly.lon=lon
    
        
    badmag=np.where(uly.bt < -10000)
    uly.bt[badmag]=np.nan  
    uly.bx[badmag]=np.nan  
    uly.by[badmag]=np.nan  
    uly.bz[badmag]=np.nan  
    
    badv=np.where(uly.vp < -100000)
    uly.vp[badv]=np.nan  
    
    badn=np.where(uly.np < -100000)
    uly.np[badn]=np.nan  
    
    badt=np.where(uly.tp < -100000)
    uly.tp[badt]=np.nan  
    
 
    header='Ulysses merged magnetic field and plasma data, obtained from CDAWEB. '+ \
    'Timerange: '+uly.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+uly.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
    'Units are btxyz [nT, RTN], vt [km/s], np [#/cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
            

    file=data_path+'ulysses_1990_2009_helcats.p'
    pickle.dump([uly,header], open(file, "wb"))
    
    
    print('Ulysses done')




















      






############################# HELCATS DATA into single file ###############################





  
  
  
def save_helcats_datacat(data_path,removed):  
    ''' to save all of helcats DATACAT into a single file'''
      
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
 
    win.r=winin.r/(astropy.constants.au.value/1e3)
    win.lat=winin.lat
    win.lon=winin.lon
     
    [win.x, win.y, win.z]=sphere2cart(win.r,win.lat,win.lon)
    win.lon=np.rad2deg(win.lon)   
    win.lat=np.rad2deg(win.lat)

    #**remove spikes in v
    #https://datascience.stackexchange.com/questions/27031/how-to-get-spike-values-from-a-value-sequence
    del(winin)
    
    hwin='Wind merged magnetic field and plasma data, obtained from HELCATS (A. Isavnin). '+ \
    'Timerange: '+win.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+win.time[-1].strftime("%d-%b-%Y %H:%M:%S") +\
    'Units are btxyz [nT, SCEQ], vtxyz [km/s, SCEQ], np [#/cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    
    pickle.dump([win,hwin], open(data_path+ "wind_2007_2018_helcats.p", "wb" ) )
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
     
    [sta.x, sta.y, sta.z]=sphere2cart(sta.r,sta.lat,sta.lon)
    sta.lon=np.rad2deg(sta.lon)   
    sta.lat=np.rad2deg(sta.lat)

    #**remove spikes in v
    #https://datascience.stackexchange.com/questions/27031/how-to-get-spike-values-from-a-value-sequence
    del(stain)
  
    hsta='STEREO-A merged magnetic field and plasma data, obtained from HELCATS (A. Isavnin). '+ \
    'Timerange: '+sta.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+sta.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
    'Units are btxyz [nT, SCEQ], vtxyz [km/s, SCEQ], np [#/cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    pickle.dump([sta,hsta], open(data_path+ "stereoa_2007_2015_helcats.p", "wb" ) )
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
     
    [stb.x, stb.y, stb.z]=sphere2cart(stb.r,stb.lat,stb.lon)
    stb.lon=np.rad2deg(stb.lon)   
    stb.lat=np.rad2deg(stb.lat)

    #**remove spikes in v
    #https://datascience.stbckexchange.com/questions/27031/how-to-get-spike-values-from-a-value-sequence
    del(stbin)
    hstb='STEREO-B merged magnetic field and plasma data, obtained from HELCATS (A. Isavnin). '+ \
    'Timerange: '+stb.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+stb.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
    'Units are btxyz [nT, SCEQ], vtxyz [km/s, SCEQ], np [#/cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    pickle.dump([stb,hstb], open(data_path+ "stereob_2007_2014_helcats.p", "wb" ) )
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
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64'),\
                    ('xo', 'float64'),('yo', 'float64'), ('zo', 'float64'),\
                    ('ro', 'float64'),('lato', 'float64'), ('lono', 'float64')])
                    
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
    [mes.x, mes.y, mes.z]=sphere2cart(mes.r,\
                                      mesin.mes_latitude_in_radians_heeq.astype('float64'),\
                                      mesin.mes_longitude_in_radians_heeq.astype('float64'))
    #convert to degree
    mes.lon=np.rad2deg(mesin.mes_longitude_in_radians_heeq.astype('float64'))   
    mes.lat=np.rad2deg(mesin.mes_latitude_in_radians_heeq.astype('float64'))

    
    #add orbit position after orbit insertion March 18, 2011, 01:00
    #https://naif.jpl.nasa.gov/pub/naif/pds/data/mess-e_v_h-spice-6-v1.0/messsp_1000/aareadme.htm
    #or from MES_2007to2015_RTN.sav       
    mes2in=pickle.load( open( datacat_path+"MES_2007to2015_RTN.p", "rb" ) )
    orbit_insertion=mdates.date2num(datetime.datetime(2011,3,18,1,0,0))
    before_orbit=np.where(mdates.date2num(mes.time) < orbit_insertion)
    
    mes2in.x.astype('float64')[before_orbit]=np.nan
    mes2in.y.astype('float64')[before_orbit]=np.nan
    mes2in.z.astype('float64')[before_orbit]=np.nan #for some reason this was saved as int
    
    mes.xo=mes2in.x
    mes.yo=mes2in.y
    mes.zo=mes2in.z
    
    [mes.ro, mes.lato, mes.lono]=cart2sphere(mes.xo,mes.yo,mes.zo)
    mes.lono=np.rad2deg(mes.lono)   
    mes.lato=np.rad2deg(mes.lato)
    del(mesin)
    del(mes2in)
    
    
    
    if removed == True:   
      hmes='MESSENGER magnetic field data, obtained from NASA PDS. '+ \
      'Timerange: '+mes.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+mes.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
      '. The magnetosphere is removed with a manual magnetopause crossings list (Lydia Philpott, Reka Winslow, Brian Anderson). '+ \
      'Units are btxyz [nT, SCEQ], orbital position: '+ \
      'xo/yo/zo/ro/lono/lato [km, degree, MSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
      pickle.dump([mes,hmes], open(data_path+ "messenger_2007_2015_helcats_removed.p", "wb" ) )
     
    
    if removed == False:  
       hmes='MESSENGER magnetic field data, obtained from NASA PDS. '+ \
       'Timerange: '+mes.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+mes.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
       '. The magnetosphere is removed with a manual magnetopause crossings list (Lydia Philpott, Reka Winslow, Brian Anderson). '+ \
       'Units are btxyz [nT, SCEQ], orbital position: '+ \
       'xo/yo/zo/ro/lono/lato [km, degree, MSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
       pickle.dump([mes,hmes], open(data_path+ "messenger_2007_2015_helcats.p", "wb" ) )
    
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
    [vex.x, vex.y, vex.z]=sphere2cart(vex.r,\
                                      vexin.vex_latitude_in_radians_heeq.astype('float64'),\
                                      vexin.vex_longitude_in_radians_heeq.astype('float64'))
    #convert to degree
    vex.lon=np.rad2deg(vexin.vex_longitude_in_radians_heeq)   
    vex.lat=np.rad2deg(vexin.vex_latitude_in_radians_heeq)

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
     '. The magnetosphere was removed with the ???? model. '+ \
     'Units are btxyz [nT, SCEQ], orbital position: '+ \
     'xo/yo/zo/ro/lono/lato [km, degree, VSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
     pickle.dump([vex,hvex], open(data_path+ "vex_2007_2014_helcats_removed.p", "wb" ) )
     
    
    if removed == False:  
     hvex='VEX magnetic field data, obtained from the VEX magnetometer PI T. Zhang IWF Graz, Austria. '+ \
     'Timerange: '+vex.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+vex.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
     'Units are btxyz [nT, SCEQ], orbital position: '+ \
     'xo/yo/zo/ro/lono/lato [km, degree, VSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
     pickle.dump([vex,hvex], open(data_path+ "vex_2007_2014_helcats.p", "wb" ) )
    
    print( 'convert VEX done.')
    
    
    
    

    
    #the Ulysses file has been generated by selecting the merged Ulysses data in CDAWEB
    #and then saved as one cdf 2.7 file
    print('read Ulysses from CDAWEB cdf')    
    save_ulysses_data(data_path)
    fileuly=data_path+'ulysses_1990_2009_helcats.p'
    [uly,huly]=pickle.load(open(fileuly, 'rb' ) )
  

    if removed==True: 
        pickle.dump([vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly], open(data_path+ "helcats_all_data_removed.p", "wb" ) )
        print('saved as ' +data_path+ 'helcats_all_non_removed.p')

    if removed==False: 
        pickle.dump([vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly], open(data_path+ "helcats_all_data_non_removed.p", "wb" ) )
        print('saved as ' +data_path+ 'helcats_all_data_non_removed.p')




  
def load_helcats_datacat(file):  
    ''' to load all of helcats DATACAT from a single file'''
    
    print('load all helcats DATACAT from single file: ', file)
    [vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly]=pickle.load( open(file, "rb" ) )
    print('Use vex,win,sta,stb,mes,uly to access data and position, hvex,hwin, hmes, hsta, hstb, huly for headers.')
    return [vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly]
    
    
    
    
    
    
    
    
    
    
    
    
    
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
    


























##################### OLD CODE






     


'''
#cme
t1=t1[40000:48000]
m1=m1[40000:48000,:]
#t_swe1=t_swe1[40000:48000]
#swe1=swe1[40000:48000,:]
r1=r1[40000:48000]
lon1=lon1[40000:48000]


bx1=m1[:,0]  
by1=m1[:,1]  
bz1=m1[:,2]  
bt1=np.sqrt(bx1**2+by1**2+bz1**2)

bx2=m2[:,0]  
by2=m2[:,1]  
bz2=m2[:,2]  
bt2=np.sqrt(bx2**2+by2**2+bz2**2)


v1=swe1[:,2] 
v2=swe2[:,2] 


plt.close('all')



def save_psp_data2(file):
     
    t_start = datetime.datetime(2018, 10, 14)
    t_end = datetime.datetime(2018, 12, 20)
    
    #fields
    psp_t1, psp_m1 = psp_sat.get_data_raw(t_start, t_end, "mag")
    #t1p, p1 = psp_sat.get_data_raw(t_start, t_end, "proton")
    psp_t1=parse_time(psp_t1,format='unix').datetime  
    
    #sweap
    t_swe1, swe1 = psp_sat.get_data_raw(t_start, t_end, "spc_l3i")
    t_swe1=parse_time(t_swe1,format='unix').datetime  
    

    #fields
    t_start2 = datetime.datetime(2019, 3, 1)
    t_end2 = datetime.datetime(2019, 5, 30)
    psp_t2, psp_m2 = psp_sat.get_data_raw(t_start2, t_end2, "mag")
    psp_t2=parse_time(psp_t2,format='unix').datetime  
    
    #sweap
    t_swe2, swe2 = psp_sat.get_data_raw(t_start2, t_end2, "spc_l3i")
    t_swe2=parse_time(t_swe2,format='unix').datetime  





    pos1=psp_sat.trajectory(psp_t1, frame="HEEQ")
    pos2=psp_sat.trajectory(psp_t2, frame="HEEQ")

    starttime =datetime.datetime(2018, 8,13)
    endtime = datetime.datetime(2019, 8, 31)
    psp_time = []
    while starttime < endtime:
        psp_time.append(starttime)
        starttime += datetime.timedelta(days=1/24.)
    psp_time_num=mdates.date2num(psp_time)     

    spice.furnish(spicedata.get_kernel('psp_pred'))
    psp=spice.Trajectory('SPP')
    psp.generate_positions(psp_time,'Sun','HEEQ')
    print('PSP pos')

    psp.change_units(astropy.units.AU)  
    
    [pos1_r, pos1_lat, pos1_lon]=cart2sphere(pos1[:,0],pos1[:,1],pos1[:,2])
    [pos2_r, pos2_lat, pos2_lon]=cart2sphere(pos2[:,0],pos2[:,1],pos2[:,2])


    pickle.dump([pos1_r, pos1_lat, pos1_lon, pos2_r, pos2_lat, pos2_lon,  psp_t1,psp_m1,psp_t2,psp_m2,t_swe1,swe1,t_swe2,swe2], open(file, "wb"))






sns.set_style('darkgrid')
fig1=plt.figure(figsize=[25, 12],dpi=100)
 

plt.suptitle('PSP first 2 orbits')

ax1=plt.subplot(421)
plt.plot_date(t1,bx1,'-r',label='BR',linewidth=0.5)
plt.plot_date(t1,by1,'-g',label='BT',linewidth=0.5)
plt.plot_date(t1,bz1,'-b',label='BN',linewidth=0.5)
plt.plot_date(t1,bt1,'-k',label='Btotal')
ax1.set_ylabel('B [nT]')
plt.legend(loc=2)
ax1.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d %H') )
ax1.set_xlim(t1[0],t1[-1])

ax2=plt.subplot(422)
plt.plot_date(t2,bx2,'-r',label='BR', linewidth=0.5)
plt.plot_date(t2,by2,'-g',label='BT', linewidth=0.5)
plt.plot_date(t2,bz2,'-b',label='BN', linewidth=0.5)
plt.plot_date(t2,bt2,'-k',label='Btotal')
ax2.set_ylabel('B [nT]')
plt.legend(loc=2)
ax2.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d %H') )
ax2.set_xlim(t2[0],t2[-1])



ax3=plt.subplot(423, sharex=ax1)
plt.plot_date(t_swe1,v1,'-k',label='V', linewidth=0.8)
ax3.set_ylabel('V [km/s]')
ax3.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d %H') )
ax3.set_xlim(t1[0],t1[-1])



ax4=plt.subplot(424,sharex=ax2)
plt.plot_date(t_swe2,v2,'-k',label='V', linewidth=0.8)
ax4.set_xlim(t2[0],t2[-1])
ax4.set_ylabel('V [km/s]')
ax4.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d %H') )


ax5=plt.subplot(425,sharex=ax1)
plt.plot_date(t1,r1,'-k')
ax5.set_ylabel('R [AU]')
ax5.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d') )
ax5.set_xlim(t1[0],t1[-1])


ax6=plt.subplot(426,sharex=ax2)
plt.plot_date(t2,r2,'-k')
ax6.set_ylabel('R [AU]')
ax6.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d') )
ax6.set_xlim(t2[0],t2[-1])


ax7=plt.subplot(427,sharex=ax1)
plt.plot_date(t1,np.rad2deg(lon1),'-r')
ax7.set_xlim(t1[0],t1[-1])
ax7.set_ylabel('HEEQ longitude [deg]')
ax7.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d') )


ax8=plt.subplot(428,sharex=ax2)
plt.plot_date(t2,np.rad2deg(lon2),'-r')
ax8.set_xlim(t2[0],t2[-1])
ax8.set_ylabel('HEEQ longitude [deg]')
ax8.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d') )


plt.tight_layout()



plt.show()



plt.savefig('results/psp_orbit12.jpg')

sys.exit()

#hd.convert_MAVEN_mat_to_pickle()




# use
# import importlib
# importlib.reload(cme_stats_module)
# to update module while working in command line 

#set breakpoints with pdb.set_trace



# LIST OF FUNCTIONS

# load_url_current_directory
# getpositions
# getcat
# gaussian
# dynamic_pressure
# decode_array
# time_to_num_cat
# get_omni2_data


def load_url_current_directory(filename,url):
#loads a file from any url to the current directory
#I use owncloud for the direct url links, 
#also works for dropbox when changing the last 0 to 1 in the url-> gives a direct link to files

 if not os.path.exists(filename):
  print('download file ', filename, ' from')
  print(url)
  try: 
    urllib.request.urlretrieve(url, filename)
    print('done')
  except urllib.error.URLError as e:
    print(' ', data_url,' ',e.reason)




def getpositions(filename):  
    pos=scipy.io.readsav(filename)  
    print
    print('positions file:', filename) 
    return pos


def getcat(filename):
   cat=scipy.io.readsav(filename)#, verbose='true')  
   return cat  
  
  
def gaussian(x, amp, mu, sig):
   return amp * exp(-(x-cen)**2 /wid)



def dynamic_pressure(density, speed):
   # make dynamic pressure from density and speed
   #assume pdyn is only due to protons
   #pdyn=np.zeros(len([density])) #in nano Pascals
   protonmass=1.6726219*1e-27  #kg
   pdyn=np.multiply(np.square(speed*1e3),density)*1e6*protonmass*1e9  #in nanoPascal
   return pdyn
  
def decode_array(bytearrin):
  #for decoding the strings from the IDL .sav file to a list of python strings, not bytes 
  #make list of python lists with arbitrary length
  bytearrout= ['' for x in range(len(bytearrin))]
  for i in range(0,len(bytearrin)-1):
    bytearrout[i]=bytearrin[i].decode()
  #has to be np array so to be used with numpy "where"
  bytearrout=np.array(bytearrout)
  return bytearrout  
  
def time_to_num_cat(time_in):  

  #for time conversion from catalogue .sav to numerical time
  #this for 1-minute data or lower time resolution

  #for all catalogues
  #time_in is the time in format: 2007-11-17T07:20:00 or 2007-11-17T07:20Z
  #for times help see: 
  #http://docs.sunpy.org/en/latest/guide/time.html
  #http://matplotlib.org/examples/pylab_examples/date_demo2.html
  
  j=0
  #time_str=np.empty(np.size(time_in),dtype='S19')
  time_str= ['' for x in range(len(time_in))]
  #=np.chararray(np.size(time_in),itemsize=19)
  time_num=np.zeros(np.size(time_in))
  
  for i in time_in:

    #convert from bytes (output of scipy.readsav) to string
    time_str[j]=time_in[j][0:16].decode()+':00'
    year=int(time_str[j][0:4])
    time_str[j]
    #convert time to sunpy friendly time and to matplotlibdatetime
    #only for valid times so 9999 in year is not converted
    #pdb.set_trace()
    if year < 2100:
     	  time_num[j]=mdates.date2num(sunpy.time.parse_time(time_str[j]))
    j=j+1  
    #the date format in matplotlib is e.g. 735202.67569444
    #this is time in days since 0001-01-01 UTC, plus 1.
   
    #return time_num which is already an array and convert the list of strings to an array
  return time_num, np.array(time_str)



def get_omni2_data():

  #FORMAT(2I4,I3,I5,2I3,2I4,14F6.1,F9.0,F6.1,F6.0,2F6.1,F6.3,F6.2, F9.0,F6.1,F6.0,2F6.1,F6.3,2F7.2,F6.1,I3,I4,I6,I5,F10.2,5F9.2,I3,I4,2F6.1,2I6,F5.1)
 #1963   1  0 1771 99 99 999 999 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 9999999. 999.9 9999. 999.9 999.9 9.999 99.99 9999999. 999.9 9999. 999.9 999.9 9.999 999.99 999.99 999.9  7  23    -6  119 999999.99 99999.99 99999.99 99999.99 99999.99 99999.99  0   3 999.9 999.9 99999 99999 99.9

#define variables from OMNI2 dataset
 #see http://omniweb.gsfc.nasa.gov/html/ow_data.html

 #omni2_url='ftp://nssdcftp.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
 
 #check how many rows exist in this file
 f=open('omni2_all_years.dat')
 dataset= len(f.readlines())
 #print(dataset)
 #global Variables
 spot=np.zeros(dataset) 
 btot=np.zeros(dataset) #floating points
 bx=np.zeros(dataset) #floating points
 by=np.zeros(dataset) #floating points
 bz=np.zeros(dataset) #floating points
 bzgsm=np.zeros(dataset) #floating points
 bygsm=np.zeros(dataset) #floating points

 speed=np.zeros(dataset) #floating points
 speedx=np.zeros(dataset) #floating points
 speed_phi=np.zeros(dataset) #floating points
 speed_theta=np.zeros(dataset) #floating points

 dst=np.zeros(dataset) #float
 kp=np.zeros(dataset) #float

 den=np.zeros(dataset) #float
 pdyn=np.zeros(dataset) #float
 year=np.zeros(dataset)
 day=np.zeros(dataset)
 hour=np.zeros(dataset)
 t=np.zeros(dataset) #index time
 
 
 j=0
 print('Read OMNI2 data ...')
 with open('omni2_all_years.dat') as f:
  for line in f:
   line = line.split() # to deal with blank 
   #print line #41 is Dst index, in nT
   dst[j]=line[40]
   kp[j]=line[38]
   
   if dst[j] == 99999: dst[j]=np.NaN
   #40 is sunspot number
   spot[j]=line[39]
   if spot[j] == 999: spot[j]=NaN

   #25 is bulkspeed F6.0, in km/s
   speed[j]=line[24]
   if speed[j] == 9999: speed[j]=np.NaN
 
   #get speed angles F6.1
   speed_phi[j]=line[25]
   if speed_phi[j] == 999.9: speed_phi[j]=np.NaN

   speed_theta[j]=line[26]
   if speed_theta[j] == 999.9: speed_theta[j]=np.NaN
   #convert speed to GSE x see OMNI website footnote
   speedx[j] = - speed[j] * np.cos(np.radians(speed_theta[j])) * np.cos(np.radians(speed_phi[j]))



   #9 is total B  F6.1 also fill ist 999.9, in nT
   btot[j]=line[9]
   if btot[j] == 999.9: btot[j]=np.NaN

   #GSE components from 13 to 15, so 12 to 14 index, in nT
   bx[j]=line[12]
   if bx[j] == 999.9: bx[j]=np.NaN
   by[j]=line[13]
   if by[j] == 999.9: by[j]=np.NaN
   bz[j]=line[14]
   if bz[j] == 999.9: bz[j]=np.NaN
 
   #GSM
   bygsm[j]=line[15]
   if bygsm[j] == 999.9: bygsm[j]=np.NaN
 
   bzgsm[j]=line[16]
   if bzgsm[j] == 999.9: bzgsm[j]=np.NaN 	
 
 
   #24 in file, index 23 proton density /ccm
   den[j]=line[23]
   if den[j] == 999.9: den[j]=np.NaN
 
   #29 in file, index 28 Pdyn, F6.2, fill values sind 99.99, in nPa
   pdyn[j]=line[28]
   if pdyn[j] == 99.99: pdyn[j]=np.NaN 		
 
   year[j]=line[0]
   day[j]=line[1]
   hour[j]=line[2]
   j=j+1     
   

 #convert time to matplotlib format
 #http://docs.sunpy.org/en/latest/guide/time.html
 #http://matplotlib.org/examples/pylab_examples/date_demo2.html

 times1=np.zeros(len(year)) #datetime time
 print('convert time start')
 for index in range(0,len(year)):
      #first to datetimeobject 
      timedum=datetime.datetime(int(year[index]), 1, 1) + datetime.timedelta(day[index] - 1) +datetime.timedelta(hours=hour[index])
      #then to matlibplot dateformat:
      times1[index] = mdates.date2num(timedum)
 print('convert time done')   #for time conversion

 print('all done.')
 print(j, ' datapoints')   #for reading data from OMNI file
 
 #make structured array of data
 omni_data=np.rec.array([times1,btot,bx,by,bz,bygsm,bzgsm,speed,speedx,den,pdyn,dst,kp,spot], \
 dtype=[('time','f8'),('btot','f8'),('bx','f8'),('by','f8'),('bz','f8'),\
 ('bygsm','f8'),('bzgsm','f8'),('speed','f8'),('speedx','f8'),('den','f8'),('pdyn','f8'),('dst','f8'),('kp','f8'), ('spot','f8')])
 
 return omni_data
 
 



import pickle
import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import sunpy.time




def IDL_time_to_num(time_in):  
 #convert IDL time to matplotlib datetime
 time_num=np.zeros(np.size(time_in))
 for ii in np.arange(0,np.size(time_in)):
   time_num[ii]=mdates.date2num(sunpy.time.parse_time(time_in[ii]))
   
 return time_num 
  

print ('read VEX')
#get insitu data
vex= pickle.load( open( "../catpy/DATACAT/VEX_2007to2014_SCEQ_removed.p", "rb" ) )

#time conversion
vex_time=IDL_time_to_num(vex.time)
print( 'read VEX done.')



plt.figure(1)
plt.plot_date(vex_time,vex.btot,'-k')
plt.show()

'''
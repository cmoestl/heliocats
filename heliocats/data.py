#data.py
#load data for heliocats
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy
import copy
import sunpy
import matplotlib.dates as mdates
import datetime
import urllib
import json
import os
import pdb
import scipy.io
import pickle
import sys
import cdflib
import matplotlib.pyplot as plt
import heliosat
from numba import njit
from sunpy.time import parse_time
import heliopy.data.cassini as cassinidata
import heliopy.data.helios as heliosdata
import heliopy.data.spice as spicedata
import heliopy.spice as spice
import astropy

data_path='/nas/helio/data/insitu_python/'




####################################### get new data ####################################



#see https://github.com/ajefweiss/HelioSat/blob/master/heliosat/json/spacecraft.json

def save_wind_data(file):
    
    print('start wind')
    wind_sat = heliosat.WIND()
    t_start = datetime.datetime(2018, 1, 1)
    t_end = datetime.datetime(2019, 11, 30)
    
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
    [r, lat, lon]=cart2sphere(earth.x,earth.y,earth.z)
    print('position end ')
    
    
    
    
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
        
    #p0 = np.interp(time_mat, tp_mat, pro[:,0])
    #p1 = np.interp(time_mat, tp_mat, pro[:,1])
    v = np.interp(time_mat, tp_mat, pro[:,1])
    #p3 = np.interp(time_mat, tp_mat, pro[:,3])
    #p4 = np.interp(time_mat, tp_mat, pro[:,4])

    
    #make array
    win=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('p0', float),('v', float),('p2', float),\
                ('p3', float),('p4', float),('x', float),('y', float),('z', float),\
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
    win.lat=lat
    win.lon=lon

    
    #wind.p0=p0
    #wind.p1=p1    
    win.v=v    
    #wind.p3=p3
    #wind.p4=p4
    
       

    #pickle.dump([tm,mag, tp,pro], open(file, "wb"))
    #[tm,mag, tp,pro]=pickle.load(open( "data/wind_oct2018_may2019.p", "rb" ) )  
    pickle.dump(win, open(file, "wb"))
    
    
    print('wind done')
    print()
    



 
def save_stereoa_data(file):

    print('start STA')
    sta_sat = heliosat.STA()
    t_start = datetime.datetime(2018, 1, 1)
    t_end = datetime.datetime(2019, 5, 31)
     
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    #tm, mag = sta_sat.get_data_raw(t_start, t_end, "sta_impact_l1")
    #tp, pro = sta_sat.get_data_raw(t_start, t_end, "sta_plastic_l2")

    tm, mag = sta_sat.get_data_raw(t_start, t_end, "sta_impact_beacon")
    tp, pro = sta_sat.get_data_raw(t_start, t_end, "sta_plastic_beacon")

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
    statra=spice.Trajectory('-234')
    statra.generate_positions(time,'Sun',frame)
    statra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(statra.x,statra.y,statra.z)
    print('position end ')
    
    
  
        
    
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
        
    #p0 = np.interp(time_mat, tp_mat, pro[:,0])
    #p1 = np.interp(time_mat, tp_mat, pro[:,1])
    #v = np.interp(time_mat, tp_mat, pro[:,1])
    #p3 = np.interp(time_mat, tp_mat, pro[:,3])
    #p4 = np.interp(time_mat, tp_mat, pro[:,4])

    
    #make array
    sta=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('p0', float),('v', float),('p2', float),\
                ('p3', float),('p4', float),('x', float),('y', float),('z', float),\
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
    sta.y=statra.x
    sta.z=statra.x
    
    sta.r=r
    sta.lat=lat
    sta.lon=lon


    
    #sta.p0=p0
    #sta.p1=p1    
    #sta.v=v    
    #sta.p3=p3
    #sta.p4=p4
    
       

    #pickle.dump([tm,mag, tp,pro], open(file, "wb"))
    #[tm,mag, tp,pro]=pickle.load(open( "data/sta_oct2018_may2019.p", "rb" ) )  
    pickle.dump(sta, open(file, "wb"))
    
    

    #pickle.dump([tm,mag, tp,pro], open(file, "wb"))
    #[tm,mag, tp,pro]=pickle.load(open( "data/sta_oct2018_may2019.p", "rb" ) )  
    #pickle.dump(sta, open(file, "wb"))
    
       
    
    

    #pickle.dump([tm,mag, tp, pro], open(file, "wb"))
    print('done sta')
    print()





    
def save_psp_data(path, file):
     
    print('start PSP')
     
    psp_sat = heliosat.PSP()
    t_start = datetime.datetime(2018, 10, 14,14,14, 30)
    #t_end = datetime.datetime(2018, 12, 12,23,59,30)
    t_end = datetime.datetime(2019, 5, 31,23,59,30)
    
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
    t_end = datetime.datetime(2019, 5, 31,23,59,0)
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
    tp = np.interp(time_mat, tp_mat, pro[:,1])
    vx = np.interp(time_mat, tp_mat, pro[:,2])
    vy = np.interp(time_mat, tp_mat, pro[:,3])
    vz = np.interp(time_mat, tp_mat, pro[:,4])
    
    den[setnan]=np.nan
    tp[setnan]=np.nan
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
    psp.lat=lat
    psp.lon=lon
  
    psp.vt=vt
    psp.vx=vx    
    psp.vy=vy  
    psp.vz=vz
    psp.np=den
    #https://en.wikipedia.org/wiki/Thermal_velocity convert from km/s to K
    from astropy.constants import m_p,k_B
    psp.tp=np.pi*m_p*((tp*1e3)**2)/(8*k_B) 
    
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








def convert_MAVEN_mat_to_pickle(data_path):

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

    mav.tp=mavraw['Tp'][:,0]      
    mav.np=mavraw['np'][:,0]      
    
    
    smooth=0
    if smooth >0:
    
       print('smoothing')
       #smooth with median for each orbit, take times of apogees (search with scipy)
       #***np.gradient for getting maxima and check sign reversal
       #
       #
    
    
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
    if smooth==0: 
         header='MAVEN merged magnetic field and plasma data, obtained from Toulouse. '+ \
         'Timerange: '+mav.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+mav.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
         '. The magnetosphere is removed with the Gruesbeck et al. 3D model (C.S. Wedlund). '+ \
         'Units are btxyz [nT, MSO], vtxyz [km/s, MSO], np [#/cm-3], tp[?], orbital position: '+ \
         'xo/yo/zo/ro/lono/lato [km, degree, MSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
         
         pickle.dump([mav,header], open(data_path+'maven_2014_2018_removed.p', "wb"))

    if smooth>0: 
         header='MAVEN merged magnetic field and plasma data, obtained from Toulouse. '+ \
         'Timerange: '+mav.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+mav.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
         '. The magnetosphere is removed with the Gruesbeck et al. 3D model (C.S. Wedlund). '+ \
         'Units are btxyz [nT], vtxyz [km/s], np [#/cm-3], tp[?], orbital position: '+ \
         'xo/yo/zo/ro/lono/lato [km, MSO], heliospheric position x/y/z/r/lon/lat [AU, HEEQ]'
                
         pickle.dump([mav,header], open(data_path+'maven_2014_2018_removed_smoothed.p', "wb"))

    







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






def save_omni_data(file):

    return 0









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
    







###################################### ICMECAT operations ################################


def load_helcats_icmecat_master_from_excel(file):

    print('load HELCATS ICMECAT from file:', file)
    ic=pd.read_excel(file)

    #convert times to datetime objects
    for i in np.arange(0,ic.shape[0]):    
    
        a=str(ic.icme_start_time[i]).strip() #remove leading and ending blank spaces if any
        ic.icme_start_time.loc[i]=sunpy.time.parse_time(a).datetime

        a=str(ic.mo_start_time[i]).strip() #remove leading and ending blank spaces if any
        ic.mo_start_time.loc[i]=sunpy.time.parse_time(a).datetime 

        a=str(ic.mo_end_time[i]).strip() #remove leading and ending blank spaces if any
        ic.mo_end_time.loc[i]=sunpy.time.parse_time(a).datetime 

        
        a=str(ic.icme_end_time[i]).strip() #remove leading and ending blank spaces if any
        if a!= '9999-99-99T99:99Z':
            ic.icme_end_time.loc[i]=sunpy.time.parse_time(a).datetime 
        else: ic.icme_end_time.loc[i]=np.nan

    return ic























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
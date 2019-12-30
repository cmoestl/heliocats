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
import heliosat
from numba import njit
from sunpy.time import parse_time
import heliopy.data.spice as spicedata
import heliopy.spice as spice
import astropy



###################################### HELCATS #############################################




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











######################################################### MAVEN ####################################
'''
Hi Christian,

I just did a run on the MAVEN datasets that you had with my new technique (and using the Gruesbeck model). The mat file is here:

https://oeawcloud.oeaw.ac.at/index.php/s/Hjc8BVmtQT3k0ub

Password is:      maven2019

The save command in matlab I have used is the following:
save('Data-MAVEN-SolarWind.mat', 'timeD','np','Vx','Vy','Vz','VT','Tp','Bx','By','Bz','BT','Xsc','Ysc','Zsc');

Note that Xsc, Ysc, and Zsc are given in units of Mars radius (Rp = 3389.5 km). All other ones are the original units of the cdf file you gave me. I have also gone through the data and set to NaN all negative values of density np, temperature Tp and total magnetic field BT.

It would be interesting to compare your old results with these ones, also as a double check that the data was filtered correctly, although the 3D model of Gruesbeck+ 2018 does not assume an aberration angle of 4 degrees like the polar model of Edberg+ 2008 -- which I am not sure you took into account originally when processing the data.

Cheers,

Cyril

P.S. I'm writing a compendium for all the details of the technique I used. Hopefully this should be finished soon.

Hallo Christian,

Super, happy that it looks fine! :) Yes, agreed with the median filtering, this should take care of the spikes. 
There were also many strange spikes in Xsc, Ysc and Zsc (very large values >2.9e27, probably due to an issue with the SPICE kernel, 
about 1685 points, i.e., 0.087% of the data), so I set all of these anomalous data points to NaN too (all variables including Xsc, Ysc, Zsc).

Cheers,

Cyril


Mehr anzeigen von Christian Möstl






'''



@njit
def cart2sphere(x,y,z):
    r = np.sqrt(x**2+ y**2 + z**2)            # r
    theta = np.arctan2(z,np.sqrt(x**2+ y**2))     # theta
    phi = np.arctan2(y,x)                        # phi
    return (r, theta, phi)
    



def convert_MAVEN_mat_to_pickle():

    print('load MAVEN from MAT')
    file='data/MAVEN_2014to2018_removed_cyril.mat'
    mavraw = scipy.io.loadmat(file)
    
    
    #now make recarray
    #make array 
    #******************** no 2D array
    mav=np.zeros([len(mavraw['Bx']),len(mavraw)],dtype=[('time',object),('bx', float),('by', float),('bz', float),('bt', float),('tp', float),('np', float),('vt', float),('vx', float),('vy', float),('vz', float),('xsc', float),('ysc', float),('zsc', float),('r', float),('lat', float),('lon', float)])   
    #convert to recarray
    mav = mav.view(np.recarray)  
    mav.time=mdates.num2date(mavraw['timeD'])      
    mav.bx=mavraw['Bx']      
    mav.by=mavraw['By']      
    mav.bz=mavraw['Bz']      
    mav.bt=mavraw['BT']      
    
    mav.tp=mavraw['Tp']      
    mav.np=mavraw['np']      
    
    
    
    #add position
    
    
    
    #smooth with median for each orbit, take times of apogees (search with scipy)
    
    print('save MAVEN as pickle')
    pickle.dump(mav, open("data/MAVEN_2014to2018_removed_cyril_2.p", "wb"))

    



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
    win=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),('bz', float),('bt', float),('p0', float),('v', float),('p2', float),('p3', float),('p4', float),('x', float),('y', float),('z', float),('r', float),('lat', float),('lon', float)])   
       
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
    sta=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),('bz', float),('bt', float),('p0', float),('v', float),('p2', float),('p3', float),('p4', float),('x', float),('y', float),('z', float),('r', float),('lat', float),('lon', float)])   
       
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





    
def save_psp_data(file):
     
    print('start PSP')
     
    psp_sat = heliosat.PSP()
    t_start = datetime.datetime(2018, 10, 15)
    t_end = datetime.datetime(2019, 5, 31)
    
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
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
        
    p0 = np.interp(time_mat, tp_mat, pro[:,0])
    p1 = np.interp(time_mat, tp_mat, pro[:,1])
    v = np.interp(time_mat, tp_mat, pro[:,2])
    p3 = np.interp(time_mat, tp_mat, pro[:,3])
    p4 = np.interp(time_mat, tp_mat, pro[:,4])

    
    #make array
    psp=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),('bz', float),('bt', float),('p0', float),('p1', float),('v', float),('p3', float),('p4', float),('x', float),('y', float),('z', float),('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    psp = psp.view(np.recarray)  

    #fill with data
    psp.time=time
    psp.bx=bx
    psp.by=by
    psp.bz=bz 
    psp.bz=bt
    
    psp.x=psptra.x
    psp.y=psptra.x
    psp.z=psptra.x
    
    psp.r=r
    psp.lat=lat
    psp.lon=lon

    
    
    psp.p0=p0
    psp.p1=p1    
    psp.v=v    
    psp.p3=p3
    psp.p4=p4
    
       

    #pickle.dump([tm,mag, tp,pro], open(file, "wb"))
    #[tm,mag, tp,pro]=pickle.load(open( "data/psp_oct2018_may2019.p", "rb" ) )  
    pickle.dump(psp, open(file, "wb"))

    print('done psp')
    print()

  


def save_ulysses_data():

   
    print('read Ulysses data from cdf and convert to pickle')   

    #load cdf
    ulycdf = cdflib.CDF('data/ulysses_1990_2009_CDAWEB.cdf') 
    #check variables
    #ulycdf.cdf_info()

    #time conversion to datetime      
    time=ulycdf.varget('Epoch')
    t=parse_time(time,format='cdf_epoch').datetime  
    
    
    #cut so that it starts with available position on Oct 6 1990
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
    
    
    uly=np.zeros(len(t),dtype=[('time',object),('bx', float),('by', float),('bz', float),('bt', float),('np', float),('vp', float),('tp', float),('x', float),('y', float),('z', float),('r', float),('lat', float),('lon', float)])   
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
    

    file='data/ulysses.p'
    pickle.dump(uly, open(file, "wb"))
    
    print('Ulysses done')

    return 0
  
  
  
def save_helcats_datacat():  
    ''' to save all of helcats DATACAT into a single file'''
      
    print('save all helcats DATACAT into single file')

    datacat_path='/nas/helio/data/DATACAT/'

    print( 'read MESSENGER')
    #get insitu data
    mes= pickle.load( open( datacat_path+"MES_2007to2015_SCEQ_removed.p", "rb" ) )
    #time conversion
    mes_time=parse_time(mes.time,format='utime').datetime
    #replace mes.time with datetime object
    mes=mes.astype([(('time', 'TIME'), 'object'), (('btot', 'BTOT'), '>f8'), (('bx', 'BX'), '>f8'), (('by', 'BY'), '>f8'), (('bz', 'BZ'), '>f8'), (('mes_radius_in_km_heeq', 'MES_RADIUS_IN_KM_HEEQ'), '>f8'), (('mes_latitude_in_radians_heeq', 'MES_LATITUDE_IN_RADIANS_HEEQ'), '>f8'), (('mes_longitude_in_radians_heeq', 'MES_LONGITUDE_IN_RADIANS_HEEQ'), '>f8')])       
    mes.time=mes_time
    print( 'read MESSENGER done.')

    print ('read VEX')
    vex= pickle.load( open(datacat_path+ "VEX_2007to2014_SCEQ_removed.p", "rb" ) )
    vex_time=parse_time(vex.time,format='utime').datetime
    vex=vex.astype([(('time', 'TIME'), 'object'), (('btot', 'BTOT'), '>f8'), (('bx', 'BX'), '>f8'), (('by', 'BY'), '>f8'), (('bz', 'BZ'), '>f8'), (('vex_radius_in_km_heeq', 'VEX_RADIUS_IN_KM_HEEQ'), '>f8'), (('vex_latitude_in_radians_heeq', 'VEX_LATITUDE_IN_RADIANS_HEEQ'), '>f8'), (('vex_longitude_in_radians_heeq', 'VEX_LONGITUDE_IN_RADIANS_HEEQ'), '>f8')])       
    vex.time=vex_time
    print( 'read VEX done.')

    print( 'read Wind')
    win= pickle.load( open(datacat_path+ "WIND_2007to2018_HEEQ.p", "rb" ) )
    win_time=parse_time(win.time,format='utime').datetime
    win=win.astype([(('time', 'TIME'),'object'), (('btot', 'BTOT'), '>f8'), (('bx', 'BX'), '>f8'), (('by', 'BY'), '>f8'), (('bz', 'BZ'), '>f8'), (('vtot', 'VTOT'), '>f8'), (('vx', 'VX'), '>f8'), (('vy', 'VY'), '>f8'), (('vz', 'VZ'), '>f8'), (('temperature', 'TEMPERATURE'), '>f8'), (('density', 'DENSITY'), '>f8'), (('win_radius_in_km_heeq', 'WIN_RADIUS_IN_KM_HEEQ'), '>f8'), (('win_latitude_in_radians_heeq', 'WIN_LATITUDE_IN_RADIANS_HEEQ'), '>f8'), (('win_longitude_in_radians_heeq', 'WIN_LONGITUDE_IN_RADIANS_HEEQ'), '>f8')])       
    win.time=win_time
    print( 'read Wind done.')

    print( 'read STEREO-A')
    sta= pickle.load( open(datacat_path+ "STA_2007to2015_SCEQ.p", "rb" ) )
    sta_time=parse_time(sta.time,format='utime').datetime
    sta=sta.astype([(('time', 'TIME'),'object'), (('btot', 'BTOT'), '>f8'), (('bx', 'BX'), '>f8'), (('by', 'BY'), '>f8'), (('bz', 'BZ'), '>f8'), (('vtot', 'VTOT'), '>f8'), (('vx', 'VX'), '>f8'), (('vy', 'VY'), '>f8'), (('vz', 'VZ'), '>f8'), (('temperature', 'TEMPERATURE'), '>f8'), (('density', 'DENSITY'), '>f8'), (('sta_radius_in_km_heeq', 'STA_RADIUS_IN_KM_HEEQ'), '>f8'), (('sta_latitude_in_radians_heeq', 'STA_LATITUDE_IN_RADIANS_HEEQ'), '>f8'), (('sta_longitude_in_radians_heeq', 'STA_LONGITUDE_IN_RADIANS_HEEQ'), '>f8')])       
    sta.time=sta_time
    print( 'read STA done.')

    print( 'read STEREO-B')
    stb= pickle.load( open(datacat_path+ "STB_2007to2014_SCEQ.p", "rb" ) )
    stb_time=parse_time(stb.time,format='utime').datetime
    stb=stb.astype([(('time', 'TIME'),'object'), (('btot', 'BTOT'), '>f8'), (('bx', 'BX'), '>f8'), (('by', 'BY'), '>f8'), (('bz', 'BZ'), '>f8'), (('vtot', 'VTOT'), '>f8'), (('vx', 'VX'), '>f8'), (('vy', 'VY'), '>f8'), (('vb', 'VZ'), '>f8'), (('temperature', 'TEMPERATURE'), '>f8'), (('density', 'DENSITY'), '>f8'), (('stb_radius_in_km_heeq', 'STB_RADIUS_IN_KM_HEEQ'), '>f8'), (('stb_latitude_in_radians_heeq', 'STB_LATITUDE_IN_RADIANS_HEEQ'), '>f8'), (('stb_longitude_in_radians_heeq', 'STB_LONGITUDE_IN_RADIANS_HEEQ'), '>f8')])       
    stb.time=stb_time
    print( 'read STB done.')
    
    
    print( 'read Ulysses from CDAWEB cdf')    
    hd.save_ulysses_data()
    fileuly='data/ulysses.p'
    uly=pickle.load(open(fileuly, 'rb' ) )


    

    pickle.dump([vex,win,mes,sta,stb,uly], open(datacat_path+ "helcats_all_data.p", "wb" ) )





  
def load_helcats_datacat(file):  
    ''' to load all of helcats DATACAT from a single file'''
    
    print('load all helcats DATACAT from single file: ', file)
    [vex,win,mes,sta,stb]=pickle.load( open(file, "rb" ) )
    print('use vex,win,sta,stb,mes to access data and position')
    return [vex,win,mes,sta,stb]
    
     


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
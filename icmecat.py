#!/usr/bin/env python
# coding: utf-8

# ## icmecat
# 
# makes the interplanetary coronal mass ejection catalog ICMECAT
# 
# available at https://helioforecast.space/icmecat
# 
# **Author**: C. Möstl, IWF Graz, Austria
# twitter @chrisoutofspace, part of https://github.com/cmoestl/heliocats
# 
# **latest release: version 2.0, 2020 June 3, updated 2020 September 28 doi: 10.6084/m9.figshare.6356420
# 
# Install a specific conda environment to run this code, see readme at https://github.com/cmoestl/heliocats
# 
# **Adding a new ICME event:** 
# - edit the file icmecat/HELCATS_ICMECAT_v20_master.xlsx to add 3 times, the id and spacecraft name
# - delete the file for the respective spacecraft under icmecat/indices_icmecat/
# - run section 3 in this notebook or script.
# 
# Convert this notebook to a script with jupyter nbconvert --to script icmecat.ipynb (automatically done in first cell)
# 
# see also https://github.com/helioforecast/Papers/blob/master/Moestl2020_PSP_rate/icme_rate.ipynb
# 
# 
# 

# In[1]:


import numpy as np
import scipy.io
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import  DateFormatter
from datetime import timedelta
import seaborn as sns
import datetime
import astropy
import astropy.constants as const
from sunpy.time import parse_time
import time
import pickle
import sys
import os
import urllib
import json
import importlib
import pandas as pd
import copy
import openpyxl
import h5py


import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames



from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging

from heliocats import data as hd
importlib.reload(hd) #reload again while debugging

from heliocats import cats as hc
importlib.reload(hc) #reload again while debugging

from heliocats import stats as hs
importlib.reload(hs) #reload again while debugging

#where the in situ data files are located is read 
#from config.py 
import config
importlib.reload(config)
from config import data_path
from config import data_path_ML


########### make directories first time if not there

resdir='results'
if os.path.isdir(resdir) == False: os.mkdir(resdir)

datadir='data'
if os.path.isdir(datadir) == False: os.mkdir(datadir)

indexdir='icmecat/indices_icmecat' 
if os.path.isdir(indexdir) == False: os.mkdir(indexdir) 

catdir='icmecat'
if os.path.isdir(catdir) == False: os.mkdir(catdir)

icplotsdir='icmecat/plots_icmecat/' 
if os.path.isdir(icplotsdir) == False: os.mkdir(icplotsdir) 

#Convert this notebook to a script with jupyter nbconvert --to script icmecat.ipynb
os.system('jupyter nbconvert --to script icmecat.ipynb')    

import warnings
warnings.filterwarnings('ignore')


# ## (0) process in situ data into similar format

# In[2]:


# make data
# from heliocats import data as hd
# importlib.reload(hd) #reload again while debugging

# # ################################# PSP


################################## USE THIS ################################
# load PSP data from server on linux command line onto leo server
# go to heliosat directory /nas/helio/data/heliosat/data/psp_fields_l2
# wget -nc "ftps://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_1min/2019/*.cdf"
# wget -nc "ftps://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_1min/2020/*.cdf"
# psp_spc_l3
# wget -nc "ftps://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spc/l3/l3i/2019/*.cdf"
# wget -nc "ftps://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spc/l3/l3i/2020/*.cdf"
############################################################################


# print('load PSP data') #from heliosat, converted to SCEQ similar to STEREO-A/B

#filepsp='psp_2018_2020_rtn.p'
#hd.save_psp_data(data_path,filepsp, sceq=False)   
#[psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) )  
#print('done')

#filepsp='psp_2018_2020_sceq.p'
#hd.save_psp_data(data_path,filepsp, sceq=True)   
#[psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) )  
# plt.xlim(parse_time('2007-08-15').plot_date,parse_time('2007-08-15 12:00').plot_date)


# ################################# Wind

# filewin="wind_2018_2019_gse.p" 
# start=datetime.datetime(2018, 1, 1)
# end=datetime.datetime(2020, 1, 1)
# hd.save_wind_data(data_path,filewin,start,end,heeq=False)

# filewin="wind_2018_2019_heeq.p" 
# start=datetime.datetime(2018, 1, 1)
# #end=datetime.datetime(2019, 12, 31)
# end=datetime.datetime.utcnow()
# hd.save_wind_data(data_path,filewin,start,end,heeq=True)

# filewin1="wind_2007_2018_heeq_helcats.p" 
# [win1,hwin1]=pickle.load(open(data_path+filewin1, "rb" ) )  
# filewin2="wind_2018_2019_heeq.p" 
# [win2,hwin2]=pickle.load(open(data_path+filewin2, "rb" ) )  
# filewin3="wind_2018_2019_gse.p" 
# [win3,hwin3]=pickle.load(open(data_path+filewin3, "rb" ) )  



################################### STEREO-A
#filesta="stereoa_2019_2020_sceq_beacon.p" 
#start=datetime.datetime(2019, 1, 1)
#end=datetime.datetime.utcnow()
#hd.save_stereoa_beacon_data(data_path,filesta,start,end,sceq=True)
#[sta,hsta]=pickle.load(open(data_path+filesta, "rb" ) ) 





########################## SAVE MSL rad data into recarray as pickle
#hd.save_msl_rad()



# In[3]:


############################# make Ulysses files
#hd.save_ulysses_data(data_path)


############################# make STEREO-A data files

#filesta_all='stereoa_2007_2019_rtn.p'
#hd.save_all_stereoa_science_data(data_path, filesta_all,sceq=False)
#[sa1,hsa1]=pickle.load(open(data_path+filesta_all, "rb" ) )  

#filesta_all='stereoa_2007_2019_sceq.p'
#hd.save_all_stereoa_science_data(data_path, filesta_all,sceq=True)
#[sa2,hsa2]=pickle.load(open(data_path+filesta_all, "rb" ) )  

# plt.plot(sta.time,sta.by,'-g',linewidth=5)
# plt.plot(sa2.time,sa2.by,'-k')
# plt.plot(sa1.time,sa1.by,'-b')
# plt.plot(sta.time,sta.lat,'-r')

# plt.xlim(parse_time('2007-08-15').plot_date,parse_time('2007-08-15 12:00').plot_date)
# plt.ylim(-5,4)

# #merge STEREO-A old and new data    
# #make array
# sta=np.zeros(np.size(sta1.time)+np.size(sta2.time),dtype=[('time',object),('bx', float),('by', float),\
#             ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
#             ('x', float),('y', float),('z', float),\
#             ('r', float),('lat', float),('lon', float)])   

# #convert to recarray
# sta = sta.view(np.recarray)  
# #add merged variables
# sta.time=np.hstack((sta1.time,sta2.time))
# sta.bx=np.hstack((sta1.bx,sta2.bx))
# sta.by=np.hstack((sta1.by,sta2.by))
# sta.bz=np.hstack((sta1.bz,sta2.bz))
# sta.bt=np.hstack((sta1.bt,sta2.bt))
# sta.vt=np.hstack((sta1.vt,sta2.vt))
# sta.np=np.hstack((sta1.np,sta2.np))
# sta.tp=np.hstack((sta1.tp,sta2.tp))
# sta.x=np.hstack((sta1.x,sta2.x))
# sta.y=np.hstack((sta1.y,sta2.y))
# sta.z=np.hstack((sta1.z,sta2.z))
# sta.r=np.hstack((sta1.r,sta2.r))
# sta.lon=np.hstack((sta1.lon,sta2.lon))
# sta.lat=np.hstack((sta1.lat,sta2.lat))


############################## make STEREO-B data files
# STEREO-B

# filestb_all='stereob_2007_2014_rtn.p'
# hd.save_all_stereob_science_data(data_path, filestb_all,sceq=False)
# [sb1,hsb1]=pickle.load(open(data_path+filestb_all, "rb" ) )  

# filestb_all='stereob_2007_2014_sceq.p'
# hd.save_all_stereob_science_data(data_path, filestb_all,sceq=True)
# [sb2,hsb2]=pickle.load(open(data_path+filestb_all, "rb" ) )  

# plt.plot(stb.time,stb.by,'-g',linewidth=5)
# plt.plot(sb2.time,sb2.by,'-k')
# plt.plot(sb1.time,sb1.by,'-b')
# #plt.plot(stb.time,stb.lat,'-r')

# plt.xlim(parse_time('2007-08-15').plot_date,parse_time('2007-08-15 12:00').plot_date)
# plt.ylim(-5,4)



############################## make Wind data files

#filewin="wind_2018_2019_gse.p" 
#for updating data
#start=datetime.datetime(2018, 1, 1)
#end=datetime.datetime(2020, 1, 1)
#hd.save_wind_data(data_path,filewin,start,end,heeq=False)


########## check GSE to HEEQ conversion

#save=1

#if save > 0:

#    filewin="wind_2018_2020_sept_heeq.p" 
#    start=datetime.datetime(2018, 1 , 1)
#    end=datetime.datetime(2020, 8, 31)

    #filewin="wind_heeq_test.p" 
    #start=datetime.datetime(2019, 8, 1)
    #end=datetime.datetime(2019, 10, 1)


    #end=datetime.datetime.utcnow()
#    hd.save_wind_data(data_path,filewin,start,end,heeq=True)

    
    #filewin="wind_2018_2020_sept_gse.p" 
    #start=datetime.datetime(2018, 1, 1)
    #end=datetime.datetime(2020, 7, 31)

    #filewin="wind_gse_test.p" 
    #start=datetime.datetime(2019, 8, 1)
    #end=datetime.datetime(2019, 10, 1)


    #end=datetime.datetime.utcnow()
    #hd.save_wind_data(data_path,filewin,start,end,heeq=False)



#filewin="wind_2018_now_gse.p" 
#[wing,hwing]=pickle.load(open(data_path+filewin, "rb" ) )  

#filewin="wind_2018_2020_sept_heeq.p" 
#[winh,hwinh]=pickle.load(open(data_path+filewin, "rb" ) )  

    
#plt.plot(wing.time,wing.bx,'-k',label='gse')
#plt.plot(winh.time,winh.bx,'-g',label='heeq')
#plt.plot(sc.time,ang[:,1],'-r',label='theta')
#plt.plot(sc.time,ang[:,2],'-b',label='lambda - omega mod 360')

#plt.legend()



#check solar orbiter conversions
#from heliocats import data as hd
#importlib.reload(hd)

#print('load Solar Orbiter RTN')
#filesolo='solo_april2020.p'
#solo_rtn=pickle.load(open(data_path+filesolo, "rb" ) )   
#solo_sceq=hd.convert_RTN_to_SCEQ(solo_rtn,'SolO')
#solo_heeq=hd.convert_RTN_to_HEEQ(solo_rtn,'SolO')

#fig=plt.figure(3,figsize=(15,10),dpi=70)

#plt.plot(solo_rtn.time,solo_rtn.bx,'-k')
#plt.plot(solo_heeq.time,solo_heeq.bx,'-g')






######################## ---------------------------------Be


'''
from heliocats import data as hd
importlib.reload(hd) #reload again while debugging

print('load Bepi Colombo HEE')
filebepi='bepi_2020_march_sept_hee.p'
bepi_hee=pickle.load(open(data_path+filebepi, "rb" ) )      
#convert other coordinate systems - SCEQ for ICMECAT in 2 steps
bepi_heeq=hd.convert_HEE_to_HEEQ(bepi_hee)

filemag=data_path+'bepi_2020_march_sept_heeq.p'
pickle.dump(bepi_heeq, open(filemag, "wb"))


#final in SCEQ
bepi=hd.convert_HEEQ_to_SCEQ(bepi_heeq)
#may add RTN for Bepi in the future
filemag=data_path+'bepi_2020_march_sept_sceq.p'
pickle.dump(bepi, open(filemag, "wb"))



print('load Solar Orbiter RTN')
filesolo='solo_april2020.p'
solo_rtn=pickle.load(open(data_path+filesolo, "rb" ) ) 

filemag=data_path+'solo_2020_april_rtn.p'
pickle.dump(solo_rtn, open(filemag, "wb"))

print('convert SolO to SCEQ')
#convert other coordinate systems - SCEQ for ICMECAT
solo=hd.convert_RTN_to_SCEQ(solo_rtn,'SolO')
filemag=data_path+'solo_2020_april_sceq.p'
pickle.dump(solo, open(filemag, "wb"))


#if needed for ICME prediction from < 1 AU
solo_heeq=hd.convert_RTN_to_HEEQ(solo_rtn,'SolO')
solo_gse=hd.convert_HEEQ_to_GSE(solo_heeq)
filemag=data_path+'solo_2020_april_gse.p'
pickle.dump(solo, open(filemag, "wb"))
'''


print('done')


# ## (1) load data from HELCATS, or made with HelioSat and heliocats.data

# In[4]:


load_data=1

if load_data > 0:
            
    print('load Ulysses RTN') #made with heliocats.data.save_ulysses_data
    fileuly='ulysses_1990_2009_rtn.p'
    [uly,huly]=pickle.load(open(data_path+fileuly, "rb" ) )     
 
    print('load VEX data (Venus magnetosphere removed) SCEQ') #legacy from HELCATS project in SCEQ, removed magnetosphere
    filevex='vex_2007_2014_sceq_removed.p'
    [vex,hvex]=pickle.load(open(data_path+filevex, 'rb' ) )

    print('load MESSENGER data (Mercury magnetosphere removed) SCEQ') #legacy from HELCATS project in SCEQ, removed magnetosphere
    filemes='messenger_2007_2015_sceq_removed.p'
    [mes,hmes]=pickle.load(open(data_path+filemes, 'rb' ) )
 
    print('load STEREO-B data SCEQ') #yearly magplasma files from stereo science center, conversion to SCEQ 
    filestb='stereob_2007_2014_sceq.p'
    [stb,hstb]=pickle.load(open(data_path+filestb, "rb" ) )      
 

    ########### CURRENT ACTIVE SPACECRAFT    
    
    print('load Bepi Colombo SCEQ')
    filebepi='bepi_2020_march_sept_sceq.p'
    bepi=pickle.load(open(data_path+filebepi, "rb" ) )      
   
    
    print('load Solar Orbiter SCEQ')
    filesolo='solo_2020_april_sceq.p'
    solo=pickle.load(open(data_path+filesolo, "rb" ) ) 
   

           
    print('load MAVEN data MSO') 
    #filemav='maven_2014_2018.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )
    filemav='maven_2014_2018_removed.p'
    [mavr,hmavr]=pickle.load(open(data_path+filemav, 'rb' ) )    
    #removed magnetosphere by C. Simon Wedlund, 1 data point per orbit, MSO
    filemav='maven_2014_2018_removed_smoothed.p'
    [mav,hmav]=pickle.load(open(data_path+filemav, 'rb' ) )
    
        
    #use hd.save_msl_rad() first to convert data doseE_sol_filter_2019.dat to pickle file
    print('load MSL RAD')
    #MSL RAD
    rad=hd.load_msl_rad()#, rad.time,rad.dose_sol
    
    
    print('load PSP data SCEQ') #from heliosat, converted to SCEQ similar to STEREO-A/B
    filepsp='psp_2018_2020_sceq.p'
    [psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) ) 
    
    
    print('load and merge STEREO-A data SCEQ') #yearly magplasma files from stereo science center, conversion to SCEQ 
    filesta1='stereoa_2007_2019_sceq.p'
    [sta1,hsta1]=pickle.load(open(data_path+filesta1, "rb" ) )  
    sta1=sta1[np.where(sta1.time < parse_time('2019-Sep-01 00:00').datetime)[0]]

    #beacon data
    #filesta2="stereoa_2019_2020_sceq_beacon.p"
    filesta2='stereoa_2019_2020_sept_sceq_beacon.p'
    [sta2,hsta2]=pickle.load(open(data_path+filesta2, "rb" ) )  
    sta2=sta2[np.where(sta2.time >= parse_time('2019-Sep-01 00:00').datetime)[0]]

    #make array
    sta=np.zeros(np.size(sta1.time)+np.size(sta2.time),dtype=[('time',object),('bx', float),('by', float),                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),                ('x', float),('y', float),('z', float),                ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    sta = sta.view(np.recarray)  
    sta.time=np.hstack((sta1.time,sta2.time))
    sta.bx=np.hstack((sta1.bx,sta2.bx))
    sta.by=np.hstack((sta1.by,sta2.by))
    sta.bz=np.hstack((sta1.bz,sta2.bz))
    sta.bt=np.hstack((sta1.bt,sta2.bt))
    sta.vt=np.hstack((sta1.vt,sta2.vt))
    sta.np=np.hstack((sta1.np,sta2.np))
    sta.tp=np.hstack((sta1.tp,sta2.tp))
    sta.x=np.hstack((sta1.x,sta2.x))
    sta.y=np.hstack((sta1.y,sta2.y))
    sta.z=np.hstack((sta1.z,sta2.z))
    sta.r=np.hstack((sta1.r,sta2.r))
    sta.lon=np.hstack((sta1.lon,sta2.lon))
    sta.lat=np.hstack((sta1.lat,sta2.lat))
    print('STA Merging done')

    
    
    print('load and merge Wind data HEEQ') 
    #from HELCATS HEEQ until 2018 1 1 + new self-processed data with heliosat and hd.save_wind_data
    filewin="wind_2007_2018_heeq_helcats.p" 
    [win1,hwin1]=pickle.load(open(data_path+filewin, "rb" ) )  
    
    #or use: filewin2="wind_2018_now_heeq.p" 
    #filewin2="wind_2018_2019_heeq.p" 
    filewin2="wind_2018_2020_sept_heeq.p" 
    [win2,hwin2]=pickle.load(open(data_path+filewin2, "rb" ) )  

    #merge Wind old and new data 
    #cut off HELCATS data at end of 2017, win2 begins exactly after this
    win1=win1[np.where(win1.time < parse_time('2018-Jan-01 00:00').datetime)[0]]
    #make array
    win=np.zeros(np.size(win1.time)+np.size(win2.time),dtype=[('time',object),('bx', float),('by', float),                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),                ('x', float),('y', float),('z', float),                ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    win = win.view(np.recarray)  
    win.time=np.hstack((win1.time,win2.time))
    win.bx=np.hstack((win1.bx,win2.bx))
    win.by=np.hstack((win1.by,win2.by))
    win.bz=np.hstack((win1.bz,win2.bz))
    win.bt=np.hstack((win1.bt,win2.bt))
    win.vt=np.hstack((win1.vt,win2.vt))
    win.np=np.hstack((win1.np,win2.np))
    win.tp=np.hstack((win1.tp,win2.tp))
    win.x=np.hstack((win1.x,win2.x))
    win.y=np.hstack((win1.y,win2.y))
    win.z=np.hstack((win1.z,win2.z))
    win.r=np.hstack((win1.r,win2.r))
    win.lon=np.hstack((win1.lon,win2.lon))
    win.lat=np.hstack((win1.lat,win2.lat))

    print('Wind merging done')
    

    
#LOAD HELCATS catalogs

#HIGEOCAT
higeocat=hc.load_higeocat_vot('data/HCME_WP3_V06.vot')
higeocat_time=parse_time(higeocat['Date']).datetime    

#load arrcat as pandas dataframe
file='arrcat/HELCATS_ARRCAT_v20_pandas.p'
[ac_pandas,h]=pickle.load( open(file, 'rb')) 
    
              
print()
    
print()       
print('time ranges of the in situ data: ')    
print()
print('active spacecraft:')

print('Solar Orbiter        ',str(solo.time[0])[0:10],str(solo.time[-1])[0:10])
print('Bepi Colombo         ',str(bepi.time[0])[0:10],str(bepi.time[-1])[0:10])
print('Parker Solar Probe   ',str(psp.time[0])[0:10],str(psp.time[-1])[0:10])
print('Wind                 ',str(win.time[0])[0:10],str(win.time[-1])[0:10])
print('STEREO-A             ',str(sta.time[0])[0:10],str(sta.time[-1])[0:10])
print('MAVEN                ',str(mav.time[0])[0:10],str(mav.time[-1])[0:10])
print('MSL/RAD              ',str(rad.time[0])[0:10],str(rad.time[-1])[0:10])
print()
print('missions finished:')
print('VEX                  ',str(vex.time[0])[0:10],str(vex.time[-1])[0:10])
print('MESSENGER            ',str(mes.time[0])[0:10],str(mes.time[-1])[0:10])
print('STEREO-B             ',str(stb.time[0])[0:10],str(stb.time[-1])[0:10])
print('Ulysses              ',str(uly.time[0])[0:10],str(uly.time[-1])[0:10])
print()
print('catalogs:')
print()
print('HELCATS HIGeoCAT     ',str(higeocat_time[0])[0:10],str(higeocat_time[-1])[0:10])
print('HELCATS ARRCAT       ',np.sort(ac_pandas.sse_launch_time)[0][0:10],np.sort(ac_pandas.sse_launch_time)[-1][0:10])



print('done')


# ### 1a save data as numpy structured arrays for machine learning if needed

# In[18]:


# save data as numpy structured arrays for machine learning
data_to_numpy=0


if data_to_numpy > 0:  

    print('convert data to numpy structured arrays suitable for machine learning')
    

    print('STEREO-B')
    #STEREO-B
    #make extra header with variable attributes
    hstb_att='dtype=[(time [matplotlib format], < f8], (bt [nT], <f8), (bx, [nT, SCEQ], <f8), (by  [nT, SCEQ], <f8),    (bz, SCEQ  [nT], <f8), (vt  [km/s], <f8), (np [ccm -3], <f8), (tp [K], <f8), (x [AU, HEEQ], <f8), (y [AU, HEEQ], <f8),    (z [AU, HEEQ], <f8), (r, <f8), (lat [deg, HEEQ], <f8), (lon [deg, HEEQ], <f8 )]'
    stb_nd=hd.recarray_to_numpy_array(stb)
    #change dtype everywhere to float
    stb_nd=stb_nd.astype(dtype=[('time', '<f8'), ('bx', '<f8'), ('by', '<f8'), ('bz', '<f8'), ('bt', '<f8'),                                 ('vt', '<f8'), ('np', '<f8'), ('tp', '<f8'), ('x', '<f8'), ('y', '<f8'),                                 ('z', '<f8'), ('r', '<f8'), ('lat', '<f8'), ('lon', '<f8')])
    pickle.dump([stb_nd,hstb_att], open(data_path_ML+ "stereob_2007_2014_sceq_ndarray.p", "wb" ) )
    
     
    
    print('Ulysses')
    huly_att='dtype=[(time [matplotlib format]), (bt [nT], <f8), (bx  [nT, RTN], <f8), (by  [nT, RTN], <f8),     (bz  [nT, RTN], <f8), (vt  [km/s], <f8), (np [ccm -3], <f8), (tp [K], <f8), (x [AU, HEEQ], <f8),    (y [AU, HEEQ], <f8), (z [AU, HEEQ], <f8), (r [AU], <f8), (lat [deg, HEEQ], <f8), (lon [deg, HEEQ], <f8 )]'
    uly_nd=hd.recarray_to_numpy_array(uly)
    uly_nd=uly_nd.astype(dtype=[('time', '<f8'), ('bx', '<f8'), ('by', '<f8'), ('bz', '<f8'), ('bt', '<f8'),                                 ('vt', '<f8'), ('np', '<f8'), ('tp', '<f8'), ('x', '<f8'), ('y', '<f8'),                                 ('z', '<f8'), ('r', '<f8'), ('lat', '<f8'), ('lon', '<f8')])
    pickle.dump([uly_nd,huly_att], open(data_path_ML+ "ulysses_1990_2009_rtn_ndarray.p", "wb" ) )
       
    

    print('VEX')
    #header no plasma, with planetary orbit
    hatt6='dtype=[(time [matplotlib format]), (bt [nT], <f8), (bx  [nT, SCEQ], <f8), (by  [nT, SCEQ], <f8),     (bz  [nT, SCEQ], <f8), (x [AU, HEEQ], <f8), (y [AU, HEEQ], <f8), (z [AU, HEEQ], <f8), (r [AU], <f8),     (lat [deg, HEEQ], <f8), (lon [deg, HEEQ]), (xo [km, VSO], <f8), (yo [km, VSO], <f8), (zo [km, VSO], <f8),     (ro [km], <f8), (lato [deg, VSO], <f8), (lono [deg, VSO], <f8)]'
    vex_nd=hd.recarray_to_numpy_array(vex)
    vex_nd=vex_nd.astype(dtype=[('time', '<f8'), ('bx', '<f8'), ('by', '<f8'), ('bz', '<f8'), ('bt', '<f8'),                                ('x', '<f8'), ('y', '<f8'),  ('z', '<f8'), ('r', '<f8'), ('lat', '<f8'), ('lon', '<f8'),                                ('xo', '<f8'), ('yo', '<f8'),  ('zo', '<f8'), ('ro', '<f8'), ('lato', '<f8'), ('lono', '<f8')])    
    pickle.dump([vex_nd,hatt6], open(data_path_ML+ "vex_2007_2014_sceq_removed_ndarray.p", "wb" ) )


    print('MESSENGER') #no plasma, no planetary orbit
    hatt7='dtype=[(time [matplotlib format]), (bt [nT], <f8), (bx  [nT, SCEQ], <f8), (by  [nT, SCEQ], <f8),     (bz  [nT, SCEQ], <f8), (x [AU, HEEQ], <f8), (y [AU, HEEQ], <f8), (z [AU, HEEQ], <f8), (r [AU], <f8),     (lat [deg, HEEQ], <f8), (lon [deg, HEEQ])'
    mes_nd=hd.recarray_to_numpy_array(mes)
    mes_nd=mes_nd.astype(dtype=[('time', '<f8'), ('bx', '<f8'), ('by', '<f8'), ('bz', '<f8'), ('bt', '<f8'),                                ('x', '<f8'), ('y', '<f8'),  ('z', '<f8'), ('r', '<f8'), ('lat', '<f8'), ('lon', '<f8')])    
    pickle.dump([mes_nd,hatt7], open(data_path_ML+ "mes_2007_2015_sceq_removed_ndarray.p", "wb" ) ) 

    
    
    
    ###################################### UPDATED
    
         
    
    print('STEREO-A')
    hsta_att='dtype=[(time [matplotlib format], < f8], (bt [nT], <f8), (bx, [nT, SCEQ], <f8), (by  [nT, SCEQ], <f8),    (bz, SCEQ  [nT], <f8), (vt  [km/s], <f8), (np [ccm -3], <f8), (tp [K], <f8), (x [AU, HEEQ], <f8), (y [AU, HEEQ], <f8),    (z [AU, HEEQ], <f8), (r, <f8), (lat [deg, HEEQ], <f8), (lon [deg, HEEQ], <f8 )]'
    sta_nd=hd.recarray_to_numpy_array(sta)
    sta_nd=sta_nd.astype(dtype=[('time', '<f8'), ('bx', '<f8'), ('by', '<f8'), ('bz', '<f8'), ('bt', '<f8'),                                 ('vt', '<f8'), ('np', '<f8'), ('tp', '<f8'), ('x', '<f8'), ('y', '<f8'),                                 ('z', '<f8'), ('r', '<f8'), ('lat', '<f8'), ('lon', '<f8')])
    pickle.dump([sta_nd,hsta_att], open(data_path_ML+ "stereoa_2007_2020_sceq_ndarray.p", "wb" ) )
    
    
  
    print('Wind')
    hwind_att='dtype=[((time [matplotlib format]), (bt [nT], <f8), (bx  [nT, HEEQ], <f8), (by  [nT, HEEQ], <f8),                         (bz  [nT, HEEQ], <f8), (vt  [km/s], <f8), (np [ccm -3], <f8),                         (tp [K], <f8), (x [AU, HEEQ], <f8), (y [AU, HEEQ], <f8), (z [AU, HEEQ], <f8), (r [AU], <f8),                        (lat [deg, HEEQ], <f8), (lon [deg, HEEQ],<f8 )]'
    win_nd=hd.recarray_to_numpy_array(win)
    win_nd=win_nd.astype(dtype=[('time', '<f8'), ('bx', '<f8'), ('by', '<f8'), ('bz', '<f8'), ('bt', '<f8'),                                 ('vt', '<f8'), ('np', '<f8'), ('tp', '<f8'), ('x', '<f8'), ('y', '<f8'),                                 ('z', '<f8'), ('r', '<f8'), ('lat', '<f8'), ('lon', '<f8')])
    pickle.dump([win_nd,hwind_att], open(data_path_ML+ "wind_2007_2020_heeq_ndarray.p", "wb" ) )

   
    
    print('PSP')  
    hpsp_att='dtype=[(time, matplotlib), (bt [nT], <f8), (bx, SCEQ  [nT], <f8), (by  [nT, SCEQ], <f8),                        (bz, SCEQ  [nT], <f8), (vt  [km/s], <f8),(vx  [km/s, SCEQ], <f8)(vy  [km/s, SCEQ], <f8),                        (vz  [km/s, SCEQ], <f8), (np [ccm -3], <f8), (tp [K], <f8), (x [AU, HEEQ], <f8), (y [AU, HEEQ], <f8),                        (z [AU, HEEQ], <f8), (r, <f8), (lat [deg, HEEQ], <f8), (lon [deg, HEEQ], <f8])'

    psp_nd=hd.recarray_to_numpy_array(psp)
    psp_nd=psp_nd.astype(dtype=[('time', '<f8'), ('bx', '<f8'), ('by', '<f8'), ('bz', '<f8'), ('bt', '<f8'),                                 ('vt', '<f8'),('vx', '<f8'),('vy', '<f8'),('vz', '<f8'), ('np', '<f8'), ('tp', '<f8'),                                ('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('r', '<f8'), ('lat', '<f8'), ('lon', '<f8')])
    pickle.dump([psp_nd,hpsp_att], open(data_path_ML+ "psp_2018_2020_sceq_ndarray.p", "wb" ) ) 

 

    print('done')


# ## (2) measure new events 

# In[5]:


#for measuring new events use these functions from heliocats.plot 

from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging

from heliocats import data as hd
importlib.reload(hd) #reload again while debugging

plt.close('all')
#works in jupyter notebooks

#works in scripts
#matplotlib.use('qt5agg')  
#%matplotlib 
#plt.ion()


#%matplotlib

###################### make huge overview plot after Aug 2018


#Solar Orbiter
##hp.plot_insitu_measure(solo, '2020-Apr-15','2020-Apr-29', 'SolO', 'results/plots_icmecat/')


#Bepi Colombo
##hp.plot_insitu_measure(bepi, '2018-Oct-30','2020-Apr-29', 'PSP', 'results/plots_icmecat/')



#PSP
#hp.plot_insitu_measure(psp, '2018-Oct-30','2020-Apr-29', 'PSP', 'results/plots_icmecat/')

#or use plot_insitu_measure_mag

#for plotting single events
#hp.plot_insitu(psp, ic.icme,'2018-Nov-15', 'PSP', icplotsdir)

#STEREO-A
#hp.plot_insitu_measure(sta, '2018-Jan-01 12:00','2018-Feb-01 12:00', 'STEREO-A', 'results/')

#Wind
#hp.plot_insitu_measure(win, '2019-Jan-29','2020-Aug-28', 'Wind', 'results/')

#hp.plot_insitu_measure(win, '2011-Feb-10','2011-Feb-25', 'Wind', 'results/')



#STEREO-A
#hp.plot_insitu_measure(sta, '2019-Nov-01 12:00','2020-Aug-01 12:00', 'STEREO-A', 'results/')



################### MAVEN needs special attention

# ########### LOAD other data sets

# # 1 load RAD data from Jingnan Guo
# rad=hd.load_msl_rad()   

# # 2 load HELCATS ARRCAT as numpy recarray
# file='arrcat/HELCATS_ARRCAT_v20_pandas.p'
# [arrcat,arrcat_header]=pickle.load( open(file, 'rb'))   

# # 3 load WSA speed file at Mars from Martin Reiss 
# #save wsa hux file if new
# #file='wsa_hux_mars_aug2014_jan2018.txt'
# #hd.save_wsa_hux(file)
# #get mars wsa hux #Reiss et al., 2019 (DOI: https://doi.org/10.3847/1538-4365/aaf8b3
# w1=hd.load_mars_wsa_hux()

# #load maven SIR catalog Huang et al. 2019 https://doi.org/10.3847/1538-4357/ab25e9
# msir=hd.load_maven_sir_huang()

# #mavr are the removed but not smoothed data

# #use for extra window
# %matplotlib


# #plot for measuring overview
# plt.figure(11)
# plt.rcParams["figure.figsize"] = (16,8)
# plt.plot_date(w1.time,w1.vt,'-k')
# plt.plot_date(rad.time,rad.dose_sol,'-b')
# plt.plot_date(mav.time,mav.vt,'-r')
# #plt.plot_date(mavr.time,mavr.vt,'-r')

# #https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JA023495%4010.1002/%28ISSN%292169-9402.MAVEN2

# #zoomed plot
# hp.plot_insitu_measure_maven(mav, '2016-Jan-1','2016-Feb-01', 'MAVEN', icplotsdir,rad,arrcat,w1,msir)


# ## (3) make ICMECAT 

# In[37]:


print('data loaded')


#public
ic=hc.load_helcats_icmecat_master_from_excel('icmecat/HELCATS_ICMECAT_v20_master.xlsx')


#non public catalog
#ic2=hc.load_helcats_icmecat_master_from_excel('icmecat/HELCATS_ICMECAT_v20_master_internal.xlsx')



####### 3a get indices for all spacecraft

wini=np.where(ic.sc_insitu == 'Wind')[:][0] 
stai=np.where(ic.sc_insitu == 'STEREO-A')[:][0]    
stbi=np.where(ic.sc_insitu == 'STEREO-B')[:][0]    
vexi=np.where(ic.sc_insitu == 'VEX')[:][0]  
mesi=np.where(ic.sc_insitu == 'MESSENGER')[:][0]   
ulyi=np.where(ic.sc_insitu == 'ULYSSES')[:][0]    
mavi=np.where(ic.sc_insitu == 'MAVEN')[:][0]    
pspi=np.where(ic.sc_insitu == 'PSP')[:][0]    

####### 3b get parameters for all spacecraft one after another


#os.system('rm icmecat/indices_icmecat/ICMECAT_indices_Wind.p')

#misssions to be updated
ic=hc.get_cat_parameters(psp,pspi,ic,'PSP')
ic=hc.get_cat_parameters(win,wini,ic,'Wind')
ic=hc.get_cat_parameters(sta,stai,ic,'STEREO-A')
ic=hc.get_cat_parameters(mav,mavi,ic,'MAVEN')

#finished missions
ic=hc.get_cat_parameters(stb,stbi,ic,'STEREO-B')
ic=hc.get_cat_parameters(vex,vexi,ic,'VEX')
ic=hc.get_cat_parameters(mes,mesi,ic,'MESSENGER')
ic=hc.get_cat_parameters(uly,ulyi,ic,'ULYSSES')


# ###### 3c make all plots if wanted
matplotlib.use('Agg')

#from heliocats import plot as hp
#importlib.reload(hp) #reload again while debugging

#missions to be updated
#hp.plot_icmecat_events(sta,stai,ic,'STEREO-A',icplotsdir)
#hp.plot_icmecat_events(psp,pspi,ic,'PSP',icplotsdir)
#hp.plot_icmecat_events(win,wini,ic,'Wind',icplotsdir)
#hp.plot_icmecat_events(mav,mavi,ic,'MAVEN',icplotsdir)

#finished missions
# hp.plot_icmecat_events(stb,stbi,ic,'STEREO-B',icplotsdir)
# hp.plot_icmecat_events(vex,vexi,ic,'VEX',icplotsdir)
# hp.plot_icmecat_events(mes,mesi,ic,'MESSENGER',icplotsdir)
# hp.plot_icmecat_events(uly,ulyi,ic,'ULYSSES',icplotsdir)
print('done')

get_ipython().run_line_magic('matplotlib', 'inline')


# ### (4) save ICMECAT 

# ### 4a save header

# In[39]:


#save header and parameters as text file and prepare for html website
header='ICME CATALOGUE v2.0 \n\nThis is the HELCATS interplanetary coronal mass ejection (ICME) catalog, based on in situ magnetic field and bulk plasma observations in the heliosphere. \n\nThis is version 2.0, released 2020-06-03, updated 2020-09-28, doi: 10.6084/m9.figshare.6356420 \n\nThe catalog is available as a python pandas dataframe (pickle), python numpy structured array (pickle), json, csv, xlsx, txt, hdf5, at \nhttps://helioforecast.space/icmecat and https://www.helcats-fp7.eu/catalogues/wp4_icmecat.html \n\nNumber of events in ICMECAT: '+str(len(ic))+' \nICME observatories: Parker Solar Probe (PSP), Wind, STEREO-A, MAVEN, STEREO-B, Venus Express (VEX), MESSENGER, Ulysses.   \nTime range: January 2007 - September 2020. \n \nAuthors: Christian Moestl, Andreas J. Weiss, Rachel L. Bailey, Martin A. Reiss, Space Research Institute, Austrian Academy of Sciences, Graz, Austria.\nContributors: Peter Boakes, Alexey Isavnin, Emilia Kilpua, David Stansby, Reka Winslow, Brian Anderson, Lydia Philpott, Vratislav Krupar, Jonathan Eastwood, Simon Good, Lan Jian, Teresa Nieves-Chinchilla, Cyril Simon Wedlund, Jingnan Guo, Johan von Forstner, Mateja Dumbovic, Benoit Lavraud.  \n\nRules: If results are produced with this catalog for peer-reviewed scientific publications, please contact christian.moestl@oeaw.ac.at for possible co-authorship. \n\nThis catalog has been made by getting the 3 times of each ICME (shock or disturbance begin, magnetic obstacle start and end) from the individual catalogs below, and then calculating all parameters again consistently from the data by us. \nThe in situ data that were used for the catalog, with a size of 8 GB in total, including extra data files with magnetic field components in RTN coordinates that are not used for producing the catalog, can be downloaded in python pickle format as recarrays from https://doi.org/10.6084/m9.figshare.11973693 \nThe python code for producing this catalog is available as part of https://github.com/cmoestl/heliocats \n\nEach event has unique identifier - the icmecat_id - which has a tag in it that indicates from which catalog the ICME times were taken: \n\nWind:       Nieves-Chinchilla et al. (2018), tag: NASA. \nSTEREO-A:   Jian et al. (2018), tag: JIAN. \nSTEREO-B:   Jian et al. (2018), tag: JIAN. \nVEX:        Good et al. (2018), tag: SGOOD \nMESSENGER:  Good et al. (2018), Winslow et al. (2018), tags: SGOOD, WINSLOW. \nMAVEN:      Made by us according to the method in the comments, tag: MOESTL.\nUlysses:    Added by us, tag: MOESTL. \nPSP:        Added by us, tag: MOESTL. \n\nWe have also added extra events at VEX, MESSENGER, Wind and STEREO-A (all tagged with MOESTL in icmecat_id).\n\nReferences: \nNieves-Chinchilla, T. et al. (2018),  https://doi.org/10.1007/s11207-018-1247-z \n                                      https://wind.nasa.gov/ICME_catalog/ICME_catalog_viewer.php \nJian, L. et al. (2018), https://doi.org/10.3847/1538-4357/aab189 \n                        https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level3/ \nGood, S. et al. (2018) https://doi.org/10.1007/s11207-015-0828-3 \nWinslow, R. et al. (2015), https://doi.org/10.1002/2015JA021200 \n\n\nComments: \n\n- Spacecraft positions are given in Heliocentric Earth Equatorial Coordinates (HEEQ) coordinates. \n\n- The coordinate system for all magnetic field components is SCEQ, except for Wind (HEEQ, which is the equivalent for SCEQ for Earth) and Ulysses (RTN, because of high latitude positions) and MAVEN (MSO). \n        Definition of SpaceCraft Equatorial Coordinates (SCEQ): \n        Z is the solar rotation axis. \n        Y is the cross product of Z and R, with R being the vector that points from the Sun to the spacecraft.\n        X completes the right handed triad (and points away from the Sun). \nThis system is thus like HEEQ but centered on the respective in situ spacecraft, so the SCEQ X and Y \nbase vectors are rotated by the HEEQ longitude of the in situ spacecraft from HEEQ X and Y.\nThe Y vector is similar to the T vector in an RTN system for each spacecraft, but the X and Z vectors \nare rotated around Y compared to an RTN system. The differences between RTN and SCEQ for spacecraft within \na few degrees of the solar equatorial plane are very small (within a few 0.1 nT usually).\nWe choose SCEQ for the advantage that a comparison of the data of multipoint CME events \nand comparisons to simulations are simpler because there is always a similar reference plane (the solar equatorial plane).\n\n- Venus Express and MESSENGER do not have plasma parameters available. \n- If there is no sheath or density pileup region, so the ICME starts immediately with a magnetic obstacle, the icme_start_time is similar to mo_start_time.\n\n- At MESSENGER and VEX, for events cataloged by Simon Good, icme_start_time has been added by V. Krupar (Imperial College) and C. Möstl (IWF Graz). \n\n- For the calculation of the parameters at MESSENGER during the orbit around Mercury, all data points inside the bowshock of Mercury have been removed, according to a list thankfully provided to us by by R. Winslow, UNH, B. Anderson, APL, and Lydia Philpott, UBC. \n\n- Calculation of the magnetic obstacle parameters at VEX is done after approximate removal of the induced magnetosphere, with a modified equation as in Zhang et al. 2008 (doi: 10.1016/j.pss.2007.09.012), with a constant of 3.5 instead of 2.14/2.364,in order to account for a larger bowshock distance during solar maximum than studied in the Zhang et al. (2008) paper. \n\n- For MAVEN, all data inside the bow shock were removed with the model from Gruesbeck et al. (2018, doi:10.1029/2018JA025366) by C. Simon Wedlund (IWF Graz, Austria). From the remaining data, the median for each orbit is taken as 1 data point, resulting in a solar wind dataset at Mars with 4.5 hour time resolution. This is a far lower resolution than available at 1 AU, so the identification of ICMEs in the MAVEN data is not as straightforward as for data sets with higher time resolution.\n\n- The identification of ICMEs for MAVEN was done as follows: (1) We looked for typical profiles for ICMEs in magnetic field (B) and speed (V). (2) B needs to be elevated in ICMEs accompanied by a flat or declining profile in V. An elevation in B accompanied by a later gradual rise in V is a strong indicator for a high speed stream (HSS). (3) The discrimination between HSS and ICMEs was further strengthened with the help of a stream interaction region list for MAVEN by Huang et al. (2019, doi: 10.3847/1538-4357/ab25e9). (4) We additionally plotted the predicted CME arrivals with STEREO/HI (ARRCATv2.0), the dose rate measured on the Mars surface by MSL/RAD (e.g. Guo et al. (2018, doi: 10.1051/0004-6361/201732087)) and the speed of the WSA/HUX model for the background solar wind (Reiss et al. 2019, doi: 10.3847/1538-4365/aaf8b3). General guidance on space weather observed by MAVEN can be found in Lee et al. (2018, doi:10.1002/2016JA023495).\n\n\n\n'


parameters='Parameters:\n00: icmecat_id: The unique identifier for the observed ICME. unit: string. \n01: sc insitu: The name of the in situ observing spacecraft. unit: string. \n02: icme_start_time: Shock arrival or density enhancement time, can be similar to mo_start_time. unit: UTC. \n03: mo_start_time: Start time of the magnetic obstacle (MO), including flux ropes, flux-rope-like, and ejecta signatures. unit: UTC. \n04: mo_end_time: End time of the magnetic obstacle. unit: UTC. \n05: mo_sc_heliodistance: Heliocentric distance of the spacecraft at mo_start_time. unit: AU.\n06: mo_sc_long_heeq: Heliospheric longitude of the spacecraft at mo_start_time, range [-180,180]. unit: degree (HEEQ).\n07: mo_sc_lat_heeq: Heliospheric latitude of the spacecraft at mo_start_time, range [-90,90]. unit: degree (HEEQ).\n08: icme_duration: Duration of the interval between icme_start_time and mo_endtime. unit: hours.\n09: icme_bmax: Maximum total magnetic field in the full icme interval (icme_start_time to mo_end_time). unit: nT.\n10: icme_bmean: Mean total magnetic field during the full icme interval (icme_start_time to mo_end_time). unit: nT.\n11: icme_bstd: Standard deviation of the total magnetic field from icme_start_time to mo_end_time. unit: nT.\n12: icme_speed_mean: Mean proton speed from icme_start_time to mo_end_time. unit: km/s.\n13: icme_speed_std: Standard deviation of proton speed from icme_start_time to mo_end_time. unit: km/s.\n14: mo_duration: Duration of the interval between mo_start_time and mo_endtime. unit: hours.\n15: mo_bmax: Maximum total magnetic field in the magnetic obstacle interval (mo_start_time to mo_end_time). unit: nT.\n16: mo_bmean: Mean total magnetic field in the magnetic obstacle. unit: nT.\n17: mo_bstd: Standard deviation of the total magnetic field in the magnetic obstacle. unit: nT.\n18: mo_bzmean: Mean magnetic field Bz component in the magnetic obstacle. unit: nT.\n19: mo_bzmin: Minimum magnetic field Bz component in the magnetic obstacle. unit: nT.\n20: mo_bzstd: Standard deviation of the magnetic field Bz component in the magnetic obstacle. unit: nT.\n21: mo_bymean: Mean magnetic field By component in the magnetic obstacle. unit: nT.\n22: mo_bystd: Standard deviation of the magnetic field By component in the magnetic obstacle. unit: nT.\n23: mo_speed_mean: Mean proton speed from mo_start_time to mo_end_time. unit: km/s.\n24: mo_speed_std: Standard deviation of proton speed from mo_start_time to mo_end_time. unit: km/s.\n25: mo_expansion_speed: Difference between proton speed at mo_start_time to proton speed at mo_end_time. unit: km/s.\n26: mo_pdyn_mean: Mean proton dynamic pressure from mo_start_time to mo_start_time. unit: nPa.\n27: mo_pdyn_std: Standard deviation of proton dynamic pressure from mo_start_time to mo_start_time. unit: nPa.\n28: mo_density_mean: Mean proton density from mo_start_time to mo_start_time. unit: cm^-3.\n29: mo_density_std: Standard deviation of proton density from mo_start_time to mo_start_time. unit: cm^-3.\n30: mo_temperature_mean: Mean proton temperature from mo_start_time to mo_start_time. unit: K.\n31: mo_temperature_std: Standard deviation of proton temperature from mo_start_time to mo_end_time. unit: K.\n32: sheath_speed_mean: Mean proton speed from icme_start_time to mo_start_time, NaN if these times are similar. unit: km/s.\n33: sheath_speed_std: Standard deviation of proton speed from icme_start_time to mo_start_time, NaN if these times are similar. unit: km/s.\n34: sheath_density_mean: Mean proton density from icme_start_time to mo_start_time, NaN if these times are similar. unit: cm^-3.\n35: sheath_density_std: Standard deviation of proton density from icme_start_time to mo_start_time, NaN if these times are similar. unit: cm^-3.\n36: sheath_pdyn_mean: Mean proton dynamic pressure, from icme_start_time to mo_start_time, NaN if these times are similar. unit: nPa.\n37: sheath_pdyn_std: Standard deviation of proton dynamic pressure, from icme_start_time to mo_start_time, NaN if these times are similar. unit: nPa.\n\n\n'



print(header)
print(parameters)


#make header file
file='icmecat/HELCATS_ICMECAT_v20_header.txt'
with open(file, "w") as text_file:
    text_file.write(header)
    text_file.write(parameters)
print()    
print('header saved as '+file)
print()    

#Convert to html regarding line breaks, paragraph beginning and spaces
header_spaces=header.replace(" ", "&nbsp;")
header_html= "<p>" +header_spaces.replace('\n', '<br>')+ "</p>" 
parameters_spaces=parameters.replace(" ", "&nbsp;")
parameters_html= "<p>" +parameters.replace('\n', '<br>')+ "</p>"
print('header converted to HTML')
print()    
print()    


# ### 4b save into different formats

# In[40]:


########## python formats

# save ICMECAT as pandas dataframe with times as datetime objects as pickle
file='icmecat/HELCATS_ICMECAT_v20_pandas.p'
pickle.dump([ic,header,parameters], open(file, 'wb'))
print('ICMECAT saved as '+file)


# save ICMECAT as numpy array with times as matplotlib datetime as pickle
ic_num=copy.deepcopy(ic) 
ic_num.icme_start_time=parse_time(ic_num.icme_start_time).plot_date
ic_num.mo_start_time=parse_time(ic_num.mo_start_time).plot_date
ic_num.mo_end_time=parse_time(ic_num.mo_end_time).plot_date
#convert to recarray
ic_num_rec=ic_num.to_records()
#create structured array
dtype1=[('index','i8'),('icmecat_id', '<U30'),('sc_insitu', '<U20')] +[(i, '<f8') for i in ic.keys()[2:len(ic.keys())]]
ic_num_struct=np.array(ic_num_rec,dtype=dtype1)



file='icmecat/HELCATS_ICMECAT_v20_numpy.p'
pickle.dump([ic_num,ic_num_struct,header,parameters], open(file, 'wb'))
print('ICMECAT saved as '+file)





################ save to different formats

#copy pandas dataframe first to change time format consistent with HELCATS
ic_copy=copy.deepcopy(ic)  
ic_copy.icme_start_time=parse_time(ic.icme_start_time).isot
ic_copy.mo_start_time=parse_time(ic.mo_start_time).isot
ic_copy.mo_end_time=parse_time(ic.mo_end_time).isot

#change time format
for i in np.arange(len(ic)):

    dum=ic_copy.icme_start_time[i] 
    ic_copy.at[i,'icme_start_time']=dum[0:16]+'Z'
     
    dum=ic_copy.mo_start_time[i] 
    ic_copy.at[i,'mo_start_time']=dum[0:16]+'Z'
     
    dum=ic_copy.mo_end_time[i] 
    ic_copy.at[i,'mo_end_time']=dum[0:16]+'Z'


#save as Excel
file='icmecat/HELCATS_ICMECAT_v20.xlsx'
ic_copy.to_excel(file,sheet_name='ICMECATv2.0')
print('ICMECAT saved as '+file)

#save as json
file='icmecat/HELCATS_ICMECAT_v20.json'
ic_copy.to_json(file)
print('ICMECAT saved as '+file)

#save as csv
file='icmecat/HELCATS_ICMECAT_v20.csv'
ic_copy.to_csv(file)
print('ICMECAT saved as '+file)


#save as txt
file='icmecat/HELCATS_ICMECAT_v20.txt'
np.savetxt(file, ic_copy.values.astype(str), fmt='%s' )
print('ICMECAT saved as '+file)








#########################


#########save into hdf5 format , use S for strings http://docs.h5py.org/en/stable/strings.html#what-about-numpy-s-u-type
dtype2=[('index','i8'),('icmecat_id', 'S30'),('sc_insitu', 'S20')] +[(i, '<f8') for i in ic.keys()[2:len(ic.keys())]]
ich5=np.array(ic_num_rec,dtype=dtype2)
file='icmecat/HELCATS_ICMECAT_v20.h5'
f=h5py.File(file,mode='w')
f["icmecat"]= ich5
#add attributes
#************************
#***********************

print('ICMECAT saved as '+file)
f.close()

#reading h5py files http://docs.h5py.org/en/latest/quick.html
#fr = h5py.File('icmecat/HELCATS_ICMECAT_v20.h5', 'r')
#list(fr.keys())
#ich5=fr['icmecat']
#ich5['mo_bstd']
#ich5.dtype
#fr.close()
##################


#save as .npy without pickle
file='icmecat/HELCATS_ICMECAT_v20_numpy.npy'
np.save(file,ich5, allow_pickle=False)
print('ICMECAT saved as '+file)

#for loading do:
#icnpy=np.load(file)
#decode strings:
#icnpy['icmecat_id'][0].decode()






############ other formats

#copy pandas dataframe first to change time format 
ic_copy2=copy.deepcopy(ic)  
ic_copy2.icme_start_time=parse_time(ic.icme_start_time).iso
ic_copy2.mo_start_time=parse_time(ic.mo_start_time).iso
ic_copy2.mo_end_time=parse_time(ic.mo_end_time).iso

#change time format
for i in np.arange(len(ic)):

    dum=ic_copy2.icme_start_time[i] 
    ic_copy2.at[i,'icme_start_time']=dum[0:16]
     
    dum=ic_copy2.mo_start_time[i] 
    ic_copy2.at[i,'mo_start_time']=dum[0:16]
     
    dum=ic_copy2.mo_end_time[i] 
    ic_copy2.at[i,'mo_end_time']=dum[0:16]


#save as json for webpage with different time format
file='icmecat/HELCATS_ICMECAT_v20_isot.json'
ic_copy2.to_json(file)
print('ICMECAT saved as '+file)


#save as html no header
file='icmecat/HELCATS_ICMECAT_v20_simple.html'
ic_copy.to_html(file)
print('ICMECAT saved as '+file)


############ save as html file with header
#save as html
file='icmecat/HELCATS_ICMECAT_v20.html'
#ic.to_html(file,justify='center')

#ichtml='{% extends "_base.html" %} \n \n {% block content %} \n \n \n '
ichtml = header_html
ichtml += parameters_html
ichtml += ic_copy.to_html()
#ichtml +='\n \n {% endblock %}'


with open(file,'w') as f:
    f.write(ichtml)
    f.close()
    
print('ICMECAT saved as '+file)    


# ## 4c load ICMECAT pickle files

# In[9]:


#load icmecat as pandas dataframe
file='icmecat/HELCATS_ICMECAT_v20_pandas.p'
[ic_pandas,h,p]=pickle.load( open(file, 'rb'))   

#load icmecat as numpy array
file='icmecat/HELCATS_ICMECAT_v20_numpy.p'
[ic_nprec,ic_np,h,p]=pickle.load( open(file, 'rb'))   


# In[10]:


print(ic_pandas.keys())


# In[11]:


ic_nprec


# In[12]:


ic_nprec.icmecat_id


# In[13]:


ic=ic_pandas

ic_mo_start_time_num=parse_time(ic.mo_start_time).plot_date

#get indices for each target
imes=np.where(ic.sc_insitu=='MESSENGER')[0]
ivex=np.where(ic.sc_insitu=='VEX')[0]
iwin=np.where(ic.sc_insitu=='Wind')[0]
imav=np.where(ic.sc_insitu=='MAVEN')[0]

ista=np.where(ic.sc_insitu=='STEREO-A')[0]
istb=np.where(ic.sc_insitu=='STEREO-B')[0]
ipsp=np.where(ic.sc_insitu=='PSP')[0]
#soloi=np.where(ac.target_name=='SolarOrbiter')[0]
#bepii=np.where(ac.target_name=='BepiColombo')[0]
iuly=np.where(ic.sc_insitu=='Ulysses')[0]




sns.set_context("talk")     
sns.set_style('darkgrid')

###############################################################################
fig=plt.figure(3,figsize=(20,14),dpi=200)



##############################################################################
ax2 = plt.subplot(221,projection='polar') 

#markersize
ms=25
#alpha
al=0.7


plt.title('ICMECAT events [HEEQ longitude, AU]')

ax2.scatter(np.radians(ic.mo_sc_long_heeq[imes]),ic.mo_sc_heliodistance	[imes],s=ms,c='dimgrey', alpha=al)
ax2.scatter(np.radians(ic.mo_sc_long_heeq[ivex]),ic.mo_sc_heliodistance	[ivex],s=ms,c='orange', alpha=al)
ax2.scatter(np.radians(ic.mo_sc_long_heeq[iwin]),ic.mo_sc_heliodistance	[iwin],s=ms,c='mediumseagreen', alpha=al)
ax2.scatter(np.radians(ic.mo_sc_long_heeq[imav]),ic.mo_sc_heliodistance	[imav],s=ms,c='orangered', alpha=al)

ax2.scatter(np.radians(ic.mo_sc_long_heeq[ista]),ic.mo_sc_heliodistance[ista],s=ms,c='red', alpha=al)
ax2.scatter(np.radians(ic.mo_sc_long_heeq[istb]),ic.mo_sc_heliodistance[istb],s=ms,c='blue', alpha=al)

ax2.scatter(np.radians(ic.mo_sc_long_heeq[ipsp]),ic.mo_sc_heliodistance[ipsp],s=ms,c='black', alpha=al)
#ax.scatter(np.radians(ac.target_heeq_lon[soloi]),ac.target_distance[soloi],s=ms,c='green', alpha=al)
#ax.scatter(np.radians(ac.target_heeq_lon[bepii]),ac.target_distance[bepii],s=ms,c='violet', alpha=al)
ax2.scatter(np.radians(ic.mo_sc_long_heeq[iuly]),ic.mo_sc_heliodistance[iuly],s=ms,c='brown', alpha=al)
ax2.set_ylim([0,1.8])
ax2.tick_params(labelsize=12,zorder=3)



#########################################################################
ax1=plt.subplot(222)
plt.title('ICMECAT event times')

#markersize
ms=6
#alpha
al=0.7

ax1.plot_date(ic_mo_start_time_num[imes],ic.mo_sc_heliodistance[imes],'ok',c='dimgrey', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[ivex],ic.mo_sc_heliodistance[ivex],'ok',c='orange', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[iwin],ic.mo_sc_heliodistance[iwin],'ok',c='mediumseagreen', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[imav],ic.mo_sc_heliodistance[imav],'ok',c='orangered', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[ista],ic.mo_sc_heliodistance[ista],'ok',c='red', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[istb],ic.mo_sc_heliodistance[istb],'ok',c='blue', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[ipsp],ic.mo_sc_heliodistance[ipsp],'ok',c='black', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[iuly],ic.mo_sc_heliodistance[iuly],'ok',c='brown', alpha=al,ms=ms)


ax1.set_ylabel('Heliocentric distance [AU]')
ax1.set_xlabel('year')
ax1.set_ylim([0,1.8])



##############################################################################
ax3=plt.subplot(223)
plt.title('ICMECAT mean magnetic field in the magnetic obstacle')
ax3.set_xlabel('Heliocentric distance [AU]')
ax3.set_ylabel('<B> [nT]')


ax3.plot(ic.mo_sc_heliodistance[imes],ic.mo_bmean[imes],'o',c='dimgrey', alpha=al,ms=ms,label='MESSENGER')
ax3.plot(ic.mo_sc_heliodistance[ivex],ic.mo_bmean[ivex],'o',c='orange', alpha=al,ms=ms,label='Venus Express')
ax3.plot(ic.mo_sc_heliodistance[iwin],ic.mo_bmean[iwin],'o',c='mediumseagreen', alpha=al,ms=ms,label='Wind at L1')
ax3.plot(ic.mo_sc_heliodistance[imav],ic.mo_bmean[imav],'o',c='orangered', alpha=al,ms=ms, label='MAVEN')
ax3.plot(ic.mo_sc_heliodistance[ista],ic.mo_bmean[ista],'o',c='red', alpha=al,ms=ms, label='STEREO-A')
ax3.plot(ic.mo_sc_heliodistance[istb],ic.mo_bmean[istb],'o',c='blue', alpha=al,ms=ms, label='STEREO-B')
ax3.plot(ic.mo_sc_heliodistance[ipsp],ic.mo_bmean[ipsp],'o',c='black', alpha=al,ms=ms, label='Parker Solar Probe')
ax3.plot(ic.mo_sc_heliodistance[iuly],ic.mo_bmean[iuly],'o',c='brown', alpha=al,ms=ms, label='Ulysses')
ax3.legend(loc=1)
ax3.set_xlim([0,1.8])


plt.tight_layout()
plt.savefig('icmecat/icmecat_overview.png', dpi=100,bbox_inches='tight')


# In[ ]:





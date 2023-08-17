#!/usr/bin/env python
# coding: utf-8

# ## icmecat
# 
# Makes the interplanetary coronal mass ejection catalog ICMECAT, available at https://helioforecast.space/icmecat.
# 
# latest release: version 2.1, 2021 November 29, updated 2023 August TBD
# 
# **Authors**: Christian Möstl, Eva Weiler, Emma E. Davies, Austrian Space Weather Office, Geosphere Austria
# 
# This catalog was built over several years, and we want to thank these persons for contributions along the way:
# 
# **Contributors**: 
# Andreas J. Weiss, Rachel L. Bailey, Martin A. Reiss, Tarik Mohammad Salman, Peter Boakes, Alexey Isavnin, Emilia Kilpua, David Stansby, Reka Winslow, Brian Anderson, Lydia Philpott, Vratislav Krupar, Jonathan Eastwood, Benoit Lavraud, Simon Good, Lan Jian, Teresa Nieves-Chinchilla, Cyril Simon Wedlund, Jingnan Guo, Johan von Forstner, Mateja Dumbovic. 
# 
# https://twitter.com/chrisoutofspace <br /> https://mastodon.social/@chrisoutofspace
# 
# This notebook is part of the heliocats package https://github.com/cmoestl/heliocats
# 
# **Cite this catalog with the doi 10.6084/m9.figshare.6356420**  <br /> https://10.6084/m9.figshare.6356420 <br />
# If this catalog is used for results that are published in peer-reviewed international journals, please contact chris.moestl@outlook.com for possible co-authorship. 
# 
# Install the helio4 conda environment to run this code, see readme at https://github.com/cmoestl/heliocats
# 
# **Adding a new ICME event:** 
# - use the notebook data_update_web_science.ipynb in this package to create pickle files for new science data. The current data can be found on figshare.
# - use measure.ipynb to manually derive the 3 times for each ICME event
# - manually edit the file icmecat/HELCATS_ICMECAT_v21_master.xlsx to add 3 times for each event, the event id and spacecraft name
# - set the switch to create_indices greater 0 and the indices will be redone for the new events so the script quickly loads the info where the ICMEs are in the data files
# - for a new release, set the the last_update variable to the current date
# 
# 
# **KNOWN ISSUES**
# 
# - none

# In[1]:


last_update='2023-August-TBD'

debug_mode=1

#redo positions file
make_positions=1
#red indices file
create_indices=0

#define number of processes for plotting
used=8 
#which plots to make
solo_plots=1
bepi_plots=0
wind_plots=0
psp_plots=0
sta_plots=0



import os
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
import multiprocessing as mp
import time
import pickle
import sys
import urllib
import json
import importlib
import pandas as pd
import copy
import openpyxl
import h5py


import plotly.graph_objects as go
from plotly.offline import iplot, init_notebook_mode


import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames

from heliocats import plot as hp
from heliocats import data as hd
from heliocats import cats as hc
from heliocats import stats as hs

if debug_mode > 0: 
    importlib.reload(hp) #reload again while debugging
    importlib.reload(hd) #reload again while debugging
    importlib.reload(hc) #reload again while debugging
    importlib.reload(hs) #reload again while debugging

#where the in situ data files are located is read, this catalog is produced only locally so far
import config_local
from config_local import data_path

print(data_path)


########### make directories first time if not there

resdir='results'
datadir='data'
indexdir='icmecat/indices_icmecat' 
catdir='icmecat'
icplotsdir='icmecat/plots_icmecat/' 

if os.path.isdir(resdir) == False: os.mkdir(resdir)
if os.path.isdir(datadir) == False: os.mkdir(datadir)
if os.path.isdir(indexdir) == False: os.mkdir(indexdir) 
if os.path.isdir(catdir) == False: os.mkdir(catdir)
if os.path.isdir(icplotsdir) == False: os.mkdir(icplotsdir) 



import warnings
# Ignore all warnings
warnings.filterwarnings("ignore")

t0all = time.time()

print('  ')
print('switches')
print('Debug mode',debug_mode)
print('Make positions',make_positions)
print('Create indices',create_indices)



#Convert this notebook to a script 
os.system('jupyter nbconvert --to script icmecat.ipynb')    


# 
# 
# ### Load positions file

# In[2]:


# the positions file is generated with positions.ipynb, and the positions from messenger, ulysses and stereob are taken from an older file
# these are merged here into a single file that can be found on the figshare repository

if make_positions > 0:

    t0 = time.time()
    #get uly stb and messenger from this file 
    [dum1, dum2, dum3, dum4, p_stb, p_mes, p_uly, dum5, dum6, dum7, dum8,dum9,dum10,dum11,dum12,dumframe]=pickle.load( open( 'results/positions/positions_for_stb_uly_mes_HEEQ_1hour.p', "rb" ) )
    
    #change old matplotlib date to datetime objects, and make angles in degrees    
    p_stb_new = np.zeros(np.size(p_stb),dtype=[('time',object),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
    p_stb_new = p_stb_new.view(np.recarray)        
    p_stb_new.time = p_stb.time + mdates.date2num(np.datetime64('0000-12-31'))
    p_stb_new.x=p_stb.x
    p_stb_new.y=p_stb.y
    p_stb_new.z=p_stb.z    
    p_stb_new.r=p_stb.r   
    #p_stb_new.lat=np.rad2deg(p_stb.lat)
    #p_stb_new.lon=np.rad2deg(p_stb.lon)
    p_stb_new.lat=p_stb.lat
    p_stb_new.lon=p_stb.lon

    
            
    #change old matplotlib date to datetime objects, and make angles in degrees    
    p_uly_new = np.zeros(np.size(p_uly),dtype=[('time',object),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
    p_uly_new = p_uly_new.view(np.recarray)        
    p_uly_new.time = p_uly.time + mdates.date2num(np.datetime64('0000-12-31'))
    p_uly_new.x=p_uly.x
    p_uly_new.y=p_uly.y
    p_uly_new.z=p_uly.z    
    p_uly_new.r=p_uly.r   
    #p_uly_new.lat=np.rad2deg(p_uly.lat)
    #p_uly_new.lon=np.rad2deg(p_uly.lon)
    p_uly_new.lat=p_uly.lat
    p_uly_new.lon=p_uly.lon
    
    
        
    #change old matplotlib date to datetime objects, and make angles in degrees    
    p_mes_new = np.zeros(np.size(p_mes),dtype=[('time',object),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
    p_mes_new = p_mes_new.view(np.recarray)  
    #p_mes_new.time = mdates.num2date(p_mes.time + mdates.date2num(np.datetime64('0000-12-31')))
    p_mes_new.time = p_mes.time + mdates.date2num(np.datetime64('0000-12-31'))
    
    #TBD interpolate messenger to hourly resolution
    #t_start=p_mes.time[0]
    #t_end=p_mes.time[-1]    
    #p_mes_new.time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int (((t_end - t_start).days+1)*60*24))]  
    #p_mes_new.time = np.interp(p_mes., bepi_raw.time, bepi_raw.bx)
    


    p_mes_new.x=p_mes.x
    p_mes_new.y=p_mes.y
    p_mes_new.z=p_mes.z    
    p_mes_new.r=p_mes.r   
    #p_mes_new.lat=np.rad2deg(p_mes.lat)
    #p_mes_new.lon=np.rad2deg(p_mes.lon)
    p_mes_new.lat=p_mes.lat
    p_mes_new.lon=p_mes.lon
    
    
    
    
    #get new spacecraft positions from a pickle file that has been done with positions.ipynb file    
    [p_psp, p_solo, p_sta, p_bepi, p_l1, p_earth, p_mercury, p_venus, p_mars, p_jupiter, p_saturn, p_uranus, p_neptune]=pickle.load( open( 'results/positions/positions_psp_solo_sta_bepi_wind_planets_HEEQ_1hour.p', "rb" ) )

    p_psp.time=mdates.date2num(p_psp.time)
    p_solo.time=mdates.date2num(p_solo.time)    
    p_sta.time=mdates.date2num(p_sta.time)
    p_bepi.time=mdates.date2num(p_bepi.time)
    p_l1.time=mdates.date2num(p_l1.time)
    p_earth.time=mdates.date2num(p_earth.time)
    p_mercury.time=mdates.date2num(p_mercury.time)
    p_venus.time=mdates.date2num(p_venus.time)
    p_mars.time=mdates.date2num(p_mars.time)    

    #all angels in radians
    
    p_psp.lat=np.deg2rad(p_psp.lat)
    p_psp.lon=np.deg2rad(p_psp.lon)

    p_solo.lat=np.deg2rad(p_solo.lat)
    p_solo.lon=np.deg2rad(p_solo.lon)

    p_sta.lat=np.deg2rad(p_sta.lat)
    p_sta.lon=np.deg2rad(p_sta.lon)
    
    p_bepi.lat=np.deg2rad(p_bepi.lat)
    p_bepi.lon=np.deg2rad(p_bepi.lon)

    p_l1.lat=np.deg2rad(p_l1.lat)
    p_l1.lon=np.deg2rad(p_l1.lon)
    
    p_earth.lat=np.deg2rad(p_earth.lat)
    p_earth.lon=np.deg2rad(p_earth.lon)

    p_mercury.lat=np.deg2rad(p_mercury.lat)
    p_mercury.lon=np.deg2rad(p_mercury.lon)

    p_venus.lat=np.deg2rad(p_venus.lat)
    p_venus.lon=np.deg2rad(p_venus.lon)
        
    p_mars.lat=np.deg2rad(p_mars.lat)
    p_mars.lon=np.deg2rad(p_mars.lon)


    pos = np.array([p_psp, p_solo, p_sta, p_bepi, p_l1, p_stb_new, p_uly_new, p_mes_new, p_earth, p_mercury, p_venus, p_mars, p_jupiter, p_saturn, p_uranus, p_neptune])
    ## !!!! add name variable for each spacecraft and then read the whole array at once
    #pos = [p_psp, p_solo, p_sta, p_bepi, p_l1, p_stb_new, p_uly_new, p_mes_new, p_earth, p_mercury, p_venus, p_mars, p_jupiter, p_saturn, p_uranus, p_neptune]
    pickle.dump(pos, open('results/positions/positions_HEEQ_1hr_new.p', 'wb'))
    
    t1 = time.time()
    print('making positions takes', int(np.round(t1-t0,0)), 'seconds')
    print('merged positions file made')
    

pos=pickle.load( open( 'results/positions/positions_HEEQ_1hr_new.p', "rb" ) )

print('positions file loaded')

print('!!!!!!!!!!! fix bug in Earth and L1 latitude in positions.ipynb')


# ## (1) load data 

# In[3]:


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
    
    #use pickle5 to read
    #print('load Juno data ') #Emma Davies https://figshare.com/articles/dataset/Juno_Cruise_Phase_Magnetometer_and_Position_Data/19517257
    #juno_df = pd.read_pickle(data_path+'juno_2011_2016_rtn.pkl') 

    

    ########### CURRENT ACTIVE SPACECRAFT    

    ############### convert MAVEN from Cyril's MAT file to pickle
    
    #from heliocats import data as hd
    #importlib.reload(hd) #reload again while debugging
    
    #file_input=data_path+'input/Data-MAVEN-MAG_SolarWind_102014-012021.mat'        
    #filename=data_path+'maven_2014_2021_removed_no_plasma.p'
    #hd.convert_MAVEN_mat_removed(file_input,filename)

    #filemav=data_path+'maven_2014_2021_removed_no_plasma.p'   
    #filename=data_path+'maven_2014_2021_removed_smoothed_no_plasma.p'   
    #hd.MAVEN_smooth_orbit(filemav,filename)
    
    
    print('load MAVEN data MSO') 
    #filemav='maven_2014_2018.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )

    #combined plasma and mag
    filemav='maven_2014_2018_removed.p'
    [mavr,hmavr]=pickle.load(open(data_path+filemav, 'rb' ) )    
    filemav='maven_2014_2018_removed_smoothed.p'
    [mav,hmav]=pickle.load(open(data_path+filemav, 'rb' ) )

    #only mag
    filemav='maven_2014_2021_removed_no_plasma.p'
    [mav2,hmav2]=pickle.load(open(data_path+filemav, 'rb' ) )    

    filemav='maven_2014_2021_removed_smoothed_no_plasma.p'
    [mavr2,hmavr2]=pickle.load(open(data_path+filemav, 'rb' ) )    

    
    #removed magnetosphere by C. Simon Wedlund, 1 data point per orbit, MSO
    #filemav='maven_2014_2021_removed_smoothed.p'
    #[mav,hmav]=pickle.load(open(data_path+filemav, 'rb' ) )
    
        
    #use hd.save_msl_rad() first to convert data doseE_sol_filter_2019.dat to pickle file
    print('load MSL RAD')
    #MSL RAD
    rad=hd.load_msl_rad(data_path)#, rad.time,rad.dose_sol
    

    ################################ Bepi Colombo
    print('load Bepi Colombo RTN')
    filebepi='bepi_ib_2019_now_rtn.p'
    [bepi,bepiheader]=pickle.load(open(data_path+filebepi, "rb" ) )  
  
   
    ########### STA
    
    print('load and merge STEREO-A data') #yearly magplasma files from stereo science center, conversion to SCEQ 
    filesta1='stereoa_2007_now_rtn.p'
    [sta1,hsta1]=pickle.load(open(data_path+filesta1, "rb" ) )  
    
    #beacon data
    filesta2='stereoa_beacon_last_300days_now.p'
    
    [sta2,hsta2]=pickle.load(open(data_path+filesta2, "rb" ) )  
    #cutoff with end of science data
    sta2=sta2[np.where(sta2.time >= parse_time('2023-Jan-01 00:00').datetime)[0]]

    #make array
    sta=np.zeros(np.size(sta1.time)+np.size(sta2.time),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   

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
    
    del sta1
    del sta2
    
    print('STA Merging done')

   
    ##################### Wind

    filewin="wind_1995_now_rtn.p" 
    [win,hwin]=pickle.load(open(data_path+filewin, "rb" ) )  
   

    ################### SolO
    print('load Solar Orbiter RTN')
    filesolo='solo_2020_now_rtn.p'
    [solo,hsolo]=pickle.load(open(data_path+filesolo, "rb" ) )    
    

    ################### PSP
    print('load PSP data RTN')
    filepsp='psp_2018_now_rtn.p'
    [psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) ) 
        
   

    
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
#print('Juno cruise phase    ',str(uly.time[0])[0:10],str(uly.time[-1])[0:10])


print()
print('catalogs:')
print()
print('HELCATS HIGeoCAT     ',str(higeocat_time[0])[0:10],str(higeocat_time[-1])[0:10])
print('HELCATS ARRCAT       ',np.sort(ac_pandas.sse_launch_time)[0][0:10],np.sort(ac_pandas.sse_launch_time)[-1][0:10])



print('done')

t1 = time.time()
print('loading data takes', int(np.round(t1-t0,0)), 'seconds')


# ## (2) measure new events with measure.ipynb

# In[4]:


#for adding Juno
#junocat='data/HELIO4CAST_ICMECAT_juno_davies_chris.csv'

#jc=pd.read_csv(junocat)

#convert all times to datetime objects
#for i in np.arange(0,ic.shape[0]):    
#remove leading and ending blank spaces if any and write datetime object into dataframe
#        ic.at[i,'icme_start_time']= parse_time(str(ic.icme_start_time[i]).strip()).datetime 
#        ic.at[i,'mo_start_time']=parse_time(str(ic.mo_start_time[i]).strip()).datetime
#        ic.at[i,'mo_end_time']=parse_time(str(ic.mo_end_time[i]).strip()).datetime


# ## (3) make ICMECAT 

# In[32]:


if debug_mode > 0: 
    importlib.reload(hc) 
    importlib.reload(hp)

#master file from 1995 onwards
ic=hc.load_helcats_icmecat_master_from_excel('icmecat/HELIO4CAST_ICMECAT_v21_master.xlsx')


####### 3a get indices for all spacecraft from excel cat
stbi=np.where(ic.sc_insitu == 'STEREO-B')[:][0]    
vexi=np.where(ic.sc_insitu == 'VEX')[:][0]  
mesi=np.where(ic.sc_insitu == 'MESSENGER')[:][0]   
ulyi=np.where(ic.sc_insitu == 'ULYSSES')[:][0]    
mavi=np.where(ic.sc_insitu == 'MAVEN')[:][0]    
stai=np.where(ic.sc_insitu == 'STEREO-A')[:][0]    
wini=np.where(ic.sc_insitu == 'Wind')[:][0] 
pspi=np.where(ic.sc_insitu == 'PSP')[:][0]    
soli=np.where(ic.sc_insitu == 'SolarOrbiter')[:][0]    
beci=np.where(ic.sc_insitu == 'BepiColombo')[:][0]    


if create_indices > 0: 

    print(' ')
    print('switch create_indices on - make indices for all events')
    
    t0 = time.time()

    hc.create_icme_indices(solo,soli,ic,'SolarOrbiter')
    hc.create_icme_indices(win,wini,ic,'Wind')
    hc.create_icme_indices(psp,pspi,ic,'PSP')
    hc.create_icme_indices(sta,stai,ic,'STEREO-A')
    hc.create_icme_indices(bepi,beci,ic,'BepiColombo')
    
    t1 = time.time()
    print('takes', np.round(t1-t0,2), 'seconds')
    print('done ')
    print(' ')

####### 3b get parameters for all spacecraft one after another

#missions to be updated

ic=hc.get_cat_parameters(solo,soli,ic,'SolarOrbiter')
ic=hc.get_cat_parameters(win,wini,ic,'Wind') 
ic=hc.get_cat_parameters(psp,pspi,ic,'PSP')
ic=hc.get_cat_parameters(sta,stai,ic,'STEREO-A')
ic=hc.get_cat_parameters(bepi,beci,ic,'BepiColombo')


ic=hc.get_cat_parameters(mav,mavi,ic,'MAVEN')
#finished missions
ic=hc.get_cat_parameters(stb,stbi,ic,'STEREO-B')
ic=hc.get_cat_parameters(vex,vexi,ic,'VEX')
ic=hc.get_cat_parameters(mes,mesi,ic,'MESSENGER')
ic=hc.get_cat_parameters(uly,ulyi,ic,'ULYSSES')

print('done')


# In[33]:


###### 3c make all plots if wanted
if debug_mode > 0: 
    importlib.reload(hc) 
    importlib.reload(hp)

matplotlib.use('Agg')
t0plot = time.time()


#### check for system type
#server
if sys.platform == 'linux': 
    print('system is linux')
    
#mac
if sys.platform =='darwin':  
    print('system is mac')
    print(os.system('sysctl -n hw.physicalcpu'))
    print(os.system('sysctl -n hw.logicalcpu'))
    
print()

#https://docs.python.org/3/library/multiprocessing.html
            
##### make plots with multiprocessing

def plot_icmecat_events_multi(i):

    #sc,sci,ic,name needs to change for each event and set before this function is called
    #global: pos, data_path, ic, icplotsdir,data_path,pos
    
    outer=60*36 #36h
    dayslim=1.5 #36h

    if name=='BepiColombo':
        hp.plot_insitu_icmecat_mag(sc[icme_start_ind[i]-outer:mo_end_ind[i]+outer],\
                             ic.icme_start_time[sci[i]]-datetime.timedelta(dayslim), \
                             ic.mo_end_time[sci[i]]+datetime.timedelta(dayslim),name, icplotsdir,ic,sci[i],pos)
    else:    
        hp.plot_insitu_icmecat_mag_plasma(sc[icme_start_ind[i]-outer:mo_end_ind[i]+outer],\
                             ic.icme_start_time[sci[i]]-datetime.timedelta(dayslim), \
                             ic.mo_end_time[sci[i]]+datetime.timedelta(dayslim),name, icplotsdir,ic,sci[i],pos)
       
    


#missions to be updated, only make one at a time otherwise variables are changed!

print('Using multiprocessing, nr of cores available: ',mp.cpu_count(), '. Nr of processes used: ',used)

print('start plots')

if solo_plots > 0:
    sc=solo; sci=soli; name='SolarOrbiter'
    counter=np.arange(0,len(sci))
    fileind='icmecat/indices_icmecat/ICMECAT_indices_'+name+'.p'
    [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))  

    #define pool using fork and number of processes
    pool=mp.get_context('fork').Pool(processes=used)
    # Map the worker function onto the parameters    
    t0 = time.time()
    pool.map(plot_icmecat_events_multi, counter) #or use apply_async?,imap
    pool.close()
    pool.join()     
    t1 = time.time()
    multi_time=np.round(t1-t0,2)
    print('SolO done')
    print('plotting takes', np.round(multi_time,2), 'seconds')   
    print(' ')

if bepi_plots > 0:
    sc=bepi; sci=beci; name='BepiColombo'
    counter=np.arange(0,len(sci))
    fileind='icmecat/indices_icmecat/ICMECAT_indices_'+name+'.p'
    [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))  

    #define pool using fork and number of processes
    pool=mp.get_context('fork').Pool(processes=used)
    # Map the worker function onto the parameters    
    t0 = time.time()
    pool.map(plot_icmecat_events_multi, counter) #or use apply_async?,imap
    pool.close()
    pool.join()     
    t1 = time.time()
    multi_time=np.round(t1-t0,2)
    print('Bepi done')
    print('plotting takes', np.round(multi_time,2), 'seconds')   
    print(' ')

    
if wind_plots > 0:
    sc=win; sci=wini; name='Wind'
    counter=np.arange(0,len(sci))
    fileind='icmecat/indices_icmecat/ICMECAT_indices_'+name+'.p'
    [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))  
    
    #define pool using fork and number of processes
    pool=mp.get_context('fork').Pool(processes=used)
    # Map the worker function onto the parameters    
    t0 = time.time()
    pool.map(plot_icmecat_events_multi, counter) #or use apply_async?,imap
    pool.close()
    pool.join()     
    t1 = time.time()
    multi_time=np.round(t1-t0,2)
    print('Wind done')
    print('plotting takes', np.round(multi_time,2), 'seconds')   
    print(' ')

if psp_plots > 0:
    sc=psp; sci=pspi; name='PSP'
    counter=np.arange(0,len(sci))
    fileind='icmecat/indices_icmecat/ICMECAT_indices_'+name+'.p'
    [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))  
    
    #define pool using fork and number of processes
    pool=mp.get_context('fork').Pool(processes=used)
    # Map the worker function onto the parameters    
    t0 = time.time()
    pool.map(plot_icmecat_events_multi, counter) #or use apply_async?,imap
    pool.close()
    pool.join()     
    t1 = time.time()
    multi_time=np.round(t1-t0,2)
    print('PSP done')
    print('plotting takes', np.round(multi_time,2), 'seconds')   
    print(' ')


if sta_plots > 0:
    sc=sta; sci=stai; name='STEREO-A'
    counter=np.arange(0,len(sci))
    fileind='icmecat/indices_icmecat/ICMECAT_indices_'+name+'.p'
    [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))  

    #define pool using fork and number of processes
    pool=mp.get_context('fork').Pool(processes=used)
    # Map the worker function onto the parameters    
    t0 = time.time()
    pool.map(plot_icmecat_events_multi, counter) #or use apply_async?,imap
    pool.close()
    pool.join()     
    t1 = time.time()
    multi_time=np.round(t1-t0,2)
    print('STEREO-A done')
    print('plotting takes', np.round(multi_time,2), 'seconds')   
    print(' ')

print(' ')
print('----------')
t1plot = time.time()
all_plot_time=np.round(t1plot-t0plot,2)
print('all plotting takes', np.round(multi_time/60,2), 'minutes')  


#mag only
#hp.plot_icmecat_events(bepi,beci,ic,'BepiColombo',icplotsdir,data_path,pos)
#
#mag and plasma
#hp.plot_icmecat_events(solo,soli,ic,'SolarOrbiter',icplotsdir,data_path,pos)
#hp.plot_icmecat_events(sta,stai,ic,'STEREO-A',icplotsdir,data_path,pos)
#hp.plot_icmecat_events(psp,pspi,ic,'PSP',icplotsdir,data_path,pos)
#hp.plot_icmecat_events(win,wini,ic,'Wind',icplotsdir,data_path,pos)

#with single processing for the older missions:
#hp.plot_icmecat_events(mav,mavi,ic,'MAVEN',icplotsdir)
#finished missions
#mag only
#hp.plot_icmecat_events(vex,vexi,ic,'VEX',icplotsdir,data_path,pos)
#hp.plot_icmecat_events(mes,mesi,ic,'MESSENGER',icplotsdir)
#mag and plasma
#hp.plot_icmecat_events(stb,stbi,ic,'STEREO-B',icplotsdir)
#hp.plot_icmecat_events(uly,ulyi,ic,'ULYSSES',icplotsdir)
print('done')


get_ipython().run_line_magic('matplotlib', 'inline')


# ### (4) save ICMECAT 

# ### 4a save header

# In[18]:


######## sort ICMECAT by date

ic = ic.sort_values(by='icme_start_time',ascending=False)
ic = ic.reset_index(drop=True)

#save header and parameters as text file and prepare for html website
header='ICME CATALOGUE v2.1 \n\n\
This is the HELIO4CAST interplanetary coronal mass ejection (ICME) catalog, based on in situ magnetic field and bulk plasma observations in the heliosphere. \n\n\
This is version 2.1, released 2021-November-29, updated '+last_update+', doi: 10.6084/m9.figshare.6356420 \n\n\
Rules: If results are produced with this catalog for peer-reviewed scientific publications, please contact chris.moestl@outlook.com for possible co-authorship. \n\
References for this catalog are Moestl et al. 2017 (https://doi.org/10.1002/2017SW001614) and Moestl et al. 2020 (https://doi.org/10.3847/1538-4357/abb9a1).  \n\n\
The catalog is available as a python pandas dataframe (pickle), python numpy structured array (pickle), json, csv, xlsx, txt, hdf5, at \n\
https://helioforecast.space/icmecat \n\n\
Number of events in ICMECAT: '+str(len(ic))+' \n\
ICME observatories: Solar Orbiter, Parker Solar Probe (PSP), BepiColombo, Wind, STEREO-A, MAVEN, STEREO-B, Venus Express (VEX), MESSENGER, Ulysses.   \n\
Time range: January 1995 - July 2023. \n \n\
Authors: Christian Moestl, E. Weiler, Emma E. Davies, Andreas J. Weiss, Rachel L. Bailey, Martin A. Reiss, \n\
Austrian Space Weather Office, GeoSphere Austria, Graz, Austria; NASA CCMC, USA.\n\
Contributors: Tarik Mohammad Salman, Peter Boakes, Alexey Isavnin, Emilia Kilpua, David Stansby, Reka Winslow, Brian Anderson, Lydia Philpott, \n\
Vratislav Krupar, Jonathan Eastwood, Simon Good, Lan Jian, Teresa Nieves-Chinchilla, Cyril Simon Wedlund, Jingnan Guo, \n\
Johan von Forstner, Mateja Dumbovic, Benoit Lavraud.  \n\n\
This catalog has been made by getting the 3 times of each ICME (shock or disturbance begin, magnetic obstacle start and end) \n\
from the individual catalogs below, and then calculating all parameters again consistently from the data by us. \n\
The selection criteria are relatively simple: the event needs to have a clearly organized large-scale magnetic structure,\n \
which manifests itself through (1) an elevated total magnetic field and (2) a rotation in the field observed through the field components. \n\n\
The in situ data that were used for the catalog, with a size of around 10 GB in total, including extra data files with magnetic field components \n\
in RTN coordinates that are not used for producing the catalog, can be downloaded in python pickle format as recarrays from \
https://doi.org/10.6084/m9.figshare.11973693 \n\n\
The python source code for producing this catalog is available at https://github.com/cmoestl/heliocats icmecat.ipynb\n\n\
Each event has unique identifier - the icmecat_id - which has a tag in it that indicates from which catalog the ICME times were taken: \n\n\
Wind:         Nieves-Chinchilla et al. (2018), tag: NASA. \n\
STEREO-A:     Jian et al. (2018), tag: JIAN. \n\
STEREO-B:     Jian et al. (2018), tag: JIAN. \n\
VEX:          Good et al. (2018), tag: SGOOD \n\
MESSENGER:    Good et al. (2018), Winslow et al. (2018), tags: SGOOD, WINSLOW. \n\
MAVEN:        Made by us according to the method in the comments, tag: MOESTL.\n\
Ulysses:      Added by us, tag: MOESTL. \n\
SolarOrbiter: Added by us, tag: MOESTL. \n\
BepiColomo:   Added by us, tag: MOESTL. \n\
PSP:          Added by us, tag: MOESTL. \n\n\
We have also added extra events at VEX, MESSENGER, Wind and STEREO-A (all tagged with MOESTL or MOESTL_SALMAN in icmecat_id).\n\n\
References: \n\
Nieves-Chinchilla, T. et al. (2018),  https://doi.org/10.1007/s11207-018-1247-z \n\
                                      https://wind.nasa.gov/ICME_catalog/ICME_catalog_viewer.php \n\
Jian, L. et al. (2018), https://doi.org/10.3847/1538-4357/aab189 \n\
                        https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level3/ \n\
Good, S. et al. (2018) https://doi.org/10.1007/s11207-015-0828-3 \n\
Winslow, R. et al. (2015), https://doi.org/10.1002/2015JA021200 \n\n\n\
Comments: \n\n\
- Spacecraft positions are given in Heliocentric Earth Equatorial Coordinates (HEEQ) coordinates. \n\n\
- The coordinate system for all magnetic field components is SCEQ, except for Wind (HEEQ, which is the equivalent for SCEQ for Earth) and Ulysses (RTN, because of high latitude positions) and MAVEN (MSO). \n\
        Definition of SpaceCraft Equatorial Coordinates (SCEQ): \n\
        Z is the solar rotation axis. \n\
        Y is the cross product of Z and R, with R being the vector that points from the Sun to the spacecraft.\n\
        X completes the right handed triad (and points away from the Sun). \n\
This system is thus like HEEQ but centered on the respective in situ spacecraft, so the SCEQ X and Y \n\
base vectors are rotated by the HEEQ longitude of the in situ spacecraft from HEEQ X and Y.\n\
The Y vector is similar to the T vector in an RTN system for each spacecraft, but the X and Z vectors \n\
are rotated around Y compared to an RTN system. The differences between RTN and SCEQ for spacecraft within \n\
a few degrees of the solar equatorial plane are very small (within a few 0.1 nT usually).\n\
We choose SCEQ for the advantage that a comparison of the data of multipoint CME events \n\
and comparisons to simulations are simpler because there is always a similar reference plane (the solar equatorial plane).\n\n\
- Venus Express and MESSENGER do not have plasma parameters available. For events that do not have plasma parameters at Solar Orbiter or PSP, they may become available in the future.\n\n\
- If there is no sheath or density pileup region, so the ICME starts immediately with a magnetic obstacle, the icme_start_time is similar to mo_start_time.\n\n\
- At MESSENGER and VEX, for events cataloged by Simon Good, icme_start_time has been added by V. Krupar (Imperial College) and C. Moestl (IWF Graz). \n\n\
- For the calculation of the parameters at MESSENGER during the orbit around Mercury, all data points inside the bowshock of Mercury have been removed, \n\
according to a list thankfully provided to us by by R. Winslow, UNH, B. Anderson, APL, and Lydia Philpott, UBC. \n\n\
- Calculation of the magnetic obstacle parameters at VEX is done after approximate removal of the induced magnetosphere, with a modified equation \
as in \n\
Zhang et al. 2008 (doi: 10.1016/j.pss.2007.09.012), with a constant of 3.5 instead of 2.14/2.364, \n\
in order to account for a larger bowshock distance during solar maximum than studied in the Zhang et al. (2008) paper. \n\n\
- For MAVEN, all data inside the bow shock were removed with the model from Gruesbeck et al. (2018, doi:10.1029/2018JA025366) by C. Simon \
Wedlund (IWF Graz, Austria). \n\
From the remaining data, the median for each orbit is taken as 1 data point, resulting in a solar wind \n\
dataset at Mars with 4.5 hour time resolution. This is a far lower resolution than available at 1 AU, so the identification \n\
of ICMEs in the MAVEN data is not as straightforward as for data sets with higher time resolution.\n\n\
- The identification of ICMEs for MAVEN was done as follows: (1) We looked for typical profiles for ICMEs in magnetic field (B) and speed (V).\n\
(2) B needs to be elevated in ICMEs accompanied by a flat or declining profile in V. An elevation in B accompanied by a later gradual rise \n\
in V is a strong indicator for a high speed stream (HSS). (3) The discrimination between HSS and ICMEs was further strengthened with \n\
the help of a stream interaction region list for MAVEN by Huang et al. (2019, doi: 10.3847/1538-4357/ab25e9). \n\
(4) We additionally plotted the predicted CME arrivals with STEREO/HI (ARRCATv2.0), the dose rate measured on the Mars surface \n\
by MSL/RAD (e.g. Guo et al. (2018, doi: 10.1051/0004-6361/201732087)) and the speed of the WSA/HUX model for the \n\
background solar wind (Reiss et al. 2019, doi: 10.3847/1538-4365/aaf8b3). \n\
General guidance on space weather observed by MAVEN can be found in Lee et al. (2018, doi:10.1002/2016JA023495).\n\n\n\n'


parameters='Parameters:\n\
00: icmecat_id: The unique identifier for the observed ICME. unit: string. \n\
01: sc insitu: The name of the in situ observing spacecraft. unit: string. \n\
02: icme_start_time: Shock arrival or density enhancement time, can be similar to mo_start_time. unit: UTC. \n\
03: mo_start_time: Start time of the magnetic obstacle (MO), including flux ropes, flux-rope-like, and ejecta signatures. unit: UTC. \n\
04: mo_end_time: End time of the magnetic obstacle. unit: UTC. \n\
05: mo_sc_heliodistance: Heliocentric distance of the spacecraft at mo_start_time. unit: AU.\n\
06: mo_sc_long_heeq: Heliospheric longitude of the spacecraft at mo_start_time, range [-180,180]. unit: degree (HEEQ).\n\
07: mo_sc_lat_heeq: Heliospheric latitude of the spacecraft at mo_start_time, range [-90,90]. unit: degree (HEEQ).\n\
08: icme_duration: Duration of the interval between icme_start_time and mo_endtime. unit: hours.\n\
09: icme_bmax: Maximum total magnetic field in the full icme interval (icme_start_time to mo_end_time). unit: nT.\n\
10: icme_bmean: Mean total magnetic field during the full icme interval (icme_start_time to mo_end_time). unit: nT.\n\
11: icme_bstd: Standard deviation of the total magnetic field from icme_start_time to mo_end_time. unit: nT.\n\
12: icme_speed_mean: Mean proton speed from icme_start_time to mo_end_time. unit: km/s.\n\
13: icme_speed_std: Standard deviation of proton speed from icme_start_time to mo_end_time. unit: km/s.\n\
14: mo_duration: Duration of the interval between mo_start_time and mo_endtime. unit: hours.\n\
15: mo_bmax: Maximum total magnetic field in the magnetic obstacle interval (mo_start_time to mo_end_time). unit: nT.\n\
16: mo_bmean: Mean total magnetic field in the magnetic obstacle. unit: nT.\n\
17: mo_bstd: Standard deviation of the total magnetic field in the magnetic obstacle. unit: nT.\n\
18: mo_bzmean: Mean magnetic field Bz component in the magnetic obstacle. unit: nT.\n\
19: mo_bzmin: Minimum magnetic field Bz component in the magnetic obstacle. unit: nT.\n\
20: mo_bzstd: Standard deviation of the magnetic field Bz component in the magnetic obstacle. unit: nT.\n\
21: mo_bymean: Mean magnetic field By component in the magnetic obstacle. unit: nT.\n\
22: mo_bystd: Standard deviation of the magnetic field By component in the magnetic obstacle. unit: nT.\n\
23: mo_speed_mean: Mean proton speed from mo_start_time to mo_end_time. unit: km/s.\n\
24: mo_speed_std: Standard deviation of proton speed from mo_start_time to mo_end_time. unit: km/s.\n\
25: mo_expansion_speed: Difference between proton speed at mo_start_time to proton speed at mo_end_time. unit: km/s.\n\
26: mo_pdyn_mean: Mean proton dynamic pressure from mo_start_time to mo_start_time. unit: nPa.\n\
27: mo_pdyn_std: Standard deviation of proton dynamic pressure from mo_start_time to mo_start_time. unit: nPa.\n\
28: mo_density_mean: Mean proton density from mo_start_time to mo_start_time. unit: cm^-3.\n\
29: mo_density_std: Standard deviation of proton density from mo_start_time to mo_start_time. unit: cm^-3.\n\
30: mo_temperature_mean: Mean proton temperature from mo_start_time to mo_start_time. unit: K.\n\
31: mo_temperature_std: Standard deviation of proton temperature from mo_start_time to mo_end_time. unit: K.\n\
32: sheath_speed_mean: Mean proton speed from icme_start_time to mo_start_time, NaN if these times are similar. unit: km/s.\n\
33: sheath_speed_std: Standard deviation of proton speed from icme_start_time to mo_start_time, NaN if these times are similar. unit: km/s.\n\
34: sheath_density_mean: Mean proton density from icme_start_time to mo_start_time, NaN if these times are similar. unit: cm^-3.\n\
35: sheath_density_std: Standard deviation of proton density from icme_start_time to mo_start_time, NaN if these times are similar. unit: cm^-3.\n\
36: sheath_pdyn_mean: Mean proton dynamic pressure, from icme_start_time to mo_start_time, NaN if these times are similar. unit: nPa.\n\
37: sheath_pdyn_std: Standard deviation of proton dynamic pressure, from icme_start_time to mo_start_time, NaN if these times are similar. unit: nPa.\n\n\n'



print(header)
print(parameters)


#make header file
file='icmecat/HELIO4CAST_ICMECAT_v21_header.txt'
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

# In[19]:


########## python formats

# save ICMECAT as pandas dataframe with times as datetime objects as pickle
file='icmecat/HELIO4CAST_ICMECAT_v21_pandas.p'
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



file='icmecat/HELIO4CAST_ICMECAT_v21_numpy.p'
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
file='icmecat/HELIO4CAST_ICMECAT_v21.xlsx'
ic_copy.to_excel(file,sheet_name='ICMECATv2.0')
print('ICMECAT saved as '+file)

#save as json
file='icmecat/HELIO4CAST_ICMECAT_v21.json'
ic_copy.to_json(file)
print('ICMECAT saved as '+file)

#save as csv
file='icmecat/HELIO4CAST_ICMECAT_v21.csv'
ic_copy.to_csv(file)
print('ICMECAT saved as '+file)


#save as txt
file='icmecat/HELIO4CAST_ICMECAT_v21.txt'
np.savetxt(file, ic_copy.values.astype(str), fmt='%s' )
print('ICMECAT saved as '+file)








#########################


#########save into hdf5 format , use S for strings http://docs.h5py.org/en/stable/strings.html#what-about-numpy-s-u-type
dtype2=[('index','i8'),('icmecat_id', 'S30'),('sc_insitu', 'S20')] +[(i, '<f8') for i in ic.keys()[2:len(ic.keys())]]
ich5=np.array(ic_num_rec,dtype=dtype2)
file='icmecat/HELIO4CAST_ICMECAT_v21.h5'
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
file='icmecat/HELIO4CAST_ICMECAT_v21_numpy.npy'
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
file='icmecat/HELIO4CAST_ICMECAT_v21_isot.json'
ic_copy2.to_json(file)
print('ICMECAT saved as '+file)


#save as html no header
file='icmecat/HELIO4CAST_ICMECAT_v21_simple.html'
ic_copy.to_html(file)
print('ICMECAT saved as '+file)


############ save as html file with header
#save as html
file='icmecat/HELIO4CAST_ICMECAT_v21.html'
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

# In[20]:


#load icmecat as pandas dataframe
file='icmecat/HELIO4CAST_ICMECAT_v21_pandas.p'
[ic_pandas,h,p]=pickle.load( open(file, 'rb'))   

#load icmecat as numpy array
file='icmecat/HELIO4CAST_ICMECAT_v21_numpy.p'
[ic_nprec,ic_np,h,p]=pickle.load( open(file, 'rb'))   


# In[21]:


print(ic_pandas.keys())



# In[22]:


ic_pandas


# In[23]:


#
ic_nprec


# In[24]:


ic_nprec.icmecat_id


# ## 5 plots

# In[25]:


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
isol=np.where(ic.sc_insitu=='SolarOrbiter')[0]
ibep=np.where(ic.sc_insitu=='BepiColombo')[0]
iuly=np.where(ic.sc_insitu=='ULYSSES')[0]


sns.set_context("talk")     
sns.set_style('darkgrid')

###############################################################################
fig=plt.figure(3,figsize=(18,7),dpi=200)

#########################################################################
ax1=plt.subplot(121)
plt.title('ICMECAT event times and distance (from 2007 onwards)')

#markersize
ms=6
#alpha
al=0.8

ax1.plot_date(ic_mo_start_time_num[iuly],ic.mo_sc_heliodistance[iuly],'o',c='brown', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[imes],ic.mo_sc_heliodistance[imes],'o',c='coral', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[ivex],ic.mo_sc_heliodistance[ivex],'o',c='orange', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[istb],ic.mo_sc_heliodistance[istb],'o',c='royalblue', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[iwin],ic.mo_sc_heliodistance[iwin],'o',c='mediumseagreen', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[imav],ic.mo_sc_heliodistance[imav],'s',c='orangered', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[ista],ic.mo_sc_heliodistance[ista],'o',c='red', alpha=al,ms=ms)

ax1.plot_date(ic_mo_start_time_num[ipsp],ic.mo_sc_heliodistance[ipsp],'o',c='black', alpha=al,ms=ms)
ax1.plot_date(ic_mo_start_time_num[isol],ic.mo_sc_heliodistance[isol],'o',c='black',markerfacecolor='white', alpha=1.0,ms=ms)
ax1.plot_date(ic_mo_start_time_num[ibep],ic.mo_sc_heliodistance[ibep],'s',c='darkblue',markerfacecolor='lightgrey', alpha=al,ms=ms)



ax1.set_ylabel('Heliocentric distance $R$ [AU]')
ax1.set_xlabel('Year')
ax1.set_ylim([0,1.7])

years = mdates.YearLocator()   # every year
ax1.xaxis.set_major_locator(years)
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)

ax1.tick_params(axis="x", labelsize=12)
ax1.tick_params(axis="y", labelsize=12)

ax1.set_xlim([datetime.datetime(2007,1,1),datetime.datetime.utcnow()])


ax1.set_yticks(np.arange(0,2,0.1))
ax1.set_ylim([0,1.7])



##############################################################################
ax3=plt.subplot(122)
plt.title('ICMECAT mean magnetic field in the magnetic obstacle')
ax3.set_xlabel('Heliocentric distance $R$ [AU]')
ax3.set_ylabel('$<B>$ [nT]')



ax3.plot(ic.mo_sc_heliodistance[istb],ic.mo_bmean[istb],'o',c='royalblue', alpha=al,ms=ms, label='STEREO-B')
ax3.plot(ic.mo_sc_heliodistance[imes],ic.mo_bmean[imes],'o',c='coral', alpha=al,ms=ms,label='MESSENGER')
ax3.plot(ic.mo_sc_heliodistance[ivex],ic.mo_bmean[ivex],'o',c='orange', alpha=al,ms=ms,label='Venus Express')
ax3.plot(ic.mo_sc_heliodistance[iuly],ic.mo_bmean[iuly],'o',c='brown', alpha=al,ms=ms, label='Ulysses')
ax3.plot(ic.mo_sc_heliodistance[imav],ic.mo_bmean[imav],'s',c='orangered', alpha=al,ms=ms, label='MAVEN')



ax3.plot(ic.mo_sc_heliodistance[ista],ic.mo_bmean[ista],'o',c='red', alpha=al,ms=ms, label='STEREO-A')
ax3.plot(ic.mo_sc_heliodistance[iwin],ic.mo_bmean[iwin],'o',c='mediumseagreen', alpha=al,ms=ms,label='Wind at L1')
ax3.plot(ic.mo_sc_heliodistance[ibep],ic.mo_bmean[ibep],'s',c='darkblue',markerfacecolor='lightgrey', alpha=al,ms=ms,label='Bepi Colombo')
ax3.plot(ic.mo_sc_heliodistance[ipsp],ic.mo_bmean[ipsp],'o',c='black', alpha=al,ms=ms, label='Parker Solar Probe',zorder=3)
ax3.plot(ic.mo_sc_heliodistance[isol],ic.mo_bmean[isol],'o',c='black', markerfacecolor='white',alpha=al,ms=ms, label='Solar Orbiter',zorder=3)


ax3.legend(loc=1,fontsize=15)

ax3.set_xticks(np.arange(0,2,0.1))
ax3.tick_params(axis="x", labelsize=12)
ax3.set_xlim([0,1.7])

ax3.set_yscale('log')

#ax3.set_ylim([0,np.max(ic.mo_bmean)+50])
#ax3.set_yticks(np.arange(0,1000,10))
#ax3.set_ylim([0,1000])

#ax3.tick_params(axis="y", labelsize=12)




plt.tight_layout()
plt.savefig('icmecat/icmecat_overview.png', dpi=150,bbox_inches='tight')


# In[26]:


#markersize
ms=25
ms2=7
#alpha
al=0.7

sns.set_context("talk")     
sns.set_style('darkgrid')

###############################################################################
fig=plt.figure(3,figsize=(9,9),dpi=200)

#########################################################################
ax2=plt.subplot(111,projection='polar')

plt.title('ICMECAT events [HEEQ longitude, AU]')

ax2.scatter(np.radians(ic.mo_sc_long_heeq[iuly]),ic.mo_sc_heliodistance[iuly],s=ms,c='brown', alpha=al)
ax2.scatter(np.radians(ic.mo_sc_long_heeq[imes]),ic.mo_sc_heliodistance	[imes],s=ms,c='coral', alpha=0.6)
ax2.scatter(np.radians(ic.mo_sc_long_heeq[ivex]),ic.mo_sc_heliodistance	[ivex],s=ms,c='orange', alpha=al)
ax2.scatter(np.radians(ic.mo_sc_long_heeq[iwin]),ic.mo_sc_heliodistance	[iwin],s=ms,c='mediumseagreen', alpha=al)
ax2.scatter(np.radians(ic.mo_sc_long_heeq[istb]),ic.mo_sc_heliodistance[istb],s=ms,c='royalblue', alpha=al)

ax2.scatter(np.radians(ic.mo_sc_long_heeq[imav]),ic.mo_sc_heliodistance	[imav],s=ms,c='orangered', alpha=al)
ax2.scatter(np.radians(ic.mo_sc_long_heeq[ista]),ic.mo_sc_heliodistance[ista],s=ms,c='red', alpha=al)

ax2.plot(np.radians(ic.mo_sc_long_heeq[ipsp]),ic.mo_sc_heliodistance[ipsp],'s',markersize=ms2, c='black', alpha=1.0)
ax2.plot(np.radians(ic.mo_sc_long_heeq[isol]),ic.mo_sc_heliodistance[isol],'s',markersize=ms2, c='black',markerfacecolor='white', alpha=al)
ax2.plot(np.radians(ic.mo_sc_long_heeq[ibep]),ic.mo_sc_heliodistance[ibep],'s',markersize=ms2, c='darkblue',markerfacecolor='lightgrey', alpha=al)



ax2.set_ylim([0,1.8])
ax2.tick_params(labelsize=12,zorder=3)

plt.tight_layout()
plt.savefig('icmecat/icmecat_overview2.png', dpi=100,bbox_inches='tight')


# ## Parameter distribution plots

# In[27]:


#make distribution plots
plt.figure(10,figsize=(21,12),dpi=200)



ax2=plt.subplot(231)
sns.histplot(ic.mo_bmax, label='magnetic obstacle $max(B_t)$',color='mediumseagreen',kde=True, bins=np.arange(0,100,2.5), stat='probability')
sns.histplot(ic.mo_bmean, label='magnetic obstacle $<B_t>$',color='coral',alpha=0.5,kde=True, bins=np.arange(0,100,2.5),stat='probability')
plt.legend(loc=1)
ax2.set_xlabel('magnetic field $B$ [nT]')
ax2.set_xlim(0,100)


ax3=plt.subplot(232)
sns.histplot(ic.mo_bzmin,label='MO $min(B_z)$',color='mediumseagreen',kde=True, bins=np.arange(-50,50,2.5),stat='probability')
sns.histplot(ic.mo_bzmean,label='MO $<B_z>$',color='coral',alpha=0.5,kde=True, bins=np.arange(-50,50,2.5),stat='probability')
plt.legend(loc=1)
ax3.set_xlabel('magnetic field $B$ [nT]')
plt.xlim(-60,60)


ax4=plt.subplot(233)
sns.histplot(ic.mo_speed_mean,label='MO $<V>$',color='mediumseagreen',kde=True,stat='probability')
sns.histplot(ic.mo_expansion_speed,label='MO $V_{exp}$',color='coral',alpha=0.5,kde=True,stat='probability')
plt.legend(loc=1)
ax4.set_xlabel('plasma bulk speed $V$ [km/s]')


ax5=plt.subplot(234)
sns.histplot(ic.mo_sc_heliodistance,color='coral',bins=np.arange(0,1.8,0.1))
ax5.set_xlim(0,1.8)
ax5.set_xlabel('Heliocentric distance $R$[AU]')


ax6=plt.subplot(235)
sns.histplot(ic.icme_duration,color='coral',alpha=0.5,bins=np.arange(0,np.max(ic.icme_duration)+5,3),kde=True,stat='count',label='ICME duration')
sns.histplot(ic.mo_duration,color='mediumseagreen',bins=np.arange(0,np.max(ic.mo_duration)+5,3),kde=True, stat='count',label='Magnetic obstacle duration')
ax6.set_xlabel('duration [hours]')
ax6.set_xlim(0,np.max(ic.icme_duration)+5)
plt.legend(loc=1)

ax7=plt.subplot(236)
#sns.histplot(ic.mo_duration,color='coral',bins=np.arange(0,np.max(ic.mo_duration)+5,3),kde=True, stat='count')
#sheath duration
sns.histplot(ic.icme_duration-ic.mo_duration,color='mediumseagreen',bins=np.arange(0,np.max(ic.icme_duration-ic.mo_duration)+5,2),kde=True, stat='count')
ax7.set_xlabel('sheath duration [hours]')
ax7.set_xlim(0,np.max(ic.icme_duration-ic.mo_duration)+5)



plt.tight_layout()
plt.savefig('icmecat/icmecat_overview3.png', dpi=150,bbox_inches='tight')


# In[29]:


t1all = time.time()
print('the full ICMECAT takes', np.round((t1all-t0all)/60,2), 'minutes')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





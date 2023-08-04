#!/usr/bin/env python
# coding: utf-8

# ## data updates with high frequency for the helioforecast.space website
# 
# Main author: C. Möstl, Austrian Space Weather Office, GeoSphere Austria
# 
# https://github.com/cmoestl/heliocats
# 
# uses environment 'envs/env_helio4.yml'

# In[1]:


# https://github.com/cmoestl/heliocats  data_update_web_hf.py

# for updating data with high frequency on the servers

import pickle
import importlib
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
import sys
import numpy as np
import datetime
import scipy.signal
import urllib
import json
import os   
import time
import h5py
import pytz
import copy
import cdflib
import sunpy

from heliocats import data as hd
from heliocats import plot as hp

##### check for system type
#server
if sys.platform == 'linux': 
    print('system is linux')
    matplotlib.use('Agg') 
#mac
if sys.platform =='darwin':  
    print('system is mac')
    #for testing
    get_ipython().run_line_magic('matplotlib', 'inline')
    #matplotlib.use('Agg') 



################################################ CHECK  ##############################################

#make sure to convert the current notebook to a script
os.system('jupyter nbconvert --to script data_update_web_hf.ipynb')   

#switches
debug_mode=0

#download and save as pickle
get_noaa=1
get_stereoa=1

####################################################################################################################

#test execution times
t0all = time.time()


# ### Configure paths depending on server or local machine
# 

# In[2]:


if sys.platform == 'linux': 
    
    from config_server import data_path
    from config_server import noaa_path
    from config_server import wind_path
    from config_server import stereoa_path
    from config_server import data_path_ml
    
if sys.platform =='darwin':  

    from config_local import data_path
    from config_local import noaa_path
    from config_local import wind_path
    from config_local import stereoa_path
    from config_local import data_path_ml

print(' ')
print('------ PATHS ')

print(data_path)
print(noaa_path)
print(wind_path)
print(stereoa_path)
print(data_path_ml)


plot_path=data_path+'plots/'
position_path=data_path+'plots_positions/'
sun_path=data_path+'plots_sun/'

print(plot_path)
print(position_path)



########### make directories first time
if os.path.isdir(plot_path) == False: os.mkdir(plot_path)
if os.path.isdir(plot_path+'omni2') == False: os.mkdir(plot_path+'omni2')
if os.path.isdir(plot_path+'wind') == False: os.mkdir(plot_path+'wind')
if os.path.isdir(plot_path+'stereoa') == False: os.mkdir(plot_path+'stereoa')
if os.path.isdir(plot_path+'combined') == False: os.mkdir(plot_path+'combined')


if os.path.isdir(position_path) == False: os.mkdir(position_path)
if os.path.isdir(sun_path) == False: os.mkdir(sun_path)
if os.path.isdir(noaa_path) == False: os.mkdir(noaa_path)
if os.path.isdir(data_path_ml) == False: os.mkdir(data_path_ml)


# ### positions and SDO plot

# In[3]:


#if debug_mode > 0: 
importlib.reload(hd) 
importlib.reload(hp) 


t0 = time.time()
print(' ')
print('------ POSITIONS ')

# spacecraft positions image
hp.plot_positions(datetime.datetime.utcnow(),position_path, 'HEEQ',now=True)

print(' ')
print('------ SDO realtime images ')
# get current SDO images 
hd.get_sdo_realtime_image(sun_path)

t1 = time.time()
print()
print('Positions and SDO images takes', np.round(t1-t0,2), 'seconds')


# In[ ]:





# In[ ]:





# ### NOAA real time solar wind 

# In[4]:


if debug_mode > 0: 
    from heliocats import data as hd
    importlib.reload(hd) 

    from heliocats import plot as hp
    importlib.reload(hp) 

print(' ')
print('------ NOAA real time solar wind data ')

#test execution times
t0 = time.time()

#define file name
filenoaa='noaa_rtsw_last_30files_now.p'

if get_noaa > 0:

    ######## DOWNLOAD NOAA data 
    
    print('download NOAA real time solar wind plasma and mag and dst')
    datestr=str(datetime.datetime.utcnow().strftime("%Y-%m-%d"))
    print(datestr+' UTC')

    plasma='http://services.swpc.noaa.gov/products/solar-wind/plasma-7-day.json'
    mag='http://services.swpc.noaa.gov/products/solar-wind/mag-7-day.json'
    dst='https://services.swpc.noaa.gov/products/kyoto-dst.json'

    try: 
        urllib.request.urlretrieve(plasma, noaa_path+'plasma/plasma-7-day_'+datestr+'.json')
        print(noaa_path+'plasma/plasma-7-day_'+datestr+'.json')
    except urllib.error.URLError as e:
        print(' ', plasma,' ',e.reason)

    try: 
        urllib.request.urlretrieve(mag, noaa_path+'mag/mag-7-day_'+datestr+'.json')
        print(noaa_path+'mag/mag-7-day_'+datestr+'.json')
    except urllib.error.URLError as e:
        print(' ', mag,' ',e.reason)
    
    try: 
        urllib.request.urlretrieve(dst, noaa_path+'dst/dst-7-day_'+datestr+'.json')
        print(noaa_path+'dst/dst-7-day_'+datestr+'.json')
    except urllib.error.URLError as e:
        print(' ', dst,' ',e.reason)
 
    print('NOAA RTSW download complete')
        
    ######## SAVE NOAA DATA AS PICKLE
    # last parameter gives a cutoff, so only the latest N files are taken for the NOAA data pickle file
    
    ############# !!!!!!!!! add Dst in this function as extra parameter to pickle file
    
    hd.save_noaa_rtsw_data(data_path,noaa_path,filenoaa,30)
    print('NOAA RTSW saved as pickle file complete')
    

else:
    print('NOAA data NOT downloaded and saved as pickle, turn on switch')  

    
[noaa,hnoaa]=pickle.load(open(data_path+filenoaa, "rb" ) ) 

 
#plot
start=datetime.datetime.utcnow() - datetime.timedelta(days=10)
end=datetime.datetime.utcnow() 
hp.plot_insitu_update_noaa_rtsw(noaa, start, end,'NOAA_RTSW',plot_path+'noaa/',now=True)


t1 = time.time()

print()
print('NOAA download latest file, save as pickle last 30 files and plotting takes', np.round(t1-t0,2), 'seconds')


    


# ### STEREO-A beacon data

# In[5]:


if debug_mode > 0: 
    importlib.reload(hd) 
    importlib.reload(hp) 

#test execution times
t0 = time.time()

print(' ')
print('------ STEREO-A beacon data ')

#define filename
file_sta_beacon='stereoa_beacon_last_30days_now.p'   

if get_stereoa > 0:

    print(' ')
    print('--- download STEREO-A beacon data last 10 days ')
    hd.stereoa_download_beacon(start=datetime.datetime.utcnow()-datetime.timedelta(days=10),end=datetime.datetime.utcnow(),stereoa_path=stereoa_path)   
    print(' ')

    print('--- process STEREO-A beacon data to pickle last 30 days') 
    #save pickle file
    hd.save_stereoa_beacon_data(data_path,stereoa_path,file_sta_beacon,datetime.datetime.utcnow()-datetime.timedelta(days=30),datetime.datetime.utcnow(),coord='RTN' )   
    
#load pickle    
[sta,hsta]=pickle.load(open(data_path+file_sta_beacon, "rb" ) )  

#plot
start=datetime.datetime.utcnow() - datetime.timedelta(days=10)
end=datetime.datetime.utcnow() 
hp.plot_insitu_update_stereoa_beacon(sta, start, end,'STEREO-A_beacon',plot_path+'stereoa/',now=True)
    
    
    
t1 = time.time()

print()
print('STEREO-A beacon downloading last 30 days, save as pickle last 30 days and plotting takes', np.round(t1-t0,2), 'seconds')


# ## Combined plot STEREO-A NOAA RTSW

# In[20]:


if debug_mode > 0:     
    importlib.reload(hd) 
    importlib.reload(hp) 

start=datetime.datetime.utcnow() - datetime.timedelta(days=3)
end=datetime.datetime.utcnow() 
hp.plot_insitu_update_stereoa_noaa(noaa, sta, start, end,'NOAA_RTSW_STEREO-A_beacon',plot_path+'combined/',now=True)
    


# In[7]:


t1all = time.time()

print(' ')
print(' ')
print(' ')
print('------------------')
print('Runtime for full high frequency data update:', np.round((t1all-t0all),2), 'seconds')
print('--------------------------------------------------------------------------------------')


# In[ ]:





# In[ ]:





# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# ## daily data updates for the helioforecast.space website of science data
# 
# Main author: C. Möstl, Austrian Space Weather Office, GeoSphere Austria
# 
# https://github.com/cmoestl/heliocats
# 
# uses environment 'envs/env_helio4.yml'
# 

# In[20]:


# https://github.com/cmoestl/heliocats  data_update_web.py

# for updating data every day on the servers

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
import pandas as pd

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
    #%matplotlib inline     
    matplotlib.use('Agg') 



################################################ CHECK  ##############################################
os.system('jupyter nbconvert --to script data_update_web_science.ipynb')   

#switches
debug_mode=0


#switches
get_wind=1
get_bepi=1
get_stereoa=0
get_solo=0
get_psp=0

####################################################################################################################


#test execution times
t0all = time.time()


# ### Configure paths depending on server or local machine
# 

# In[21]:


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


# ### Wind data

# In[18]:


if debug_mode > 0: 
    importlib.reload(hd) 
    importlib.reload(hp) 
    
print(' ')
print('------ download Wind data ')

#download data for current year only    
if get_wind > 0:
    hd.wind_download_ascii(start_year=datetime.datetime.utcnow().year, wind_path=wind_path) 
else:
    print('Wind data NOT downloaded, turn on switch')  

    
    
    
#filewin="wind_2018_now_heeq.p" 
#start=datetime.datetime(2022, 12, 1)
#start=datetime.datetime(2022, 12, 1)
#end=datetime.datetime.utcnow()

#hd.save_wind_data_ascii(data_path,filewin,start,end,coord='HEEQ')
#[win,winh]=pickle.load(open(data_path+filewin, "rb"))
#start=win.time[-1]-datetime.timedelta(days=365)
#end=datetime.datetime.utcnow()         
#hp.plot_insitu_update(win, start, end,'Wind',plot_path+'wind/',now=True)


#on linux
#print('download new Wind data files without overwriting existing files')

#on mac for testing
#wind_data_path='/Users/chris/python/data/wind/wind_mfi_k0'
#os.system('curl -nc --directory-prefix='+wind_data_path+' "ftps://spdf.gsfc.nasa.gov/pub/data/wind/mfi/mfi_k0/2020/*.cdf"')



#wind_data_path='/perm/aswo/data/wind/wind_mfi_k0'
#print(wind_data_path)
#os.system('wget -nc --directory-prefix='+wind_data_path+' "ftps://spdf.gsfc.nasa.gov/pub/data/wind/mfi/mfi_k0/2020/*.cdf"')
#wind_data_path='/nas/helio/data/heliosat/data/wind_swe_h1'
#os.system('wget -nc --directory-prefix='+wind_data_path+' "ftps://spdf.gsfc.nasa.gov/pub/data/wind/swe/swe_h1/2020/*.cdf"')

#filewin="wind_2018_now_gse.p" 
##xstart=datetime.datetime(2018, 1, 1)
#end=datetime.datetime.utcnow()
#if get_new_data: hd.save_wind_data(data_path,filewin,start,end,heeq=False)
#[win,hwin]=pickle.load(open(data_path+filewin, "rb" ) )  

#filewin="wind_2018_now_heeq.p" 
#start=datetime.datetime(2018, 1, 1)
#end=datetime.datetime.utcnow()
#if get_new_data: hd.save_wind_data(data_path,filewin,start,end,heeq=True)

#start=win.time[-1]-datetime.timedelta(days=100)
#end=datetime.datetime.utcnow()         
#hp.plot_insitu_update(win, start, end,'Wind',plot_path,now=True)


# ### BepiColombo

# In[10]:


if debug_mode > 0: 
    importlib.reload(hd) 
    importlib.reload(hp) 

print(' ')
print('------ download BepiColombo data ')

#if get_bepi >0: 
#    hd.download_stereoa_science_merged()


# In[10]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ### STEREO-A science data

# In[19]:


if debug_mode > 0: 
    importlib.reload(hd) 
    importlib.reload(hp) 

print(' ')
print('------ download STEREO-A science data ')

if get_stereoa >0: 
    hd.download_stereoa_science_merged()


# ### Parker Solar Probe
# 

# In[5]:


if debug_mode > 0: 
    importlib.reload(hd) 
    importlib.reload(hp) 

print(' ')
print('------ download PSP data ')

#if get_psp >0: 
#    hd.download_stereoa_science_merged()


################################# PSP data update

################################## USE THIS ################################
# old
# load PSP data from server on linux command line onto leo server
# go to heliosat directory /nas/helio/data/heliosat/data/psp_fields_l2
# wget -nc "ftps://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_1min/2019/*.cdf"
# wget -nc "ftps://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_1min/2020/*.cdf"
#  /nas/helio/data/heliosat/data/psp_spc_l3
# wget -nc "ftps://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spc/l3/l3i/2019/*.cdf"
# wget -nc "ftps://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spc/l3/l3i/2020/*.cdf"
#new
#use beginning of icmecat.ipynb
############################################################################

# print('load PSP data') #from heliosat, converted to SCEQ similar to STEREO-A/B
#set time there
#change date in hd.save_psp_data
#filepsp='psp_2018_2020_nov_rtn.p'
#hd.save_psp_data(data_path,filepsp, sceq=False)   

#filepsp='psp_2018_2020_sceq.p'
#hd.save_psp_data(data_path,filepsp, sceq=True)   z


# ### Solar Orbiter
# 

# In[15]:


if debug_mode > 0: 
    importlib.reload(hd) 
    importlib.reload(hp) 

print(' ')
print('------ download Solar Orbiter data ')

#if get_solo >0: 
#    hd.download_stereoa_science_merged()


# ### BepiColombo

# In[10]:


if debug_mode > 0: 
    importlib.reload(hd) 
    importlib.reload(hp) 

print(' ')
print('------ download BepiColombo data ')

#if get_bepi >0: 
#    hd.download_stereoa_science_merged()


# In[9]:


#sys.exit()


# In[ ]:





# In[16]:


############### write header file for daily updates
#text = open('/nas/helio/data/insitu_python/data_update_headers.txt', 'w')
#text.write('Contains headers for the data files which are updated in real time.'+'\n \n')
#text.write('File creation date:  '+datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M") +' \n \n')


#text.write('NOAA real time solar wind: '+filenoaa+'\n \n'+ hnoaa+' \n \n')
#text.write('load with: >> [noaa,hnoaa]=pickle.load(open("'+data_path+filenoaa+'", "rb"))') 
#text.write(' \n \n \n \n')

#text.write('STEREO-A beacon: '+filesta_sceq+'\n \n'+ hsta+' \n \n')
#text.write('load with: >> [sta,hsta]=pickle.load(open("'+data_path+filesta+'", "rb"))') 
#text.write(' \n \n \n \n')

#text.write('Wind: '+filewin+'\n \n'+ hwin+' \n \n')
#text.write('load with: >> [win,hwin]=pickle.load(open("'+data_path+filewin+'", "rb" ))') 
#text.write(' \n \n \n \n')


#text.write('OMNI2: '+fileomni+'\n \n'+ ho+' \n \n')
#text.write('load with: >> [o,ho]=pickle.load(open("'+data_path+fileomni+'", "rb" ))') 
#text.write(' \n \n \n \n')

#text.close()


# In[ ]:





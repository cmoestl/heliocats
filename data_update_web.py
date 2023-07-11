#!/usr/bin/env python
# coding: utf-8

# ## data updates for the helioforecast.space website
# 
# Main author: C. Möstl, Austrian Space Weather Office, GeoSphere Austria
# 
# https://github.com/cmoestl/heliocats
# 
# uses environment 'envs/env_helio4.yml'

# In[1]:


# https://github.com/cmoestl/heliocats  data_update_web.py

# for updating data every day on the servers

# MIT LICENSE
# Copyright 2020-2023, Christian Moestl 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this 
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify, 
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
# permit persons to whom the Software is furnished to do so, subject to the following 
# conditions:
# The above copyright notice and this permission notice shall be included in all copies 
# or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

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
import h5py
import pytz

#import 
from heliocats import data as hd
importlib.reload(hd) #reload again while debugging

from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging

#for server
matplotlib.use('Agg')
#for mac
#%matplotlib inline

#Convert this notebook to a script with:
import os
os.system('jupyter nbconvert --to script data_update_web.ipynb')   


# ### Configure paths
# 

# In[2]:


from config import data_path

print(data_path)

plot_path=data_path+'plots/'
position_path=data_path+'plots_positions/'
sun_path=data_path+'plots_sun/'

print(plot_path)
print(position_path)

from config import noaa_path
print(noaa_path)

from config import data_path_ml
print(data_path_ml)


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


# spacecraft positions image
hp.plot_positions(datetime.datetime.utcnow(),position_path, 'HEEQ',now=True)

# get current SDO images 
hd.get_sdo_realtime_image(sun_path)


# ### OMNI2 data
# 

# In[12]:


get_omni=1

# OMNI2
fileomni="omni_1963_now.p"
if get_omni: hd.save_omni_data(data_path,fileomni)
[o,ho]=pickle.load(open(data_path+fileomni, "rb" ) )  

start=datetime.datetime.utcnow() - datetime.timedelta(days=365)
end=datetime.datetime.utcnow() 
hp.plot_insitu_update(o, start, end,'OMNI2',plot_path+'omni2/',now=True)


# ### NOAA real time solar wind

# In[17]:


get_noaa=0


if get_noaa > 0:
    print('download NOAA real time solar wind plasma and mag')
    datestr=str(datetime.datetime.utcnow().strftime("%Y_%b_%d_%H_%M"))
    print(datestr+' UTC')

    plasma='http://services.swpc.noaa.gov/products/solar-wind/plasma-7-day.json'
    mag='http://services.swpc.noaa.gov/products/solar-wind/mag-7-day.json'

    try: urllib.request.urlretrieve(plasma, noaa_path+'plasma-7-day_'+datestr+'.json')
    except urllib.error.URLError as e:
      print(' ', plasma,' ',e.reason)

    try: urllib.request.urlretrieve(mag, noaa_path+'mag-7-day_'+datestr+'.json')
    except urllib.error.URLError as e:
      print(' ', mag,' ',e.reason)
  
    print('NOAA download complete')

save_noaa=1
    
filenoaa='noaa_rtsw_jan_2023_now.p'
if save_noaa > 0: hd.save_noaa_rtsw_data(data_path,noaa_path,filenoaa)
[noaa,hnoaa]=pickle.load(open(data_path+filenoaa, "rb" ) ) 

    


# In[15]:


sys.exit()


# ### Wind data

# In[10]:


#from heliocats import data as hd
#importlib.reload(hd) #reload again while debugging

#from heliocats import plot as hp
#importlib.reload(hp) #reload again while debugging

#get_wind=1

#filewin="wind_2018_now_heeq.p" 
#start=datetime.datetime(2022, 12, 1)
#start=datetime.datetime(2022, 12, 1)
#end=datetime.datetime.utcnow()

###!!use urllib.request.urlretrieve(plasma, noaa_path+'plasma-7-day_'+datestr+'.json')

#if get_wind: 
    #hd.wind_download_ascii()    
    #hd.save_wind_data_ascii(data_path,filewin,start,end,coord='HEEQ')

#[win,winh]=pickle.load(open(data_path+filewin, "rb"))
    
#start=win.time[-1]-datetime.timedelta(days=365)
#end=datetime.datetime.utcnow()         
#hp.plot_insitu_update(win, start, end,'Wind',plot_path+'wind/',now=True)











# ### STEREO-A beacon data

# In[ ]:





# In[ ]:





# In[ ]:


sys.exit()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[19]:












# In[58]:


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


# In[ ]:


filewin="wind_2018_now_gse.p" 
start=datetime.datetime(2018, 1, 1)
end=datetime.datetime.utcnow()
if get_new_data: hd.save_wind_data(data_path,filewin,start,end,heeq=False)
[win,hwin]=pickle.load(open(data_path+filewin, "rb" ) )  

filewin="wind_2018_now_heeq.p" 
start=datetime.datetime(2018, 1, 1)
end=datetime.datetime.utcnow()
if get_new_data: hd.save_wind_data(data_path,filewin,start,end,heeq=True)

start=win.time[-1]-datetime.timedelta(days=100)
end=datetime.datetime.utcnow()         
hp.plot_insitu_update(win, start, end,'Wind',plot_path,now=True)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


sys.exit()


# In[ ]:


############### write header file for daily updates
text = open('/nas/helio/data/insitu_python/data_update_headers.txt', 'w')
text.write('Contains headers for the data files which are updated in real time.'+'\n \n')
text.write('File creation date:  '+datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M") +' \n \n')


text.write('NOAA real time solar wind: '+filenoaa+'\n \n'+ hnoaa+' \n \n')
text.write('load with: >> [noaa,hnoaa]=pickle.load(open("'+data_path+filenoaa+'", "rb"))') 
text.write(' \n \n \n \n')

text.write('STEREO-A beacon: '+filesta_sceq+'\n \n'+ hsta+' \n \n')
text.write('load with: >> [sta,hsta]=pickle.load(open("'+data_path+filesta+'", "rb"))') 
text.write(' \n \n \n \n')

text.write('Wind: '+filewin+'\n \n'+ hwin+' \n \n')
text.write('load with: >> [win,hwin]=pickle.load(open("'+data_path+filewin+'", "rb" ))') 
text.write(' \n \n \n \n')


text.write('OMNI2: '+fileomni+'\n \n'+ ho+' \n \n')
text.write('load with: >> [o,ho]=pickle.load(open("'+data_path+fileomni+'", "rb" ))') 
text.write(' \n \n \n \n')

text.close()


# In[4]:


#for easier debugging - do not download and process data but do everything else
#get_new_data=True


#~/miniconda/envs/helio/bin/python /home/cmoestl/pycode/heliocats/sc_positions_insitu_orbit1.py
#~/miniconda/envs/helio/bin/python /home/cmoestl/pycode/heliocats/sc_positions_insitu_orbit2.py


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
#hd.save_psp_data(data_path,filepsp, sceq=True)   



############################################ STEREO - A science data update


# 
# 
# -----------------------------
# data files:
# 
# sftp.put(path+'Wind_now.png')  # upload file to public/ on remote
# sftp.put(path+'STEREO-A_beacon_14_days_now.png')
# sftp.put(path+'OMNI2_now.png')
# 
# 
# sftp.put(path+'NOAA_RTSW_PREDSTORM_55days_now.png')
# sftp.put(path+'NOAA_RTSW_PREDSTORM_14days_now.png')
# sftp.put(path+'NOAA_RTSW_PREDSTORM_3days_now.png')
# sftp.put(path+'OMNI2_and_NOAA_RTSW_now.png')
# 
# 
# ----------------------------

# In[ ]:





# In[ ]:





# In[ ]:


################## NOAA real time
################ load jsons

print('download NOAA real time solar wind plasma and mag')
datestr=str(datetime.datetime.utcnow().strftime("%Y_%b_%d_%H_%M"))
print(datestr+' UTC')

plasma='http://services.swpc.noaa.gov/products/solar-wind/plasma-7-day.json'
mag='http://services.swpc.noaa.gov/products/solar-wind/mag-7-day.json'

try: urllib.request.urlretrieve(plasma, noaa_path+'plasma-7-day_'+datestr+'.json')
except urllib.error.URLError as e:
  print(' ', plasma,' ',e.reason)

try: urllib.request.urlretrieve(mag, noaa_path+'mag-7-day_'+datestr+'.json')
except urllib.error.URLError as e:
  print(' ', mag,' ',e.reason)
  
print()


# my version:
#filenoaa='noaa_rtsw_jan_2020_now.p'
#if get_new_data: hd.save_noaa_rtsw_data(data_path,noaa_path,filenoaa)
#[noaa,hnoaa]=pickle.load(open(data_path+filenoaa, "rb" ) ) 

#################### with Rachel Bailey's data from predstorm
noaapath='/nas/helio/realcode/real/predstorm/data/rtsw_min_last100days.h5'
filenoaa='noaa_rtsw_last100days_now.p'
if get_new_data: hd.save_noaa_rtsw_data_predstorm(data_path,noaapath,filenoaa)        
[noaa,hnoaa]=pickle.load(open(data_path+filenoaa, "rb" ) ) 
    
              
############### make 3 plots

#for ICMEs
start=noaa.time[-1]-datetime.timedelta(days=3)
end=datetime.datetime.utcnow() #noaa.time[-1]     
hp.plot_insitu_update(noaa, start, end,'NOAA_RTSW_PREDSTORM_3days',plot_path,now=True)

start=noaa.time[-1]-datetime.timedelta(days=14)
end=datetime.datetime.utcnow() #noaa.time[-1]     
hp.plot_insitu_update(noaa, start, end,'NOAA_RTSW_PREDSTORM_14days',plot_path,now=True)

start=noaa.time[-1]-datetime.timedelta(days=55)
end=datetime.datetime.utcnow() #noaa.time[-1]     
hp.plot_insitu_update(noaa, start, end,'NOAA_RTSW_PREDSTORM_55days',plot_path,now=True)

print()

#make quick new plot for ICMEs
# hd.save_noaa_rtsw_data(data_path,noaa_path,filenoaa)
# start=noaa.time[-1]-datetime.timedelta(days=3)
# end=datetime.datetime.utcnow() #noaa.time[-1]     
# hp.plot_insitu(noaa, start, end,'NOAA_RTSW','/home/cmoestl/pycode/heliocats')





####################################### OMNI2
fileomni="omni_1963_now.p"
if get_new_data: hd.save_omni_data(data_path,fileomni)
[o,ho]=pickle.load(open(data_path+fileomni, "rb" ) )  

start=datetime.datetime.utcnow() - datetime.timedelta(days=365)
end=datetime.datetime.utcnow() 
hp.plot_insitu_update(o, start, end,'OMNI2',plot_path,now=True)


########################## add NOAA RTSW to OMNI data and make combined plot

#get index of last OMNI data
last_omni_index=np.where(np.isfinite(o.bt) == True)[0][-1]
#get time for this index
last_omni_time=o.time[last_omni_index]
#add utc timezone awareness
last_omni_time_utc=last_omni_time.astimezone(tz=datetime.timezone.utc)
#get index where omni ends in noaa data
noaa_omni_end=np.where(noaa.time>last_omni_time_utc)[0][0]


#length of NOAA data
size_noaa=np.size(noaa)-noaa_omni_end

combi_omni_noaa=np.zeros(last_omni_index+size_noaa,dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   

#convert to recarray
combi_omni_noaa = combi_omni_noaa.view(np.recarray)  
#stack omni and noaa data
combi_omni_noaa.time=np.hstack((o.time[0:last_omni_index],noaa.time[noaa_omni_end-1:-1]))
combi_omni_noaa.bx=np.hstack((o.bx[0:last_omni_index],noaa.bx[noaa_omni_end-1:-1]))
combi_omni_noaa.by=np.hstack((o.by[0:last_omni_index],noaa.by[noaa_omni_end-1:-1]))
combi_omni_noaa.bz=np.hstack((o.bz[0:last_omni_index],noaa.bz[noaa_omni_end-1:-1]))
combi_omni_noaa.bt=np.hstack((o.bt[0:last_omni_index],noaa.bt[noaa_omni_end-1:-1]))
combi_omni_noaa.vt=np.hstack((o.vt[0:last_omni_index],noaa.vt[noaa_omni_end-1:-1]))
combi_omni_noaa.np=np.hstack((o.np[0:last_omni_index],noaa.np[noaa_omni_end-1:-1]))
combi_omni_noaa.tp=np.hstack((o.tp[0:last_omni_index],noaa.tp[noaa_omni_end-1:-1]))

print('Merging NOAA OMNI done')

start=datetime.datetime.utcnow() -datetime.timedelta(days=365)
end=datetime.datetime.utcnow() 
hp.plot_insitu_update(combi_omni_noaa, start, end,'OMNI2_and_NOAA_RTSW',plot_path,now=True)




#for making quick plots
# plot_path='/home/cmoestl/pycode/heliocats/'
# fileomni="omni_1963_now.p"
# [o,ho]=pickle.load(open(data_path+fileomni, "rb" ) )  
# start=datetime.datetime.utcnow() -datetime.timedelta(days=365)
# end=datetime.datetime.utcnow() 
# hp.plot_insitu_update(o, start, end,'OMNI2',plot_path,now=True)





################################### STEREO-A

from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging



#start=datetime.datetime(2020, 8, 1)
#end=datetime.datetime.utcnow()
#filesta_sceq="stereoa_2020_august_now_sceq_beacon.p" 
#filesta_rtn="stereoa_2020_august_now_rtn_beacon.p" 

#if get_new_data: 
#    hd.save_stereoa_beacon_data(data_path,filesta_sceq,start,end,sceq=True)
#    hd.save_stereoa_beacon_data(data_path,filesta_rtn,start,end,sceq=False)


start=datetime.datetime(2020, 1, 1)
end=datetime.datetime.utcnow()
filesta_sceq="stereoa_2020_now_sceq_beacon.p" 
filesta_rtn="stereoa_2020_now_rtn_beacon.p" 

if get_new_data: 
    hd.save_stereoa_beacon_data(data_path,filesta_sceq,start,end,sceq=True)
    hd.save_stereoa_beacon_data(data_path,filesta_rtn,start,end,sceq=False)

[sta,hsta]=pickle.load(open(data_path+filesta_rtn, "rb" ) ) 

start=sta.time[-1]-datetime.timedelta(days=14)
end=datetime.datetime.utcnow()     
hp.plot_insitu_update_stereoa_beacon(sta, start, end,'STEREO-A_beacon_14_days',plot_path,now=True)




####################################### Wind

print('download Wind 2020 files without overwriting existing files into ',data_path)
wind_data_path='/nas/helio/data/heliosat/data/wind_mfi_k0'
os.system('wget -nc --directory-prefix='+wind_data_path+' "ftps://spdf.gsfc.nasa.gov/pub/data/wind/mfi/mfi_k0/2020/*.cdf"')
wind_data_path='/nas/helio/data/heliosat/data/wind_swe_h1'
os.system('wget -nc --directory-prefix='+wind_data_path+' "ftps://spdf.gsfc.nasa.gov/pub/data/wind/swe/swe_h1/2020/*.cdf"')


filewin="wind_2018_now_gse.p" 
start=datetime.datetime(2018, 1, 1)
end=datetime.datetime.utcnow()
if get_new_data: hd.save_wind_data(data_path,filewin,start,end,heeq=False)
[win,hwin]=pickle.load(open(data_path+filewin, "rb" ) )  

filewin="wind_2018_now_heeq.p" 
start=datetime.datetime(2018, 1, 1)
end=datetime.datetime.utcnow()
if get_new_data: hd.save_wind_data(data_path,filewin,start,end,heeq=True)

start=win.time[-1]-datetime.timedelta(days=100)
end=datetime.datetime.utcnow()         
hp.plot_insitu_update(win, start, end,'Wind',plot_path,now=True)



############### write header file for daily updates
text = open('/nas/helio/data/insitu_python/data_update_headers.txt', 'w')
text.write('Contains headers for the data files which are updated in real time.'+'\n \n')
text.write('File creation date:  '+datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M") +' \n \n')


text.write('NOAA real time solar wind: '+filenoaa+'\n \n'+ hnoaa+' \n \n')
text.write('load with: >> [noaa,hnoaa]=pickle.load(open("'+data_path+filenoaa+'", "rb"))') 
text.write(' \n \n \n \n')

text.write('STEREO-A beacon: '+filesta_sceq+'\n \n'+ hsta+' \n \n')
text.write('load with: >> [sta,hsta]=pickle.load(open("'+data_path+filesta+'", "rb"))') 
text.write(' \n \n \n \n')

text.write('Wind: '+filewin+'\n \n'+ hwin+' \n \n')
text.write('load with: >> [win,hwin]=pickle.load(open("'+data_path+filewin+'", "rb" ))') 
text.write(' \n \n \n \n')


text.write('OMNI2: '+fileomni+'\n \n'+ ho+' \n \n')
text.write('load with: >> [o,ho]=pickle.load(open("'+data_path+fileomni+'", "rb" ))') 
text.write(' \n \n \n \n')

text.close()


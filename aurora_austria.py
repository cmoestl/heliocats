#!/usr/bin/env python
# coding: utf-8

# ## Aurora days in Austria
# 
# astropy for sunset and sunrise in central point of Austria:
# Bad Aussee
# 47,6964° N , 13,3458° E
# 
# - get Dst from OMNI
# count days when Dst < - 200 nT and time between sunset+ 1 hour and sunrise -1 hour
# 
# - uses environment helio5, see /envs/env_helio5.yml in the heliocats package
# 

# In[3]:


from astropy.coordinates import EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import get_sun
from datetime import datetime, timedelta

#import astroplan #sunset sunrise can also be done with astroplan

import numpy as np
import pickle
import sys
import os


#from heliocats import data as hd
#from heliocats import plot as hp

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

if sys.platform == 'linux': 
    
    from config_server import data_path
    
if sys.platform =='darwin':  

    from config_local import data_path
################################################ CHECK  ##############################################

#make sure to convert the current notebook to a script
os.system('jupyter nbconvert --to script aurora_austria.ipynb')   


####################################################################################################################

#test execution times
#t0all = time.time()

fileomni="omni_1963_now.p"
[o,ho]=pickle.load(open(data_path+fileomni, "rb" ) )  
start=datetime.utcnow() - timedelta(days=365)
end=datetime.utcnow() 

#hp.plot_insitu_update(o, start, end,'OMNI2',plot_path+'omni2/',now=True)#compare with cloud statistics in Austria ?


# In[94]:


# Define the geographic location (latitude, longitude, and elevation)
latitude = 47.70  # 
longitude = -13.34  # 
elevation = 100  # Elevation in meters (optional)

# Create an EarthLocation object
location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=elevation*u.m)
print(location)

#make hourly times from 1963 onwards

time = Time(datetime(1983,3,14,7))
print(time)

#gets the altitude of the Sun at our location
sun_altaz = get_sun(time).transform_to(AltAz(obstime=time, location=location))

#night is when Sun < 0 altitude?

print(sun_altaz)

#get all these hours


# In[95]:


# Get the altitude and azimuth of the Sun
sun_altaz = get_sun(time).transform_to(AltAz(obstime=time, location=location))

print(sun_altaz.alt.value)


# In[91]:


#compare to Dst values

#compare to cloud statistics for a given day in Austria (ask meteorologists)


#make report, put on social media


# In[ ]:





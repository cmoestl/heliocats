#!/usr/bin/env python
# coding: utf-8

# ## arrcat
# 
# makes the HELCATS HI ARRIVAL catalog
# 
# Author: C. Möstl, IWF Graz, Austria
# twitter @chrisoutofspace, part of https://github.com/cmoestl/heliocats
# 
# last update April 2020
# 
# Install a specific conda environment to run this code, see readme at https://github.com/cmoestl/heliocats
# 
# 
# Convert this notebook to a script with jupyter nbconvert --to script arrcat.ipynb
# 
# **current status: work in progress**
# 
# features to be added: 
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
import heliopy

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

catdir='arrcat'
if os.path.isdir(catdir) == False: os.mkdir(catdir)

icplotsdir='arrcat/plots_arrcat/' 
if os.path.isdir(icplotsdir) == False: os.mkdir(icplotsdir) 

#Convert this notebook to a script with jupyter nbconvert --to script icmecat.ipynb
os.system('jupyter nbconvert --to script arrcat.ipynb')    


# ## Make HI SSEF30 arrival catalog

# In[2]:


from heliocats import cats as hc
importlib.reload(hc) #reload again while debugging 


#https://www.helcats-fp7.eu/

#LOAD HELCATS HIGeoCAT
url_higeocat='https://www.helcats-fp7.eu/catalogues/data/HCME_WP3_V06.vot'

try: urllib.request.urlretrieve(url_higeocat,'data/HCME_WP3_V06.vot')
except urllib.error.URLError as e:
    print('higeocat not loaded')

higeocat=hc.load_higeocat_vot('data/HCME_WP3_V06.vot')
higeocat_time=parse_time(higeocat['Date']).datetime    
higeocat_t0=parse_time(higeocat['SSE Launch']).datetime   #backprojected launch time


#define empty pandas dataframe for arrival catalog with column names

ac = pd.DataFrame([], columns = ['id', 'sc','sse_launch_time','target_name','sse_heeq_long',                                 'target_delta','sse_speed','target_speed','target_arrival_time',                                 'target_distance','target_heeq_lat','target_heeq_lon',                                 'target_pa','pa_fit','pa_n','pa_s','pa_center'])



ac=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, 'PSP')
ac=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac,'Solo')
ac=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac,'Bepi')
   


#Make arrival catalog from HIGEOCAT
ac=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, 'STA')
ac=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, 'STB')
ac=hc.make_arrival_catalog_insitu_ssef30(higeocat,ac, 'Earth')



ac=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac,'Mercury')
ac=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac,'Venus')
ac=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac,'Mars')


ac = ac.reset_index(drop=True)

# #ac.reset_index(drop=True,inplace=True)

ac


# ## save ARRCAT
# 

# #### save header

# In[101]:


#save header and parameters as text file and prepare for html website
header='ARRIVAL CATALOGUE v2.0 \n\nThis ARRival CATalog (ARRCAT) models the arrivals of CMEs tracked in the STEREO heliospheric imagers EU HELCATS project (2014-2017). \nIt lists predicted arrivals of solar coronal mass ejections at various spacecraft and planets with the STEREO heliospheric imager \ninstruments, between April 2007 - April 2020, based on the HIGeoCAT catalog of CMEs established at RAL Space, UK (Harrison, Davies, Barnes). \nThis is version 2.0, released 2020-**-**. DOI: 10.6084/m9.figshare.6356420 \n\nBased on HIGeoCAT version ***\nThe catalog is available as  python pandas dataframe (pickle), python numpy structured array (pickle), json, csv, xlsx, txt, hdf5, at \nhttps://helioforecast.space/arrcat \nNumber of events in ARRCAT: '+str(len(ac))+'\nTargets: Earth L1, STEREO-A, STEREO-B,  Solar Orbiter, Parker Solar Probe, Bepi Colombo, Venus, Mercury, Mars.\n\nParameters: \n    1: id: From HIGeoCAT, the unique identifier for the observed CME.\n    2: sc: From HIGeoCAT, the HI observing STEREO spacecraft, (A=Ahead or B=Behind)\n    3: sse_launch_time: From HIGeoCAT, launch time of the CME on the Sun, unit: UTC.\n    4: target_name: Name of in situ target.\n    5: sse_heeq_long: From HIGeoCAT, the HEEQ longitude of the CME apex propagation direction, unit: degree.\n    6: target_detla: Difference in HEEQ longitude between central CME direction and target location, positive values: spacecraft is west of CME apex. unit: degree.\n    7: sse_speed: From HIGeoCAT, speed of CME apex, unit: km/s.\n    8: target_speed: CME arrival speed at target location, corrected for SSE shape. unit: km/s.\n    9: target_arrival_time: CME arrival time at target location, corrected for SSE shape. unit: UTC.\n    10: target_distance: Target distance from Sun, at CME launch time. unit: AU.\n    11: target_heeq_lat: Target latitude in HEEQ, at CME launch time. unit: degree.\n    12: target_heeq_lon: Target longitude in HEEQ, at CME launch time. unit: degree.\n    13: target_pa: PA of target from HI observing STEREO spacecraft, unit: degree.\n    14: pa_fit: From HICAT, PA along which time-elongation profile is extracted, unit: degree.\n    15: pa_n: From HICAT, northern position angle of CME, unit: degree.\n    16: pa_s: From HICAT, southernmost position angle of CME, unit: degree.\n    17: pa_center: average of pa_n and pa_s, unit: degree.\n\nNotes:\n\n- We have applied modified method from Möstl & Davies (2013, Solar Physics) for calculating speeds \nand arrival times of the CMEs modeled with SSEF30 to all CMEs in the HELCATS HIGeoCAT catalog \n(see website helcats-fp7.eu, and Möstl et al. 2014, ApJ, for more details). \nIf the SSEF30 circle hits a spacecraft or planet, an entry in ARRCAT is produced. \n\n- A new iterative method was used ***\nReferences \nMöstl et al. (2017),  https://doi.org/************* \n'



print(header)

#make header file
file='icmecat/HELCATS_ICMECAT_v20_header.txt'
with open(file, "w") as text_file:
    text_file.write(header)
print()    
print('header saved as '+file)
print()    

#Convert to html regarding line breaks, paragraph beginning and spaces
header_spaces=header.replace(" ", "&nbsp;")
header_html= "<p>" +header_spaces.replace('\n', '<br>')+ "</p>" 
print('header converted to HTML')
print()    
print()    


# ### 4b save into different formats

# In[109]:


########## python formats

# save ICMECAT as pandas dataframe with times as datetime objects as pickle
file='arrcat/HELCATS_ARRCAT_v20_pandas.p'
pickle.dump([ac,header], open(file, 'wb'))
print('ARRCAT saved as '+file)

# # save ICMECAT as numpy array with times as matplotlib datetime as pickle
# ac_num=copy.deepcopy(ac) 
# ac_num.icme_start_time=parse_time(ac_num.sse_launch_time).plot_date
# ac_num.mo_start_time=parse_time(ac_num.target_arrival_time).plot_date
# #convert to recarray
# ac_num_rec=ac_num.to_records()
# #create structured array ******************
# dtype1=[('index','i8'),('icmecat_id', '<U30'),('sc_insitu', '<U20')] +[(i, '<f8') for i in ac.keys()[2:len(ac.keys())]]
# ac_num_struct=np.array(ic_num_rec,dtype=dtype1)

# file='arrcat/HELCATS_ARRCAT_v20_numpy.p'
# pickle.dump([ac_num,ic_num_struct,header], open(file, 'wb'))
# print('ARRCAT saved as '+file)



################ save to different formats

#copy pandas dataframe first to change time format consistent with HELCATS
ac_copy=copy.deepcopy(ac)  
ac_copy.sse_launch_time=parse_time(ac.sse_launch_time).isot
ac_copy.target_arrival_time=parse_time(ac.target_arrival_time).isot


#change time format
for i in np.arange(len(ac)):

    dum=ac_copy.sse_launch_time[i] 
    ac_copy.at[i,'sse_launch_time']=dum[0:16]+'Z'
     
    dum=ac_copy.target_arrival_time[i] 
    ac_copy.at[i,'target_arrival_time']=dum[0:16]+'Z'


#save as Excel
file='arrcat/HELCATS_ARRCAT_v20.xlsx'
ac_copy.to_excel(file,sheet_name='acMECATv2.0')
print('ARRCAT saved as '+file)

#save as json
file='arrcat/HELCATS_ARRCAT_v20.json'
ac_copy.to_json(file)
print('ARRCAT saved as '+file)

#save as csv
file='arrcat/HELCATS_ARRCAT_v20.csv'
ac_copy.to_csv(file)
print('ARRCAT saved as '+file)


#save as txt
file='arrcat/HELCATS_ARRCAT_v20.txt'
np.savetxt(file, ac_copy.values.astype(str), fmt='%s' )
print('ARRCAT saved as '+file)



sys.exit()








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
ic_copy2.sse_launch_time=parse_time(ac.sse_launch_time).iso
ic_copy2.mo_start_time=parse_time(ac.target_arrival_time).iso

#change time format
for i in np.arange(len(ic)):

    dum=ic_copy2.icme_start_time[i] 
    ic_copy2.at[i,'icme_start_time']=dum[0:16]
     
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

# In[ ]:


#load icmecat as pandas dataframe
file='icmecat/HELCATS_ICMECAT_v20_pandas.p'
[ic_pandas,h,p]=pickle.load( open(file, 'rb'))   

#load icmecat as numpy array
file='icmecat/HELCATS_ICMECAT_v20_numpy.p'
[ic_nprec,ic_np,h,p]=pickle.load( open(file, 'rb'))   


# In[ ]:


ic_pandas
ic_pandas.keys()


# In[ ]:


ic_nprec


# In[ ]:


ic_nprec


# In[ ]:


ic_nprec.icmecat_id


# In[ ]:





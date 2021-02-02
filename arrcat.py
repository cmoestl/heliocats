#!/usr/bin/env python
# coding: utf-8

# ## arrcat
# 
# makes the HELCATS HI ARRIVAL catalog
# 
# Author: C. Möstl, IWF Graz, Austria
# twitter @chrisoutofspace, part of https://github.com/cmoestl/heliocats
# 
# current version ARRCAT 2.0, released 2020 May 13, updated 2021 January 12
# 
# Install a specific conda environment to run this code, see readme at https://github.com/cmoestl/heliocats
# 
# Convert this notebook to a script with "jupyter nbconvert --to script arrcat.ipynb", automatically done in first cell
# 
# https://figshare.com/articles/dataset/Arrival_catalog_of_coronal_mass_ejections_observed_with_STEREO_Heliospheric_Imagers_HELCATS_ARRCAT_2_0_/12271292
# 

# In[ ]:


last_update='2021 January 12'


# In[ ]:


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


# ## 1 Make HI SSEF30 arrival catalog ARRCAT

# In[ ]:


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

column_list=['id', 'sc','target_name','sse_launch_time','target_arrival_time','target_arrival_time_err',                                 'target_distance','target_heeq_lon','target_heeq_lat',                                 'target_delta','sse_heeq_lon','sse_heeq_lat',                                 'sse_speed','sse_speed_err','target_speed','target_speed_err',                                 'pa_fit','pa_n','pa_s','pa_center']
             
#pandas dataframe for current version with iteration in calculating arrival time
ac = pd.DataFrame([], columns = column_list)

#pandas dataframe for old method of arrival time prediction without iteration
ac_old = pd.DataFrame([], columns = column_list)

#Make arrival catalog from HIGEOCAT for each spacecraft

[ac,ac_old]=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, ac_old, 'PSP',column_list)
[ac,ac_old]=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, ac_old, 'Solo',column_list)
[ac,ac_old]=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, ac_old, 'Bepi',column_list)
   
[ac,ac_old]=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, ac_old, 'STA',column_list)
[ac,ac_old]=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, ac_old, 'STB',column_list)
[ac,ac_old]=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, ac_old, 'Earth',column_list)

[ac,ac_old]=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, ac_old, 'Ulysses',column_list)
[ac,ac_old]=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, ac_old, 'Mercury',column_list)
[ac,ac_old]=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, ac_old, 'Venus',column_list)
[ac,ac_old]=hc.make_arrival_catalog_insitu_ssef30(higeocat, ac, ac_old, 'Mars',column_list)



ac = ac.sort_values(by='target_arrival_time',ascending=False)
ac = ac.reset_index(drop=True)

ac_old = ac_old.sort_values(by='target_arrival_time',ascending=False)
ac_old = ac_old.reset_index(drop=True)




ac


# ## 2 save ARRCAT
# 

# #### save header

# In[ ]:


#save header and parameters as text file and prepare for html website
header='ARRIVAL CATALOGUE 2.0 \n\nIn this ARRival CATalog (ARRCAT), the arrivals of solar coronal mass ejections that were \ntracked in the STEREO heliospheric imagers in the HELCATS project are calculated.\nARRCAT lists modeled arrivals of CMEs at various spacecraft and planets, \nbased on the HIGeoCAT catalog of CMEs established at RAL Space, UK (D. Barnes, J. A. Davies, R. A. Harrison). \n\nThis is version 2.0, released 2020-May-13, updated 2020-Jan-12. https://doi.org/10.6084/m9.figshare.12271292\n\nIt is based on this HIGeoCAT version: '+url_higeocat+' using the SSEF30 model.\n\nThe catalog is available as python pandas dataframe (pickle), \npython numpy arrays (pickle, as recarray and structured array), \nnpy (numpy, no pickle), json, csv, xlsx, txt, hdf5, html at \nhttps://helioforecast.space/arrcat and https://www.helcats-fp7.eu\n\nNumber of events in ARRCAT: '+str(len(ac))+'\nTargets: Earth_L1, STEREO-A, STEREO-B, Solar Orbiter, Parker Solar Probe (PSP), Bepi Colombo, Ulysses, Venus, Mercury, Mars.\n\nAuthors: Christian Moestl, D. Barnes, R. A. Harrison, J. A. Davies, Andreas J. Weiss, David Stansby;\nSpace Research Institute (IWF), Austrian Academy of Sciences, Graz, Austria; RAL Space, UK; Mullard Space Science Laboratory, UK.\n\nRules: If results are produced with this catalog for peer-reviewed scientific publications,\nplease contact christian.moestl@oeaw.ac.at for co-authorships.\n\nParameters \n\n    0: id: From HIGeoCAT, the unique identifier for the observed CME.\n    1: sc: From HIGeoCAT, the HI observing STEREO spacecraft, (A=Ahead or B=Behind)\n    2: target_name: Name of in situ target.\n    3: sse_launch_time: From HIGeoCAT, launch time of the CME at Sun center, unit: UTC.\n    4: target_arrival_time: CME arrival time at target location, corrected for SSE shape (Moestl and Davies, 2013). unit: UTC.\n    5: target_arrival_time_err: Error of CME arrival time at target location, solely based on error in parameter sse_speed_err. unit: hours.\n    6: target_distance: Target distance from Sun, at CME arrival time. unit: AU.\n    7: target_heeq_lon: Target longitude in HEEQ, at CME arrival time. unit: degree.\n    8: target_heeq_lat: Target latitude in HEEQ, at CME arrival time. unit: degree.\n    9: target_delta: Difference in HEEQ longitude between central CME direction and target location at CME arrival time, \n       positive values: spacecraft is west of CME apex. unit: degree.\n    10: sse_heeq_lon: From HIGeoCAT, the HEEQ longitude of the CME apex propagation direction, unit: degree.\n    11: sse_heeq_lat: From HIGeoCAT, the HEEQ latitude of the CME apex propagation direction, unit: degree.\n    12: sse_speed: From HIGeoCAT, speed of the CME apex, unit: km/s.\n    13: sse_speed_err: From HIGeoCAT, error in the speed of the CME apex, unit: km/s.\n    14: target_speed: CME arrival speed at target location, corrected for SSE shape (Moestl and Davies, 2013). unit: km/s.\n    15: target_speed_err: Error in CME arrival speed at target location, corrected for SSE shape (Moestl and Davies, 2013). unit: km/s.\n    16: pa_fit: From HIGeoCAT, PA along which time-elongation profile is extracted, unit: degree.\n    17: pa_n: From HIGeoCAT, northern position angle of CME, unit: degree.\n    18: pa_s: From HIGeoCAT, southernmost position angle of CME, unit: degree.\n    19: pa_center: average of pa_n and pa_s, unit: degree.\n\nComments\n\n    - We have modified the calculation of CME arrival time here by a new iterative method compared to Moestl et al. (2014, 2017). \n      In the first iteration, the sse_launch_time is used to calculate the target HEEQ position. \n      In subsequent three iterations, we use the arrival time from the previous iteration (instead of sse_launch time) \n      to calculate the target HEEQ position. \n      In this way we avoid an error of a few degrees in the arrival target location (see also Moestl et al. 2017). \n    - The analytical formulas for calculating the speeds and arrival times of the CMEs modeled with SSEF30, \n      corrected for the SSEF30 circular shape, can be found in Moestl & Davies (2013). \n\n\nReferences \n\nMoestl & Davies (2013) https://doi.org/10.1007/s11207-012-9978-8 arxiv: https://arxiv.org/abs/1202.1299\nMoestl et al. (2014)   https://doi.org/10.1088/0004-637X/787/2/119 (open access) \nMoestl et al. (2017)   https://doi.org/10.1002/2017SW001614 (open access)\n'

print(header)

#make header file
file='arrcat/HELCATS_ARRCAT_v20_header.txt'
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

print('-----------------------')
print('ARRCAT sorted by time gives:')
np.sort(ac.target_arrival_time)


# #### save into different formats

# In[ ]:


########## python formats

# save ARRCAT as pandas dataframe with times as datetime objects as pickle
file='arrcat/HELCATS_ARRCAT_v20_pandas.p'
pickle.dump([ac,header], open(file, 'wb'))
print('ARRCAT saved as '+file)

# save ARRCAT as numpy array with times as matplotlib datetime as pickle
ac_num=copy.deepcopy(ac) 
ac_num.sse_launch_time=parse_time(ac_num.sse_launch_time).plot_date
ac_num.target_arrival_time=parse_time(ac_num.target_arrival_time).plot_date
#convert to recarray
ac_num_rec=ac_num.to_records()
#create structured array 
dtype1=[('index','i8'),('id', '<U30'),('sc', '<U20'),('target_name', '<U20')] +[(i, '<f8') for i in ac.keys()[3:len(ac.keys())]]
ac_num_struct=np.array(ac_num_rec,dtype=dtype1)

file='arrcat/HELCATS_ARRCAT_v20_numpy.p'
pickle.dump([ac_num,ac_num_struct,header], open(file, 'wb'))
print('ARRCAT saved as '+file)




################ save to different formats

#copy pandas dataframe first to change time format consistent with HELCATS
ac_copy=copy.deepcopy(ac)  
#change time format and add Z
ac_copy.sse_launch_time=parse_time(ac.sse_launch_time).isot
ac_copy.target_arrival_time=parse_time(ac.target_arrival_time).isot

for i in np.arange(len(ac)):

    dum=ac_copy.sse_launch_time[i] 
    ac_copy.at[i,'sse_launch_time']=dum[0:16]+'Z'
     
    dum=ac_copy.target_arrival_time[i] 
    ac_copy.at[i,'target_arrival_time']=dum[0:16]+'Z'


#save as Excel
file='arrcat/HELCATS_ARRCAT_v20.xlsx'
ac_copy.to_excel(file,sheet_name='ARRCATv2.0')
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



#**************************########### save old version for comparison

#copy pandas dataframe first to change time format consistent with HELCATS
ac_copy_old=copy.deepcopy(ac_old)  
#change time format and add Z
ac_copy_old.sse_launch_time=parse_time(ac_old.sse_launch_time).isot
ac_copy_old.target_arrival_time=parse_time(ac_old.target_arrival_time).isot

for i in np.arange(len(ac_old)):

    dum=ac_copy_old.sse_launch_time[i] 
    ac_copy_old.at[i,'sse_launch_time']=dum[0:16]+'Z'
     
    dum=ac_copy_old.target_arrival_time[i] 
    ac_copy_old.at[i,'target_arrival_time']=dum[0:16]+'Z'


#save as Excel
file='arrcat/HELCATS_ARRCAT_v10_old.xlsx'
ac_old.to_excel(file,sheet_name='ARRCATv2.0')
print('ARRCAT saved as '+file)





#########################

#########save into hdf5 format , use S for strings http://docs.h5py.org/en/stable/strings.html#what-about-numpy-s-u-type
dtype2=[('index','i8'),('id', 'S30'),('sc', 'S20'),('target_name', 'S20')] +[(i, '<f8') for i in ac.keys()[3:len(ac.keys())]]
ach5=np.array(ac_num_rec,dtype=dtype2)
file='arrcat/HELCATS_ARRCAT_v20.h5'
f=h5py.File(file,mode='w')
f["arrcat"]= ach5

#add attributes
#************************
#***********************

print('ARRCAT saved as '+file)
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
file='arrcat/HELCATS_ARRCAT_v20_numpy.npy'
np.save(file,ach5, allow_pickle=False)
print('ARRCAT saved as '+file)

#for loading do:
#icnpy=np.load(file)
#decode strings:
#icnpy['icmecat_id'][0].decode()



############ other formats

#copy pandas dataframe first to change time format 
ac_copy2=copy.deepcopy(ac)  
ac_copy2.sse_launch_time=parse_time(ac.sse_launch_time).iso
ac_copy2.target_arrival_time=parse_time(ac.target_arrival_time).iso

#change time format
for i in np.arange(len(ac)):

    dum=ac_copy2.sse_launch_time[i] 
    ac_copy2.at[i,'sse_launch_time']=dum[0:16]
     
    dum=ac_copy2.target_arrival_time[i] 
    ac_copy2.at[i,'target_arrival_time']=dum[0:16]


#save as json for webpage with different time format
file='arrcat/HELCATS_ARRCAT_v20_iso.json'
ac_copy2.to_json(file)
print('ARRCAT saved as '+file)

#save as html no header
file='arrcat/HELCATS_ARRCAT_v20_simple.html'
ac_copy.to_html(file)
print('ARRCAT saved as '+file)

############ save as html file with header
#save as html
file='arrcat/HELCATS_ARRCAT_v20.html'
#ic.to_html(file,justify='center')

#ichtml='{% extends "_base.html" %} \n \n {% block content %} \n \n \n '
achtml = header_html
achtml += ac_copy.to_html()
#ichtml +='\n \n {% endblock %}'


with open(file,'w') as f:
    f.write(achtml)
    f.close()
    
print('ARRCAT saved as '+file)    


# ## 3 load ARRCAT examples

# In[ ]:


#load icmecat as pandas dataframe
file='arrcat/HELCATS_ARRCAT_v20_pandas.p'
[ac_pandas,h]=pickle.load( open(file, 'rb'))   

#load icmecat as numpy array
file='arrcat/HELCATS_ARRCAT_v20_numpy.p'
[ac_rec,ac_struct,h]=pickle.load( open(file, 'rb'))   


# In[ ]:


ac_pandas
ac_pandas.keys()


# In[ ]:


ac


# In[ ]:


ac_rec.id
ac_rec.target_name[5]


# In[ ]:


ac_struct


# In[ ]:


ac_struct['id']
ac_struct['id'][100]

#comparison old and new method for arrival time prediction
deltata=(parse_time(ac.target_arrival_time[0:100]).plot_date-parse_time(ac_old.target_arrival_time).plot_date[0:100])*24
print(deltata)


# ### plot directions and targets

# In[ ]:



plt.rcParams["figure.figsize"] = (19,11)

sns.set_context('talk')
sns.set_style('darkgrid')



fig=plt.figure(1, figsize=(18,10), dpi=100) 
ax = plt.subplot(121,projection='polar') 

plt.title('ARRCAT CME arrivals at targets [HEEQ longitude]')

#get indices for each target
merci=np.where(ac.target_name=='Mercury')[0]
venusi=np.where(ac.target_name=='Venus')[0]
earthi=np.where(ac.target_name=='Earth_L1')[0]
marsi=np.where(ac.target_name=='Mars')[0]

stai=np.where(ac.target_name=='STEREO-A')[0]
stbi=np.where(ac.target_name=='STEREO-B')[0]
pspi=np.where(ac.target_name=='PSP')[0]
soloi=np.where(ac.target_name=='SolarOrbiter')[0]
bepii=np.where(ac.target_name=='BepiColombo')[0]
ulyi=np.where(ac.target_name=='Ulysses')[0]

#markersize
ms=15
#alpha
al=0.6

ax.scatter(np.radians(ac.target_heeq_lon[merci]),ac.target_distance[merci],s=ms,c='dimgrey', alpha=al)
ax.scatter(np.radians(ac.target_heeq_lon[venusi]),ac.target_distance[venusi],s=ms,c='orange', alpha=al)
ax.scatter(np.radians(ac.target_heeq_lon[earthi]),ac.target_distance[earthi],s=ms,c='mediumseagreen', alpha=al)
ax.scatter(np.radians(ac.target_heeq_lon[marsi]),ac.target_distance[marsi],s=ms,c='orangered', alpha=al)

ax.scatter(np.radians(ac.target_heeq_lon[stai]),ac.target_distance[stai],s=ms,c='red', alpha=al)
ax.scatter(np.radians(ac.target_heeq_lon[stbi]),ac.target_distance[stbi],s=ms,c='blue', alpha=al)

ax.scatter(np.radians(ac.target_heeq_lon[pspi]),ac.target_distance[pspi],s=ms,c='black', alpha=al)
ax.scatter(np.radians(ac.target_heeq_lon[soloi]),ac.target_distance[soloi],s=ms,c='green', alpha=al)
ax.scatter(np.radians(ac.target_heeq_lon[bepii]),ac.target_distance[bepii],s=ms,c='violet', alpha=al)
ax.scatter(np.radians(ac.target_heeq_lon[ulyi]),ac.target_distance[ulyi],s=ms,c='brown', alpha=al)

plt.ylim([0,np.max(ac.target_distance)+0.03])



ax = plt.subplot(122,projection='polar')
plt.title('HIGeoCAT SSE speed')
hiai=np.where(ac.sc=='A')[0]
hibi=np.where(ac.sc=='B')[0]

ax.scatter(np.radians(ac.sse_heeq_lon[hiai].astype(float)),ac.sse_speed[hiai],s=ms,c='red', alpha=al)
ax.scatter(np.radians(ac.sse_heeq_lon[hibi].astype(float)),ac.sse_speed[hibi],s=ms,c='blue', alpha=al)

plt.ylim([0,2000])

plt.tight_layout()

plotfile='arrcat/plots_arrcat/arrcat_targets.png'
plt.savefig(plotfile,dpi=100)
print('saved as ',plotfile)


# ### plot error distributions

# In[ ]:


fig=plt.figure(2, figsize=(18,8), dpi=100)

ax = plt.subplot(121) 
sns.distplot(ac.sse_speed_err,bins=200,kde=False,rug=True)
plt.xlim(0,200)
plt.ylabel('events')
plt.xlabel('sse_speed_err [km/s]')
plt.title('HIGeoCAT SSE speed error distribution')
print('mean sse_speed error',np.mean(ac.sse_speed_err).astype(int),'+/-',np.std(ac.sse_speed_err).astype(int),' km/s' )

ax2 = plt.subplot(122) 
sns.distplot(ac.target_arrival_time_err,bins=800,kde=False,rug=True)
plt.xlim(0,70)
plt.title('ARRCAT target arrival time error distribution')
plt.ylabel('events')
plt.xlabel('target_arrival_time_err [hours]')
print('mean target arrival time error from sse_speed error',np.mean(ac.target_arrival_time_err).astype(int),'+/-',np.std(ac.target_arrival_time_err).astype(int),'hours' )

plt.tight_layout()

plotfile='arrcat/plots_arrcat/arrcat_errors_dist.png'
plt.savefig(plotfile,dpi=100)
print('saved as ',plotfile)


# ## HIGeoCAT CME rates plot
# 
# 
# 
# 
# 

# In[ ]:


'''
(histmerc1, bin_edgesmerc) = np.histogram(parse_time(ac.sse_launch_time[merci]).jd, yearly_bin_edges)
(histvenus1, bin_edgesvenus) = np.histogram(parse_time(ac.sse_launch_time[venusi]).jd, yearly_bin_edges)
(histearth1, bin_edgesearth) = np.histogram(parse_time(ac.sse_launch_time[earthi]).jd, yearly_bin_edges)
(histmars1, bin_edgesmars) = np.histogram(parse_time(ac.sse_launch_time[marsi]).jd, yearly_bin_edges)
(histstb1, bin_edgesstb) = np.histogram(parse_time(ac.sse_launch_time[stbi]).jd, yearly_bin_edges)
(histsta1, bin_edgessta) = np.histogram(parse_time(ac.sse_launch_time[stai]).jd, yearly_bin_edges)
(histpsp1, bin_edgespsp) = np.histogram(parse_time(ac.sse_launch_time[pspi]).jd, yearly_bin_edges)
(histbepi1, bin_edgesbepi) = np.histogram(parse_time(ac.sse_launch_time[bepi]).jd, yearly_bin_edges)
(histsolo1, bin_edgessolo) = np.histogram(parse_time(ac.sse_launch_time[soloi]).jd, yearly_bin_edges)


# these are for arrcat
hiai=np.where(ac.sc=='A')[0]
hibi=np.where(ac.sc=='B')[0]


(hist_hia, binedges) = np.histogram(parse_time(ac.sse_launch_time[hiai]).plot_date, yearly_bin_edges)
(hist_hib, binedges) = np.histogram(parse_time(ac.sse_launch_time[hibi]).plot_date, yearly_bin_edges)


'''



# In[ ]:


last_year=2021
#define dates of January 1 from 2007 to 2020
years_jan_1_str=[str(i)+'-01-01' for i in np.arange(2007,last_year+1) ] 
yearly_bin_edges=parse_time(years_jan_1_str).plot_date

#
hiai=np.where(higeocat['SC']==b'A')[0]
hibi=np.where(higeocat['SC']==b'B')[0]


hia_t0=higeocat['SSE Launch'][hiai]

hib_t0=higeocat['SSE Launch'][hibi]


(hist_hia, binedges) = np.histogram(parse_time(hia_t0).plot_date, yearly_bin_edges)
(hist_hib, binedges) = np.histogram(parse_time(hib_t0).plot_date, yearly_bin_edges)

#convert back to datetime and make shape consistent
binedges1=binedges[:-1]+mdates.date2num(np.datetime64('0000-12-31'))

#change matplotlib time before plotting
yearly_bin_edges2=yearly_bin_edges + mdates.date2num(np.datetime64('0000-12-31'))


#---------------------------------------------------------------
sns.set_context("talk")     
#sns.set_style('darkgrid')
#sns.set_style('whitegrid',{'grid.linestyle': '--'})

sns.set_style("ticks",{'grid.linestyle': '--'})
fsize=15

fig=plt.figure(1,figsize=(14,8),dpi=80)



ax1 = plt.subplot(111) 
ax1.set_title('Yearly CME rate observed by the STEREO heliospheric imagers',fontsize=15)

ax1.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)

binweite=365/3
alp=1.0
ax1.bar(binedges1+binweite,hist_hia, width=binweite,color='tomato', alpha=alp,label='HI on STEREO-A')
ax1.bar(binedges1+binweite*2,hist_hib, width=binweite,color='mediumslateblue', alpha=alp,label='HI on STEREO-B')

ax1.set_xlim(yearly_bin_edges2[0],yearly_bin_edges2[-1])
ax1.legend(loc=2,fontsize=14)
ax1.grid(alpha=0.5)
ax1.set_xticks(yearly_bin_edges2) 
ax1.set_ylabel('CMEs observed per year')

ax1.text(yearly_bin_edges2[11],165,'last update: '+last_update,fontsize=15,zorder=2,horizontalalignment='center')


plotfile='arrcat/plots_arrcat/higeocat_rate.png'
plt.savefig(plotfile,dpi=100)
print('saved as ',plotfile)


# In[ ]:





# In[ ]:





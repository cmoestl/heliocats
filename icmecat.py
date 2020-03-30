#!/usr/bin/env python
# coding: utf-8

# ## icmecat
# 
# makes the ICMECATv2.0
# 
# Author: C. Moestl, IWF Graz, Austria
# twitter @chrisoutofspace, part of https://github.com/cmoestl/heliocats
# last update March 2020
# 
# Install a specific conda environment to run this code, see https://github.com/cmoestl/heliocats
# 
# **current status: work in progress**
# 
# * to do: despike sta stb wind all, new STA, Wind and PSP converted to SCEQ components 
# 
# * Convert to script with 
# 
#     jupyter nbconvert --to script icmecat.ipynb
#     
#     
# * Adding a new event: edit the file icmecat/HELCATS_ICMECAT_v20_master.xlsx to add 3 times, the id and spacecraft name, 
# delete the file for the respective spacecraft under data/indices_icmecat, and run this notebook or script.
# 

# In[22]:


import numpy as np
import scipy.io
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import  DateFormatter
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


from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging

from heliocats import data as hd
importlib.reload(hd) #reload again while debugging

from heliocats import cats as hc
importlib.reload(hc) #reload again while debugging

from heliocats import stats as hs
importlib.reload(hs) #reload again while debugging

#where the 6 in situ data files are located is read from input.py
#as data_path=....
from input import data_path


########### make directories first time if not there

resdir='results'
if os.path.isdir(resdir) == False: os.mkdir(resdir)

datadir='data'
if os.path.isdir(datadir) == False: os.mkdir(datadir)

indexdir='data/indices_icmecat' 
if os.path.isdir(indexdir) == False: os.mkdir(indexdir) 


catdir='icmecat'
if os.path.isdir(catdir) == False: os.mkdir(catdir)

icplotsdir='icmecat/plots_icmecat/' 
if os.path.isdir(icplotsdir) == False: os.mkdir(icplotsdir) 


# ### (1) load new data with HelioSat and heliocats.data

# In[2]:



load_data=1

if load_data > 0:

 print('load new Wind, STEREO-A, MAVEN, Parker Solar Probe data')

 # MAVEN
 #filemav='maven_2014_2018.p'
 #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )

 #filemav='maven_2014_2018_removed.p'
 #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )
 
 filemav='maven_2014_2018_removed_smoothed.p'
 [mav,hmav]=pickle.load(open(data_path+filemav, 'rb' ) )

 # Wind
 filewin="wind_2018_2019.p" 
 #for updating data
 #start=datetime.datetime(2018, 1, 1)
 #end=datetime.datetime.utcnow()
 #hd.save_wind_data(data_path,filewin,start,end)
 [win2,hwin2]=pickle.load(open(data_path+filewin, "rb" ) )  

 # STEREO-A    
 filesta2='sta_2018_2019_beacon.p'
 #start=datetime.datetime(2018, 1, 1)
 #end=datetime.datetime(2019, 12, 31)
 #hd.save_stereoa_beacon_data(data_path,filesta,start,end)

 [sta2,hsta2]=pickle.load(open(data_path+filesta2, "rb" ) )  

 # Parker Solar Probe
 filepsp='psp_2018_2019.p'
 [psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) )  


 # ADD BepiColombo  
 
 
 # ADD Solar Orbiter
 
 
 # Ulysses is currently taken from the full 
 # helcats data below, but a file is available on figshare


 # get data file from helcats with headers
 [vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly]=hd.load_helcats_datacat(data_path+'helcats_all_data_removed.p') 
 print('done')


# ## (2) measure new events 

# In[5]:


#for measuring new events use this function from heliocats.plot 
#plt.close('all')
#works in jupyter notebooks

#works in scripts
#matplotlib.use('qt5agg')  
#plt.ion()
#hp.plot_insitu_measure(psp, '2018-Nov-10','2018-Nov-15', 'PSP', 'results/plots_icmecat/')
#%matplotlib 
#hp.plot_insitu_measure(win, '2010-Oct-29','2010-Oct-31', 'Wind', 'results/')


#for plotting single events
#hp.plot_insitu(psp, ic.icme,'2018-Nov-15', 'PSP', icplotsdir)

#-----------------------

#read HIGEOCAT from https://www.helcats-fp7.eu/catalogues/wp3_cat.html
#https://docs.astropy.org/en/stable/io/votable/
from astropy.io.votable import parse_single_table
table = parse_single_table('data/HCME_WP3_V06.vot')
data = table.array
a=table.array['HM HEEQ Long'][10]
print(a)



#    "columns" : [ "ID", "Date [UTC]", "SC", "L-N", "PA-N [deg]", "L-S", "PA-S [deg]", "Quality" 
#       , "PA-fit [deg]"
#       , "FP Speed [kms-1]", "FP Speed Err [kms-1]", "FP Phi [deg]", "FP Phi Err [deg]","FP HEEQ Long [deg]",  "FP HEEQ Lat [deg]",  "FP Carr Long [deg]", "FP Launch [UTC]"
#       , "SSE Speed [kms-1]", "SSE Speed Err [kms-1]", "SSE Phi [deg]", "SSE Phi Err [deg]", "SSE HEEQ Long [deg]", "SSE HEEQ Lat [deg]",  "SSE Carr Long [deg]","SSE Launch [UTC]"
#       , "HM Speed [kms-1]", "HM Speed Err [kms-1]", "HM Phi [deg]", "HM Phi Err [deg]", "HM HEEQ Long [deg]", "HM HEEQ Lat [deg]", "HM Carr Long [deg]", "HM Launch [UTC]"
#  ],




# ### (3) make ICMECAT 

# In[6]:


print('data loaded')
ic=hc.load_helcats_icmecat_master_from_excel('icmecat/HELCATS_ICMECAT_v20_master.xlsx')

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

ic=hc.get_cat_parameters(win,wini,ic,'Wind')
ic=hc.get_cat_parameters(sta,stai,ic,'STEREO-A')
ic=hc.get_cat_parameters(stb,stbi,ic,'STEREO-B')
ic=hc.get_cat_parameters(vex,vexi,ic,'VEX')
ic=hc.get_cat_parameters(mes,mesi,ic,'MESSENGER')
ic=hc.get_cat_parameters(uly,ulyi,ic,'ULYSSES')
ic=hc.get_cat_parameters(mav,mavi,ic,'MAVEN')
ic=hc.get_cat_parameters(psp,pspi,ic,'PSP')


####### 3c make all plots if wanted
# matplotlib.use('Agg')
# hp.plot_icmecat_events(win,wini,ic,'Wind',icplotsdir)
# hp.plot_icmecat_events(sta,stai,ic,'STEREO-A',icplotsdir)
# hp.plot_icmecat_events(stb,stbi,ic,'STEREO-B',icplotsdir)
# hp.plot_icmecat_events(vex,vexi,ic,'VEX',icplotsdir)
# hp.plot_icmecat_events(mes,mesi,ic,'MESSENGER',icplotsdir)
# hp.plot_icmecat_events(uly,ulyi,ic,'ULYSSES',icplotsdir)
# hp.plot_icmecat_events(mav,mavi,ic,'MAVEN',icplotsdir)
# hp.plot_icmecat_events(psp,pspi,ic,'PSP',icplotsdir)


# ### (4) save ICMECAT 

# ### 4a save header

# In[146]:


#save header and parameters as text file and prepare for html website
print('header:')
header='ICME CATALOGUE v2.0 \n\nThis is the HELCATS interplanetary coronal mass ejection (ICME) catalog, based on in situ magnetic field and bulk plasma observations in the heliosphere. \n\nThis is version 2.0, released in April 2020 with a major update to the original version 1.0, originally a product of EU HELCATS project (2014-2017). \n\nReleased 2020-**-**. DOI: 10.6084/m9.figshare.4588315.v2 *** \nThe catalog is available as .p (python pickle), .xlsx, .csv, .html, .json at https://helioforecast.space/catalogs \n\nNumber of events in ICMECAT: '+str(len(ic))+' \nICME observatories: Wind, STEREO-A, STEREO-B, Venus Express (VEX), MESSENGER, Ulysses, MAVEN, Parker Solar Probe (PSP).   \nTime range: January 2007 - December 2018. \n \nAuthors: Christian Moestl, Andreas Weiss, Space Research Institute, Austrian Academy of Sciences, Graz, Austria.\nContributors: Peter Boakes, Alexey Isavnin, Emilia Kilpua, Reka Winslow, Brian Anderson, Lydia Philpott, Vratislav Krupar, Jonathan Eastwood, Simon Good, Lan Jian, Teresa Nieves-Chinchilla, Cyril Simon Wedlund, Jingnan Guo, Mateja Dumbovic, Benoit Lavraud.  \n\nThis catalog has been made by getting the 3 times of each ICME (shock or disturbance begin, magnetic obstacle start and end) from the individual catalogs below, and then calculating all parameters again consistently from the data by us. \nEach icmecat_id has a tag in it that indicates from which catalog the ICME times were taken. \n\nWind:       Nieves-Chinchilla et al. (2018), tag: NASA. \nSTEREO-A:   Jian et al. (2018), tag: JIAN. \nSTEREO-B:   Jian et al. (2018), tag: JIAN. \nVEX:        Good et al. (2018), tag: SGOOD \nMESSENGER:  Good et al. (2018), Winslow et al. (2018), tags: SGOOD, WINSLOW. \nMAVEN:      Möstl et al. (2020, in prep.), tag: MOESTL.\nUlysses:    Added by us, tag: MOESTL. \nPSP:        Added by us, tag: MOESTL. \nWe have also added extra events at Ulysses, STEREO-A, PSP, MESSENGER, VEX, Wind and MAVEN (all tagged with MOESTL in icmecat_id).\n\nThe in situ data (4 GB) can be downloaded in python pickle format as recarrays from https://doi.org/10.6084/m9.figshare.11973693 \nThe python code for producing this catalog is available as part of https://github.com/cmoestl \n\nReferences: \nNieves-Chinchilla, T. et al. (2018), https://doi.org/10.1007/s11207-018-1247-z  https://wind.nasa.gov/fullcatalogue.php \nJian, L. et al. (2018), https://doi.org/10.3847/1538-4357/aab189 https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level3/ \nGood, S. et al. (2018) https://doi.org/10.1007/s11207-015-0828-3 \nWinslow, R. et al. (2015), https://doi.org/10.1002/2015JA021200 \nMöstl, C. et al. (2020) in preparation \n\n\nComments: \n- Spacecraft positions are given in Heliocentric Earth Equatorial Coordinates (HEEQ) coordinates. \n- Coordinate system for all magnetic field components is SCEQ except Ulysses and Parker Solar Probe (RTN), MAVEN (MSO) and Wind (HEEQ). \n        Definition of SpaceCraft Equatorial Coordinates (SCEQ): \n        Z is the solar rotation axis. \n        X points from the Sun to the spacecraft, projected in the solar equatorial plane. \n        Y completes the right handed triad and points to solar west. \n        This system is thus centered on the respective in situ spacecraft. \n        The solar equatorial plane as the reference plane is similar for all spacecraft.\n- Venus Express and MESSENGER do not have plasma parameters available. \n- If there is no sheath region, so the ICME starts immediately with a magnetic obstacle, the icme_start_time is similar to mo_start_time.\n- At MESSENGER and VEX, for events cataloged by Simon Good, icme_start_time has been added by V. Krupar (Imperial College) and C. Möstl (IWF Graz). \n- For the calculation of the parameters at MESSENGER during the orbit around Mercury, all data points inside the bowshock of Mercury have been removed, according to a list thankfully provided to us by by R. Winslow, UNH, B. Anderson, APL, and Lydia Philpott, UBC. \n- Calculation of the magnetic obstacle parameters at VEX is done after approximate removal of the induced magnetosphere, with a modified equation \nin Zhang et al. 2008 (doi: 10.1016/j.pss.2007.09.012), with a constant of 3.5 instead of 2.14/2.364,in order to account for a larger bowshock distance during solar maximum than studied in this paper. \n- For MAVEN, all data inside the bow shock were removed with the model from ****************** (by C. Simon Wedlund). From the remaining data, the median for each orbit is taken as 1 data point, resulting in a solar wind dataset at Mars with 4.5 hour time resolution.\n- The identification of ICMEs for MAVEN is a mixture of methods using data from MSL/RAD, MAVEN and STEREO/HI (see Möstl et al. 2020, in prep.). \n\n\n\n'



parameters='Parameters:\n01: icmecat_id: The unique identifier for the observed ICME. unit: string. \n02: icme_start_time: Shock arrival or density enhancement time, can be similar to mo_start_time. unit: UTC. \n03: mo_start_time: Start time of the magnetic obstacle (MO), including flux ropes, flux-rope-like, and ejecta signatures. unit: UTC. \n04: mo_end_time: End time of the magnetic obstacle. unit: UTC. \n05: mo_sc_heliodistance: Heliocentric distance of the spacecraft at mo_start_time. unit: AU.\n06: mo_sc_long_heeq: Heliospheric longitude of the spacecraft at mo_start_time, range [-180,180]. unit: degree (HEEQ).\n07: mo_sc_lat_heeq: Heliospheric latitude of the spacecraft at mo_start_time, range [-90,90]. unit: degree (HEEQ).\n08: icme_duration: Duration of the interval between icme_start_time and mo_endtime. unit: hours.\n09: icme_bmax: Maximum total magnetic field in the full icme interval (icme_start_time to mo_end_time). unit: nT.\n10: icme_bmean: Mean total magnetic field during the full icme interval (icme_start_time to mo_end_time). unit: nT.\n11: icme_bstd: Standard deviation of the total magnetic field from icme_start_time to mo_end_time. unit: nT.\n12: icme_speed_mean: Mean proton speed from icme_start_time to mo_end_time. unit: km/s.\n13: icme_speed_std: Standard deviation of proton speed from icme_start_time to mo_end_time. unit: km/s.\n14: mo_duration: Duration of the interval between mo_start_time and mo_endtime. unit: hours.\n15: mo_bmax: Maximum total magnetic field in the magnetic obstacle interval (mo_start_time to mo_end_time). unit: nT.\n16: mo_bmean: Mean total magnetic field in the magnetic obstacle. unit: nT.\n17: mo_bstd: Standard deviation of the total magnetic field in the magnetic obstacle. unit: nT.\n18: mo_bzmean: Mean magnetic field Bz component in the magnetic obstacle. unit: nT.\n19: mo_bzmin: Minimum magnetic field Bz component in the magnetic obstacle. unit: nT.\n20: mo_bzstd: Standard deviation of the magnetic field Bz component in the magnetic obstacle. unit: nT.\n21: mo_bymean: Mean magnetic field By component in the magnetic obstacle. unit: nT.\n22: mo_bystd: Standard deviation of the magnetic field By component in the magnetic obstacle. unit: nT.\n23: mo_speed_mean: Mean proton speed from mo_start_time to mo_end_time. unit: km/s.\n24: mo_speed_std: Standard deviation of proton speed from mo_start_time to mo_end_time. unit: km/s.\n25: mo_expansion_speed: Difference between proton speed at mo_start_time to proton speed at mo_end_time. unit: km/s.\n26: mo_pdyn_mean: Mean proton dynamic pressure from mo_start_time to mo_start_time. unit: nPa.\n27: mo_pdyn_std: Standard deviation of proton dynamic pressure from mo_start_time to mo_start_time. unit: nPa.\n28: mo_density_mean: Mean proton density from mo_start_time to mo_start_time. unit: cm^-3.\n29: mo_density_std: Standard deviation of proton density from mo_start_time to mo_start_time. unit: cm^-3.\n30: mo_temperature_mean: Mean proton temperature from mo_start_time to mo_start_time. unit: K.\n31: mo_temperature_std: Standard deviation of proton temperature from mo_start_time to mo_end_time. unit: K.\n32: sheath_speed_mean: Mean proton speed from icme_start_time to mo_start_time, NaN if these times are similar. unit: km/s.\n33: sheath_speed_std: Standard deviation of proton speed from icme_start_time to mo_start_time, NaN if these times are similar. unit: km/s.\n34: sheath_density_mean: Mean proton density from icme_start_time to mo_start_time, NaN if these times are similar. unit: cm^-3.\n35: sheath_density_std: Standard deviation of proton density from icme_start_time to mo_start_time, NaN if these times are similar. unit: cm^-3.\n36: sheath_pdyn_mean: Mean proton dynamic pressure, from icme_start_time to mo_start_time, NaN if these times are similar. unit: nPa.\n37: sheath_pdyn_std: Standard deviation of proton dynamic pressure, from icme_start_time to mo_start_time, NaN if these times are similar. unit: nPa.\n\n\n'


print('-------------------------')
print(header)
print(parameters)
print('-------------------------')


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

# In[174]:


########## python formats

# save ICMECAT as pandas dataframe with times as datetime objects as pickle
file='icmecat/HELCATS_ICMECAT_v20_pandas.p'
pickle.dump(ic, open(file, 'wb'))
print('ICMECAT saved as '+file)


# save ICMECAT as numpy array with times as matplotlib datetime as pickle
ic_num=copy.deepcopy(ic) 
ic_num.icme_start_time=parse_time(ic_num.icme_start_time).plot_date
ic_num.mo_start_time=parse_time(ic_num.mo_start_time).plot_date
ic_num.mo_end_time=parse_time(ic_num.mo_end_time).plot_date
ic_num=ic_num.to_numpy()
file='icmecat/HELCATS_ICMECAT_v20_numpy.p'
pickle.dump(ic_num, open(file, 'wb'))
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

#save as hdf needs pip install tables
#file='icmecat/HELCATS_ICMECAT_v20.hdf'
#ic.to_hdf(file,key='icmecat')

#save as txt
file='icmecat/HELCATS_ICMECAT_v20.txt'
np.savetxt(file, ic_copy.values.astype(str), fmt='%s' )

print('ICMECAT saved as '+file)


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

# In[170]:


#load icmecat as pandas dataframe
file='icmecat/HELCATS_ICMECAT_v20_pandas.p'
ic_pandas=pickle.load( open(file, 'rb'))   

#load icmecat as numpy array
file='icmecat/HELCATS_ICMECAT_v20_numpy.p'
ic_np=pickle.load( open(file, 'rb'))   
print(ic_np)


# In[171]:


ic_pandas


# In[ ]:





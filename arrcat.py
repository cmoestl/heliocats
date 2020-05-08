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

# In[14]:


from heliocats import cats as hc
importlib.reload(hc) #reload again while debugging 


#LOAD HELCATS HIGeoCAT

url_higeocat='https://www.helcats-fp7.eu/catalogues/data/HCME_WP3_V06.vot'

try: urllib.request.urlretrieve(url_higeocat,'data/HCME_WP3_V06.vot')
except urllib.error.URLError as e:
    print('higeocat not loaded')


higeocat=hc.load_higeocat_vot('data/HCME_WP3_V06.vot')
higeocat_time=parse_time(higeocat['Date']).datetime    
higeocat_t0=parse_time(higeocat['SSE Launch']).datetime   #backprojected launch time
#Make arrival catalog from HIGEOCAT
arrcat_sta=hc.make_arrival_catalog_insitu_ssef30(higeocat, 'STA')

arrcat_stb=hc.make_arrival_catalog_insitu_ssef30(higeocat, 'STB')



arrcat_psp=hc.make_arrival_catalog_insitu_ssef30(higeocat, 'PSP')
arrcat_solo=hc.make_arrival_catalog_insitu_ssef30(higeocat, 'Solo')
arrcat_bepi=hc.make_arrival_catalog_insitu_ssef30(higeocat, 'Bepi')

    

arrcat_mercury=hc.make_arrival_catalog_insitu_ssef30(higeocat, 'Mercury')
arrcat_venus=hc.make_arrival_catalog_insitu_ssef30(higeocat, 'Venus')
arrcat_earth=hc.make_arrival_catalog_insitu_ssef30(higeocat, 'Earth')
arrcat_mars=hc.make_arrival_catalog_insitu_ssef30(higeocat, 'Mars')

#https://www.helcats-fp7.eu/catalogues/data/HELCATS_ARRCAT_v6.txt
#**TO DO: make catalog
#    HCME_A__20070520_01  A   2007-05-19T17:54Z   EARTH_L1     -13     13.0     444      420     2007-05-23T20:52Z       1.0017      -2.17      -0.00       97     120     110     130     120
# - Earth -> L1!!
# add STEREO-A, STEREO-B, keine cruise phases
    


# ## save ARRCAT
# 

# #### save header

# In[ ]:



This ARRival CATalog (ARRCAT) is a product of working packages 3 and 4 in the EU HELCATS project (2014-2017). 

It lists predicted arrivals of solar coronal mass ejections at various spacecraft and planets with the STEREO heliospheric imager instruments, between April 2007 - September 2014.

HEADER FILE FOR: HELCATS_ARRCAT_v6.txt, HELCATS_ARRCAT_v6.sav

the .sav file can be read directly in IDL ("restore" function) and python ("scipy.io.readsav").
 
AUTHORS: Christian Moestl, Peter Boakes, University of Graz, Austria; SRI, Austrian Academy of Sciences, Graz, Austria. Based on the HIGeoCAT catalog of CMEs established at RAL Space, UK (Harrison, Davies, Barnes).
 
FILE CREATION DATE: Wed Sep 21 11:27:24 2016
 
INPUT FILES: HCME_WP3_V03.json
 
Number of events in ARRCAT: 1995

Targets: EARTH-L1, STEREO-A, STEREO-B, VENUS, MESSENGER, MARS, SATURN, ULYSSES, MSL, MAVEN, ROSETTA

 VARIABLES: 
	1: ID: From HICAT, the unique identifier for the observed CME.
	2: SC: From HICAT, the HI observing STEREO spacecraft, (A=Ahead or B=Behind)
	3: SSE_LAUNCH: From HICAT, launch time of the CME on the Sun, unit: UTC.
	4: TARGET_NAME: Name of in situ target.
	5: SSE_HEEQ_LONG: From HICAT, the HEEQ longitude of the CME apex propagation direction, unit: degree.
	6: TARGET_DELTA: Difference in HEEQ longitude between central CME direction and target location, positive values: spacecraft is west of CME apex. unit: degree.
	7: SSE_SPEED: From HICAT, speed of CME apex, unit: km/s.
	8: TARGET_SPEED: CME arrival speed at target location, corrected for SSE shape. unit: km/s.
	9: TARGET_ARRIVAL: CME arrival time at target location, corrected for SSE shape. unit: UTC.
	10: TARGET_DISTANCE: Target distance from Sun, at CME launch time. unit: AU.
	11: TARGET_HEEQ_LAT: Target latitude in HEEQ, at CME launch time. unit: degree.
	12: TARGET_HEEQ_LONG: Target longitude in HEEQ, at CME launch time. unit: degree.
	13: TARGET_PA: PA of target from HI observing STEREO spacecraft, unit: degree.
	14: PA_FIT: From HICAT, PA along which time-elongation profile is extracted, unit: degree.
	15: PA_N: From HICAT, northern position angle of CME, unit: degree.
	16: PA_S: From HICAT, southernmost position angle of CME, unit: degree.
	17: PA_CENTER: average of pa_n and pa_s, unit: degree.

Notes:

1. We have applied the method from Möstl & Davies (2013, Solar Physics) for calculating speeds and arrival times of the CMEs modeled with SSEF30 to all CMEs in the HELCATS HIGeoCAT catalog (see website helcats-fp7.eu, and Möstl et al. 2014, ApJ, for more details). If the SSEF30 circle hits a spacecraft or planet, an entry in ARRCAT is produced.

2. The position of Venus Express is assumed equal to the location of Venus. Arrivals at Ulysses are calculated only around its last ecliptic pass in August 2007. For Rosetta, no arrivals are calculated during its deep space hibernation from 2011 June 8 to 2014 January 20. For MESSENGER, MSL and MAVEN ARRCAT covers both the cruise and orbit phases of those missions. 





















#save header and parameters as text file and prepare for html website
header='ICME CATALOGUE v2.0 \n\nThis is the HELCATS interplanetary coronal mass ejection (ICME) catalog, based on in situ magnetic field and bulk plasma observations in the heliosphere. \n\nThis is version 2.0, released 2020-**-**. DOI: 10.6084/m9.figshare.6356420 \n\nThe catalog is available as  python pandas dataframe (pickle), python numpy structured array (pickle), json, csv, xlsx, txt, hdf5, at \nhttps://helioforecast.space/icmecat \nhttps://www.helcats-fp7.eu/catalogues/wp4_icmecat.html \n\nNumber of events in ICMECAT: '+str(len(ic))+' \nICME observatories: Parker Solar Probe (PSP), Wind, STEREO-A, MAVEN, STEREO-B, Venus Express (VEX), MESSENGER, Ulysses.   \nTime range: January 2007 - December 2019. \n \nAuthors: Christian Moestl, Andreas Weiss, Space Research Institute, Austrian Academy of Sciences, Graz, Austria.\nContributors: Peter Boakes, Alexey Isavnin, Emilia Kilpua, Reka Winslow, Brian Anderson, Lydia Philpott, Vratislav Krupar, Jonathan Eastwood, Simon Good, Lan Jian, Teresa Nieves-Chinchilla, Cyril Simon Wedlund, Jingnan Guo, Mateja Dumbovic, Benoit Lavraud.  \n\nRules: If results are produced with this catalog for peer-reviewed scientific publications, please contact christian.moestl@oeaw.ac.at for possible co-authorship. \n\nThis catalog has been made by getting the 3 times of each ICME (shock or disturbance begin, magnetic obstacle start and end) from the individual catalogs below, and then calculating all parameters again consistently from the data by us. \nThe in situ data that were used for the catalog, with a size of 8 GB in total, including extra data files with magnetic field components in RTN coordinates that are not used for producing the catalog, can be downloaded in python pickle format as recarrays from https://doi.org/10.6084/m9.figshare.11973693 \nThe python code for producing this catalog is available as part of https://github.com/cmoestl/heliocats \n\nEach icmecat_id has a tag in it that indicates from which catalog the ICME times were taken: \n\nWind:       Nieves-Chinchilla et al. (2018), tag: NASA. \nSTEREO-A:   Jian et al. (2018), tag: JIAN. \nSTEREO-B:   Jian et al. (2018), tag: JIAN. \nVEX:        Good et al. (2018), tag: SGOOD \nMESSENGER:  Good et al. (2018), Winslow et al. (2018), tags: SGOOD, WINSLOW. \nMAVEN:      Möstl et al. (2020, in prep.), tag: MOESTL.\nUlysses:    Added by us, tag: MOESTL. \nPSP:        Added by us, tag: MOESTL. \n\nWe have also added extra events at VEX, MESSENGER, Wind and STEREO-A (all tagged with MOESTL in icmecat_id).\n\nReferences \nNieves-Chinchilla, T. et al. (2018),  https://doi.org/10.1007/s11207-018-1247-z \n                                      https://wind.nasa.gov/fullcatalogue.php \nJian, L. et al. (2018), https://doi.org/10.3847/1538-4357/aab189 \n                        https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level3/ \nGood, S. et al. (2018) https://doi.org/10.1007/s11207-015-0828-3 \nWinslow, R. et al. (2015), https://doi.org/10.1002/2015JA021200 \nMöstl, C. et al. (2020) in preparation \n\n\nComments: \n- Spacecraft positions are given in Heliocentric Earth Equatorial Coordinates (HEEQ) coordinates. \n- The coordinate system for all magnetic field components is SCEQ, except for Wind (HEEQ, which is the equivalent for SCEQ for Earth) and Ulysses (RTN, because of high latitude positions) and MAVEN (MSO). \n        Definition of SpaceCraft Equatorial Coordinates (SCEQ): \n        Z is the solar rotation axis. \n        Y is the cross product of Z and R, with R being the vector that points from the Sun to the spacecraft.\n        X completes the right handed triad (and points away from the Sun). \nThis system is thus like HEEQ but centered on the respective in situ spacecraft, so the SCEQ X and Y \nbase vectors are rotated by the HEEQ longitude of the in situ spacecraft from HEEQ X and Y.\nThe Y vector is similar to the T vector in an RTN system for each spacecraft, but the X and Z vectors \nare rotated around Y compared to an RTN system. The differences between RTN and SCEQ for spacecraft within \na few degrees of the solar equatorial plane are very small (within a few 0.1 nT usually).\nWe choose SCEQ because it has the advantage that a comparison between multipoint CME events \nand for comparison to simulations there is always a similar reference plane (the solar equatorial plane).\n\n- Venus Express and MESSENGER do not have plasma parameters available. \n- If there is no sheath or density pileup region, so the ICME starts immediately with a magnetic obstacle, the icme_start_time is similar to mo_start_time.\n- At MESSENGER and VEX, for events cataloged by Simon Good, icme_start_time has been added by V. Krupar (Imperial College) and C. Möstl (IWF Graz). \n- For the calculation of the parameters at MESSENGER during the orbit around Mercury, all data points inside the bowshock of Mercury have been removed, according to a list thankfully provided to us by by R. Winslow, UNH, B. Anderson, APL, and Lydia Philpott, UBC. \n- Calculation of the magnetic obstacle parameters at VEX is done after approximate removal of the induced magnetosphere, with a modified equation \nin Zhang et al. 2008 (doi: 10.1016/j.pss.2007.09.012), with a constant of 3.5 instead of 2.14/2.364,in order to account for a larger bowshock distance during solar maximum than studied in this paper. \n- For MAVEN, all data inside the bow shock were removed with the model from Gruesbeck et al. (2018, https://doi.org/10.1029/2018JA025366) by C. Simon Wedlund (IWF Graz, Austria). From the remaining data, the median for each orbit is taken as 1 data point, resulting in a solar wind dataset at Mars with 4.5 hour time resolution. The identification of ICMEs for MAVEN is a mixture of methods using data from MSL/RAD, MAVEN and STEREO/HI (see Möstl et al. 2020, in prep.).\n\n\n\n'


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

# In[ ]:


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





#!/usr/bin/env python
# coding: utf-8

# ## sircat
# 
# makes a catalog of solar wind stream interaction regions (SIRs) for the Wind and STEREO spacecraft since 2007
# 
# Authors: [C. Möstl](https://www.iwf.oeaw.ac.at/en/user-site/christian-moestl/) (twitter @chrisoutofspace), A. J. Weiss, R. L. Bailey, IWF Graz, Austria; Lan Jian, NASA, USA.
# 
# 
# **current status: work in progress**
# 
# If you want to use parts of this code for generating results for peer-reviewed scientific publications, please contact me per email (christian.moestl@oeaw.ac.at) or via https://twitter.com/chrisoutofspace for co-authorships.
# 
# 
# part of https://github.com/cmoestl/heliocats, last update April 2020
# 
# ### Installation 
# Install a specific conda environment to run this code, see readme at https://github.com/cmoestl/heliocats
# 
# 
# ### Updates
# 
# Adding a new SIR event: add times in section 2 (to be done) before the master file sircat/HELCATS_SIRCAT_v10_master.xlsx is produced. Then delete the file for the respective spacecraft under sircat/indices_sircat, and run this notebook or script.
# 
# 
# ### Data sources
# 
# **STEREO SIR list**: Lan Jian, https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level3/
# published in:
# 
# L. K. Jian, J. G. Luhmann, C. T. Russell, A. B. Galvin, Solar-Terrestrial Relations Observatory (STEREO) Observations of Stream Interaction Regions in 2007-2016: Relationship with Heliospheric Current Sheets, Solar Cycle Variations, and Dual Observations, Solar Phys., 294, 31, doi: 10.1007/s11207-019-1416-8, 2019.
# 
# The STEREO SIR list currently includes events from 2007 January 1 to 2018 July 31.
# 
# 
# **Wind SIR list** 
# ? https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018JA026396
# 
# 
# 
# 
# **Review on SIRs** by Ian G. Richardson: https://link.springer.com/article/10.1007/s41116-017-0011-z
# 
# 
# 
# Convert this notebook to a script with "jupyter nbconvert --to script sircat.ipynb"
# 
# 

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

indexdir='sircat/indices_sircat' 
if os.path.isdir(indexdir) == False: os.mkdir(indexdir) 

catdir='sircat'
if os.path.isdir(catdir) == False: os.mkdir(catdir)

scplotsdir='sircat/plots_sircat/' 
if os.path.isdir(scplotsdir) == False: os.mkdir(scplotsdir) 

#Convert this notebook to a script with jupyter nbconvert --to script icmecat.ipynb
os.system('jupyter nbconvert --to script sircat.ipynb')    


#data is updated via the icmecat.ipynb notebook
    


# ## (1) load data from HELCATS, or made with HelioSat and heliocats.data

# In[3]:


load_data=1

if load_data > 0:
    
        
    #print('load Ulysses RTN') #made with heliocats.data.save_ulysses_data
    #fileuly='ulysses_1990_2009_rtn.p'
    #[uly,huly]=pickle.load(open(data_path+fileuly, "rb" ) )      
    
    print('load STEREO-B data SCEQ') #yearly magplasma files from stereo science center, conversion to SCEQ 
    filestb='stereob_2007_2014_sceq.p'
    [stb,hstb]=pickle.load(open(data_path+filestb, "rb" ) )      
 

    ########### CURRENT ACTIVE SPACECRAFT    

    
    # ADD BepiColombo  
    
    
    # ADD Solar Orbiter
    
       
    #print('load MAVEN data MSO') #removed magnetosphere by C. Simon Wedlund, 1 data point per orbit, MSO 
    #filemav='maven_2014_2018.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )
    #filemav='maven_2014_2018_removed.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )    
    #filemav='maven_2014_2018_removed_smoothed.p'
    #[mav,hmav]=pickle.load(open(data_path+filemav, 'rb' ) )
    
    #print('load MSL RAD')
    #MSL RAD
    #rad=hd.load_msl_rad()#, rad.time,rad.dose_sol

    
    #print('load PSP data SCEQ') #from heliosat, converted to SCEQ similar to STEREO-A/B
    #filepsp='psp_2018_2019_sceq.p'
    #[psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) ) 
    
    
    print('load and merge STEREO-A data SCEQ') #yearly magplasma files from stereo science center, conversion to SCEQ 
    filesta1='stereoa_2007_2019_sceq.p'
    [sta1,hsta1]=pickle.load(open(data_path+filesta1, "rb" ) )  
    sta1=sta1[np.where(sta1.time < parse_time('2019-Sep-01 00:00').datetime)[0]]

    #beacon data
    filesta2="stereoa_2019_2020_sceq_beacon.p"
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
    filewin2="wind_2018_2019_heeq.p" 
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
    
    

#HIGEOCAT
higeocat=hc.load_higeocat()
higeocat_time=parse_time(higeocat['Date']).datetime    
    

     
         
print()
    
print()       
print('time ranges of the in situ data: ')    
print()
print('active spacecraft:')
print('Wind                 ',str(win.time[0])[0:10],str(win.time[-1])[0:10])
print('STEREO-A             ',str(sta.time[0])[0:10],str(sta.time[-1])[0:10])
#print('Parker Solar Probe   ',str(psp.time[0])[0:10],str(psp.time[-1])[0:10])
#print('MAVEN                ',str(mav.time[0])[0:10],str(mav.time[-1])[0:10])
#print('MSL/RAD              ',str(rad.time[0])[0:10],str(rad.time[-1])[0:10])
print()
print('missions finished:')
#print('VEX                  ',str(vex.time[0])[0:10],str(vex.time[-1])[0:10])
#print('MESSENGER            ',str(mes.time[0])[0:10],str(mes.time[-1])[0:10])
print('STEREO-B             ',str(stb.time[0])[0:10],str(stb.time[-1])[0:10])
#print('Ulysses              ',str(uly.time[0])[0:10],str(uly.time[-1])[0:10])
print()
# print('catalogs:')
# print()
# print('HELCATS HIGeoCAT     ',str(higeocat_time[0])[0:10],str(higeocat_time[-1])[0:10])



print('done')


# ## (2) make SIRCAT masterfile from STEREO and Wind catalogs

# In[5]:


###################### read raw STEREO SIR catalog

file='sircat/STEREO_Level3_SIR_data.xlsx'
print('load Jian STEREO catalog from excel file:', file)
sraw=pd.read_excel(file)

print('Events in STEREO SIR cat:', sraw.shape[0])

sc=sraw.loc[:,'spacecraft']

year_start=sraw.loc[:,'year_start']
stime=sraw.loc[:,'start_time']

year_end=sraw.loc[:,'year_end']
etime=sraw.loc[:,'end_time']

year_pt=sraw.loc[:,'year_pt']
ptime=sraw.loc[:,'pt_time']

print('convert times to datetime objects')

#make list for all basic times, ids etc. for master file
rows_list = []


for i in np.arange(0,sraw.shape[0]):
    
    s=stime[i]    
    y=year_start[i]
    doy=int(s[0:3])
    hour=int(s[-5:-3])
    minute=int(s[-2:])
    #print(y,doy,hour, min)
    sir_start_time=datetime.datetime(y,1,1)+timedelta(days=doy-1)+timedelta(hours=hour)+timedelta(minutes=minute)

    e=etime[i]    
    y=year_end[i]
    doy=int(e[0:3])
    hour=int(e[-5:-3])
    minute=int(e[-2:])
    #print(y,doy,hour, min)
    sir_end_time=datetime.datetime(y,1,1)+timedelta(days=doy-1)+timedelta(hours=hour)+timedelta(minutes=minute)
    
    #make id for event    
    id_time=parse_time(sir_start_time).isot
    if sc[i]=='A': sc_idstring='SIR_STEREO_A_JIAN_'
    if sc[i]=='B': sc_idstring='SIR_STEREO_B_JIAN_'

    if sc[i]=='A': sc_string='STEREO-A'
    if sc[i]=='B': sc_string='STEREO-B'
        
    sircat_id=sc_idstring+id_time[0:4]+id_time[5:7]+id_time[8:10]+'_01'
    
    #put all data for this event in a list
    list1 = [sircat_id,sc_string,parse_time(sir_start_time).isot,parse_time(sir_end_time).isot,np.nan,np.nan,np.nan,             np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
    
    #print(list1)    
    #append to full list
    rows_list.append(list1)
    



########################## read raw Wind catalog

#********** TO DO similar to STEREO-A
    
    

    
    

    
###################################  add new events
#for measuring new events use this function from heliocats.plot 
#plt.close('all')
#works in jupyter notebooks

#works in scripts
#matplotlib.use('qt5agg')  
#plt.ion()

#STEREO-A
#hp.plot_insitu_measure(sta, '2018-Jan-01 12:00','2018-Feb-01 12:00', 'STEREO-A', 'results/')

#Wind
#hp.plot_insitu_measure(win, '2019-Jan-29','2019-Feb-28', 'Wind', 'results/')

    
    
    
################ make pandas data frame for master file
        
parameters =['sircat_id','sc_insitu','sir_start_time','sir_end_time','hss_end_time',             'sir_vmax','sir_vmean','sir_vstd','sir_btmax','sir_btmean','sir_bstd','hss_vmax',             ' hss_vmean','hss_vstd','hss_btmax','hss_btmean','hss_bstd']


master=pd.DataFrame(rows_list,columns=parameters)
master   
    
   



#save master file as Excel
file='sircat/HELCATS_SIRCAT_v10_master.xlsx'
master.to_excel(file,sheet_name='SIRCATv1.0')
print()
print('SIRCAT master saved as '+file)


# ## (3) make SIRCAT 

# In[39]:


#load master file
scat=hc.load_helcats_sircat_master_from_excel('sircat/HELCATS_SIRCAT_v10_master.xlsx')
scat

####### 3a get indices for all spacecraft
stai=np.where(ic.sc_insitu == 'STEREO-A')[:][0]    
stbi=np.where(ic.sc_insitu == 'STEREO-B')[:][0]    
#wini=np.where(ic.sc_insitu == 'Wind')[:][0] 


####### 3b get parameters for all spacecraft one after another

#remove indices for new events
# os.system('rm sirecat/indices_icmecat/ICMECAT_indices_STEREO-A.p')

#scat=hc.get_sircat_parameters(sta,stai,scat,'STEREO-A')
#scat=hc.get_sircat_parameters(stb,stbi,scat,'STEREO-B')
#scat=hc.get_sircat_parameters(win,wini,scat,'Wind')


# ###### 3c make all plots if wanted
# matplotlib.use('Agg')
# hp.plot_icmecat_events(sta,stai,ic,'STEREO-A',icplotsdir)
# hp.plot_icmecat_events(sta,stai,ic,'STEREO-A',icplotsdir)


print('done')


# ### (4) save SIRCAT 

# ### 4a save header

# In[76]:


#save header and parameters as text file and prepare for html website
header='SIR CATALOGUE v1.0 \n\nThis is the HELCATS stream interaction region (SRI) catalog, based on in situ magnetic field and bulk plasma observations in the heliosphere. \n\nThis is version 1.0, released 2020-**-**. DOI: 10.6084/m9.**** \n\nThe catalog is available as  python pandas dataframe (pickle), python numpy structured array (pickle), json, csv, xlsx, txt, hdf5, at \nhttps://helioforecast.space/sircat \nhttps://www.helcats-fp7.eu/catalogues/wp4_sircat.html \n\nNumber of events in SIRCAT: '+str(len(scat))+' \nICME observatories: Wind, STEREO-A, STEREO-B.   \nTime range: January 2007 - July 2018. \n\nAuthors: Christian Moestl, Andreas Weiss, R. L. Bailey, Space Research Institute, Austrian Academy of Sciences, Graz, Austria. Lan Jian, NASA, USA. \n\nRules: If results are produced with this catalog for peer-reviewed scientific publications, please contact christian.moestl@oeaw.ac.at for possible co-authorship. \n\nThis catalog has been made by getting the start and end times of each SIR from the individual catalogs below, and then calculating all parameters again consistently from the data by us. \nThe in situ data that were used for the catalog, with a size of 8 GB in total, including extra data files with magnetic field components in RTN coordinates that are not used for producing the catalog and many other in situ data files, can be downloaded in python pickle format as recarrays from https://doi.org/10.6084/m9.figshare.11973693 \nThe python code for producing this catalog is available as part of https://github.com/cmoestl/heliocats \n\nEach icmecat_id has a tag in it that indicates from which catalog the ICME times were taken: \n\nWind:       ****************** \nSTEREO-A:   Jian et al. (2019), tag: JIAN. \nSTEREO-B:   Jian et al. (2019), tag: JIAN. \n\nReferences \nJian, L. et al. (2019), doi: 10.1007/s11207-019-1416-8, 2019. \n\nComments: \n- Spacecraft positions are given in Heliocentric Earth Equatorial Coordinates (HEEQ) coordinates. \n- The coordinate system for all magnetic field components is SCEQ, except for Wind (HEEQ, which is the equivalent for SCEQ for Earth). \n        Definition of SpaceCraft Equatorial Coordinates (SCEQ): \n        Z is the solar rotation axis. \n        Y is the cross product of Z and R, with R being the vector that points from the Sun to the spacecraft.\n        X completes the right handed triad (and points away from the Sun). \nThis system is thus like HEEQ but centered on the respective in situ spacecraft, so the SCEQ X and Y \nbase vectors are rotated by the HEEQ longitude of the in situ spacecraft from HEEQ X and Y.\nThe Y vector is similar to the T vector in an RTN system for each spacecraft, but the X and Z vectors \nare rotated around Y compared to an RTN system. The differences between RTN and SCEQ for spacecraft within \na few degrees of the solar equatorial plane are very small (within a few 0.1 nT usually).\nWe choose SCEQ because it has the advantage that a comparison between multipoint CME events \nand for comparison to simulations there is always a similar reference plane (the solar equatorial plane). \n\n '     


                       
parameters='Parameters: \n00: icmecat_id: The unique identifier for the observed ICME. unit: string. \n01: sc insitu: The name of the in situ observing spacecraft. unit: string. \n02: icme_start_time: Shock arrival or density enhancement time, can be similar to mo_start_time. unit: UTC. \n03: mo_start_time: Start time of the magnetic obstacle (MO), including flux ropes, flux-rope-like, and ejecta signatures. unit: UTC. \n04: mo_end_time: End time of the magnetic obstacle. unit: UTC. \n05: mo_sc_heliodistance: Heliocentric distance of the spacecraft at mo_start_time. unit: AU.\n06: mo_sc_long_heeq: Heliospheric longitude of the spacecraft at mo_start_time, range [-180,180]. unit: degree (HEEQ).\n07: mo_sc_lat_heeq: Heliospheric latitude of the spacecraft at mo_start_time, range [-90,90]. unit: degree (HEEQ).\n08: icme_duration: Duration of the interval between icme_start_time and mo_endtime. unit: hours.\n09: icme_bmax: Maximum total magnetic field in the full icme interval (icme_start_time to mo_end_time). unit: nT.\n10: icme_bmean: Mean total magnetic field during the full icme interval (icme_start_time to mo_end_time). unit: nT.\n11: icme_bstd: Standard deviation of the total magnetic field from icme_start_time to mo_end_time. unit: nT.\n12: icme_speed_mean: Mean proton speed from icme_start_time to mo_end_time. unit: km/s.\n13: icme_speed_std: Standard deviation of proton speed from icme_start_time to mo_end_time. unit: km/s.\n14: mo_duration: Duration of the interval between mo_start_time and mo_endtime. unit: hours.\n15: mo_bmax: Maximum total magnetic field in the magnetic obstacle interval (mo_start_time to mo_end_time). unit: nT.\n16: mo_bmean: Mean total magnetic field in the magnetic obstacle. unit: nT.\n17: mo_bstd: Standard deviation of the total magnetic field in the magnetic obstacle. unit: nT.\n18: mo_bzmean: Mean magnetic field Bz component in the magnetic obstacle. unit: nT.\n19: mo_bzmin: Minimum magnetic field Bz component in the magnetic obstacle. unit: nT.\n20: mo_bzstd: Standard deviation of the magnetic field Bz component in the magnetic obstacle. unit: nT.\n21: mo_bymean: Mean magnetic field By component in the magnetic obstacle. unit: nT.\n22: mo_bystd: Standard deviation of the magnetic field By component in the magnetic obstacle. unit: nT.\n23: mo_speed_mean: Mean proton speed from mo_start_time to mo_end_time. unit: km/s.\n24: mo_speed_std: Standard deviation of proton speed from mo_start_time to mo_end_time. unit: km/s.\n25: mo_expansion_speed: Difference between proton speed at mo_start_time to proton speed at mo_end_time. unit: km/s.\n26: mo_pdyn_mean: Mean proton dynamic pressure from mo_start_time to mo_start_time. unit: nPa.\n27: mo_pdyn_std: Standard deviation of proton dynamic pressure from mo_start_time to mo_start_time. unit: nPa.\n28: mo_density_mean: Mean proton density from mo_start_time to mo_start_time. unit: cm^-3.\n29: mo_density_std: Standard deviation of proton density from mo_start_time to mo_start_time. unit: cm^-3.\n30: mo_temperature_mean: Mean proton temperature from mo_start_time to mo_start_time. unit: K.\n31: mo_temperature_std: Standard deviation of proton temperature from mo_start_time to mo_end_time. unit: K.\n32: sheath_speed_mean: Mean proton speed from icme_start_time to mo_start_time, NaN if these times are similar. unit: km/s.\n33: sheath_speed_std: Standard deviation of proton speed from icme_start_time to mo_start_time, NaN if these times are similar. unit: km/s.\n34: sheath_density_mean: Mean proton density from icme_start_time to mo_start_time, NaN if these times are similar. unit: cm^-3.\n35: sheath_density_std: Standard deviation of proton density from icme_start_time to mo_start_time, NaN if these times are similar. unit: cm^-3.\n36: sheath_pdyn_mean: Mean proton dynamic pressure, from icme_start_time to mo_start_time, NaN if these times are similar. unit: nPa.\n37: sheath_pdyn_std: Standard deviation of proton dynamic pressure, from icme_start_time to mo_start_time, NaN if these times are similar. unit: nPa.\n\n\n'



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

# In[79]:


########## python formats

# save ICMECAT as pandas dataframe with times as datetime objects as pickle
file='icmecat/HELCATS_SIRCAT_v20_pandas.p'
pickle.dump([scat,header,parameters], open(file, 'wb'))
print('SIRCAT saved as '+file)



#load icmecat as pandas dataframe
file='icmecat/HELCATS_ICMECAT_v20_pandas.p'
[scat_pandas,h,p]=pickle.load( open(file, 'rb'))   
scat.keys()
scat



sys.exit()


# In[70]:


# save SIRCAT as numpy array with times as matplotlib datetime as pickle
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

# In[10]:


#load icmecat as pandas dataframe
file='icmecat/HELCATS_ICMECAT_v20_pandas.p'
[ic_pandas,h,p]=pickle.load( open(file, 'rb'))   

#load icmecat as numpy array
file='icmecat/HELCATS_ICMECAT_v20_numpy.p'
[ic_nprec,ic_np,h,p]=pickle.load( open(file, 'rb'))   


# In[11]:


ic_pandas
ic_pandas.keys()


# In[12]:


ic_nprec


# In[13]:


ic_nprec


# In[14]:


ic_nprec.icmecat_id


# In[ ]:





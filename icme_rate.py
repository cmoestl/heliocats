#!/usr/bin/env python
# coding: utf-8

# # ICME rate update plots
# 
# adapted from
# https://github.com/helioforecast/Papers/tree/master/Moestl2020_PSP_rate
# makes a prediction of the ICME rate in solar cycle 25
# 
# Main author: C. Moestl, IWF Graz, Austria; twitter @chrisoutofspace; https://github.com/cmoestl
# 
# ssn is automatically loaded from http://www.sidc.be/silso/DATA/SN_d_tot_V2.0.csv
# 
# plots are saved to '/nas/helio/data/insitu_python/icme_rate_cycle_update'
# 
# Convert this notebook to a script with:
# 
# import os
# 
# os.system('jupyter nbconvert --to script icme_rate.ipynb')    
# 
# 
# 
# 
# ---
# 
# **MIT LICENSE**
# 
# Copyright 2020-2022, Christian Moestl
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this 
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify, 
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
# permit persons to whom the Software is furnished to do so, subject to the following 
# conditions:
# 
# The above copyright notice and this permission notice shall be included in all copies 
# or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# In[1]:


#real time updates: icme_rate.py

#import matplotlib

#for server runs
#matplotlib.use('Agg')

##############!!!!for notebook use

#%matplotlib inline
#outputdirectory='results/icme_rate_cycle_update'





#############!!!!!!!!!!!!!!!!!!
#for automatic updates

outputdirectory='/nas/helio/data/insitu_python/icme_rate_cycle_update'


#Convert this notebook to a script with:

#import os
#os.system('jupyter nbconvert --to script icme_rate.ipynb')    


# In[2]:


from scipy import stats
import scipy.io
from matplotlib import cm
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import datetime
from datetime import timedelta
import astropy.constants as const
from sunpy.time import parse_time
import sunpy.time
import time
import pickle
import seaborn as sns
import os
import urllib
import json
import warnings
import importlib
import heliopy.spice as spice
import heliopy.data.spice as spicedata
import astropy
import copy


#our own package
from heliocats import stats as hs
from heliocats import data as hd



#where the 6 in situ data files are located is read from input.py
#as data_path=....
from config import data_path
#reload again while debugging

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()



#%matplotlib inline
#matplotlib.use('Qt5Agg')
#matplotlib.use('Agg')
#warnings.filterwarnings('ignore') # some numpy mean-of-empty-slice runtime warnings

########### make directories first time
resdir='results'
if os.path.isdir(resdir) == False: os.mkdir(resdir)

datadir='data'
if os.path.isdir(datadir) == False: os.mkdir(datadir)

if os.path.isdir(outputdirectory) == False: os.mkdir(outputdirectory)
    
    
#animdirectory='results/plots_rate/anim'
#if os.path.isdir(animdirectory) == False: os.mkdir(animdirectory)
    
#animdirectory2='results/plots_rate/anim2'
#if os.path.isdir(animdirectory2) == False: os.mkdir(animdirectory2)
    
plt.rcParams["figure.figsize"] = (15,8)
print('done')


# ## 1 Settings and load data

# In[3]:


plt.close('all')

print('icme_rate main program.')
print('Christian Moestl et al., IWF Graz, Austria')

#constants: 
#solar radiusx
Rs_in_AU=float(const.R_sun/const.au)
#define AU in km
AU_in_km=const.au.value/1e3

#set for loading
load_data=1
get_new_sunspots=1
get_new_sunspots_ms=1
get_new_sunspots_mean=1


if load_data > 0:
    
    print('load data (takes a minute or so)')
    print('')
    
    #####################
    print('get RC ICME list')    
    #download richardson and cane list
    rc_url='http://www.srl.caltech.edu/ACE/ASC/DATA/level3/icmetable2.htm'

    try: urllib.request.urlretrieve(rc_url,data_path+'rc_list.htm')
    except urllib.error.URLError as e:
        print('Failed downloading ', rc_url,' ',e)

    #read RC list into pandas dataframe    
    rc_dfull=pd.read_html(data_path+'rc_list.htm')
    rc_df=rc_dfull[0]
  

    ##################
    print('get sunspot number from SIDC')    
    #get daily sunspot number from SIDC
    #http://www.sidc.be/silso/datafiles
    #parameters
    #http://www.sidc.be/silso/infosndtot
    #daily sunspot number

    if get_new_sunspots==1:
        ssn=pd.read_csv('http://www.sidc.be/silso/DATA/SN_d_tot_V2.0.csv',sep=';',names=['year','month','day','year2','spot','stand','obs','prov'])     
        
   
        ssn_time=np.zeros(len(ssn))
        for k in np.arange(len(ssn)):
            ssn_time[k]=parse_time(str(ssn.year[k])+'-'+str(ssn.month[k])+'-'+str(ssn.day[k])).plot_date

        print('time convert done')
        ssn.insert(0,'time',ssn_time)
        ssn.spot.loc[np.where(ssn.spot< 0)[0]]=np.nan    
        ssn.stand.loc[np.where(ssn.stand< 0)[0]]=np.nan    
        fileout='ssn.p'
        pickle.dump(ssn, open(data_path+fileout, "wb"))
        

        #also get preliminary data for current month for plotting    
        #does not work
        #ssn_prelim_raw=pd.csv('http://www.sidc.be/silso/DATA/EISN/EISN_current.csv',sep=',',names=['year','month','day','year2','spot','stand','stat calc','stat avail'])       
        
        
        #download manually
        ssn_prelim_url='http://www.sidc.be/silso/DATA/EISN/EISN_current.csv'
        try: urllib.request.urlretrieve(ssn_prelim_url,'data/EISN_current.csv')
        except urllib.error.URLError as e:
            print('Failed downloading ', ssn_prelim_url,' ',e)
        
        ssn_prelim_raw = np.loadtxt('data/EISN_current.csv', delimiter=',',usecols=(0,1,2,4,5))
        ssn_p_int=ssn_prelim_raw.astype(int)
        ssn_p=pd.DataFrame(ssn_p_int,columns=['year','month','day','spot','stand'])        
        
        ssn_p_time=np.zeros(len(ssn_p))
        for k in np.arange(len(ssn_p)):
            ssn_p_time[k]=parse_time(str(ssn_p.year[k])+'-'+str(ssn_p.month[k])+'-'+str(ssn_p.day[k])).plot_date
        
        ssn_p.insert(0,'time',ssn_p_time)
        ssn_p.spot.loc[np.where(ssn_p.spot< 0)[0]]=np.nan    
        ssn_p.stand.loc[np.where(ssn_p.stand< 0)[0]]=np.nan    
    

        fileout='ssn_prelim.p'
        pickle.dump(ssn_p, open(data_path+fileout, "wb"))



        

    if get_new_sunspots_ms==1:

        #ssn_ms=pd.read_csv('data/SN_ms_tot_V2.0.csv',sep=';')
        ssn_ms=pd.read_csv('http://www.sidc.be/silso/DATA/SN_ms_tot_V2.0.csv',sep=';',names=['year','month','year2','spot','stand','obs','check'])       


        ssn_ms_time=np.zeros(len(ssn_ms))

        for k in np.arange(len(ssn_ms)):
            ssn_ms_time[k]=parse_time(str(ssn_ms.year[k])+'-'+str(ssn_ms.month[k])+'-01').plot_date
            #print(mdates.num2date(ssn_ms_time[k]))
            #print(ssn_ms.spot[k])

        print('time convert done')
        ssn_ms.insert(0,'time',ssn_ms_time)
        ssn_ms.spot.loc[np.where(ssn_ms.spot< 0)[0]]=np.nan    
    
        fileout='ssn_13ms.p'
        pickle.dump(ssn_ms, open(data_path+fileout, "wb"))
        
    
    if get_new_sunspots_mean==1:

       
        ssn_m=pd.read_csv('http://www.sidc.be/silso/DATA/SN_m_tot_V2.0.csv',sep=';',names=['year','month','year2','spot','stand','obs','check'])       

        ssn_m_time=np.zeros(len(ssn_m))

        for k in np.arange(len(ssn_m)):
            ssn_m_time[k]=parse_time(str(ssn_m.year[k])+'-'+str(ssn_m.month[k])+'-01').plot_date

        print('time convert done')
        ssn_m.insert(0,'time',ssn_m_time)
        ssn_m.spot.loc[np.where(ssn_m.spot< 0)[0]]=np.nan    
    
        fileout='ssn_m.p'
        pickle.dump(ssn_m, open(data_path+fileout, "wb"))

        

    file='ssn.p'  
    ssn=pickle.load(open(data_path+file, "rb"))
    
    #make 13 month running mean
    runmean_months=13.0
    ssn_mean_13=hs.running_mean(ssn.spot,int(np.rint(30.42*runmean_months+1)))
    ssn_std_13=hs.running_mean(ssn.stand,int(np.rint(30.42*runmean_months+1)))    

    ssn.insert(1,'spot_mean_13',ssn_mean_13)    
    ssn.insert(2,'spot_std_13',ssn_std_13)    
    
    print('SIDC sunspots done')

    
    ################## Spacecraft

    #completed missions
    
    
    filevex='vex_2007_2014_sceq_removed.p'
    [vex,hvex]=pickle.load(open(data_path+filevex, 'rb' ) )

    filevex='vex_2007_2014_sceq.p'
    [vexnon,hvexnon]=pickle.load(open(data_path+filevex, 'rb' ) )

    filemes='messenger_2007_2015_sceq_removed.p'
    [mes,hmes]=pickle.load(open(data_path+filemes, 'rb' ) )

    filemes='messenger_2007_2015_sceq.p'
    [mesnon,hmesnon]=pickle.load(open(data_path+filemes, 'rb' ) )

    filestb='stereob_2007_2014_sceq.p'
    [stb,hstb]=pickle.load(open(data_path+filestb, "rb" ) )      


    fileuly='ulysses_1990_2009_rtn.p'
    [uly,huly]=pickle.load(open(data_path+fileuly, "rb" ) ) 
    
    
    
    
    ##active mission
    
    filemav='maven_2014_2018_removed_smoothed.p'
    [mav,hmav]=pickle.load(open(data_path+filemav, 'rb' ) )
    
    

    print('load and merge Wind data HEEQ') 
    #from HELCATS HEEQ until 2018 1 1 + new self-processed data with heliosat and hd.save_wind_data
    filewin="wind_2007_2018_heeq_helcats.p" 
    [win1,hwin1]=pickle.load(open(data_path+filewin, "rb" ) )  
    
    #or use: filewin2="wind_2018_now_heeq.p" 
    filewin2="wind_2018_now_heeq.p" 
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

    
    ########### STA
    
    print('load and merge STEREO-A data SCEQ') #yearly magplasma files from stereo science center, conversion to SCEQ 
    filesta1='stereoa_2007_2020_sceq.p'
    sta1=pickle.load(open(data_path+filesta1, "rb" ) )  
    
    #beacon data
    #filesta2="stereoa_2019_2020_sceq_beacon.p"
    #filesta2='stereoa_2019_2020_sept_sceq_beacon.p'
    #filesta2='stereoa_2019_now_sceq_beacon.p'
    #filesta2="stereoa_2020_august_november_sceq_beacon.p" 
    filesta2='stereoa_2020_now_sceq_beacon.p'
    
    [sta2,hsta2]=pickle.load(open(data_path+filesta2, "rb" ) )  
    #cutoff with end of science data
    sta2=sta2[np.where(sta2.time >= parse_time('2020-Aug-01 00:00').datetime)[0]]

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

       

    print('load PSP data SCEQ') #from heliosat, converted to SCEQ similar to STEREO-A/B
    filepsp='psp_2018_2021_sceq.p'
    [psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) )     
    
        
    ##############################################
    print('load Bepi Colombo SCEQ')
    filebepi='bepi_2019_2021_sceq.p'
    bepi=pickle.load(open(data_path+filebepi, "rb" ) )      
   
    ##############################################
    print('load Solar Orbiter SCEQ')
    filesolo='solo_2020_april_2021_sep_sceq.p'
    solo=pickle.load(open(data_path+filesolo, "rb" ) )    
    #set all plasma data to NaN
    solo.vt=np.nan
    solo.np=np.nan
    solo.tp=np.nan
    

    fileomni='omni_1963_2020.p'
    [omni,homni]=pickle.load(open(data_path+fileomni, "rb" ) ) 

    print('load all data done')
    


# In[4]:



########### load ICMECAT v2.0, made with icmecat.py or ipynb
file='icmecat/HELIO4CAST_ICMECAT_v21_pandas.p'
print()
print('loaded ', file)
print()
print('Keys (parameters) in this pandas data frame are:')

[ic,h,p]=pickle.load(open(file, "rb" ) )  
print(ic.keys())
print()

################### get indices of events for each spacecraft

mercury_orbit_insertion_time= parse_time('2011-03-18').datetime

#spacecraft near the 4 terrestrial planets
#get indices for Mercury after orbit insertion in March 2011
merci=np.where(np.logical_and(ic.sc_insitu =='MESSENGER', ic.icme_start_time > mercury_orbit_insertion_time))[0]
vexi=np.where(ic.sc_insitu == 'VEX')[:][0]  
wini=np.where(ic.sc_insitu == 'Wind')[:][0] 
mavi=np.where(ic.sc_insitu == 'MAVEN')[:][0]    

#other spacecraft
#all MESSENGER events including cruise phase
mesi=np.where(ic.sc_insitu == 'MESSENGER')[:][0]   
stbi=np.where(ic.sc_insitu == 'STEREO-B')[:][0]    
ulyi=np.where(ic.sc_insitu == 'ULYSSES')[:][0]   


pspi=np.where(ic.sc_insitu == 'PSP')[:][0]    
soli=np.where(ic.sc_insitu == 'SolarOrbiter')[:][0]    
beci=np.where(ic.sc_insitu == 'BepiColombo')[:][0]    
stai=np.where(ic.sc_insitu == 'STEREO-A')[:][0]    

############### set limits of solar minimum, rising/declining phase and solar maximum

# minimim maximum times as given by
#http://www.sidc.be/silso/cyclesmm
#24    2008   12    2.2   2014 04   116.4  

solarmin=parse_time('2008-12-01').datetime
minstart=solarmin-datetime.timedelta(days=366*1.5)
minend=solarmin+datetime.timedelta(days=365)
minstart_num=parse_time(minstart).plot_date
minend_num=parse_time(minend).plot_date

solarmax=parse_time('2014-04-01').datetime
maxstart=solarmax-datetime.timedelta(days=365*3)
maxend=solarmax+datetime.timedelta(days=365/2)
maxstart_num=parse_time(maxstart).plot_date
maxend_num=parse_time(maxend).plot_date


#rising phase not used
# risestart=parse_time('2010-01-01').datetime
# riseend=parse_time('2011-06-30').datetime
# risestart_num=parse_time('2010-01-01').plot_date
# riseend_num=parse_time('2011-06-30').plot_date

# declstart=parse_time('2015-01-01').datetime
# declend=parse_time('2018-12-31').datetime
# declstart_num=parse_time('2015-01-01').plot_date
# declend_num=parse_time('2018-12-31').plot_date


############### extract events by limits of solar minimum and  maximum
iall_min=np.where(np.logical_and(ic.icme_start_time > minstart,ic.icme_start_time < minend))[0]
#iall_rise=np.where(np.logical_and(ic.icme_start_time > risestart,ic.icme_start_time < riseend))[0]
iall_max=np.where(np.logical_and(ic.icme_start_time > maxstart,ic.icme_start_time < maxend))[0]

wini_min=iall_min[np.where(ic.sc_insitu[iall_min]=='Wind')]
#wini_rise=iall_rise[np.where(ic.sc_insitu[iall_rise]=='Wind')]
wini_max=iall_max[np.where(ic.sc_insitu[iall_max]=='Wind')]

pspi_min=iall_min[np.where(ic.sc_insitu[iall_min]=='PSP')]
#wini_rise=iall_rise[np.where(ic.sc_insitu[iall_rise]=='Wind')]
pspi_max=iall_max[np.where(ic.sc_insitu[iall_max]=='PSP')]


vexi_min=iall_min[np.where(ic.sc_insitu[iall_min]=='VEX')]
#vexi_rise=iall_rise[np.where(ic.sc_insitu[iall_rise]=='VEX')]
vexi_max=iall_max[np.where(ic.sc_insitu[iall_max]=='VEX')]

mesi_min=iall_min[np.where(ic.sc_insitu[iall_min]=='MESSENGER')]
#mesi_rise=iall_rise[np.where(ic.sc_insitu[iall_rise]=='MESSENGER')]
mesi_max=iall_max[np.where(ic.sc_insitu[iall_max]=='MESSENGER')]

stai_min=iall_min[np.where(ic.sc_insitu[iall_min]=='STEREO-A')]
#stai_rise=iall_rise[np.where(ic.sc_insitu[iall_rise]=='STEREO-A')]
stai_max=iall_max[np.where(ic.sc_insitu[iall_max]=='STEREO-A')]

stbi_min=iall_min[np.where(ic.sc_insitu[iall_min]=='STEREO-B')]
#stbi_rise=iall_rise[np.where(ic.sc_insitu[iall_rise]=='STEREO-B')]
stbi_max=iall_max[np.where(ic.sc_insitu[iall_max]=='STEREO-B')]

# select the events at Mercury extra after orbit insertion, note that no events available for solar minimum
merci_min=iall_min[np.where(np.logical_and(ic.sc_insitu[iall_min] =='MESSENGER',ic.icme_start_time[iall_min] > parse_time('2011-03-18').datetime))]
#merci_rise=iall_rise[np.where(np.logical_and(ic.sc_insitu[iall_rise] =='MESSENGER',ic.icme_start_time[iall_rise] > parse_time('2011-03-18').datetime))]
merci_max=iall_max[np.where(np.logical_and(ic.sc_insitu[iall_max] =='MESSENGER',ic.icme_start_time[iall_max] > parse_time('2011-03-18').datetime))]

print(len(ic))
print('done')


# In[5]:


ic


# ## 2  ICME rate for solar cycles 23/24 from the Heliophysics System Observatory (ICMECAT and Richardson and Cane)

# ### Check data days available each year for each planet or spacecraft

# In[6]:


######################## make bin for each year for yearly histograms
#define dates of January 1 from 2007 to end year

last_year=2022 #2022 means last date is 2021 Dec 31

years_jan_1_str=[str(i)+'-01-01' for i in np.arange(2007,last_year) ] 
yearly_start_times=parse_time(years_jan_1_str).datetime
yearly_start_times_num=parse_time(years_jan_1_str).plot_date

#same for July 1 as middle of the year
years_jul_1_str=[str(i)+'-07-01' for i in np.arange(2007,last_year) ] 
yearly_mid_times=parse_time(years_jul_1_str).datetime
yearly_mid_times_num=parse_time(years_jul_1_str).plot_date

#same for december 31
years_dec_31_str=[str(i)+'-12-31' for i in np.arange(2007,last_year) ] 
yearly_end_times=parse_time(years_dec_31_str).datetime
yearly_end_times_num=parse_time(years_dec_31_str).plot_date



########### define arrays for total data days and fill with nan
total_data_days_yearly_win=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_win.fill(np.nan)


total_data_days_yearly_psp=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_psp.fill(np.nan)

total_data_days_yearly_solo=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_solo.fill(np.nan)

total_data_days_yearly_bepi=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_bepi.fill(np.nan)

total_data_days_yearly_sta=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_sta.fill(np.nan)


total_data_days_yearly_stb=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_stb.fill(np.nan)

total_data_days_yearly_mes=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_mes.fill(np.nan)

total_data_days_yearly_vex=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_vex.fill(np.nan)

total_data_days_yearly_mav=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_mav.fill(np.nan)


######################## go through each year and search for available data
#time is available for all dates, so there are no NaNs in time, thus need to search for all not NaNs in Btotal variable

for i in range(np.size(yearly_mid_times)):

    print(yearly_start_times[i])

    #get indices of Wind time for the current year
    thisyear=np.where(np.logical_and((win.time > yearly_start_times[i]),(win.time < yearly_end_times[i])))[0]
    #get np.size of available data for each year
    datas=np.size(np.where(np.isnan(win.bt[thisyear])==False))
    #wind is  in 1 minute resolution until 31 Dec 2017, from 1 Jan 2018 its 2 minutes
    min_in_days=1/(60*24)
    if i > 10: min_in_days=1/(60*24)  
    #calculate available days from number of datapoints (each 1 minute) 
    #divided by number of minutes in 1 days
    #this should only be the case if data is available this year, otherwise set to NaN
    if datas > 0: total_data_days_yearly_win[i]=datas*min_in_days
   
    #manual override because Wind data for 2018 and 2019 are heavily despiked

    total_data_days_yearly_win[-4]=360
    total_data_days_yearly_win[-3]=360
    total_data_days_yearly_win[-2]=360
    total_data_days_yearly_win[-1]=180

    

    #all other data is in 1 min resolution      
    min_in_days=1/(60*24)
    
    
    #for PSP
    thisyear=np.where(np.logical_and((psp.time > yearly_start_times[i]),(psp.time < yearly_end_times[i])))[0]
    datas=np.size(np.where(np.isnan(psp.bt[thisyear])==False))
    if datas >0: total_data_days_yearly_psp[i]=datas*min_in_days
        
    #for Bepi
    thisyear=np.where(np.logical_and((bepi.time > yearly_start_times[i]),(bepi.time < yearly_end_times[i])))[0]
    datas=np.size(np.where(np.isnan(bepi.bt[thisyear])==False))
    if datas >0: total_data_days_yearly_bepi[i]=datas*min_in_days

    #for solo
    thisyear=np.where(np.logical_and((solo.time > yearly_start_times[i]),(solo.time < yearly_end_times[i])))[0]
    datas=np.size(np.where(np.isnan(solo.bt[thisyear])==False))
    if datas >0: total_data_days_yearly_solo[i]=datas*min_in_days

        
    #same for STEREO-A
    thisyear=np.where(np.logical_and((sta.time > yearly_start_times[i]),(sta.time < yearly_end_times[i])))[0]
    datas=np.size(np.where(np.isnan(sta.bt[thisyear])==False))
    if datas >0: total_data_days_yearly_sta[i]=datas*min_in_days

    #same for STEREO-B
    thisyear=np.where(np.logical_and((stb.time > yearly_start_times[i]),(stb.time < yearly_end_times[i])))[0]
    datas=np.size(np.where(np.isnan(stb.bt[thisyear])==False))
    if datas >0: total_data_days_yearly_stb[i]=datas*min_in_days

    #same for MESSENGER
    thisyear=np.where(np.logical_and((mesnon.time > yearly_start_times[i]),(mesnon.time < yearly_end_times[i])))
    datas=np.size(np.where(np.isnan(mes.bt[thisyear])==False))
    if datas >0: total_data_days_yearly_mes[i]=datas*min_in_days

    #same for Mercury alone with non-removed dataset
    #start with 2011
    #   if i == 4: 
    #    thisyear=np.where(np.logical_and((mesnon.time > mercury_orbit_insertion_time),(mesnon.time < yearly_end_times[i])))[0]
    #    datas=np.size(np.where(np.isnan(mesnon.bt[thisyear])==False))
    #    if datas >0: total_data_days_yearly_merc[i]=datas*min_in_days
    #   #2012 onwards 
    #   if i > 4: 
    #    thisyear=np.where(np.logical_and((mesnon.time > yearly_start_times[i]),(mesnon.time < yearly_end_times[i])))
    #    datas=np.size(np.where(np.isnan(mesnon.bt[thisyear])==False))
    #    if datas >0: total_data_days_yearly_merc[i]=datas*min_in_days
      

    #same for VEX
    thisyear=np.where(np.logical_and((vexnon.time > yearly_start_times[i]),(vexnon.time < yearly_end_times[i])))[0]
    datas=np.size(np.where(np.isnan(vexnon.bt[thisyear])==False))
    if datas >0: total_data_days_yearly_vex[i]=datas*min_in_days

    #for MAVEN different time resolution
    thisyear=np.where(np.logical_and((mav.time > yearly_start_times[i]),(mav.time < yearly_end_times[i])))[0]
    datas=np.size(np.where(np.isnan(mav.bt[thisyear])==False))
    datas_ind=np.where(np.isnan(mav.bt[thisyear])==False)
    #sum all time intervals for existing data points, but avoid counting gaps where diff is > 1 orbit (0.25 days)
    alldiff=np.diff(parse_time(mav.time[datas_ind]).plot_date)
    smalldiff_ind=np.where(alldiff <0.25)  
    if datas >0: total_data_days_yearly_mav[i]=np.sum(alldiff[smalldiff_ind])

        
        
print('Data days each year:')

print()
print('MESSENGER')
print(np.round(total_data_days_yearly_mes,1))
print()
print('VEX at Venus')
print(np.round(total_data_days_yearly_vex,1))
print()
print('STB')
print(np.round(total_data_days_yearly_stb,1))
print()
print('MAVEN')
print(np.round(total_data_days_yearly_mav,1))

print()

print()


print('Wind')
print(np.round(total_data_days_yearly_win,1))

print('STA')
print(np.round(total_data_days_yearly_sta,1))

print('PSP')
print(np.round(total_data_days_yearly_psp,1))


print('Bepi')
print(np.round(total_data_days_yearly_bepi,1))

print('Solar Orbiter')
print(np.round(total_data_days_yearly_solo,1))


print()
print('done')


# ### get yearly ICME rates at each spacecraft

# In[7]:


#define dates of January 1 from 2007 to 2022
years_jan_1_str_plus1=[str(i)+'-01-01' for i in np.arange(2007,last_year+1) ] 
yearly_bin_edges=parse_time(years_jan_1_str_plus1).plot_date
#bin width in days         
binweite=365/8

(histmes1, bin_edgesmes) = np.histogram(parse_time(ic.icme_start_time[mesi]).plot_date, yearly_bin_edges)
(histvex1, bin_edgesvex) = np.histogram(parse_time(ic.icme_start_time[vexi]).plot_date, yearly_bin_edges)
(histmav1, bin_edgesmav) = np.histogram(parse_time(ic.icme_start_time[mavi]).plot_date, yearly_bin_edges)
(histstb1, bin_edgesstb) = np.histogram(parse_time(ic.icme_start_time[stbi]).plot_date, yearly_bin_edges)

(histwin1, bin_edgeswin) = np.histogram(parse_time(ic.icme_start_time[wini]).plot_date, yearly_bin_edges)
(histsta1, bin_edgessta) = np.histogram(parse_time(ic.icme_start_time[stai]).plot_date, yearly_bin_edges)
(histpsp1, bin_edgespsp) = np.histogram(parse_time(ic.icme_start_time[pspi]).plot_date, yearly_bin_edges)
(histsolo1, bin_edgessolo) = np.histogram(parse_time(ic.icme_start_time[soli]).plot_date, yearly_bin_edges)
(histbepi1, bin_edgesbepi) = np.histogram(parse_time(ic.icme_start_time[beci]).plot_date, yearly_bin_edges)


binedges=bin_edgeswin



#normalize each dataset for data gaps, so correcting ICME rate for actual data availability
#note that for VEX and MESSENGER this was done with the non-removed datasets (vexnon, mesnon)
histvex=np.round(histvex1/total_data_days_yearly_vex*365.24,1)
histmes=np.round(histmes1/total_data_days_yearly_mes*365.24,1)

#ok for these spacecraft as continously in the solar wind and the MAVEN data set is made without orbit gaps
histsta=np.round(histsta1/total_data_days_yearly_sta*365.24,1)


#STA beacon data used in 2019 - set manually
histsta[-3]=13
#STA beacon data used in 2020 - set manually
histsta[-2]=10*4/3



histmav=np.round(histmav1/total_data_days_yearly_mav*365.24,1)
histmav[7]=np.nan #not enough data for 2014

histstb=np.round(histstb1/total_data_days_yearly_stb*365.24,1)
histwin=np.round(histwin1/total_data_days_yearly_win*365.24,1)
histpsp=np.round(histpsp1/total_data_days_yearly_psp*365.24,1)
histsolo=np.round(histsolo1/total_data_days_yearly_solo*365.24,1)
histbepi=np.round(histbepi1/total_data_days_yearly_bepi*365.24,1)




print('corrected ICME rates for years')
print(yearly_mid_times)
print('MESSENGER',histmes)
print('VEX',histvex)
print('Wind',histwin)
print('PSP',histpsp)
print('STA',histsta)
print('STB',histstb)
print('MAVEN',histmav)


sns.set_context("talk")     
sns.set_style('darkgrid')

plt.figure(11,figsize=(12,6),dpi=60)

plt.ylabel('ICMEs per year, ICMECAT')
plt.plot(yearly_mid_times,histmes,'-',label='MESSENGER')
plt.plot(yearly_mid_times,histvex,'-',label='VEX')
plt.plot(yearly_mid_times,histwin,'-',label='Wind')
plt.plot(yearly_mid_times,histsta,'-',label='STA')
plt.plot(yearly_mid_times,histstb,'-',label='STB')
plt.plot(yearly_mid_times,histmav,'-',label='MAVEN')
plt.plot(yearly_mid_times,histpsp,'-',label='PSP')
plt.plot(yearly_mid_times,histsolo,'-',label='SolarOrbiter')
plt.plot(yearly_mid_times,histbepi,'-',label='BepiColombo')



plt.legend(loc=1,fontsize=10)




################### calculate general parameters


print()
print()
print('calculate ICME rate matrix std, mean, max, min')

#arrange icmecat rate data so each row contains std, mean, max, min

icrate=pd.DataFrame(np.zeros([len(yearly_mid_times),6]), columns=['year','std1', 'median1','mean1', 'max1', 'min1'] )


for i in np.arange(0,len(yearly_mid_times)):
    
    icrate.at[i,'year']=yearly_mid_times_num[i]
    icrate.at[i,'median1']=np.round(np.nanmedian([histwin[i],histvex[i],histsta[i],histstb[i],histmes[i],histmav[i],histpsp[i]]),1)
    icrate.at[i,'mean1']=np.round(np.nanmean([histwin[i],histvex[i],histsta[i],histstb[i],histmes[i],histmav[i],histpsp[i]]),1)
    icrate.at[i,'std1']=np.round(np.nanstd([histwin[i],histvex[i],histsta[i],histstb[i],histmes[i],histmav[i],histpsp[i]]),1)
    icrate.at[i,'max1']=np.nanmax([histwin[i],histvex[i],histsta[i],histstb[i],histmes[i],histmav[i],histpsp[i]])
    icrate.at[i,'min1']=np.nanmin([histwin[i],histvex[i],histsta[i],histstb[i],histmes[i],histmav[i],histpsp[i]])

#change matplotlib time before plotting
icrate_year2=icrate.year + mdates.date2num(np.datetime64('0000-12-31'))
    
plt.plot([icrate_year2,icrate_year2],[icrate.mean1-icrate.std1,icrate.mean1+icrate.std1],'-k',lw=0.5)
plt.plot(icrate_year2,icrate.mean1,'ok',markerfacecolor='white')

# icrate=pd.DataFrame(np.zeros([len(histvex)*6,3]), columns=['year','rate','sc'] )
# #write all icme rates into this array
# icrate.at[0:12,'rate']=histmes
# icrate.at[0:12,'year']=yearly_start_times_num
# icrate.at[0:12,'sc']='MESSENGER'

# icrate.at[13:25,'rate']=histvex
# icrate.at[13:25,'year']=yearly_start_times_num
# icrate.at[13:25,'sc']='VEX'

# icrate.at[26:38,'rate']=histwin
# icrate.at[26:38,'year']=yearly_start_times_num
# icrate.at[26:38,'sc']='Wind'

# icrate.at[39:51,'rate']=histvex
# icrate.at[39:51,'year']=yearly_start_times_num
# icrate.at[39:51,'sc']='STA'

# icrate.at[52:64,'rate']=histvex
# icrate.at[52:64,'year']=yearly_start_times_num
# icrate.at[52:64,'sc']='STB'

# icrate.at[65:77,'rate']=histvex
# icrate.at[65:77,'year']=yearly_start_times_num
# icrate.at[65:77,'sc']='MAVEN'

# sns.boxplot(x='year',y='rate',data=icrate)



icrate


# ### get Richardson and Cane ICME rate for comparison

# In[8]:


#convert times in dataframe from richardson and cane list to numpy array
r1=np.array(rc_df['Disturbance Y/M/D (UT) (a)'])

#to get ICME rate, go through all rows
rc_year=np.zeros(len(r1))

#extract string and check whether its a viable float and non nan:
for p in np.arange(0,len(r1)):    
    rc_yearstr=str(r1[p,0])    
    if hs.is_float(rc_yearstr[0:4]):
        if np.isfinite(float(rc_yearstr[0:4])):
            rc_year[p]=float(rc_yearstr[0:4]) #rc_year contains all ICME

rc_year.sort() 
rc_icme_per_year=np.trim_zeros(rc_year)
#print(rc_year)


#plot check whats in this array

sns.set_style('darkgrid')
fig=plt.figure(12,figsize=(12,5),dpi=80)
ax11=sns.distplot(rc_icme_per_year,bins=24,kde=None)
plt.ylabel('ICMEs per year, RC list')

#count all full years from 1996-2021
bins_years=2022-1996

#get yearly ICME rate (use range to get correct numbers)
rc_rate_values=np.histogram(rc_icme_per_year,bins=bins_years,range=(1996,2022))[0]
rc_rate_time=np.histogram(rc_icme_per_year,bins=bins_years,range=(1996,2022))[1][0:-1]

print(rc_rate_values)
print(rc_rate_time)

years_jul_1_str_rc=[str(i)+'-07-01' for i in np.arange(1996,2021) ] 
yearly_mid_times_rc=parse_time(years_jul_1_str_rc).datetime
yearly_mid_times_num_rc=parse_time(years_jul_1_str_rc).plot_date
print(yearly_mid_times_rc)
#plt.figure(2)
#plt.plot(yearly_mid_times_rc,rc_rate_values)


# ### **Figure 1** plot ICME frequency cycle 24

# In[9]:


sns.set_context("talk")     
#sns.set_style('whitegrid',{'grid.linestyle': '--'})

sns.set_style("ticks",{'grid.linestyle': '--'})
fsize=15

fig=plt.figure(1,figsize=(12,10),dpi=80)

######################## Fig 1a - sc positions during ICMEs vs time



ax1 = plt.subplot(211) 
msize=5
plt.plot_date(ic.icme_start_time[mesi],ic.mo_sc_heliodistance[mesi],fmt='o',color='coral',markersize=msize,label='MESSENGER')
plt.plot_date(ic.icme_start_time[vexi],ic.mo_sc_heliodistance[vexi],fmt='o',color='orange',markersize=msize,label='VEX')
plt.plot_date(ic.icme_start_time[wini],ic.mo_sc_heliodistance[wini],fmt='o',color='mediumseagreen',markersize=msize,label='Wind')
plt.plot_date(ic.icme_start_time[stai],ic.mo_sc_heliodistance[stai],fmt='o',color='red',markersize=msize,label='STEREO-A')
plt.plot_date(ic.icme_start_time[stbi],ic.mo_sc_heliodistance[stbi],fmt='o',color='royalblue',markersize=msize,label='STEREO-B')
plt.plot_date(ic.icme_start_time[mavi],ic.mo_sc_heliodistance[mavi],fmt='o',color='steelblue',markersize=msize,label='MAVEN')
plt.plot_date(ic.icme_start_time[pspi],ic.mo_sc_heliodistance[pspi],fmt='o',color='black',markersize=msize,label='PSP')

plt.plot_date(ic.icme_start_time[soli],ic.mo_sc_heliodistance[soli],'o',c='black',markerfacecolor='white', markersize=msize,label='Solar Orbiter')
plt.plot_date(ic.icme_start_time[beci],ic.mo_sc_heliodistance[beci],'s',c='darkblue',markerfacecolor='lightgrey', markersize=msize, label='BepiColombo')


plt.legend(loc='upper left',fontsize=10)

plt.ylabel('Heliocentric distance R [AU]',fontsize=fsize)
plt.xticks(yearly_start_times,fontsize=fsize) 

#change matplotlib time before plotting
yearly_bin_edges2=yearly_bin_edges + mdates.date2num(np.datetime64('0000-12-31'))

plt.xlim(yearly_bin_edges2[0],yearly_bin_edges2[-1])
plt.ylim([0,1.7])
plt.yticks(np.arange(0,1.9,0.2),fontsize=fsize) 

ax1.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)

#grid for icme rate
for i in np.arange(0,2.0,0.2):
    ax1.plot([datetime.datetime(2007,1,1),datetime.datetime(2024,1,1)],np.zeros(2)+i,linestyle='--',color='grey',alpha=0.5,lw=0.8,zorder=0)
    
for i in np.arange(0,14):    
    ax1.plot([yearly_start_times[i],yearly_start_times[i]],[0,2],linestyle='--',color='grey',alpha=0.5,lw=0.8,zorder=0)


    
    

#################### Fig 1b
sns.set_style("ticks",{'grid.linestyle': '--'})

ax2 = plt.subplot(212) 

ax3=ax2.twinx()

#change matplotlib time before plotting
#ssn_time2=ssn.time + mdates.date2num(np.datetime64('0000-12-31'))

ax3.plot(ssn_ms.time,ssn_ms.spot,'-k',alpha=0.5,linewidth=1.5,label='monthly smoothed sunspot number',zorder=0)
ax3.set_ylabel('Sunspot number SIDC')
ax3.set_ylim(0,155)
ax3.legend(loc=1,fontsize=10)

#grid for icme rate
for i in np.arange(0,50,10):
    ax2.plot([datetime.datetime(2007,1,1),datetime.datetime(2023,1,1)],np.zeros(2)+i,linestyle='--',color='k',alpha=0.4,lw=0.8,zorder=0)

binweite=int(np.round(360/8))
bin_edges=bin_edgeswin[:-1]

#change matplotlib time before plotting
bin_edges2=bin_edges + mdates.date2num(np.datetime64('0000-12-31'))

alp=0.8
ax2.bar(bin_edges2+5+binweite,histmes, width=binweite,color='coral', alpha=alp,label='MESSENGER')
ax2.bar(bin_edges2+5+binweite*2,histvex, width=binweite,color='orange', alpha=alp,label='VEX')
ax2.bar(bin_edges2+5+binweite*3,histstb, width=binweite,color='royalblue', alpha=alp,label='STEREO-B')

ax2.bar(bin_edges2+5+binweite*2,histbepi, width=binweite,color='lightgrey', alpha=alp,label='BepiColombo')
ax2.bar(bin_edges2+5+binweite*3,histsolo, width=binweite,color='white', edgecolor='black', alpha=alp,label='Solar Orbiter')
ax2.bar(bin_edges2+5+binweite*4,histwin, width=binweite,color='mediumseagreen', alpha=alp,label='Wind')
ax2.bar(bin_edges2+5+binweite*5,histsta, width=binweite,color='red', alpha=alp,label='STEREO-A')
ax2.bar(bin_edges2+5+binweite*6,histpsp, width=binweite,color='black', alpha=alp,label='PSP')
ax2.bar(bin_edges2+5+binweite*7,histmav, width=binweite,color='steelblue', alpha=alp,label='MAVEN')



#ax2.bar(bin_edgeswin[:-1]+5+binweite*3,rc_rate_values[-14:-1], width=binweite,color='darkgreen', alpha=0.8,label='Wind')




#ax2.boxplot(histmes)
#RC values
#ax2.plot(bin_edgeswin[:-1]+5+binweite*3,rc_rate_values[-14:-1],'ok',markerfacecolor='white',marker='o',markersize=5,label='Earth RC list')


#mean and standard deviation and min max
ax2.plot([icrate_year2,icrate_year2],[icrate.mean1-icrate.std1,icrate.mean1+icrate.std1],'-k',lw=1.1)

#ax2.plot([icrate.year,icrate.year],[icrate.min1,icrate.max1],'--k',lw=1.1)
ax2.plot(icrate_year2,icrate.mean1,'ok',markerfacecolor='white',label='yearly ICME rate mean',zorder=3)
ax2.plot([icrate_year2[1],icrate_year2[1]],[icrate.mean1[1]-icrate.std1[1],icrate.mean1[1]+icrate.std1[1]],'--k',lw=1.1,label='yearly ICME rate std')

ax2.set_ylim(0,48)
ax2.set_xlim(yearly_bin_edges2[0],yearly_bin_edges2[-1])
ax2.legend(loc=2,fontsize=10)

fsize=15
ax2.set_ylabel('normalized ICME rate per year',fontsize=fsize)
#ax2.set_yticks(fontsize=fsize) 

ax2.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax2.xaxis.set_major_formatter(myformat)
plt.xticks(yearly_start_times, fontsize=fsize) 
plt.xlabel('Year',fontsize=fsize)

plt.tight_layout()

#plt.annotate('(a)',[0.01,0.96],xycoords='figure fraction')
#plt.annotate('(b)',[0.01,0.47],xycoords='figure fraction')

#plt.savefig('results/cycle25_icme_rate.pdf', dpi=100)
plt.savefig(outputdirectory+'/icmecat_icme_rate.png', dpi=100)


# 
# 
# 
# 
# 

# # 3 get solar cycle results on ICME rates and sunspot numbers

# ## solar cycle 23

# In[10]:


print('cycle 23\n')


############################# times
print('times:')
#these years cover solar cycle 23
years23=np.arange(1996,2009)
print(years23)

last_year=years23[-1] 

years_jan_1_str_23=[str(i)+'-01-01' for i in np.arange(1996,last_year+1) ] 
yearly_start_times_23=parse_time(years_jan_1_str_23).datetime
yearly_start_times_num_23=parse_time(years_jan_1_str_23).plot_date

#same for July 1 as middle of the year
years_jul_1_str_23=[str(i)+'-07-01' for i in np.arange(1996,last_year+1) ] 
yearly_mid_times_23=parse_time(years_jul_1_str_23).datetime
yearly_mid_times_num_23=parse_time(years_jul_1_str_23).plot_date


print(yearly_mid_times_23)

#same for december 31
years_dec_31_str_23=[str(i)+'-12-31' for i in np.arange(1996,last_year+1) ] 
yearly_end_times_23=parse_time(years_dec_31_str_23).datetime
yearly_end_times_num_23=parse_time(years_dec_31_str_23).plot_date


# minimim maximum times as given by
#http://www.sidc.be/silso/cyclesmm
#1996   08   11.2   2001 11  180.3   12  04

solarmin23=parse_time('1996-08-01').datetime
# minstart_23=solarmin_23-datetime.timedelta(days=366*1.5)
# minend=solarmin+datetime.timedelta(days=365)
# minstart_num=parse_time(minstart).plot_date
# minend_num=parse_time(minend).plot_date

solarmax23=parse_time('2001-11-01').datetime
# maxstart=solarmax-datetime.timedelta(days=365*3)
# maxend=solarmax+datetime.timedelta(days=365/2)
# maxstart_num=parse_time(maxstart).plot_date
# maxend_num=parse_time(maxend).plot_date

print('min/max',solarmin23,solarmax23)
print()


#################### spots 
#get yearly smoothed 12 month spot rate 

spots23=np.zeros(len(years23))
counter=0
for q in years23:
    spots23[counter]=np.mean(ssn.spot[np.where(ssn.year==q)[0] ] )
    counter=counter+1

print('spots:')
print('spots yearly mean', np.rint(spots23))
print()


#################### ICME rate
#number of MFR events in Wind ICME catalog, 
#for years as in yearly_mid_times_23 but start aug 1996 and end with 
#nov 2008 #note : halloween events at ACE! Wind not?
wind_mfr_number_23=[2,8,29,11,30,21,20,6,12,16,13,7,3]
#wind_mfr_number_23_err=[2,8,29,11,30,21,20,6,12,16,13,7,3]


rc_rate23=rc_rate_values[0:13]


print('icme rate:')
print('icmes RC',rc_rate23)
print()


# ## solar cycle 24

# In[11]:


print('cycle 24\n')

#################### times
print('times:')
#these years cover solar cycle 24
years24=np.arange(2009,2020)
print(years24)

#same for July 1 as middle of the year
last_year=2020
years_jul_1_str_24=[str(i)+'-07-01' for i in np.arange(2009,last_year) ] 
yearly_mid_times_24=parse_time(years_jul_1_str_24).datetime
yearly_mid_times_num_24=parse_time(years_jul_1_str_24).plot_date
print(yearly_mid_times_24)
print()


#################### spots 
print('spots:')

#get yearly smoothed 12 month spot rate 
spots24=np.zeros(len(years24))
counter=0
for q in years24:
    spots24[counter]=np.mean(ssn.spot[np.where(ssn.year==q)[0] ] )
    counter=counter+1

print('spots yearly mean:', np.rint(spots24))
print()

print('----------')

################# ICME rates
print('ICME rate:')
print()

#2008 to 2019 are indices 12:24
rc_rate24=rc_rate_values[13:24]
print(years24)
print('icmes RC',rc_rate24)
print()

#here also from 2008 to 2019
ic_rate24=icrate[2:13].mean1.to_numpy()
ic_rate24_std=icrate[2:13].std1.to_numpy()
print('icmes ICMECAT mean',ic_rate24)
print('icmes ICMECAT std',ic_rate24_std)
icrate_years_24=parse_time(mdates.num2date(icrate_year2[2:13])).iso
print(icrate_years_24)

print()
print('ratio RC to ICMECAT:')
print(np.round(np.mean(rc_rate24/ic_rate24),2))


# ## solar cycle 25

# In[12]:


print('cycle 25\n')

#################### times from 2020 onwards
print('times:')
#these years cover solar cycle 24
years25=np.arange(2020,2022)
print(years25)

#same for July 1 as middle of the year
last_year=2022
years_jul_1_str_25=[str(i)+'-07-01' for i in np.arange(2020,last_year) ] 
yearly_mid_times_25=parse_time(years_jul_1_str_25).datetime
yearly_mid_times_num_25=parse_time(years_jul_1_str_25).plot_date
print(yearly_mid_times_25)
print()


#################### spots 
print('spots:')

#get yearly smoothed 12 month spot rate 
spots25=np.zeros(len(years25))
counter=0
for q in years25:
    spots25[counter]=np.mean(ssn.spot[np.where(ssn.year==q)[0] ] )
    counter=counter+1

print('spots yearly mean:', np.rint(spots25))
print()

print('----------')

################# ICME rates
print('ICME rate:')
print()

#2020 is index 25
rc_rate25=rc_rate_values[24:26]
print(years25)
print('icmes RC',rc_rate25)
print()

ic_rate25=icrate[13:len(icrate)].mean1.to_numpy()
ic_rate25_std=icrate[13:len(icrate)].std1.to_numpy()
print('icmes ICMECAT mean',ic_rate25)
print('icmes ICMECAT std',ic_rate25_std)
icrate_years_25=parse_time(mdates.num2date(icrate_year2[13:len(icrate)])).iso
print(icrate_years_25)

print()
#print('ratio RC to ICMECAT:')
#print(np.round(np.mean(rc_rate24/ic_rate24),2))


# ## **Figure 2** correlation SSN with ICME rate and fit
# plot SSN vs ICME rate, linear fit with confidence interval

# In[13]:


#add spots23/24 and rc_rate23/24 into 1 array for correlation
spots_corr=np.hstack([spots23,spots24])
rc_rate_corr=np.hstack([rc_rate23,rc_rate24])

#quick check with seaborn for correlation
#seaborn uses this :import statsmodels
#kind{ “scatter” | “reg” | “resid” | “kde” | “hex” }, optional
#sns.jointplot(spots_corr,rc_rate_corr,kind='reg',xlim=[0,np.max(spots_corr)+20],ylim=[0,np.max(rc_rate_corr+10)], \
#                  marginal_kws=dict(bins=5, rug=True),x_ci=95).set_axis_labels("yearly mean SSN (12 month running mean)", "ICME rate (Richardson & Cane)")


############################## make linear fit
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
print('linear fit SSN vs. ICME rate')
linfit=scipy.stats.linregress(spots_corr,rc_rate_corr)
print(linfit)
print('icme_rate[per year] = (',np.round(linfit.slope,3),'+/-',np.round(linfit.stderr,3),') * sunspot number [yearly average] + ',np.round(linfit.intercept,3))
print('inverse slope:',np.round(1/linfit.slope,2))
print('Pearson correlation coefficient:',np.round(linfit.rvalue,2))

######################### Function for conversion SSN to ICMERATE, with errors on fit confidence 1 time std

#with these results from the linear fit, make a conversion function from ssn to icme_rate
#old
#def ssn_to_rate(ssn,fitresult):    
#    rate=linfit.slope*ssn+linfit.intercept    
#    rate_low=(linfit.slope-1*linfit.stderr)*ssn+linfit.intercept    
#    rate_up=(linfit.slope+1*linfit.stderr)*ssn+linfit.intercept    
#    return rate, rate_low, rate_up



#with these results from the linear fit, make a conversion function from ssn to icme_rate
def ssn_to_rate(ssn,fitresult):    
    rate=linfit.slope*ssn+linfit.intercept    
    rate_low=(linfit.slope-1*linfit.stderr)*ssn+linfit.intercept    
    rate_up=(linfit.slope+1*linfit.stderr)*ssn+linfit.intercept    
    return rate, rate_low, rate_up



print('mean difference icme rate SC23 to predicted rate with linear fit, over full cycle:',        np.round(np.mean(ssn_to_rate(spots23,linfit)[0]-rc_rate23),2),'+/-',       np.round(np.std(ssn_to_rate(spots23,linfit)[0]-rc_rate23),1))#
print('mean difference icme rate SC24 to predicted rate with linear fit, over full cycle:',        np.round(np.mean(ssn_to_rate(spots24,linfit)[0]-rc_rate24),2),'+/-',       np.round(np.std(ssn_to_rate(spots24,linfit)[0]-rc_rate24),1))#
print()



mean_stddev=np.mean([np.std(ssn_to_rate(spots23,linfit)[0]-rc_rate23),np.std(ssn_to_rate(spots24,linfit)[0]-rc_rate24)])
print('mean stddev for both cycles gives a 1 sigma spread in the linear fit',linfit.intercept,'+/-',    mean_stddev)
print()





#with these results from the linear fit, make a conversion function from ssn to icme_rate
def ssn_to_rate(ssn,fitresult):    
    rate=linfit.slope*ssn+linfit.intercept    
    rate_low=(linfit.slope-1*linfit.stderr)*ssn+linfit.intercept    - mean_stddev
    rate_up=(linfit.slope+1*linfit.stderr)*ssn+linfit.intercept    + mean_stddev
    return rate, rate_low, rate_up



############################ plot figure
sns.set_context("talk")     
sns.set_style('darkgrid')
fsize=15
fig=plt.figure(2,figsize=(10,7),dpi=80)

plt.plot(spots23,rc_rate23,color='black',marker='o',linestyle='',markersize=10,label='solar cycle 23')
plt.plot(spots24,rc_rate24,color='black',markerfacecolor='white',marker='o',linestyle='',markersize=10,label='solar cycle 24')
plt.plot(spots25,rc_rate25,color='black',markerfacecolor='white',marker='s',linestyle='',markersize=10,label='solar cycle 25 (not fitted)')



plt.xlim(0,320)
plt.ylim(0,np.max(rc_rate_corr)+30)
plt.xlabel("yearly mean SSN (from 13 month running mean)")
plt.ylabel("ICME rate per year (Richardson & Cane)")

#errors
#ylinfit_1=(linfit.slope-linfit.stderr)*xlinfit+linfit.intercept
#ylinfit_2=(linfit.slope+linfit.stderr)*xlinfit+linfit.intercept
#plt.plot(xlinfit,ylinfit_1,'--k')
#plt.plot(xlinfit,ylinfit_2,'--k')

#https://seaborn.pydata.org/generated/seaborn.regplot.html
#sns.regplot(spots_corr,rc_rate_corr, x_ci='ci',ci=95,label=r'fit confidence 2$\mathrm{\sigma}$',truncate=False)
#sns.regplot(spots_corr,rc_rate_corr, x_ci='ci',ci=68,label=r'fit confidence 1$\mathrm{\sigma}$',truncate=False)




xlinfit=np.arange(0,350)
#1 sigma interval by using mean difference as +/- to linear fit

ylinfit_1=(linfit.slope+linfit.stderr)*xlinfit+linfit.intercept+mean_stddev
ylinfit_2=(linfit.slope-linfit.stderr)*xlinfit+linfit.intercept-mean_stddev


plt.fill_between(xlinfit,ylinfit_1,ylinfit_2,alpha=0.2,color='coral',label='fit confidence 1$\mathrm{\sigma}$')
#plt.plot(xlinfit,ylinfit_1,'--k')
#plt.plot(xlinfit,ylinfit_2,'--k')



ylinfit=linfit.slope*xlinfit+linfit.intercept
plt.plot(xlinfit,ylinfit,'-k',label='linear fit')
plt.plot(xlinfit,np.zeros(len(xlinfit))+52,'--k',alpha=0.5)
plt.plot(xlinfit,np.zeros(len(xlinfit))+26,'--k',alpha=0.5)

plt.legend(loc=2)
plt.tight_layout()
#plt.savefig(outputdirectory+'/fig2_rate_ssn.pdf', dpi=300)
plt.savefig(outputdirectory+'/fig2_rate_ssn.png', dpi=300)


# ## predictions for solar cycle 25: SSN and ICME rate
# ### 1. Mean cycle model

# In[14]:


# from heliocats import stats as hs
# importlib.reload(hs) #reload again while debugging

print('---------------------------------')
print('cycle 25')
print()
print('calculate several Hathaway function models for SSN and predict the ICME rate')
print('')
print('---------------------------------')


############# set yearly times
last_year=2033
years_jan_1_str_25=[str(i)+'-01-01' for i in np.arange(2020,last_year) ] 
yearly_start_times_25=parse_time(years_jan_1_str_25,format='iso').datetime
yearly_start_times_num_25=parse_time(years_jan_1_str_25,format='iso').plot_date

#same for July 1 as middle of the year
years_jul_1_str_25=[str(i)+'-07-01' for i in np.arange(2020,last_year) ] 
yearly_mid_times_25=parse_time(years_jul_1_str_25,format='iso').datetime
yearly_mid_times_num_25=parse_time(years_jul_1_str_25,format='iso').plot_date

#### for smooth plotting daily list of times
#t0 of solar cycler 25 is unclear at time of writing! assumed 2020 june 1
start_25=parse_time('2020-06-01',format='iso').datetime
#not end time of next solar cycle, but end time of plotting
end_25=parse_time('2033-01-01',format='iso').datetime
    
#create an array with 1 day resolution between t start and end
times_25_daily = [ start_25 + datetime.timedelta(days=n) for n in range(int ((end_25 - start_25).days))]  
times_25_daily_mat=mdates.date2num(times_25_daily) 

######################################################## 1. Mean cycle model

#hathaway 2015: t0 minus 4 months is good match
#t0 is start date of cycle
shift_t0=timedelta(days=4*30+1)

print()
print('1. mean cycle')

#mean of all cycles function: Hathaway 1994
#parameters taken from Hathaway 2015 review section 4.5 DOI 10.1007/lrsp-2015-4
#https://link.springer.com/article/10.1007/lrsp-2015-4
#but this is based on different SSN definitions?
# am=195
# amerr=50
# bm=56
# cm=0.8

print('load min max for all cycles from SILSO')

#http://sidc.oma.be/silso/DATA/Cycles/TableCyclesMiMa.txt
#Table of minima, maxima and cycle durations based on
#13-month smoothed monthly mean sunspot numbers (Version 2.0).
mima_num=np.loadtxt('data/TableCyclesMiMa.txt', skiprows=2)
mima = pd.DataFrame(mima_num)
mima.columns = ['cycle','min_year','min_month','min_sn','max_year','max_month','max_sn','dur_year','dur_months']

print()
print('Average maximum sunspot number (SILSO) of all cycles, 13 months smoothed:',np.round(mima.max_sn.mean(),1),' +/- ',np.round(mima.max_sn.std(),1))
print()

#use all parameters as in Hathaway but adjust amplitude to average for the SILSO numbers:
am=342
bm=56
cm=0.8
print('average cycle a,b,c:', am,bm,cm)



#mean of all cycles yearly ssn numbers 
spots_predict_25m=hs.hathaway(yearly_mid_times_25, start_25-shift_t0,am,bm,cm)
#and same for 1 day resolution 
spots_predict_25m_daily=hs.hathaway(times_25_daily, start_25-shift_t0,am,bm,cm)

print('t0 in SC25 PP prediction is:',start_25-shift_t0)
#time of maximum and ssn
print('max ssn',np.rint(np.max(spots_predict_25m_daily)))
print('at time',str(times_25_daily[np.argmax(spots_predict_25m_daily)])[0:11])
print()
print(yearly_mid_times_25)
print('spots yearly: ',np.rint(spots_predict_25m))

#yearly icme numbers: convert with function
icmes_predict_25m=ssn_to_rate(spots_predict_25m,linfit)[0]
#same daily resolution
icmes_predict_25m_daily=ssn_to_rate(spots_predict_25m_daily,linfit)[0]



print('icmes yearly: ',np.rint(icmes_predict_25m))

print()
print()
print('Merge error from 1. ssn prediction from fit with 2. spread in icme rate observed with ICMECAT (not used in plots later as PP19 is quite similar)')

#1. error from SSN fit yearly, low and high
icmes_predict_25m_low=ssn_to_rate(spots_predict_25m,linfit)[1]
icmes_predict_25m_high=ssn_to_rate(spots_predict_25m,linfit)[2]

#this is the range in the icme rate arising from the fit, symmetric for high and low values 
ic_rate25_m_std_fit=np.round(icmes_predict_25m-icmes_predict_25m_low,1)
print('error from SSN to ICME fit',ic_rate25_m_std_fit)

#2. error from ICME rate
#assumption icrate_std=2 for last 3 years
ic_rate25_std=np.hstack([ic_rate24_std[0:-1],np.array([2.0,2.0,2.0])])
print('spread in ICME rate from SC24, assuming 2009=> 2020',ic_rate25_std)

#add both errors as sigma_new=sqrt(sigma1^2+sigma2^2)
ic_rate_25_m_std=np.round(np.sqrt(ic_rate25_m_std_fit**2+ic_rate25_std**2),1)
print('Std in ICME rate from fit and ICMECAT range for each year:')
print(ic_rate_25_m_std)


# In[15]:


########################################################### 2. SC25 panel prediction (SC25PP)
print('---------------------------- 2. SC25 PP 2019')
print()
print('get PP25 prediction from JSON file from NOAA (2020 May 27) ')    

#download PP25 prediction from NOAA with timestamp
#sc25pp_url='https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json'

#try: urllib.request.urlretrieve(sc25pp_url,data_path+'predicted-solar-cycle_2020_may_27.json')
#except urllib.error.URLError as e:
#    print('Failed downloading ', sc25pp_url,' ',e)

pp25_df=pd.read_json('data/predicted-solar-cycle_2020_may_27.json')


#kill first few rows and start with May 2020
pp25_df=pp25_df.drop([0,1,2,3,4,5],axis=0)
pp25_df['time-tag']
pp25_df_times=parse_time(pp25_df['time-tag']).datetime
pp25_df_times_num=parse_time(pp25_df['time-tag']).plot_date
pp25_df_ssn=pp25_df['predicted_ssn']
pp25_df_ssn_high=pp25_df['high_ssn']
pp25_df_ssn_low=pp25_df['low_ssn']

#make hathaway fit
hw_param_pp25 = scipy.optimize.curve_fit(hs.hathaway_fit, pp25_df_times_num-pp25_df_times_num[0],pp25_df_ssn,                                        p0=[-300,200,60,1])    
print('Hathaway function fit parameters x0,a,b,c:',np.round(hw_param_pp25[0][0:3],1),np.round(hw_param_pp25[0][3],2))
#get t0 date for hathway function from fit x0
x0fit_in_days=int(np.rint(hw_param_pp25[0][0]))
start_25_fit=mdates.num2date(pp25_df_times_num[0]+ mdates.date2num(np.datetime64('0000-12-31')))+timedelta(days=x0fit_in_days)
print('t0 in SC25 PP prediction is:',start_25_fit)

# same for low and high
hw_param_pp25_low = scipy.optimize.curve_fit(hs.hathaway_fit, pp25_df_times_num-pp25_df_times_num[0],pp25_df_ssn_low,                                        p0=[-300,200,60,1])    
x0fit_in_days_low=int(np.rint(hw_param_pp25_low[0][0]))
start_25_fit_low=mdates.num2date(pp25_df_times_num[0]+mdates.date2num(np.datetime64('0000-12-31')))+timedelta(days=x0fit_in_days_low)
print('lower error: Hathaway function fit parameters x0,a,b,c:',np.round(hw_param_pp25_low[0][0:3],1),np.round(hw_param_pp25_low[0][3],2))

hw_param_pp25_high = scipy.optimize.curve_fit(hs.hathaway_fit, pp25_df_times_num-pp25_df_times_num[0],pp25_df_ssn_high,                                        p0=[-300,200,60,1])    
x0fit_in_days_high=int(np.rint(hw_param_pp25_high[0][0]))
start_25_fit_high=mdates.num2date(pp25_df_times_num[0]+mdates.date2num(np.datetime64('0000-12-31')))+timedelta(days=x0fit_in_days_high)
print('higher error: Hathaway function fit parameters x0,a,b,c:',np.round(hw_param_pp25_high[0][0:3],1),np.round(hw_param_pp25_high[0][3],2))
print('note that high and low error ranges are calculated here with a Hathaway function, small error to PP forecast')


spots_predict_25pp=hs.hathaway(yearly_mid_times_25,start_25_fit,hw_param_pp25[0][1],hw_param_pp25[0][2],hw_param_pp25[0][3])
spots_predict_25pp_daily=hs.hathaway(times_25_daily,start_25_fit,hw_param_pp25[0][1],hw_param_pp25[0][2],hw_param_pp25[0][3])

#fit parameters for PP19
#start_25_fit
app=hw_param_pp25[0][1]
bpp=hw_param_pp25[0][2]
cpp=hw_param_pp25[0][3]

spots_predict_25pp_low=hs.hathaway(yearly_mid_times_25,start_25_fit_low,hw_param_pp25_low[0][1],hw_param_pp25_low[0][2],hw_param_pp25_low[0][3])
spots_predict_25pp_daily_low=hs.hathaway(times_25_daily,start_25_fit_low,hw_param_pp25_low[0][1],hw_param_pp25_low[0][2],hw_param_pp25_low[0][3])
spots_predict_25pp_high=hs.hathaway(yearly_mid_times_25,start_25_fit_high,hw_param_pp25_high[0][1],hw_param_pp25_high[0][2],hw_param_pp25_high[0][3])
spots_predict_25pp_daily_high=hs.hathaway(times_25_daily,start_25_fit_high,hw_param_pp25_high[0][1],hw_param_pp25_high[0][2],hw_param_pp25_high[0][3])

plt.figure(12,figsize=(10,6),dpi=60)
plt.plot(pp25_df_times,pp25_df_ssn,'-k',markerfacecolor='white')
plt.plot(pp25_df_times,pp25_df['high_ssn'],'-k',markerfacecolor='white')
plt.plot(pp25_df_times,pp25_df['low_ssn'],'-k',markerfacecolor='white')
plt.plot_date(yearly_mid_times_25,spots_predict_25pp,'ok')
plt.plot_date(times_25_daily,spots_predict_25pp_daily,'-y')
plt.plot_date(yearly_mid_times_25,spots_predict_25pp_low,'ok')
plt.plot_date(yearly_mid_times_25,spots_predict_25pp_high,'ok')
plt.plot_date(times_25_daily,spots_predict_25pp_daily_low,'-g')
plt.plot_date(times_25_daily,spots_predict_25pp_daily_high,'-b')
plt.plot_date(times_25_daily,spots_predict_25m_daily,'-r')

plt.ylabel('SSN, PP 25 prediction')
#time of maximum and ssn
print('max ssn',np.rint(np.max(spots_predict_25pp_daily)))
print('at time',str(times_25_daily[np.argmax(spots_predict_25pp_daily)])[0:11])
print()
print('spots yearly: ',np.rint(spots_predict_25pp))


#yearly spots numbers
icmes_predict_25pp=ssn_to_rate(spots_predict_25pp,linfit)[0]
icmes_predict_25pp_low=ssn_to_rate(spots_predict_25pp_low,linfit)[0]
icmes_predict_25pp_high=ssn_to_rate(spots_predict_25pp_high,linfit)[0]

#same daily
icmes_predict_25pp_daily=ssn_to_rate(spots_predict_25pp_daily,linfit)[0]
icmes_predict_25pp_daily_low=ssn_to_rate(spots_predict_25pp_daily_low,linfit)[0]
icmes_predict_25pp_dailyhigh=ssn_to_rate(spots_predict_25pp_daily_high,linfit)[0]


print()
print('Merge error from 1. ssn prediction 2. from SSN to ICME from fit and 3. spread in icme rate observed with ICMECAT')

#1. error from SSN prediction
ic_rate25_std_pp_ssnpred=np.round(((icmes_predict_25pp_high-icmes_predict_25pp)+abs(icmes_predict_25pp_low-icmes_predict_25pp))/2,2)
print('ICME rate error from SSN prediction',ic_rate25_std_pp_ssnpred)


#2. error from fit SSN to ICME rate
icmes_predict_25_pp=ssn_to_rate(spots_predict_25pp,linfit)[0]
icmes_predict_25_pp_low=ssn_to_rate(spots_predict_25pp,linfit)[1]
icmes_predict_25_pp_high=ssn_to_rate(spots_predict_25pp,linfit)[2]

#this is the range in the icme rate arising from the fit, symmetric for high and low values 
ic_rate25_std_pp_ssnfit=np.round(icmes_predict_25_pp-icmes_predict_25_pp_low,1)
print('error from SSN to ICME fit',ic_rate25_std_pp_ssnfit)

#3. error from ICME rate spread
print('spread in ICME rate', ic_rate25_std)

print()
#add all 3 errors as sigma_new=sqrt(sigma1^2+sigma2^2)
ic_rate_25_pp_std=np.round(np.sqrt(ic_rate25_std_pp_ssnpred**2+ic_rate25_std_pp_ssnfit**2+ic_rate25_std**2),1)
print('final Std in ICME rate from SSN prediction, SSN to ICME fit and ICMECAT range for each year:')
print(ic_rate_25_pp_std)


# In[16]:


################################### SC25MC
#SC25MC prediction or MC20 (see paper)
#https://arxiv.org/abs/2006.15263


print('-------------------------- 3. SC25 MC 2020')
a=444
aerr68=48 #MC20: 204-254 68, 153 305 95
aerr95=147 #153 305 95, +/- 76 95
b=60
c=0.8
print('a,b,c:', a,b,c)
print('range for a:',a-aerr68,a+aerr68)

print('start of sc25 here in MC20',start_25-shift_t0)

#yearly_numbers
spots_predict_25=hs.hathaway(yearly_mid_times_25, start_25-shift_t0,a,b,c)
#error ranges
spots_predict_25_lower68=hs.hathaway(yearly_mid_times_25, start_25-shift_t0,a-aerr68,b,c)
spots_predict_25_upper68=hs.hathaway(yearly_mid_times_25, start_25-shift_t0,a+aerr68,b,c)
#spots_predict_25_lower95=hs.hathaway(yearly_mid_times_25, start_25-shift_t0,a-aerr95,b,c)
#spots_predict_25_upper95=hs.hathaway(yearly_mid_times_25, start_25-shift_t0,a+aerr95,b,c)

#daily numbers
spots_predict_25_daily=hs.hathaway(times_25_daily, start_25-shift_t0,a,b,c)
#error ranges
spots_predict_25_daily_lower68=hs.hathaway(times_25_daily, start_25-shift_t0,a-aerr68,b,c)
spots_predict_25_daily_upper68=hs.hathaway(times_25_daily, start_25-shift_t0,a+aerr68,b,c)
#spots_predict_25_daily_lower95=hs.hathaway(times_25_daily, start_25-shift_t0,a-aerr95,b,c)
#spots_predict_25_daily_upper95=hs.hathaway(times_25_daily, start_25-shift_t0,a+aerr95,b,c)

#time of maximum and ssn
print('max ssn',np.rint(np.max(spots_predict_25_daily)))
print('at time',str(times_25_daily[np.argmax(spots_predict_25_daily)])[0:11])

print()
print('spots yearly',np.rint(spots_predict_25))

#
#yearly spots numbers
icmes_predict_25mc=ssn_to_rate(spots_predict_25,linfit)[0]
icmes_predict_25mc_lower68=ssn_to_rate(spots_predict_25_lower68,linfit)[0]
icmes_predict_25mc_upper68=ssn_to_rate(spots_predict_25_upper68,linfit)[0]
#icmes_predict_25mc_lower95=ssn_to_rate(spots_predict_25_lower95,linfit)[0]
#icmes_predict_25mc_upper95=ssn_to_rate(spots_predict_25_upper95,linfit)[0]

#same daily
icmes_predict_25mc_daily=ssn_to_rate(spots_predict_25_daily,linfit)[0]
icmes_predict_25mc_daily_lower68=ssn_to_rate(spots_predict_25_daily_lower68,linfit)[0]
icmes_predict_25mc_daily_upper68=ssn_to_rate(spots_predict_25_daily_upper68,linfit)[0]
#icmes_predict_25mc_daily_lower95=ssn_to_rate(spots_predict_25_daily_lower95,linfit)[0]
#icmes_predict_25mc_daily_upper95=ssn_to_rate(spots_predict_25_daily_upper95,linfit)[0]


#1. error from SSN prediction
print()
print('icmes yearly: ',np.rint(icmes_predict_25mc))
print('icmes lower68: ',np.rint(icmes_predict_25mc_lower68))
print('icmes upper68: ',np.rint(icmes_predict_25mc_upper68))

#print('icmes lower95: ',np.rint(icmes_predict_25mc_lower95))
#print('icmes upper95: ',np.rint(icmes_predict_25mc_upper95))


print()
print('Merge error from 1. ssn prediction 2. from SSN to ICME from fit and 3. spread in icme rate observed with ICMECAT')

#1. error from SSN prediction
ic_rate25_std_mc20_ssnpred=np.round(((icmes_predict_25mc_upper68-icmes_predict_25mc)+abs(icmes_predict_25mc_lower68-icmes_predict_25mc))/2,2)
print('ICME rate error from SSN prediction',ic_rate25_std_mc20_ssnpred)


#2. error from fit SSN to ICME rate
icmes_predict_25_mc20=ssn_to_rate(spots_predict_25,linfit)[0]
icmes_predict_25_mc20_low=ssn_to_rate(spots_predict_25,linfit)[1]
icmes_predict_25_mc20_high=ssn_to_rate(spots_predict_25,linfit)[2]

#this is the range in the icme rate arising from the fit, symmetric for high and low values 
ic_rate25_std_mc20_ssnfit=np.round(icmes_predict_25_mc20-icmes_predict_25_mc20_low,1)
print('error from SSN to ICME fit',ic_rate25_std_mc20_ssnfit)

#3. error from ICME rate spread
print('spread in ICME rate', ic_rate25_std)

print()
#add all 3 errors as sigma_new=sqrt(sigma1^2+sigma2^2)
ic_rate_25_mc20_std=np.round(np.sqrt(ic_rate25_std_mc20_ssnpred**2+ic_rate25_std_mc20_ssnfit**2+ic_rate25_std**2),1)
print('final Std in ICME rate from SSN prediction, SSN to ICME fit and ICMECAT range for each year:')
print(ic_rate_25_mc20_std)


# In[17]:


################################### SC25MC
#SC25MC prediction or MC20 (see paper)
#https://arxiv.org/abs/2006.15263

#90 ± 20.
#https://spaceweatherarchive.com/2022/02/25/the-termination-event-has-arrived/


print('-------------------------- 4. SC25 MC 2022 based on terminator ')
a=363
aerr68=38 #MC20: 325 -> 170 max, 
b=60#b=60
c=0.8 #c=0.8
print('a,b,c:', a,b,c)
print('range for a:',a-aerr68,a+aerr68)

#shift_t02=timedelta(days=4*30+1+12*30)

shift_t02=timedelta(days=121)


print('start of sc25 here in MC20 v2',start_25-shift_t02)

#yearly_numbers
spots_predict_25_mc2=hs.hathaway(yearly_mid_times_25, start_25-shift_t02,a,b,c)
#error ranges
spots_predict_25_lower68_mc2=hs.hathaway(yearly_mid_times_25, start_25-shift_t02,a-aerr68,b,c)
spots_predict_25_upper68_mc2=hs.hathaway(yearly_mid_times_25, start_25-shift_t02,a+aerr68,b,c)
#spots_predict_25_lower95=hs.hathaway(yearly_mid_times_25, start_25-shift_t0,a-aerr95,b,c)
#spots_predict_25_upper95=hs.hathaway(yearly_mid_times_25, start_25-shift_t0,a+aerr95,b,c)


#daily numbers
spots_predict_25_daily_mc2=hs.hathaway(times_25_daily, start_25-shift_t02,a,b,c)
#error ranges
spots_predict_25_daily_lower68_mc2=hs.hathaway(times_25_daily, start_25-shift_t02,a-aerr68,b,c)
spots_predict_25_daily_upper68_mc2=hs.hathaway(times_25_daily, start_25-shift_t02,a+aerr68,b,c)
#spots_predict_25_daily_lower95=hs.hathaway(times_25_daily, start_25-shift_t0,a-aerr95,b,c)
#spots_predict_25_daily_upper95=hs.hathaway(times_25_daily, start_25-shift_t0,a+aerr95,b,c)


plt.figure(12,figsize=(10,6),dpi=60)
plt.plot_date(times_25_daily,spots_predict_25_daily_lower68_mc2,'-g')
plt.plot_date(times_25_daily,spots_predict_25_daily_upper68_mc2,'-b')
plt.plot_date(times_25_daily,spots_predict_25_daily_mc2,'-r')


#time of maximum and ssn
print('max ssn',np.rint(np.max(spots_predict_25_daily_mc2)))
print('at time',str(times_25_daily[np.argmax(spots_predict_25_daily_mc2)])[0:11])

print()
print('spots yearly',np.rint(spots_predict_25_mc2))

#
#yearly spots numbers
icmes_predict_25mc2=ssn_to_rate(spots_predict_25_mc2,linfit)[0]
icmes_predict_25mc2_lower68=ssn_to_rate(spots_predict_25_lower68_mc2,linfit)[0]
icmes_predict_25mc2_upper68=ssn_to_rate(spots_predict_25_upper68_mc2,linfit)[0]
#icmes_predict_25mc_lower95=ssn_to_rate(spots_predict_25_lower95,linfit)[0]
#icmes_predict_25mc_upper95=ssn_to_rate(spots_predict_25_upper95,linfit)[0]

#same daily
icmes_predict_25mc2_daily=ssn_to_rate(spots_predict_25_daily_mc2,linfit)[0]
icmes_predict_25mc2_daily_lower68=ssn_to_rate(spots_predict_25_daily_lower68_mc2,linfit)[0]
icmes_predict_25mc2_daily_upper68=ssn_to_rate(spots_predict_25_daily_upper68_mc2,linfit)[0]
#icmes_predict_25mc_daily_lower95=ssn_to_rate(spots_predict_25_daily_lower95,linfit)[0]
#icmes_predict_25mc_daily_upper95=ssn_to_rate(spots_predict_25_daily_upper95,linfit)[0]


#1. error from SSN prediction
print()
print('icmes yearly: ',np.rint(icmes_predict_25mc2))
print('icmes lower68: ',np.rint(icmes_predict_25mc2_lower68))
print('icmes upper68: ',np.rint(icmes_predict_25mc2_upper68))

#print('icmes lower95: ',np.rint(icmes_predict_25mc_lower95))
#print('icmes upper95: ',np.rint(icmes_predict_25mc_upper95))


print()
print('Merge error from 1. ssn prediction 2. from SSN to ICME from fit and 3. spread in icme rate observed with ICMECAT')

#1. error from SSN prediction
ic_rate25_std_mc20_2_ssnpred=np.round(((icmes_predict_25mc2_upper68-icmes_predict_25mc2)+abs(icmes_predict_25mc2_lower68-icmes_predict_25mc2))/2,2)
print('ICME rate error from SSN prediction',ic_rate25_std_mc20_2_ssnpred)


#2. error from fit SSN to ICME rate
icmes_predict_25_mc20_2=ssn_to_rate(spots_predict_25_mc2,linfit)[0]
icmes_predict_25_mc20_2_low=ssn_to_rate(spots_predict_25_mc2,linfit)[1]
icmes_predict_25_mc20_2_high=ssn_to_rate(spots_predict_25_mc2,linfit)[2]

#this is the range in the icme rate arising from the fit, symmetric for high and low values 
ic_rate25_std_mc20_2_ssnfit=np.round(icmes_predict_25_mc20_2-icmes_predict_25_mc20_2_low,1)
print('error from SSN to ICME fit',ic_rate25_std_mc20_2_ssnfit)

#3. error from ICME rate spread
print('spread in ICME rate', ic_rate25_std)

print()
#add all 3 errors as sigma_new=sqrt(sigma1^2+sigma2^2)
ic_rate_25_mc20_2_std=np.round(np.sqrt(ic_rate25_std_mc20_2_ssnpred**2+ic_rate25_std_mc20_2_ssnfit**2+ic_rate25_std**2),1)
print('final Std in ICME rate from SSN prediction, SSN to ICME fit and ICMECAT range for each year:')
print(ic_rate_25_mc20_std)


# ## **Figure 3** ICME rate predictions

# In[18]:


sns.set_context("talk")     
sns.set_style('whitegrid')
sns.set_style("ticks",{'grid.linestyle': '--'})

fig=plt.figure(300,figsize=(13,10),dpi=70)

fsize=15
max_spot=350
max_icme=80

####################### MC20 model
ax1 = plt.subplot(212) 

################## ICR 
ax1.plot(yearly_mid_times_23,rc_rate23, color='black',marker='o',markerfacecolor='black',label='RC ICME rate SC23',linestyle='')
ax1.plot(yearly_mid_times_24,rc_rate24, color='black', marker='o',markerfacecolor='white',label='RC ICME rate SC24',linestyle='')
ax1.plot(yearly_mid_times_25[0:2],rc_rate25, color='blue', marker='o',markerfacecolor='white',label='RC ICME rate SC25',linestyle='')




ax1.plot(yearly_mid_times_24,ic_rate24, color='black', marker='s',markerfacecolor='white',label='ICMECAT ICME rate SC24',linestyle='')
ax1.plot(yearly_mid_times_25[0:2],ic_rate25, color='blue', marker='s',markerfacecolor='white',label='ICMECAT ICME rate SC25',linestyle='')

ax1.plot([icrate_year2[2:],icrate_year2[2:]],[icrate.mean1[2:]-icrate.std1[2:],icrate.mean1[2:]+icrate.std1[2:]],'-k',lw=1.1,label='')

#ax1.plot(yearly_mid_times_25,icmes_predict_25mc,color='black', marker='o',markerfacecolor='grey',label='Predicted ICME rate MC20',linestyle='')
ax1.errorbar(yearly_mid_times_25,icmes_predict_25mc,yerr=ic_rate_25_mc20_std, color='black', marker='o',markerfacecolor='grey',label='Predicted ICME rate MC20',linestyle='')

ax1.set_ylim(0,max_icme)
ax1.set_ylabel('number of ICMEs per year',fontsize=fsize)
ax1.legend(loc=2,fontsize=12)


################## SSN
ax2=ax1.twinx()
ssn_time2=ssn.time + mdates.date2num(np.datetime64('0000-12-31'))
#ax2.plot(ssn_time2,ssn.spot_mean_13,'-k',alpha=1,linewidth=1.5,label='Observed sunspot number (SIDC)')
ax2.plot(ssn_ms.time,ssn_ms.spot,'-k',alpha=1,linewidth=1.5,label='Observed sunspot number (SIDC)')



ax2.plot(times_25_daily,spots_predict_25_daily,'-r',alpha=1,linewidth=2.5,label='Predicted SSN by MC20')
ax2.plot(times_25_daily,spots_predict_25m_daily,'-g',alpha=1,linewidth=2.5,label='Predicted SSN by mean solar cycle')
#ax2.fill_between(times_25_daily, spots_predict_25_daily_lower95, spots_predict_25_daily_upper95, alpha=0.2)
ax2.fill_between(times_25_daily, spots_predict_25_daily_lower68, spots_predict_25_daily_upper68, color='red',alpha=0.2)

ax2.set_ylabel('Sunspot number ')
ax2.set_xlim(parse_time('1995-Aug-01 00:00').datetime,parse_time('2033-Jan-01 00:00').datetime)
ax2.set_ylim(0,max_spot)
years = mdates.YearLocator()   # every year
ax2.xaxis.set_minor_locator(years)
ax2.legend(loc='upper center',fontsize=12)


############################################################## PP19 panel
ax3 = plt.subplot(211) 

##########ICR axis
ax3.plot(yearly_mid_times_23,rc_rate23, color='black',marker='o',markerfacecolor='black',label='RC ICME rate SC23',linestyle='')
ax3.plot(yearly_mid_times_24,rc_rate24, color='black', marker='o',markerfacecolor='white',label='RC ICME rate SC24',linestyle='')
ax3.plot(yearly_mid_times_25[0:2],rc_rate25, color='blue', marker='o',markerfacecolor='white',label='RC ICME rate SC25',linestyle='')


ax3.plot(yearly_mid_times_24,ic_rate24, color='black', marker='s',markerfacecolor='white',label='ICMECAT ICME rate SC24',linestyle='')
ax3.plot(yearly_mid_times_25[0:2],ic_rate25, color='blue', marker='s',markerfacecolor='white',label='ICMECAT ICME rate SC25',linestyle='')


ax3.plot([icrate_year2[2:],icrate_year2[2:]],[icrate.mean1[2:]-icrate.std1[2:],icrate.mean1[2:]+icrate.std1[2:]],'-k',lw=1.1,label='')

#ax3.plot(yearly_mid_times_25,icmes_predict_25pp,color='black', marker='o',markerfacecolor='grey',label='Predicted ICME rate PP19',linestyle='')
ax3.errorbar(yearly_mid_times_25,icmes_predict_25pp,yerr=ic_rate_25_pp_std, color='black', marker='o',markerfacecolor='grey',label='Predicted ICME rate PP19',linestyle='')

ax3.set_ylim(0,max_icme)
ax3.set_ylabel('number of ICMEs per year',fontsize=fsize)
ax3.legend(loc=2,fontsize=12)

######SSN axis
ax4=ax3.twinx()
#observed SSN
#ax4.plot(ssn_time2,ssn.spot_mean_13,'-k',alpha=1,linewidth=1.5,label='Observed sunspot number (SIDC)')
ax4.plot(ssn_ms.time,ssn_ms.spot,'-k',alpha=1,linewidth=1.5,label='Observed sunspot number (SIDC)')


#PP25 prediction
ax4.plot(times_25_daily,spots_predict_25pp_daily,'-b',alpha=1,linewidth=2.5,label='Predicted SSN by PP19')
ax4.fill_between(times_25_daily,spots_predict_25pp_daily_low,spots_predict_25pp_daily_high,alpha=0.2)
#mean cycle
ax4.plot(times_25_daily,spots_predict_25m_daily,'-g',alpha=1,linewidth=2.5,label='Predicted SSN by mean solar cycle')
ax4.set_ylabel('Sunspot number ')
ax4.set_xlim(parse_time('1995-Aug-01 00:00').datetime,parse_time('2033-Jan-01 00:00').datetime)
ax4.set_ylim(0,max_spot)
years = mdates.YearLocator()   # every year
ax4.xaxis.set_minor_locator(years)
ax4.legend(loc='upper center',fontsize=12)

#grid on both panels
gridtime=[datetime.datetime(1995, 1,1), datetime.datetime(2035, 1, 1)]

for i in np.arange(0,100,10):
    ax1.plot(gridtime,[0,0]+i,linestyle='--',color='k',alpha=0.4,lw=0.8)
for i in np.arange(0,100,10):
    ax3.plot(gridtime,[0,0]+i,linestyle='--',color='k',alpha=0.4,lw=0.8)

plt.tight_layout()

#plt.annotate('(a)',[0.0,0.965],xycoords='figure fraction',weight='bold')
#plt.annotate('(b)',[0.0,0.47],xycoords='figure fraction',weight='bold')

#plt.savefig('results/plots_rate/fig3_sc25_predictions.pdf', dpi=100)
plt.savefig(outputdirectory+'/cycle25_icme_rate_predictions.png', dpi=100)


# ### with new hathaway function for cycle by McIntosh et al. shifted max 

# In[19]:


sns.set_context("talk")     
sns.set_style('whitegrid')
sns.set_style("ticks",{'grid.linestyle': '--'})

fig=plt.figure(300,figsize=(13,10),dpi=70)

fsize=15
max_spot=350
max_icme=80

####################### MC20 model
ax1 = plt.subplot(212) 

################## ICR 
ax1.plot(yearly_mid_times_23,rc_rate23, color='black',marker='o',markerfacecolor='black',label='RC ICME rate SC23',linestyle='')
ax1.plot(yearly_mid_times_24,rc_rate24, color='black', marker='o',markerfacecolor='white',label='RC ICME rate SC24',linestyle='')
ax1.plot(yearly_mid_times_24,ic_rate24, color='black', marker='s',markerfacecolor='white',label='ICMECAT ICME rate SC24',linestyle='')
ax1.plot([icrate_year2[2:],icrate_year2[2:]],[icrate.mean1[2:]-icrate.std1[2:],icrate.mean1[2:]+icrate.std1[2:]],'-k',lw=1.1,label='')

#ax1.plot(yearly_mid_times_25,icmes_predict_25mc,color='black', marker='o',markerfacecolor='grey',label='Predicted ICME rate MC20',linestyle='')
ax1.errorbar(yearly_mid_times_25,icmes_predict_25mc2,yerr=ic_rate_25_mc20_2_std, color='black', marker='o',markerfacecolor='grey',label='Predicted ICME rate MC20',linestyle='')

ax1.set_ylim(0,max_icme)
ax1.set_ylabel('number of ICMEs per year',fontsize=fsize)
ax1.legend(loc=2,fontsize=12)


################## SSN
ax2=ax1.twinx()
ax2.plot(ssn_time2,ssn.spot_mean_13,'-k',alpha=1,linewidth=1.5,label='Observed sunspot number (SIDC)')

ax2.plot(times_25_daily,spots_predict_25_daily_mc2,'-r',alpha=1,linewidth=2.5,label='Predicted SSN by MC20 shifted max')
ax2.plot(times_25_daily,spots_predict_25m_daily,'-g',alpha=1,linewidth=2.5,label='Predicted SSN by mean solar cycle')
#ax2.fill_between(times_25_daily, spots_predict_25_daily_lower95, spots_predict_25_daily_upper95, alpha=0.2)
ax2.fill_between(times_25_daily, spots_predict_25_daily_lower68_mc2, spots_predict_25_daily_upper68_mc2, color='red',alpha=0.2)

ax2.set_ylabel('Sunspot number ')
ax2.set_xlim(parse_time('1995-Aug-01 00:00').datetime,parse_time('2033-Jan-01 00:00').datetime)
ax2.set_ylim(0,max_spot)
years = mdates.YearLocator()   # every year
ax2.xaxis.set_minor_locator(years)
ax2.legend(loc='upper center',fontsize=12)


############################################################## PP19 panel
ax3 = plt.subplot(211) 

##########ICR axis
ax3.plot(yearly_mid_times_23,rc_rate23, color='black',marker='o',markerfacecolor='black',label='RC ICME rate SC23',linestyle='')
ax3.plot(yearly_mid_times_24,rc_rate24, color='black', marker='o',markerfacecolor='white',label='RC ICME rate SC24',linestyle='')
ax3.plot(yearly_mid_times_24,ic_rate24, color='black', marker='s',markerfacecolor='white',label='ICMECAT ICME rate SC24',linestyle='')
ax3.plot([icrate_year2[2:],icrate_year2[2:]],[icrate.mean1[2:]-icrate.std1[2:],icrate.mean1[2:]+icrate.std1[2:]],'-k',lw=1.1,label='')

#ax3.plot(yearly_mid_times_25,icmes_predict_25pp,color='black', marker='o',markerfacecolor='grey',label='Predicted ICME rate PP19',linestyle='')
ax3.errorbar(yearly_mid_times_25,icmes_predict_25pp,yerr=ic_rate_25_pp_std, color='black', marker='o',markerfacecolor='grey',label='Predicted ICME rate PP19',linestyle='')

ax3.set_ylim(0,max_icme)
ax3.set_ylabel('number of ICMEs per year',fontsize=fsize)
ax3.legend(loc=2,fontsize=12)

######SSN axis
ax4=ax3.twinx()
#observed SSN
ax4.plot(ssn_time2,ssn.spot_mean_13,'-k',alpha=1,linewidth=1.5,label='Observed sunspot number (SIDC)')
#PP25 prediction
ax4.plot(times_25_daily,spots_predict_25pp_daily,'-b',alpha=1,linewidth=2.5,label='Predicted SSN by PP19')
ax4.fill_between(times_25_daily,spots_predict_25pp_daily_low,spots_predict_25pp_daily_high,alpha=0.2)
#mean cycle
ax4.plot(times_25_daily,spots_predict_25m_daily,'-g',alpha=1,linewidth=2.5,label='Predicted SSN by mean solar cycle')
ax4.set_ylabel('Sunspot number ')
ax4.set_xlim(parse_time('1995-Aug-01 00:00').datetime,parse_time('2033-Jan-01 00:00').datetime)
ax4.set_ylim(0,max_spot)
years = mdates.YearLocator()   # every year
ax4.xaxis.set_minor_locator(years)
ax4.legend(loc='upper center',fontsize=12)

#grid on both panels
gridtime=[datetime.datetime(1995, 1,1), datetime.datetime(2035, 1, 1)]

for i in np.arange(0,100,10):
    ax1.plot(gridtime,[0,0]+i,linestyle='--',color='k',alpha=0.4,lw=0.8)
for i in np.arange(0,100,10):
    ax3.plot(gridtime,[0,0]+i,linestyle='--',color='k',alpha=0.4,lw=0.8)

plt.tight_layout()

#plt.annotate('(a)',[0.0,0.965],xycoords='figure fraction',weight='bold')
#plt.annotate('(b)',[0.0,0.47],xycoords='figure fraction',weight='bold')

#plt.savefig('results/plots_rate/fig3_sc25_predictions.pdf', dpi=100)
plt.savefig(outputdirectory+'/cycle25_icme_rate_predictions_shiftmax.png', dpi=100)


# In[20]:



#Extra plot for solar cycle comparison
sns.set_context('talk')
sns.set_style('darkgrid')
fig=plt.figure(30,figsize=(12,6),dpi=100)

#print('get sunspot number from SIDC')    
#get 13month smoothed sunspot number from SIDC
#http://www.sidc.be/silso/datafiles
#add year;month;year2;spot;stand;obs;check in first row to read pandas dataframe


file='ssn_m.p'  
ssn_m=pickle.load(open(data_path+file, "rb"))

fsize=15
max_spot=450

ax1 = plt.subplot(111) 
ax1.plot(ssn.time,ssn.spot,'-g',alpha=0.4,linewidth=1.0,label='Observed sunspot number (SIDC, daily)')
#ax1.plot(ssn.time,ssn.spot_mean_13,'-k',alpha=0.5,linewidth=1.5,label='Observed sunspot number (SIDC, 13 month smoothed)')
ax1.plot(ssn_m.time+15,ssn_m.spot,'-k',alpha=1,linewidth=1.5,label='Observed sunspot number (SIDC, monthly mean)')



#make pp19 begin time earlier
import copy
times_25_daily_shift=copy.deepcopy(times_25_daily)

for i in np.arange(len(times_25_daily_shift)):
    times_25_daily_shift[i]=times_25_daily_shift[i]-timedelta(days=6*30)

#PP19 prediction
ax1.plot(times_25_daily_shift,spots_predict_25pp_daily,'-b',alpha=1,linewidth=1.5,label='Prediction NOAA/NASA/ISES shifted -6 months')
ax1.fill_between(times_25_daily_shift,spots_predict_25pp_daily_low,spots_predict_25pp_daily_high,alpha=0.2)

#MC20 prediction
#ax1.plot(times_25_daily,spots_predict_25_daily,'--r',alpha=0.5,linewidth=1.0,label='McIntosh et al. (2020)')
#ax1.fill_between(times_25_daily, spots_predict_25_daily_lower68, spots_predict_25_daily_upper68, alpha=0.2)

#MC22 prediction
ax1.plot(times_25_daily,spots_predict_25_daily_mc2,'-r',alpha=1,linewidth=1.5,label='McIntosh/Leamon 2022')
ax1.fill_between(times_25_daily, spots_predict_25_daily_lower68_mc2, spots_predict_25_daily_upper68_mc2, alpha=0.1, color='red')



#mean cycle
ax1.plot(times_25_daily,spots_predict_25m_daily,'-g',alpha=1,linewidth=1.5,label='Mean solar cycle since 1750')
ax1.set_xlim(datetime.datetime(1749,1,1),datetime.datetime(2035,1,1))


#ax1.plot(times_25_daily,spots_predict_25_daily,'-k',alpha=1,linewidth=1.5,label='Solar dynamo revolution April 2021')
#ax1.fill_between(times_25_daily, spots_predict_25_daily_lower68, spots_predict_25_daily_upper68, alpha=0.2)

ax1.set_ylim(0,max_spot)
ax1.set_ylabel('Sunspot number')

plt.legend(loc='upper right',fontsize=8)
plt.tight_layout()

plt.savefig(outputdirectory+'/cycle25_prediction.png',dpi=100)

#with shorter interval


plt.legend(loc='upper right',fontsize=13)
ax1.set_xlim(datetime.datetime(1975,1,1),datetime.datetime(2033,1,1))
ax1.set_ylim(0,320)
years = mdates.YearLocator(5)   # every year
ax1.xaxis.set_major_locator(years)
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)
plt.savefig(outputdirectory+'/cycle25_prediction_short.png',dpi=100)
plt.savefig(outputdirectory+'/cycle25_prediction_short.pdf')


# In[21]:


#with shortest interval

#Extra plot for solar cycle comparison
sns.set_context('talk')
sns.set_style('ticks')
fig=plt.figure(4,figsize=(12,6),dpi=100)

#print('get sunspot number from SIDC')    
#get 13month smoothed sunspot number from SIDC
#http://www.sidc.be/silso/datafiles
#add year;month;year2;spot;stand;obs;check in first row to read pandas dataframe


file='ssn_m.p'  
ssn_m=pickle.load(open(data_path+file, "rb"))

file='ssn_prelim.p'
ssn_p=pickle.load(open(data_path+file, "rb"))

fsize=15
max_spot=420

ax1 = plt.subplot(111) 
ax1.plot(ssn.time,ssn.spot,'-g',alpha=0.5,linewidth=1.0,label='Observed sunspot number (SIDC, daily)')
ax1.plot(ssn_p.time,ssn_p.spot,color='coral',alpha=0.9,linewidth=1.0,label='Observed sunspot number (SIDC, daily, estimated)')
#ax1.plot(ssn.time,ssn.spot_mean_13,'-k',alpha=0.5,linewidth=1.5,label='Observed sunspot number (SIDC, 13 month smoothed)')

#move average to monthly midpoint
ax1.plot(ssn_m.time+15,ssn_m.spot,'-k',alpha=1,linewidth=2,label='Observed sunspot number (SIDC, monthly mean)')

    
#PP19 prediction
ax1.plot(times_25_daily,spots_predict_25pp_daily,color='deepskyblue',alpha=1,linewidth=2.5,label='Prediction Panel NOAA/NASA/ISES')
#ax1.fill_between(times_25_daily,spots_predict_25pp_daily_low,spots_predict_25pp_daily_high,alpha=0.2)

ax1.plot(times_25_daily_shift,spots_predict_25pp_daily,color='blue',alpha=1,linewidth=2.5,label='Prediction Panel NOAA/NASA/ISES shifted -6 months')
ax1.fill_between(times_25_daily_shift,spots_predict_25pp_daily_low,spots_predict_25pp_daily_high,alpha=0.2, color='dodgerblue')



#MC20
#ax1.plot(times_25_daily,spots_predict_25_daily,'--r',alpha=0.5,linewidth=2.5,label='McIntosh et al. (2020)')
#ax1.fill_between(times_25_daily, spots_predict_25_daily_lower68, spots_predict_25_daily_upper68, alpha=0.1, color='red')


#MC22


ax1.plot(times_25_daily,spots_predict_25_daily_mc2,'-r',alpha=1,linewidth=2.5,label='McIntosh/Leamon 2022')
ax1.fill_between(times_25_daily, spots_predict_25_daily_lower68_mc2, spots_predict_25_daily_upper68_mc2, alpha=0.1, color='red')




#mean cycle
ax1.plot(times_25_daily,spots_predict_25m_daily,'-g',alpha=1,linewidth=2.5,label='Mean solar cycle since 1750')
ax1.set_xlim(datetime.datetime(1749,1,1),datetime.datetime(2035,1,1))
ax1.set_ylim(0,max_spot)
ax1.set_ylabel('Sunspot number')


plt.legend(loc='upper left',fontsize=12)
plt.annotate('C. Möstl  helioforecast.space/solarcycle',xy=(0.995,0.02),xycoords='axes fraction',fontsize=9,ha='right')
ax1.set_xlim(datetime.datetime(2019,1,1),datetime.datetime(2024,6,1))
ax1.set_ylim(0,180)
months = mdates.MonthLocator()   # every year
ax1.xaxis.set_minor_locator(months)
ax1.grid(linestyle='--')
plt.tight_layout()
plt.savefig(outputdirectory+'/cycle25_prediction_focus.png',dpi=100)


# # 4 Parker Probe ICME rate prediction

# ### make PSP and Solar Orbiter position

# In[22]:


frame='HEEQ'
starttime =datetime.datetime(2018, 8,13)
endtime = datetime.datetime(2025, 8, 31)
pspt_time = []
res_in_days=1
while starttime < endtime:
    pspt_time.append(starttime)
    starttime += timedelta(days=res_in_days)
pspt_time_num=parse_time(pspt_time).plot_date
starttime =datetime.datetime(2018, 8,13)

spice.furnish(spicedata.get_kernel('psp_pred'))
pspt=spice.Trajectory('SPP')
pspt.generate_positions(pspt_time,'Sun',frame)
psp_speed=pspt.speed.value
pspt.change_units(astropy.units.AU)  
[psp_r, psp_lat, psp_lon]=hd.cart2sphere(pspt.x,pspt.y,pspt.z)
print('PSP pos')
print()

frame='HEEQ'
starttime =datetime.datetime(2020, 3,1)
endtime = datetime.datetime(2030, 9, 1)
solot_time = []
res_in_days=1
while starttime < endtime:
    solot_time.append(starttime)
    starttime += timedelta(days=res_in_days)
solot_time_num=parse_time(solot_time).plot_date     
starttime =datetime.datetime(2020, 3,1)

spice.furnish(spicedata.get_kernel('solo_2020'))
solot=spice.Trajectory('Solar Orbiter')
solot.generate_positions(solot_time,'Sun',frame)
solot.change_units(astropy.units.AU)  
[solo_r, solo_lat, solo_lon]=hd.cart2sphere(solot.x,solot.y,solot.z)
print('Solo pos')





############################################## BepiColombo

starttime =datetime.datetime(2018, 10, 21)
endtime = datetime.datetime(2025, 11, 2)
bepi_time = []
while starttime < endtime:
    bepi_time.append(starttime)
    starttime += timedelta(days=res_in_days)

spice.furnish(spicedata.get_kernel('bepi_pred'))
bepi=spice.Trajectory('BEPICOLOMBO MPO') # or BEPICOLOMBO MMO
bepi.generate_positions(bepi_time,'Sun',frame)
bepi.change_units(astropy.units.AU)  
[bepi_r, bepi_lat, bepi_lon]=hd.cart2sphere(bepi.x,bepi.y,bepi.z)

print('Bepi done')


#add mercury
planet_kernel=spicedata.get_kernel('planet_trajectories')
starttime= datetime.datetime(2025, 11, 3)
endtime = datetime.datetime(2032, 12, 31)
mercury_time = []
while starttime < endtime:
    mercury_time.append(starttime)
    starttime += timedelta(days=res_in_days)



sta_time=solot_time
spice.furnish(spicedata.get_kernel('stereo_a_pred'))
sta=spice.Trajectory('-234')  
sta.generate_positions(solot_time,'Sun',frame)  
sta.change_units(astropy.units.AU)  
[sta_r, sta_lat, sta_lon]=hd.cart2sphere(sta.x,sta.y,sta.z)    
    
    
mercury=spice.Trajectory('1')  #barycenter
mercury.generate_positions(mercury_time,'Sun',frame)  
mercury.change_units(astropy.units.AU)  
[mercury_r, mercury_lat, mercury_lon]=hd.cart2sphere(mercury.x,mercury.y,mercury.z)
print('mercury') 

#combine bepi trajectory with Mercury
bepi_time2=np.hstack([bepi_time,mercury_time])
bepi_r2=np.hstack([bepi_r,mercury_r])
bepi_lon2=np.hstack([bepi_lon,mercury_lon])
bepi_lat2=np.hstack([bepi_lat,mercury_lat])


sns.set_style('darkgrid')
plt.figure(21,dpi=70)
sns.distplot(psp_r)
sns.distplot(solo_r)
sns.distplot(bepi_r)
sns.distplot(sta_r)
plt.ylim(0,5)
plt.xlabel('AU')


# In[23]:


#get the speed in hourly resolution

frame='HEEQ'
starttime =datetime.datetime(2018, 8,13)
endtime = datetime.datetime(2025, 8, 31)
pspt_time_highres = []
res_in_days=1/24.
while starttime < endtime:
    pspt_time_highres.append(starttime)
    starttime += timedelta(days=res_in_days)
pspt_time_highres_num=parse_time(pspt_time_highres).plot_date+ mdates.date2num(np.datetime64('0000-12-31'))
starttime =datetime.datetime(2018, 8,13)

pspt_highres=spice.Trajectory('SPP')
pspt_highres.generate_positions(pspt_time_highres,'Sun',frame)
psp_highres_speed=pspt_highres.speed.value
pspt_highres.change_units(astropy.units.AU)  
[psp_highres_r, psp_highres_lat, psp_highres_lon]=hd.cart2sphere(pspt_highres.x,pspt_highres.y,pspt_highres.z)

plt.figure(22,dpi=70)
plt.title('PSP speed histogram full nominal mission')
plt.plot(psp_highres_r,psp_highres_speed,'.k')
plt.xlabel('AU')

print('psp maximum speed ',np.max(psp_highres_speed),' km/s at ',psp_highres_r[np.argmax(psp_highres_speed)], ' AU')


# In[24]:


#%matplotlib inline
sns.set_context('talk')
sns.set_style('whitegrid')

fig=plt.figure(23,figsize=(13,10),dpi=70)

ax1 = plt.subplot(211) 
ax1.plot_date(pspt_time_highres_num+mdates.date2num(np.datetime64('0000-12-31')),psp_highres_r,c='r',linestyle='-',markersize=0)

ax1.set_ylim(0,0.95)
ax1.set_xlim(datetime.datetime(2018,9,1),datetime.datetime(2025,9,1))
ax1.yaxis.set_ticks(np.arange(0,1,0.1))
ax1.xaxis.set_minor_locator(mdates.MonthLocator())
ax1.tick_params(which="both", bottom=True)
plt.ylabel('PSP distance [AU]')
ax1.grid(which='major',linestyle='--',alpha=0.7)
#ax1.set_zorder(3)

plt.title('Parker Solar Probe trajectory')

#grid
#gridtime=[datetime.datetime(2018, 1,1), datetime.datetime(2026, 1, 1)]
#for i in np.arange(0,1,0.1):
#    ax1.plot(gridtime,[0,0]+i,linestyle='--',color='grey',alpha=0.7,lw=0.8)


    
ax2=ax1.twinx()
ax2.plot_date(pspt_time_highres_num,215.03*psp_highres_r,c='r',linestyle='-',markersize=0)
ax2.set_ylabel('solar radii')
ax2.set_ylim(0,0.95*215.03)
ax2.set_xlim(datetime.datetime(2018,9,1),datetime.datetime(2025,9,1))
ax2.yaxis.set_ticks(np.arange(0,220,20))
ax2.grid(b=None)


sns.set_style('whitegrid')
ax3 = plt.subplot(212) 
ax3.plot_date(pspt_time_highres_num,psp_highres_speed,'-r',zorder=3)
plt.ylabel('PSP speed [km/s]')
#plt.xlabel('year')
ax3.set_xlim(datetime.datetime(2018,9,1),datetime.datetime(2025,9,1))
ax3.set_ylim(0,210)
ax3.xaxis.set_minor_locator(mdates.MonthLocator())
ax3.tick_params(which="both", bottom=True)
ax3.grid(which='major',linestyle='--',alpha=0.7)

#find perihelia and aphelia
dr=np.gradient(psp_highres_r)
dr_minmax=np.where(np.diff(np.sign(dr)))[0]
for i in np.arange(0,48,2):
    ax1.text(pspt_time_highres_num[dr_minmax[i]],psp_highres_r[dr_minmax[i]]-0.038,str(int(i/2)+1),fontsize=13,zorder=0,horizontalalignment='center')



#find perihelia and aphelia
for i in np.arange(0,48,2):
    ax3.text(pspt_time_highres_num[dr_minmax[i]],psp_highres_speed[dr_minmax[i]]+5,str(int(i/2)+1),fontsize=13,horizontalalignment='center')

plt.tight_layout()
plt.figtext(0.99,0.008,'C. Möstl @chrisoutofspace', fontsize=10, ha='right',color='k',style='italic')     


plt.savefig(outputdirectory+'/psp_orbits.png', dpi=100)


# In[25]:


#same thing for Solar Orbiter

frame='HEEQ'
starttime =datetime.datetime(2020, 3,1)
endtime = datetime.datetime(2029, 12, 31)
solo_time_highres = []
res_in_days=1.
while starttime < endtime:
    solo_time_highres.append(starttime)
    starttime += timedelta(days=res_in_days)
solo_time_highres_num=parse_time(solo_time_highres).plot_date+ mdates.date2num(np.datetime64('0000-12-31'))
starttime =datetime.datetime(2020, 3,1)

solo_highres=spice.Trajectory('Solar Orbiter')
solo_highres.generate_positions(solo_time_highres,'Sun',frame)
solo_highres_speed=solo_highres.speed.value
solo_highres.change_units(astropy.units.AU)  
[solo_highres_r, solo_highres_lat, solo_highres_lon]=hd.cart2sphere(solo_highres.x,solo_highres.y,solo_highres.z)

solo_highres_lat=np.degrees(solo_highres_lat)
solo_highres_lon=np.degrees(solo_highres_lon)
print('solo maximum speed ',np.max(solo_highres_speed),' km/s at ',solo_highres_r[np.argmax(solo_highres_speed)], ' AU')


#%matplotlib inline

sns.set_context('talk')
sns.set_style('whitegrid')

fig=plt.figure(24,figsize=(13,10),dpi=70)

ax1 = plt.subplot(211) 
ax1.plot_date(solo_time_highres_num,solo_highres_r,c='r',linestyle='-',markersize=0)

ax1.set_xlim(datetime.datetime(2020,3,1),datetime.datetime(2029,12,31))

ax1.yaxis.set_ticks(np.arange(0,1.3,0.1))
ax1.xaxis.set_minor_locator(mdates.MonthLocator())
ax1.tick_params(which="both", bottom=True)
plt.ylabel('Solar Orbiter distance [AU]')
ax1.grid(which='major',linestyle='--',alpha=0.7)
#ax1.set_zorder(3)

ax1.set_ylim(0,1.1)


plt.title('Solar Orbiter trajectory')

#grid
#gridtime=[datetime.datetime(2018, 1,1), datetime.datetime(2026, 1, 1)]
#for i in np.arange(0,1,0.1):
#    ax1.plot(gridtime,[0,0]+i,linestyle='--',color='grey',alpha=0.7,lw=0.8)

    
ax2=ax1.twinx()
ax2.plot_date(solo_time_highres_num,215.03*solo_highres_r,c='r',linestyle='',markersize=0)
ax2.set_ylabel('solar radii')
ax2.set_xlim(datetime.datetime(2020,3,1),datetime.datetime(2029,12,31))
ax2.yaxis.set_ticks(np.arange(0,260,20))
ax2.grid(b=None)
ax2.set_ylim(0,1.1*215.03)



sns.set_style('whitegrid')
ax3 = plt.subplot(212) 
ax3.plot_date(solo_time_highres_num,solo_highres_lat,'-r',zorder=3)
plt.ylabel('Solar Orbiter HEEQ latitude [degree]')
#plt.xlabel('year')
ax3.set_xlim(datetime.datetime(2020,3,1),datetime.datetime(2029,12,31))

ax3.set_ylim(-35,35)
ax3.xaxis.set_minor_locator(mdates.MonthLocator())
ax3.tick_params(which="both", bottom=True)
ax3.grid(which='major',linestyle='--',alpha=0.7)

#find perihelia and aphelia
dr=np.gradient(solo_highres_r)
dr_minmax=np.where(np.diff(np.sign(dr)))[0]

#plot numbers
for i in np.arange(0,39,2):
    ax1.text(solo_time_highres_num[dr_minmax[i]],solo_highres_r[dr_minmax[i]]-0.07,str(int(i/2)+1),fontsize=15,zorder=0,horizontalalignment='center')
    
    ax1.plot([solo_time_highres_num[dr_minmax[i]],solo_time_highres_num[dr_minmax[i]]],[-50,50],linestyle='-', color='grey',alpha=0.3)


    ax3.text(solo_time_highres_num[dr_minmax[i]]+40,solo_highres_lat[dr_minmax[i]]+1,str(int(i/2)+1),fontsize=15,horizontalalignment='center')    
    ax3.plot([solo_time_highres_num[dr_minmax[i]],solo_time_highres_num[dr_minmax[i]]],[-50,50],linestyle='-', color='grey',alpha=0.4)

plt.tight_layout()
plt.figtext(0.99,0.008,'C. Möstl @chrisoutofspace', fontsize=10, ha='right',color='k',style='italic')     


plt.savefig(outputdirectory+'/solo_orbits.png', dpi=100)


# ### Calculate expected number of PSP ICMEs < 0.3, 0.2, 0.1 AU 

# first calculate smooth functions for the icme rate including the derived error bars in Figure 3

# In[26]:


#fit yearly ICME rates again with hathaway function to get to daily resolution including errors


#define times a bit differently
yearly_mid_times_25_num=mdates.date2num(yearly_mid_times_25)
times_25_daily_icrange= [ yearly_mid_times_25[0] + datetime.timedelta(days=n) for n in range(int ((yearly_mid_times_25[-1] - yearly_mid_times_25[0]).days))]  
times_25_daily_icrange_num=mdates.date2num(times_25_daily_icrange)


##PP19 model

#yearly_mid_times_25,icmes_predict_25pp,yerr=ic_rate_25_pp_std, 
plt.figure(51,figsize=(10,5),dpi=70)
plt.errorbar(yearly_mid_times_25,icmes_predict_25pp,yerr=ic_rate_25_pp_std)

icratepp_low=icmes_predict_25pp-ic_rate_25_pp_std
icratepp_high=icmes_predict_25pp+ic_rate_25_pp_std
#plt.plot(yearly_mid_times_25,icratepp_low)
#plt.plot(yearly_mid_times_25,icratepp_high)

#smooth with spline
from scipy.interpolate import interp1d
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d

fpp = interp1d(yearly_mid_times_25_num, icmes_predict_25pp, kind=2)
fpp_high = interp1d(yearly_mid_times_25_num, icratepp_high, kind=2)
fpp_low = interp1d(yearly_mid_times_25_num, icratepp_low, kind=2)

plt.plot(times_25_daily_icrange_num,fpp(times_25_daily_icrange_num))
plt.plot(times_25_daily_icrange_num,fpp_high(times_25_daily_icrange_num))
plt.plot(times_25_daily_icrange_num,fpp_low(times_25_daily_icrange_num))



#MC20 model
plt.figure(52,figsize=(10,5),dpi=70)


#yearly_mid_times_25,icmes_predict_25mc,yerr=ic_rate_25_mc20_std,
plt.errorbar(yearly_mid_times_25,icmes_predict_25mc,yerr=ic_rate_25_mc20_std)

icratemc_low=icmes_predict_25mc-ic_rate_25_mc20_std
icratemc_high=icmes_predict_25mc+ic_rate_25_mc20_std
#plt.plot(yearly_mid_times_25,icratemc_low)
#plt.plot(yearly_mid_times_25,icratemc_high)

fmc = interp1d(yearly_mid_times_25_num, icmes_predict_25mc, kind=2)
fmc_high = interp1d(yearly_mid_times_25_num, icratemc_high, kind=2)
fmc_low = interp1d(yearly_mid_times_25_num, icratemc_low, kind=2)

plt.plot(times_25_daily_icrange_num,fmc(times_25_daily_icrange_num))
plt.plot(times_25_daily_icrange_num,fmc_high(times_25_daily_icrange_num))
plt.plot(times_25_daily_icrange_num,fmc_low(times_25_daily_icrange_num))


# Figure out how many ICMEs PSP sees < 0.1 AU, < 0.2 AU, < 0.3 AU for the predicted ICME rates

# In[27]:


#make position new in order to be of similar range with ICME rate spline fits
frame='HEEQ'
starttime =times_25_daily_icrange[0]
endtime = datetime.datetime(2025, 9, 1)
pspt_time = []
res_in_days=1
while starttime < endtime:
    pspt_time.append(starttime)
    starttime += timedelta(days=res_in_days)
pspt_time_num=parse_time(pspt_time).plot_date

spice.furnish(spicedata.get_kernel('psp_pred'))
pspt=spice.Trajectory('SPP')
pspt.generate_positions(pspt_time,'Sun',frame)
pspt.change_units(astropy.units.AU)  
[psp_r, psp_lat, psp_lon]=hd.cart2sphere(pspt.x,pspt.y,pspt.z)
print('PSP pos')
print()



#all positions
psp_l10=np.where(psp_r < 1.0)[0]
#positions less 0.3 AU
psp_l03=np.where(psp_r < 0.3)[0]
#positions less 0.2 AU
psp_l02=np.where(psp_r < 0.2)[0]
#positions less 0.1 AU
psp_l01=np.where(psp_r < 0.1)[0]



print('!!!!!! all results below calculated from 2020 July 1 to 2025 August 31!!!!!!!!!!!!!')
print('PSP spends ... ')
print('days < 0.3 AU:',psp_l03.size)
print('days < 0.2 AU:',psp_l02.size)
print('days < 0.1 AU:',psp_l01.size)

# pspt_time_num[psp_l03]
# psp_r[psp_l03]
# times_25_daily_num=mdates.date2num(times_25_daily) 
# ic25=icmes_predict_25mc_daily/365.24


#################### Calculate PSP ICME rate

#MC20
print('---------------')
print('MC20 prediction')
print('using spline fits to get from yearly to daily ICME rates including error bars')

#old version with hathawy function for SSN and conversion instead of spline (result similar for average, but errors not good)
#spots_predict_25_daily_mc_psp=hs.hathaway(pspt_time, start_25-shift_t0,a,b,c)
#icmes_psp_predict_mc_daily=ssn_to_rate(spots_predict_25_daily_mc_psp,linfit)[0]/365.24

icmes_psp_predict_mc_daily=fmc(times_25_daily_icrange_num)/365.24
icmes_psp_predict_mc_daily_low=fmc_low(times_25_daily_icrange_num)/365.24
icmes_psp_predict_mc_daily_high=fmc_high(times_25_daily_icrange_num)/365.24

print('ICME events in situ at PSP total mission : ', int(np.rint(np.sum(icmes_psp_predict_mc_daily[psp_l10]))))
print('range low/high: ',int(np.rint(np.sum(icmes_psp_predict_mc_daily_low[psp_l10]))),' / ',int(np.rint(np.sum(icmes_psp_predict_mc_daily_high[psp_l10]))))
print('ICME events in situ at PSP < 0.3 AU: ', int(np.rint(np.sum(icmes_psp_predict_mc_daily[psp_l03]))))
print('range low/high: ',int(np.rint(np.sum(icmes_psp_predict_mc_daily_low[psp_l03]))),' / ',int(np.rint(np.sum(icmes_psp_predict_mc_daily_high[psp_l03]))))
print('ICME events in situ at PSP < 0.2 AU: ', int(np.rint(np.sum(icmes_psp_predict_mc_daily[psp_l02]))))
print('range low/high: ',int(np.rint(np.sum(icmes_psp_predict_mc_daily_low[psp_l02]))),' / ',int(np.rint(np.sum(icmes_psp_predict_mc_daily_high[psp_l02]))))
print('ICME events in situ at PSP < 0.1 AU: ', int(np.rint(np.sum(icmes_psp_predict_mc_daily[psp_l01]))))
print('range low/high: ',int(np.rint(np.sum(icmes_psp_predict_mc_daily_low[psp_l01]))),' / ',int(np.rint(np.sum(icmes_psp_predict_mc_daily_high[psp_l01]))))


#check if correct
# plt.plot(pspt_time_num,psp_r)
# plt.plot_date(pspt_time_num[psp_l03],icmes_psp_predict_mc_daily[psp_l03],'ob')
# plt.plot_date(pspt_time_num[psp_l02],icmes_psp_predict_mc_daily[psp_l02],'ok')
# plt.plot_date(pspt_time_num[psp_l01],icmes_psp_predict_mc_daily[psp_l01],'or')




print('---------------')
print('PP19 prediction')

#old version with Hathway fit
#PP19 get SSN from Hathaway model
#spots_predict_25_daily_pp_psp=hs.hathaway(pspt_time, start_25_fit,app,bpp,cpp)
#get ICR from SSN                     
#icmes_psp_predict_pp_daily=ssn_to_rate(spots_predict_25_daily_pp_psp,linfit)[0]/365.24


#new version with spline fits
icmes_psp_predict_pp_daily=fpp(times_25_daily_icrange_num)/365.24
icmes_psp_predict_pp_daily_low=fpp_low(times_25_daily_icrange_num)/365.24
icmes_psp_predict_pp_daily_high=fpp_high(times_25_daily_icrange_num)/365.24



# plt.plot_date(pspt_time_num[psp_l03],icmes_psp_predict_pp_daily[psp_l03],'sb')
# plt.plot_date(pspt_time_num[psp_l02],icmes_psp_predict_pp_daily[psp_l02],'sk')
# plt.plot_date(pspt_time_num[psp_l01],icmes_psp_predict_pp_daily[psp_l01],'sr')

print('ICME events in situ at PSP total mission : ', int(np.rint(np.sum(icmes_psp_predict_pp_daily[psp_l10]))))
print('range low/high: ',int(np.rint(np.sum(icmes_psp_predict_pp_daily_low[psp_l10]))),' / ',int(np.rint(np.sum(icmes_psp_predict_pp_daily_high[psp_l10]))))
print('ICME events in situ at PSP < 0.3 AU: ', int(np.rint(np.sum(icmes_psp_predict_pp_daily[psp_l03]))))
print('range low/high: ',int(np.rint(np.sum(icmes_psp_predict_pp_daily_low[psp_l03]))),' / ',int(np.rint(np.sum(icmes_psp_predict_pp_daily_high[psp_l03]))))
print('ICME events in situ at PSP < 0.2 AU: ', int(np.rint(np.sum(icmes_psp_predict_pp_daily[psp_l02]))))
print('range low/high: ',int(np.rint(np.sum(icmes_psp_predict_pp_daily_low[psp_l02]))),' / ',int(np.rint(np.sum(icmes_psp_predict_pp_daily_high[psp_l02]))))
print('ICME events in situ at PSP < 0.1 AU: ', int(np.rint(np.sum(icmes_psp_predict_pp_daily[psp_l01]))))
print('range low/high: ',int(np.rint(np.sum(icmes_psp_predict_pp_daily_low[psp_l01]))),' / ',int(np.rint(np.sum(icmes_psp_predict_pp_daily_high[psp_l01]))))




#---------------------------------------------------- Solar Orbiter

solo_l05=np.where(solo_r < 0.5)[0]
solo_l04=np.where(solo_r < 0.4)[0]
solo_l03=np.where(solo_r < 0.3)[0]

print()
print()
print('-----------------------')
print('SolO spends ... ')
print('days < 0.5 AU:',solo_l05.size)
print('days < 0.4 AU:',solo_l04.size)
print('days < 0.3 AU:',solo_l03.size)


# ## **Figure 4** PSP Solar Orbiter distance and ICME rate

# In[28]:


sns.set_context("talk")     
#sns.set_style('darkgrid')
sns.set_style("ticks")
fsize=15

starttime =datetime.datetime(2020, 3,1)
endtime = datetime.datetime(2029, 12, 31)

fig=plt.figure(45,figsize=(16,12),dpi=70)
ax1 = plt.subplot(211) 
#plot R[AU]
ax1.plot_date(pspt_time,psp_r,'-k')
#grid
for i in np.arange(0.1,1.2,0.1):
    ax1.plot(pspt_time,np.zeros(len(psp_r))+i,linestyle='--',color='k',alpha=0.4,lw=0.8)

ax1.set_ylabel('PSP heliocentric distance [AU]')
ax1.set_xlim(starttime,endtime)
ax1.set_ylim(0,1.1)

#p######################## plot ICR
ax2=ax1.twinx()

#ax2.plot(yearly_mid_times_25,icmes_predict_25/365.24, marker='s',markerfacecolor='coral', alpha=0.9,label='Predicted ICME rate MC20',linestyle='')
ax2.plot(times_25_daily,icmes_predict_25mc_daily*30.25/365.24, color='coral', label='Predicted ICME rate MC20',linestyle='-')
ax2.plot(times_25_daily,icmes_predict_25pp_daily*30.25/365.24, color='magenta', label='Predicted ICME rate PP19',linestyle='-')

#error bars - 30.42/364.24 converts the yearly rate from the spline function to monthly rates like 3 lines above
#ax2.plot(times_25_daily_icrange_num,fmc_low(times_25_daily_icrange_num)*30.42/365.24,color='coral',linestyle='--')
#ax2.plot(times_25_daily_icrange_num,fmc_high(times_25_daily_icrange_num)*30.42/365.24,color='coral',linestyle='--')
ax2.fill_between(times_25_daily_icrange_num, fmc_low(times_25_daily_icrange_num)*30.42/365.248, fmc_high(times_25_daily_icrange_num)*30.42/365.24, alpha=0.2,color='coral')

#same for pp model
#ax2.plot(times_25_daily_icrange_num,fpp_low(times_25_daily_icrange_num)*30.42/365.24,color='magenta',linestyle='--')
#ax2.plot(times_25_daily_icrange_num,fpp_high(times_25_daily_icrange_num)*30.42/365.24,color='magenta',linestyle='--')
ax2.fill_between(times_25_daily_icrange_num, fpp_low(times_25_daily_icrange_num)*30.42/365.248, fpp_high(times_25_daily_icrange_num)*30.42/365.24, alpha=0.2,color='blueviolet')

ax2.set_ylabel('monthly ICME rate')
plt.legend(loc='upper center')

plotstart=datetime.datetime(2020, 7, 1)
plotend=datetime.datetime(2025, 9,1)
plt.xlim(plotstart,plotend)


################ Solar Orbiter
ax3 = plt.subplot(212) 
ax3.plot_date(solot_time,solo_r,'-k')
ax3.plot_date(bepi_time2,bepi_r2,linestyle='-',color='orange',lw=2,markersize=0)


#ax3.plot_date(solot_time,np.degrees(solo_lat)/100,'--k')
#ax3.plot(solot_time,np.zeros(len(solo_r))+0.5,linestyle='--',color='b',alpha=.6)
#ax3.plot(solot_time,np.zeros(len(solo_r))+0.4,linestyle='--',color='b',alpha=.6)
#ax3.plot(solot_time,np.zeros(len(solo_r))+0.3,linestyle='--',color='b',alpha=.6)

ax3.set_ylabel('Heliocentric distance R [AU]')
ax3.set_xlim(starttime,endtime)

for i in np.arange(0.1,1.2,0.1):
    ax3.plot(solot_time,np.zeros(len(solo_r))+i,linestyle='--',color='k',alpha=0.4,lw=0.8)

ax3.set_ylim(0,1.2)

#p######################## plot ICR
ax4=ax3.twinx()

#ax2.plot(yearly_mid_times_25,icmes_predict_25/365.24, marker='s',markerfacecolor='coral', alpha=0.9,label='Predicted ICME rate MC20',linestyle='')
ax4.plot(times_25_daily,icmes_predict_25mc_daily*30.42/365.24, color='coral', label='Predicted ICME rate MC20',linestyle='-')
ax4.plot(times_25_daily,icmes_predict_25pp_daily*30.42/365.24, color='magenta', label='Predicted ICME rate PP19',linestyle='-')

#error ranges
ax4.fill_between(times_25_daily_icrange_num, fmc_low(times_25_daily_icrange_num)*30.42/365.248, fmc_high(times_25_daily_icrange_num)*30.42/365.24, alpha=0.2,color='coral')
ax4.fill_between(times_25_daily_icrange_num, fpp_low(times_25_daily_icrange_num)*30.42/365.248, fpp_high(times_25_daily_icrange_num)*30.42/365.24, alpha=0.2,color='blueviolet')
ax4.set_ylabel('monthly ICME rate')
plt.legend(loc=1)

plotstart=datetime.datetime(2020, 7, 1)
plotend=datetime.datetime(2029, 12,31)
plt.xlim(plotstart,plotend)

plt.tight_layout()


#plt.annotate('(a)',[0.01,0.965],xycoords='figure fraction',weight='bold')
#plt.annotate('(b)',[0.01,0.475],xycoords='figure fraction',weight='bold')
plt.annotate('Parker Solar Probe',[0.08,0.95],xycoords='figure fraction',weight='bold')
plt.annotate('Solar Orbiter',[0.08,0.455],xycoords='figure fraction',weight='bold')
plt.annotate('Bepi Colombo',[0.20,0.455],xycoords='figure fraction',weight='bold',color='orange')

#plt.savefig('results/plots_rate/fig4_psp_rate.pdf', dpi=300)
plt.savefig(outputdirectory+'/cycle25_icme_rate_psp_orbiter_bepi.png', dpi=100)


# ### calculate lineups

# first homogenize spacecraft positions, load hourly position data
# 

# In[29]:


[psp, bepi, solo, sta, earth, venus, mars, mercury,frame]=       pickle.load( open( '/nas/helio/data/insitu_python/positions_psp_solo_bepi_sta_planets_HEEQ_1hour.p', "rb" ))
#combine bepi trajectory with Mercury after orbit insertion

orbit_insertion=bepi.time[-1]
moi_ind=np.where(orbit_insertion<mercury.time)[0][0]

bepi_time2=np.hstack([bepi.time,mercury.time[moi_ind:-1]])
bepi_r2=np.hstack([bepi.r,mercury.r[moi_ind:-1]])
bepi_lon2=np.hstack([bepi.lon,mercury.lon[moi_ind:-1]])
bepi_lat2=np.hstack([bepi.lat,mercury.lat[moi_ind:-1]])




print(mdates.num2date(earth.time[0]))
print(mdates.num2date(sta.time[0]))
print(mdates.num2date(psp.time[0]))
print(mdates.num2date(bepi_time2[0]))
print(mdates.num2date(solo.time[0]))


print(len(earth.lon))
print(len(sta.lon))
print(len(bepi_lon2))
print(len(psp.lon))
print(len(solo.lon))
     
     

#cut all so consistent with solo start

solo=solo[0:-2]

beg_ind=np.where(earth.time>solo.time[0])[0][0]-1
earth=earth[beg_ind:-2]

beg_ind=np.where(sta.time>solo.time[0])[0][0]-1
sta=sta[beg_ind:-2]

beg_ind=np.where(psp.time>solo.time[0])[0][0]-1
psp=psp[beg_ind:-1]

beg_ind=np.where(bepi_time2>solo.time[0])[0][0]-1
bepi_time2=bepi_time2[beg_ind:-1]
bepi_r2=bepi_r2[beg_ind:-1]
bepi_lon2=bepi_lon2[beg_ind:-1]
bepi_lat2=bepi_lat2[beg_ind:-1]



print()
print(mdates.num2date(earth.time[0]))
print(mdates.num2date(sta.time[0]))
print(mdates.num2date(bepi_time2[0]))
print(mdates.num2date(solo.time[0]))
print(mdates.num2date(psp.time[0]))

print()
print(mdates.num2date(earth.time[-1]))
print(mdates.num2date(sta.time[-1]))
print(mdates.num2date(bepi_time2[-1]))
print(mdates.num2date(solo.time[-1]))
print(mdates.num2date(psp.time[-1]))


print()

print('SolO  ',len(solo.lon))
print('Earth ',len(earth.lon))
print('STA   ',len(sta.lon))
print('Bepi  ',len(bepi_lon2))
print('PSP   ',len(psp.lon))


# ## define and calculate lineups

# ### first until end of PSP nominal, so make arrays all similar

# In[30]:


solo1=solo[0:len(psp)]
sta1=sta[0:len(psp)]
win1=earth[0:len(psp)]
#set Wind longitude to 0
win1.lon=np.zeros(len(psp))
bepi1_time=bepi_time2[0:len(psp)]
bepi1_lon=bepi_lon2[0:len(psp)]
psp1=psp

#stack longitude arrays and convert to degrees from 0 to 360 (Earth at 180)

a=np.transpose(np.vstack((win1.lon,sta1.lon,solo1.lon,bepi1_lon,psp1.lon)))
a=np.rad2deg(a)
a=a+180
sns.distplot(a)


# In[31]:


####lineup counter arrays for all spacecraft
alldiff=np.zeros((len(a),10))

lineup_counter30=0
counter30=0
all30=np.zeros((len(a)))
all30[:]=np.nan

lineup_counter20=0
counter20=0
all20=np.zeros((len(a)))
all20[:]=np.nan

lineup_counter10=0
counter10=0
all10=np.zeros((len(a)))
all10[:]=np.nan

lineup_counter5=0
counter5=0
all5=np.zeros((len(a)))
all5[:]=np.nan


####arrays for individual spacecraft
#win1.lon,sta1.lon,solo1.lon,bepi1_lon,psp1.lon

#5,10,20,30 => 4
lineup_counter_wind=0
counter_wind=0
#5,10,20,30
winddeltas=np.zeros((len(a),4))
winddeltas[:]=np.nan

lineup_counter_sta=0
counter_sta=0
stadeltas=np.zeros((len(a),4))
stadeltas[:]=np.nan

lineup_counter_solo=0
counter_solo=0
solodeltas=np.zeros((len(a),4))
solodeltas[:]=np.nan

lineup_counter_bepi=0
counter_bepi=0
bepideltas=np.zeros((len(a),4))
bepideltas[:]=np.nan

lineup_counter_psp=0
counter_psp=0
pspdeltas=np.zeros((len(a),4))
pspdeltas[:]=np.nan



#go through all hourly times
for i in np.arange(0,len(a)):
    
    #current longitude
    b=a[i]   
    
    #array of all absolute differences between all 5 spacecraft
    #order in a and b is 0 win1.lon 1 sta1.lon 2 solo1.lon 3 bepi1_lon 4 psp1.lon
    bm=abs(np.round([b[0]-b[1],b[0]-b[2],b[0]-b[3],b[0]-b[4],
    b[1]-b[2],b[1]-b[3],b[1]-b[4],
    b[2]-b[3],b[2]-b[4],
    b[3]-b[4]],1))
    alldiff[i,:]=bm
    
    
    #for all spacecraft
    lineup_ind30=np.where(np.logical_or(bm < 30,bm > 330)==True)[0]
    #sum all lineup events
    lineup_counter30=lineup_counter30+len(lineup_ind30)
    #count timestep if there is at least one lineup
    if len(lineup_ind30 >0): 
        all30[i]=1
        counter30=counter30+1

    lineup_ind20=np.where(np.logical_or(bm < 20,bm > 340)==True)[0]
    lineup_counter20=lineup_counter20+len(lineup_ind20)
    if len(lineup_ind20 >0): 
        all20[i]=1
        counter20=counter20+1

    lineup_ind10=np.where(np.logical_or(bm < 10,bm > 350)==True)[0]
    lineup_counter10=lineup_counter10+len(lineup_ind10)
    if len(lineup_ind10 >0): 
        all10[i]=1
        counter10=counter10+1

    lineup_ind5=np.where(np.logical_or(bm < 5,bm > 355)==True)[0]
    lineup_counter5=lineup_counter5+len(lineup_ind5)
    if len(lineup_ind5 >0): 
        all5[i]=1
        counter5=counter5+1
        
    #for individual sc
    
    #order in a and b is  win1.lon,sta1.lon,solo1.lon,bepi1_lon,psp1.lon
    #delta longitudes for each spacecraft to all other spacecraft
    bmwind=bm[[0,1,2,3]]
    bmsta=bm[[0,4,5,6]]
    bmsolo=bm[[1,4,7,8]]
    bmbepi=bm[[2,5,7,9]]
    bmpsp=bm[[3,6,8,9]]

    #wind
    bmnow=bmwind
    lineup_ind30=np.where(np.logical_or(bmnow < 30,bmnow > 330)==True)[0]
    lineup_ind20=np.where(np.logical_or(bmnow < 20,bmnow > 340)==True)[0]
    lineup_ind10=np.where(np.logical_or(bmnow < 10,bmnow > 350)==True)[0]
    lineup_ind5=np.where(np.logical_or(bmnow < 5,bmnow > 355)==True)[0]

    if len(lineup_ind30 >0): winddeltas[i,0]=1 
    if len(lineup_ind20 >0): winddeltas[i,1]=1 
    if len(lineup_ind10 >0): winddeltas[i,2]=1 
    if len(lineup_ind5 >0): winddeltas[i,3]=1 

    #sta
    bmnow=bmsta
    lineup_ind30=np.where(np.logical_or(bmnow < 30,bmnow > 330)==True)[0]
    lineup_ind20=np.where(np.logical_or(bmnow < 20,bmnow > 340)==True)[0]
    lineup_ind10=np.where(np.logical_or(bmnow < 10,bmnow > 350)==True)[0]
    lineup_ind5=np.where(np.logical_or(bmnow < 5,bmnow > 355)==True)[0]

    if len(lineup_ind30 >0): stadeltas[i,0]=1 
    if len(lineup_ind20 >0): stadeltas[i,1]=1 
    if len(lineup_ind10 >0): stadeltas[i,2]=1 
    if len(lineup_ind5 >0): stadeltas[i,3]=1 

    #bepi
    bmnow=bmbepi
    lineup_ind30=np.where(np.logical_or(bmnow < 30,bmnow > 330)==True)[0]
    lineup_ind20=np.where(np.logical_or(bmnow < 20,bmnow > 340)==True)[0]
    lineup_ind10=np.where(np.logical_or(bmnow < 10,bmnow > 350)==True)[0]
    lineup_ind5=np.where(np.logical_or(bmnow < 5,bmnow > 355)==True)[0]

    if len(lineup_ind30 >0): bepideltas[i,0]=1 
    if len(lineup_ind20 >0): bepideltas[i,1]=1 
    if len(lineup_ind10 >0): bepideltas[i,2]=1 
    if len(lineup_ind5 >0): bepideltas[i,3]=1 
        
    #solo
    bmnow=bmsolo
    lineup_ind30=np.where(np.logical_or(bmnow < 30,bmnow > 330)==True)[0]
    lineup_ind20=np.where(np.logical_or(bmnow < 20,bmnow > 340)==True)[0]
    lineup_ind10=np.where(np.logical_or(bmnow < 10,bmnow > 350)==True)[0]
    lineup_ind5=np.where(np.logical_or(bmnow < 5,bmnow > 355)==True)[0]

    if len(lineup_ind30 >0): solodeltas[i,0]=1 
    if len(lineup_ind20 >0): solodeltas[i,1]=1 
    if len(lineup_ind10 >0): solodeltas[i,2]=1 
    if len(lineup_ind5 >0): solodeltas[i,3]=1 

        
    #psp
    bmnow=bmpsp
    lineup_ind30=np.where(np.logical_or(bmnow < 30,bmnow > 330)==True)[0]
    lineup_ind20=np.where(np.logical_or(bmnow < 20,bmnow > 340)==True)[0]
    lineup_ind10=np.where(np.logical_or(bmnow < 10,bmnow > 350)==True)[0]
    lineup_ind5=np.where(np.logical_or(bmnow < 5,bmnow > 355)==True)[0]

    if len(lineup_ind30 >0): pspdeltas[i,0]=1 
    if len(lineup_ind20 >0): pspdeltas[i,1]=1 
    if len(lineup_ind10 >0): pspdeltas[i,2]=1 
    if len(lineup_ind5 >0): pspdeltas[i,3]=1 


    


# 
# ### lineup results

# In[32]:


#Möstl et al. 2020    #until middle 2025 total ICME
#PP19: 109
#MC20: 186
i1=109
i2=186

  

print('total data points',len(a))
print('at least 2 spacecraft within 30 deg longitude %',int(counter30/len(a)*100))
print('at least 2 spacecraft within 20 deg longitude %',int(counter20/len(a)*100))
print('at least 2 spacecraft within 10 deg longitude %',int(counter10/len(a)*100))
print('at least 2 spacecraft within  5 deg longitude %', int(counter5/len(a)*100))

print()
print('Wind')
print('delta lon 30° lineup time %:',int(np.nansum(winddeltas[:,0])/len(a)*100),'ICMEs: ',int(i1*np.nansum(winddeltas[:,0])/len(a)),int(i2*np.nansum(winddeltas[:,0])/len(a)) )
print('delta lon 20° lineup time %:',int(np.nansum(winddeltas[:,1])/len(a)*100),'ICMEs: ',int(i1*np.nansum(winddeltas[:,1])/len(a)),int(i2*np.nansum(winddeltas[:,1])/len(a)) )
print('delta lon 10° lineup time %:',int(np.nansum(winddeltas[:,2])/len(a)*100),'ICMEs: ',int(i1*np.nansum(winddeltas[:,2])/len(a)),int(i2*np.nansum(winddeltas[:,2])/len(a)) )
print('delta lon  5° lineup time %:',int(np.nansum(winddeltas[:,3])/len(a)*100),'ICMEs: ',int(i1*np.nansum(winddeltas[:,3])/len(a)),int(i2*np.nansum(winddeltas[:,3])/len(a)))




fig=plt.figure(61,figsize=(10,6),dpi=100)
plt.plot_date(win1.time,winddeltas[:,0]*30)
plt.plot_date(win1.time,winddeltas[:,1]*20)
plt.plot_date(win1.time,winddeltas[:,2]*10)
plt.plot_date(win1.time,winddeltas[:,3]*5)
plt.ylim(0,35)

print()
print('STA')
print('delta lon 30° lineup time %:',int(np.nansum(stadeltas[:,0])/len(a)*100),'ICMEs: ',int(i1*np.nansum(stadeltas[:,0])/len(a)),int(i2*np.nansum(stadeltas[:,0])/len(a)) )
print('delta lon 20° lineup time %:',int(np.nansum(stadeltas[:,1])/len(a)*100),'ICMEs: ',int(i1*np.nansum(stadeltas[:,1])/len(a)),int(i2*np.nansum(stadeltas[:,1])/len(a)) )
print('delta lon 10° lineup time %:',int(np.nansum(stadeltas[:,2])/len(a)*100),'ICMEs: ',int(i1*np.nansum(stadeltas[:,2])/len(a)),int(i2*np.nansum(stadeltas[:,2])/len(a)) )
print('delta lon  5° lineup time %:',int(np.nansum(stadeltas[:,3])/len(a)*100),'ICMEs: ',int(i1*np.nansum(stadeltas[:,3])/len(a)),int(i2*np.nansum(stadeltas[:,3])/len(a)))



print()
print('Solo')
print('delta lon 30° lineup time %:',int(np.nansum(solodeltas[:,0])/len(a)*100),'ICMEs: ',int(i1*np.nansum(solodeltas[:,0])/len(a)),int(i2*np.nansum(solodeltas[:,0])/len(a)) )
print('delta lon 20° lineup time %:',int(np.nansum(solodeltas[:,1])/len(a)*100),'ICMEs: ',int(i1*np.nansum(solodeltas[:,1])/len(a)),int(i2*np.nansum(solodeltas[:,1])/len(a)) )
print('delta lon 10° lineup time %:',int(np.nansum(solodeltas[:,2])/len(a)*100),'ICMEs: ',int(i1*np.nansum(solodeltas[:,2])/len(a)),int(i2*np.nansum(solodeltas[:,2])/len(a)) )
print('delta lon  5° lineup time %:',int(np.nansum(solodeltas[:,3])/len(a)*100),'ICMEs: ',int(i1*np.nansum(solodeltas[:,3])/len(a)),int(i2*np.nansum(solodeltas[:,3])/len(a)))



print()
print('Bepi')
print('delta lon 30° lineup time %:',int(np.nansum(bepideltas[:,0])/len(a)*100),'ICMEs: ',int(i1*np.nansum(bepideltas[:,0])/len(a)),int(i2*np.nansum(bepideltas[:,0])/len(a)) )
print('delta lon 20° lineup time %:',int(np.nansum(bepideltas[:,1])/len(a)*100),'ICMEs: ',int(i1*np.nansum(bepideltas[:,1])/len(a)),int(i2*np.nansum(bepideltas[:,1])/len(a)) )
print('delta lon 10° lineup time %:',int(np.nansum(bepideltas[:,2])/len(a)*100),'ICMEs: ',int(i1*np.nansum(bepideltas[:,2])/len(a)),int(i2*np.nansum(bepideltas[:,2])/len(a)) )
print('delta lon  5° lineup time %:',int(np.nansum(bepideltas[:,3])/len(a)*100),'ICMEs: ',int(i1*np.nansum(bepideltas[:,3])/len(a)),int(i2*np.nansum(bepideltas[:,3])/len(a)))


print()
print('PSP')
print('delta lon 30° lineup time %:',int(np.nansum(pspdeltas[:,0])/len(a)*100),'ICMEs: ',int(i1*np.nansum(pspdeltas[:,0])/len(a)),int(i2*np.nansum(pspdeltas[:,0])/len(a)) )
print('delta lon 20° lineup time %:',int(np.nansum(pspdeltas[:,1])/len(a)*100),'ICMEs: ',int(i1*np.nansum(pspdeltas[:,1])/len(a)),int(i2*np.nansum(pspdeltas[:,1])/len(a)) )
print('delta lon 10° lineup time %:',int(np.nansum(pspdeltas[:,2])/len(a)*100),'ICMEs: ',int(i1*np.nansum(pspdeltas[:,2])/len(a)),int(i2*np.nansum(pspdeltas[:,2])/len(a)) )
print('delta lon  5° lineup time %:',int(np.nansum(pspdeltas[:,3])/len(a)*100),'ICMEs: ',int(i1*np.nansum(pspdeltas[:,3])/len(a)),int(i2*np.nansum(pspdeltas[:,3])/len(a)))





#plt.plot_date(win1.time,a,'-')

fig=plt.figure(62,figsize=(10,6),dpi=100)
plt.plot_date(win1.time,alldiff,'-')  
plt.ylim(0,35)


# ### second interval from end of PSP nominal to 2029

# In[33]:


solo2=solo[len(psp):-1]
sta2=sta[len(psp):-1]
win2=earth[len(psp):-1]
bepi2_time=bepi_time2[len(psp):-1]
bepi2_lon=bepi_lon2[len(psp):-1]

#stack longitude arrays and convert to degrees
a2=np.transpose(np.vstack((win2.lon,sta2.lon,solo2.lon,bepi2_lon)))
a2=np.rad2deg(a2)

a2=a2+180
sns.distplot(a2)


# In[34]:


alldiff2=np.zeros((len(a2),6))


lineup_counter30_2=0
counter30_2=0
all30_2=np.zeros((len(a2)))
all30_2[:]=np.nan

lineup_counter20_2=0
counter20_2=0
all20_2=np.zeros((len(a2)))
all20_2[:]=np.nan

lineup_counter10_2=0
counter10_2=0
all10_2=np.zeros((len(a2)))
all10_2[:]=np.nan

lineup_counter5_2=0
counter5_2=0
all5_2=np.zeros((len(a2)))
all5_2[:]=np.nan




####arrays for individual spacecraft
#win1.lon,sta1.lon,solo1.lon,bepi1_lon,psp1.lon

#5,10,20,30 => 4
lineup_counter_wind=0
counter_wind=0
#5,10,20,30
winddeltas=np.zeros((len(a),4))
winddeltas[:]=np.nan

lineup_counter_sta=0
counter_sta=0
stadeltas=np.zeros((len(a),4))
stadeltas[:]=np.nan

lineup_counter_solo=0
counter_solo=0
solodeltas=np.zeros((len(a),4))
solodeltas[:]=np.nan

lineup_counter_bepi=0
counter_bepi=0
bepideltas=np.zeros((len(a),4))
bepideltas[:]=np.nan

lineup_counter_psp=0
counter_psp=0
pspdeltas=np.zeros((len(a),4))
pspdeltas[:]=np.nan



for i in np.arange(0,len(a2)):
    
    #current longitude
    b2=a2[i]   
    
    #array of all absolute differences between all 5 spacecraft
    #win2.lon,sta2.lon,solo2.lon,bepi2_lon
    bm2=abs(np.round([b2[0]-b2[1],b2[0]-b2[2],b2[0]-b2[3],
    b2[1]-b2[2],b2[1]-b2[3], b2[2]-b2[3]],1))
    
    alldiff2[i,:]=bm2
    
    lineup_ind30_2=np.where(np.logical_or(bm2 < 30,bm2 > 330)==True)[0]
    lineup_counter30_2=lineup_counter30_2+len(lineup_ind30_2)
    if len(lineup_ind30_2 >0): 
        all30_2[i]=1
        counter30_2=counter30_2+1

    lineup_ind20_2=np.where(np.logical_or(bm2 < 20,bm2 > 340)==True)[0]
    lineup_counter20_2=lineup_counter20_2+len(lineup_ind20_2)
    if len(lineup_ind20_2 >0): 
        all20_2[i]=1
        counter20_2=counter20_2+1

    lineup_ind10_2=np.where(np.logical_or(bm2 < 10,bm2 > 350)==True)[0]
    lineup_counter10_2=lineup_counter10_2+len(lineup_ind10_2)
    if len(lineup_ind10_2 >0): 
        all10_2[i]=1
        counter10_2=counter10_2+1

    lineup_ind5_2=np.where(np.logical_or(bm2 < 5,bm2 > 355)==True)[0]
    lineup_counter5_2=lineup_counter5_2+len(lineup_ind5_2)
    if len(lineup_ind5_2 >0): 
        all5_2[i]=1
        counter5_2=counter5_2+1
        
        
        
    #for individual sc
    
    #order in a and b is  win1.lon,sta1.lon,solo1.lon,bepi1_lon,psp1.lon
    #delta longitudes for each spacecraft to all other spacecraft
    bmwind=bm2[[0,1,2]]
    bmsta=bm2[[0,3,4]]
    bmsolo=bm2[[1,3,5]]
    bmbepi=bm2[[2,4,5]]

    #wind
    bmnow=bmwind
    lineup_ind30=np.where(np.logical_or(bmnow < 30,bmnow > 330)==True)[0]
    lineup_ind20=np.where(np.logical_or(bmnow < 20,bmnow > 340)==True)[0]
    lineup_ind10=np.where(np.logical_or(bmnow < 10,bmnow > 350)==True)[0]
    lineup_ind5=np.where(np.logical_or(bmnow < 5,bmnow > 355)==True)[0]

    if len(lineup_ind30 >0): winddeltas[i,0]=1 
    if len(lineup_ind20 >0): winddeltas[i,1]=1 
    if len(lineup_ind10 >0): winddeltas[i,2]=1 
    if len(lineup_ind5 >0): winddeltas[i,3]=1 

    #sta
    bmnow=bmsta
    lineup_ind30=np.where(np.logical_or(bmnow < 30,bmnow > 330)==True)[0]
    lineup_ind20=np.where(np.logical_or(bmnow < 20,bmnow > 340)==True)[0]
    lineup_ind10=np.where(np.logical_or(bmnow < 10,bmnow > 350)==True)[0]
    lineup_ind5=np.where(np.logical_or(bmnow < 5,bmnow > 355)==True)[0]

    if len(lineup_ind30 >0): stadeltas[i,0]=1 
    if len(lineup_ind20 >0): stadeltas[i,1]=1 
    if len(lineup_ind10 >0): stadeltas[i,2]=1 
    if len(lineup_ind5 >0): stadeltas[i,3]=1 

    #bepi
    bmnow=bmbepi
    lineup_ind30=np.where(np.logical_or(bmnow < 30,bmnow > 330)==True)[0]
    lineup_ind20=np.where(np.logical_or(bmnow < 20,bmnow > 340)==True)[0]
    lineup_ind10=np.where(np.logical_or(bmnow < 10,bmnow > 350)==True)[0]
    lineup_ind5=np.where(np.logical_or(bmnow < 5,bmnow > 355)==True)[0]

    if len(lineup_ind30 >0): bepideltas[i,0]=1 
    if len(lineup_ind20 >0): bepideltas[i,1]=1 
    if len(lineup_ind10 >0): bepideltas[i,2]=1 
    if len(lineup_ind5 >0): bepideltas[i,3]=1 
        
    #solo
    bmnow=bmsolo
    lineup_ind30=np.where(np.logical_or(bmnow < 30,bmnow > 330)==True)[0]
    lineup_ind20=np.where(np.logical_or(bmnow < 20,bmnow > 340)==True)[0]
    lineup_ind10=np.where(np.logical_or(bmnow < 10,bmnow > 350)==True)[0]
    lineup_ind5=np.where(np.logical_or(bmnow < 5,bmnow > 355)==True)[0]

    if len(lineup_ind30 >0): solodeltas[i,0]=1 
    if len(lineup_ind20 >0): solodeltas[i,1]=1 
    if len(lineup_ind10 >0): solodeltas[i,2]=1 
    if len(lineup_ind5 >0): solodeltas[i,3]=1 

        
    #psp
    bmnow=bmpsp
    lineup_ind30=np.where(np.logical_or(bmnow < 30,bmnow > 330)==True)[0]
    lineup_ind20=np.where(np.logical_or(bmnow < 20,bmnow > 340)==True)[0]
    lineup_ind10=np.where(np.logical_or(bmnow < 10,bmnow > 350)==True)[0]
    lineup_ind5=np.where(np.logical_or(bmnow < 5,bmnow > 355)==True)[0]

    if len(lineup_ind30 >0): pspdeltas[i,0]=1 
    if len(lineup_ind20 >0): pspdeltas[i,1]=1 
    if len(lineup_ind10 >0): pspdeltas[i,2]=1 
    if len(lineup_ind5 >0): pspdeltas[i,3]=1 



     


# ### Results for 2nd interval

# In[35]:


#Möstl et al. 2020    #until middle 2025 total ICME, assume middle 2025 to end of 2029 is similar, also because of declining phase
#PP19: 109
#MC20: 186
i1=109*4.5/5
i2=186*4.5/5

 

print('total data points',len(a2))
print('at least 2 spacecraft within 30 deg longitude %',int(counter30_2/len(a2)*100))
print('at least 2 spacecraft within 20 deg longitude %',int(counter20_2/len(a2)*100))   
print('at least 2 spacecraft within 10 deg longitude %',int(counter10_2/len(a2)*100))   
print('at least 2 spacecraft within 5 deg longitude %',int(counter5_2/len(a2)*100))   


#until middle 2025 total ICME
#PP19: 100
#MC20: 200


print()
print('Wind')
print('delta lon 30° lineup time %:',int(np.nansum(winddeltas[:,0])/len(a)*100),'ICMEs: ',int(i1*np.nansum(winddeltas[:,0])/len(a)),int(i2*np.nansum(winddeltas[:,0])/len(a)) )
print('delta lon 20° lineup time %:',int(np.nansum(winddeltas[:,1])/len(a)*100),'ICMEs: ',int(i1*np.nansum(winddeltas[:,1])/len(a)),int(i2*np.nansum(winddeltas[:,1])/len(a)) )
print('delta lon 10° lineup time %:',int(np.nansum(winddeltas[:,2])/len(a)*100),'ICMEs: ',int(i1*np.nansum(winddeltas[:,2])/len(a)),int(i2*np.nansum(winddeltas[:,2])/len(a)) )
print('delta lon  5° lineup time %:',int(np.nansum(winddeltas[:,3])/len(a)*100),'ICMEs: ',int(i1*np.nansum(winddeltas[:,3])/len(a)),int(i2*np.nansum(winddeltas[:,3])/len(a)))




fig=plt.figure(61,figsize=(10,6),dpi=100)
plt.plot_date(win1.time,winddeltas[:,0]*30)
plt.plot_date(win1.time,winddeltas[:,1]*20)
plt.plot_date(win1.time,winddeltas[:,2]*10)
plt.plot_date(win1.time,winddeltas[:,3]*5)
plt.ylim(0,35)

print()
print('STA')
print('delta lon 30° lineup time %:',int(np.nansum(stadeltas[:,0])/len(a)*100),'ICMEs: ',int(i1*np.nansum(stadeltas[:,0])/len(a)),int(i2*np.nansum(stadeltas[:,0])/len(a)) )
print('delta lon 20° lineup time %:',int(np.nansum(stadeltas[:,1])/len(a)*100),'ICMEs: ',int(i1*np.nansum(stadeltas[:,1])/len(a)),int(i2*np.nansum(stadeltas[:,1])/len(a)) )
print('delta lon 10° lineup time %:',int(np.nansum(stadeltas[:,2])/len(a)*100),'ICMEs: ',int(i1*np.nansum(stadeltas[:,2])/len(a)),int(i2*np.nansum(stadeltas[:,2])/len(a)) )
print('delta lon  5° lineup time %:',int(np.nansum(stadeltas[:,3])/len(a)*100),'ICMEs: ',int(i1*np.nansum(stadeltas[:,3])/len(a)),int(i2*np.nansum(stadeltas[:,3])/len(a)))



print()
print('Solo')
print('delta lon 30° lineup time %:',int(np.nansum(solodeltas[:,0])/len(a)*100),'ICMEs: ',int(i1*np.nansum(solodeltas[:,0])/len(a)),int(i2*np.nansum(solodeltas[:,0])/len(a)) )
print('delta lon 20° lineup time %:',int(np.nansum(solodeltas[:,1])/len(a)*100),'ICMEs: ',int(i1*np.nansum(solodeltas[:,1])/len(a)),int(i2*np.nansum(solodeltas[:,1])/len(a)) )
print('delta lon 10° lineup time %:',int(np.nansum(solodeltas[:,2])/len(a)*100),'ICMEs: ',int(i1*np.nansum(solodeltas[:,2])/len(a)),int(i2*np.nansum(solodeltas[:,2])/len(a)) )
print('delta lon  5° lineup time %:',int(np.nansum(solodeltas[:,3])/len(a)*100),'ICMEs: ',int(i1*np.nansum(solodeltas[:,3])/len(a)),int(i2*np.nansum(solodeltas[:,3])/len(a)))



print()
print('Bepi')
print('delta lon 30° lineup time %:',int(np.nansum(bepideltas[:,0])/len(a)*100),'ICMEs: ',int(i1*np.nansum(bepideltas[:,0])/len(a)),int(i2*np.nansum(bepideltas[:,0])/len(a)) )
print('delta lon 20° lineup time %:',int(np.nansum(bepideltas[:,1])/len(a)*100),'ICMEs: ',int(i1*np.nansum(bepideltas[:,1])/len(a)),int(i2*np.nansum(bepideltas[:,1])/len(a)) )
print('delta lon 10° lineup time %:',int(np.nansum(bepideltas[:,2])/len(a)*100),'ICMEs: ',int(i1*np.nansum(bepideltas[:,2])/len(a)),int(i2*np.nansum(bepideltas[:,2])/len(a)) )
print('delta lon  5° lineup time %:',int(np.nansum(bepideltas[:,3])/len(a)*100),'ICMEs: ',int(i1*np.nansum(bepideltas[:,3])/len(a)),int(i2*np.nansum(bepideltas[:,3])/len(a)))






plt.plot_date(win2.time,a2,'-')

#order win1.lon,sta1.lon,solo1.lon,bepi1_lon,psp1.lon    
plt.plot_date(win2.time,alldiff2,'-')  
plt.ylim(0,30)


# ## lineups and icme rate figure

# In[36]:


#times for at least 2 spacecraft in lineups



fig=plt.figure(50,figsize=(10,6),dpi=100)


### upper

################ Solar Orbiter
ax1 = plt.subplot(111) 
ax1.plot_date(win1.time,all30*30,'-r',lw=0.1)
ax1.fill_between(win1.time,all30*29,all30*31,label='At least 1 lineup with < 30°$\Delta$lon',color='r',lw=0.1,alpha=0.6)
ax1.fill_between(win1.time,all20*19,all20*21,label='At least 1 lineup with < 20°$\Delta$lon',color='b',lw=0.1,alpha=0.6)
ax1.fill_between(win1.time,all10*9,all10*11,label='At least 1 lineup with < 10°$\Delta$lon',color='y',lw=0.1,alpha=0.6)
ax1.fill_between(win1.time,all5*4,all5*6,label='At least 1 lineup with < 5°$\Delta$lon',color='g',lw=0.1,alpha=0.6)

ax1.fill_between(win2.time,all30_2*29,all30_2*31,color='r',lw=0.1,alpha=0.6)
ax1.fill_between(win2.time,all20_2*19,all20_2*21,color='b',lw=0.1,alpha=0.6)
ax1.fill_between(win2.time,all10_2*9,all10_2*11,color='y',lw=0.1,alpha=0.6)
ax1.fill_between(win2.time,all5_2*4,all5_2*6,color='g',lw=0.1,alpha=0.6)


#ax1.plot_date(win1.time,all20*20,'-b',label='At least 1 lineup with < 20°$\Delta$lon',lw=3)
#ax1.plot_date(win1.time,all10*10,'-y',label='At least 1 lineup with < 10°$\Delta$lon',lw=3)
#ax1.plot_date(win1.time,all5*5,'-g',label='At least 1 lineup with < 5°$\Delta$lon',lw=3)
#ax1.plot_date(win2.time,all30_2*30,'-r')
#ax1.plot_date(win2.time,all20_2*20,'-b')
#ax1.plot_date(win2.time,all10_2*10,'-y')
#ax1.plot_date(win2.time,all5_2*5,'-g')


ax1.set_ylabel('$\Delta$longitude [HEEQ]')

for i in np.arange(0.1,1.2,0.1):
    ax3.plot(solo.time,np.zeros(len(solo.r))+i,linestyle='--',color='k',alpha=0.4,lw=0.8)

ax1.set_ylim(0,45)
ax1.legend(loc=2,fontsize=11)

#vertical line at PSP nominal end
ax1.plot([psp1.time[-1],psp1.time[-1]],[0,50],'-k',lw=1.5)
#bepi orbit insertion
ax1.plot([bepi.time[-1],bepi.time[-1]],[0,50],'--k',lw=1.5)

#plt.annotate('PSP nominal mission end',xy=(psp1.time[-1],2),fontsize=12,ha='right')#arrowprops=dict(arrowstyle="->",
                            #connectionstyle="arc3"))
#plt.annotate('PSP nominal mission end')



#p######################## plot ICR
ax2=ax1.twinx()

#ax2.plot(yearly_mid_times_25,icmes_predict_25/365.24, marker='s',markerfacecolor='coral', alpha=0.9,label='Predicted ICME rate MC20',linestyle='')
ax2.plot(times_25_daily,icmes_predict_25mc_daily*30.42/365.24, color='coral', label='Predicted ICME rate McIntosh+ (2020)',linestyle='-')
ax2.plot(times_25_daily,icmes_predict_25pp_daily*30.42/365.24, color='magenta', label='Predicted ICME rate NOAA/NASA/ISES',linestyle='-')

#error ranges
ax2.fill_between(times_25_daily_icrange_num, fmc_low(times_25_daily_icrange_num)*30.42/365.248, fmc_high(times_25_daily_icrange_num)*30.42/365.24, alpha=0.2,color='coral')
ax2.fill_between(times_25_daily_icrange_num, fpp_low(times_25_daily_icrange_num)*30.42/365.248, fpp_high(times_25_daily_icrange_num)*30.42/365.24, alpha=0.2,color='blueviolet')
ax2.set_ylabel('monthly ICME rate')
plt.legend(loc=1,fontsize=10)

plotstart=datetime.datetime(2021, 7, 1)
plotend=datetime.datetime(2030, 1,30)
plt.xlim(plotstart,plotend)
ax2.set_ylim(0,8)


plt.tight_layout()

plt.savefig('results/icme_rate_cycle_update/cycle25_icme_rate_lineups.png', dpi=300)
plt.savefig('results/icme_rate_cycle_update/cycle25_icme_rate_lineups.pdf', dpi=100)



# ## Figure with longitude for lineups

# In[37]:


sns.set_context("talk")     
#sns.set_style('darkgrid')
sns.set_style("ticks")
fsize=15

fig=plt.figure(50,figsize=(15,13),dpi=100)


### upper

################ Solar Orbiter
ax3 = plt.subplot(311) 
ax3.plot_date(solo.time,solo.r,'-g',lw=2)
ax3.plot_date(bepi_time2,bepi_r2,linestyle='-',color='orange',lw=2,markersize=0)
ax3.plot_date(psp.time,psp.r,'-k',lw=1)

ax3.set_ylabel(' R [AU]')
ax3.set_xlim(starttime,endtime)

for i in np.arange(0.1,1.2,0.1):
    ax3.plot(solo.time,np.zeros(len(solo.r))+i,linestyle='--',color='k',alpha=0.4,lw=0.8)

ax3.set_ylim(0,1.2)

#p######################## plot ICR
ax4=ax3.twinx()

#ax2.plot(yearly_mid_times_25,icmes_predict_25/365.24, marker='s',markerfacecolor='coral', alpha=0.9,label='Predicted ICME rate MC20',linestyle='')
ax4.plot(times_25_daily,icmes_predict_25mc_daily*30.42/365.24, color='coral', label='Predicted ICME rate MC20',linestyle='-')
ax4.plot(times_25_daily,icmes_predict_25pp_daily*30.42/365.24, color='magenta', label='Predicted ICME rate PP19',linestyle='-')

#error ranges
ax4.fill_between(times_25_daily_icrange_num, fmc_low(times_25_daily_icrange_num)*30.42/365.248, fmc_high(times_25_daily_icrange_num)*30.42/365.24, alpha=0.2,color='coral')
ax4.fill_between(times_25_daily_icrange_num, fpp_low(times_25_daily_icrange_num)*30.42/365.248, fpp_high(times_25_daily_icrange_num)*30.42/365.24, alpha=0.2,color='blueviolet')
ax4.set_ylabel('monthly ICME rate')
plt.legend(loc=1,fontsize=12)

plotstart=datetime.datetime(2020, 7, 1)
plotend=datetime.datetime(2030, 1,10)
ax4.set_xlim(plotstart,plotend)


################# lower


ax1 = plt.subplot(312) 
#plot R[AU]
ax1.plot_date(bepi_time2,np.rad2deg(bepi_lon2),linestyle='-',color='orange',lw=1,markersize=0)
ax1.plot_date(psp.time,np.rad2deg(psp.lon),'-k',lw=1)
ax1.plot_date(solo.time,np.rad2deg(solo.lon),'-g',lw=1)
ax1.plot_date(sta.time,np.rad2deg(sta.lon),'-r',lw=1)


ax1.plot(solo.time,np.zeros(len(solo.r)),linestyle='-',color='b',alpha=0.8,lw=2)

ax1.set_ylabel('HEEQ longitude [°]')
ax1.set_xlim(starttime,endtime)
ax1.set_ylim(-180,180)


##############

plotstart=datetime.datetime(2020, 7, 1)
plotend=datetime.datetime(2030, 1,10)
plt.xlim(plotstart,plotend)


ax6 = plt.subplot(313) 
ax6.plot_date(bepi_time2,np.rad2deg(bepi_lat2),linestyle='-',color='orange',lw=1,markersize=0)
ax6.plot_date(psp.time,np.rad2deg(psp.lat),'-k',lw=1)
ax6.plot_date(solo.time,np.rad2deg(solo.lat),'-g',lw=2)
ax6.plot_date(sta.time,np.rad2deg(sta.lat),'-r',lw=1)


#Earth
ax6.plot(solo.time,np.zeros(len(solo.r)),linestyle='-',color='b',alpha=0.8,lw=2)

ax6.set_ylabel('HEEQ latitude [°]')
ax6.set_xlim(starttime,endtime)
ax6.set_ylim(-35,35)


#lineups

#ax6.plot_date(win1.time,all20*22,'-g')
#ax6.plot_date(win1.time,all10*20,'-y')
#ax6.plot_date(win1.time,all5*18,'-b')



plotstart=datetime.datetime(2020, 7, 1)
plotend=datetime.datetime(2030, 1,10)
plt.xlim(plotstart,plotend)

plt.grid('on')


ylevel=0.91
#plt.annotate('(a)',[0.01,0.965],xycoords='figure fraction',weight='bold')
#plt.annotate('(b)',[0.01,0.475],xycoords='figure fraction',weight='bold')
plt.annotate('Solar Orbiter',[0.09,ylevel],xycoords='figure fraction',weight='bold',color='green')
plt.annotate('Bepi Colombo',[0.22,ylevel],xycoords='figure fraction',weight='bold',color='orange')
plt.annotate('Parker Solar Probe',[0.36,ylevel],xycoords='figure fraction',weight='bold',color='black')

plt.annotate('STEREO-Ahead',[0.56,ylevel],xycoords='figure fraction',weight='bold',color='red')
plt.annotate('Earth',[0.73,ylevel],xycoords='figure fraction',weight='bold',color='blue')



#plt.tight_layout()

plt.savefig('results/icme_rate_cycle_update/cycle25_icme_rate_psp_orbiter_bepi_one_panel.png', dpi=100)
plt.savefig('results/icme_rate_cycle_update/cycle25_icme_rate_psp_orbiter_bepi_one_panel.pdf', dpi=100)


# In[38]:


#sns.set_style('darkgrid')
sns.set_style("ticks")
fsize=15

fig=plt.figure(50,figsize=(14,6),dpi=100)


### upper

################ Solar Orbiter
ax3 = plt.subplot(111) 
ax3.plot_date(solo.time,solo.r,'-g',lw=2)
ax3.plot_date(bepi_time2,bepi_r2,linestyle='-',color='orange',lw=2,markersize=0)
ax3.plot_date(psp.time,psp.r,'-k',lw=1)

ax3.set_ylabel('Heliocentric distance R [AU]')
ax3.set_xlim(starttime,endtime)

for i in np.arange(0.1,1.2,0.1):
    ax3.plot(solo.time,np.zeros(len(solo.r))+i,linestyle='--',color='k',alpha=0.4,lw=0.8)

ax3.set_ylim(0,1.2)

#p######################## plot ICR
ax4=ax3.twinx()

#ax2.plot(yearly_mid_times_25,icmes_predict_25/365.24, marker='s',markerfacecolor='coral', alpha=0.9,label='Predicted ICME rate MC20',linestyle='')
ax4.plot(times_25_daily,icmes_predict_25mc_daily*30.42/365.24, color='coral', label='Predicted ICME rate MC20',linestyle='-')
ax4.plot(times_25_daily,icmes_predict_25pp_daily*30.42/365.24, color='magenta', label='Predicted ICME rate PP19',linestyle='-')

#error ranges
ax4.fill_between(times_25_daily_icrange_num, fmc_low(times_25_daily_icrange_num)*30.42/365.248, fmc_high(times_25_daily_icrange_num)*30.42/365.24, alpha=0.2,color='coral')
ax4.fill_between(times_25_daily_icrange_num, fpp_low(times_25_daily_icrange_num)*30.42/365.248, fpp_high(times_25_daily_icrange_num)*30.42/365.24, alpha=0.2,color='blueviolet')
ax4.set_ylabel('monthly ICME rate')
plt.legend(loc=1,fontsize=15)

plotstart=datetime.datetime(2020, 7, 1)
plotend=datetime.datetime(2029, 8,30)
plt.xlim(plotstart,plotend)

plt.annotate('Parker Solar Probe',[0.09,0.885],xycoords='figure fraction',weight='bold',color='black',fontsize=14)
plt.annotate('Solar Orbiter',[0.25,0.885],xycoords='figure fraction',weight='bold',color='green',fontsize=14)
plt.annotate('BepiColombo',[0.37,0.885],xycoords='figure fraction',weight='bold',color='orange',fontsize=14)

plt.tight_layout()

plt.savefig('results/icme_rate_cycle_update/cycle25_icme_rate_psp_orbiter_bepi_b1.png', dpi=100)
plt.savefig('results/icme_rate_cycle_update/cycle25_icme_rate_psp_orbiter_bepi_b1.pdf', dpi=100)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





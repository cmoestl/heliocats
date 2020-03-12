'''
icmecat_maker.py

makes the ICMECATv2.0

Author: C. Moestl, IWF Graz, Austria
twitter @chrisoutofspace, https://github.com/cmoestl/heliocats
last update March 2020

python > 3.7, install a conda environment to run this code, see https://github.com/cmoestl/heliocats

current status:
work in progress


to do:

- despike sta stb wind all
- go through all ICMEs and extract data
- (new B and V for STA, Wind and PSP converted to SCEQ components, plasma correct for new PSP, wind, sta)


MIT LICENSE
Copyright 2020, Christian Moestl 
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, 
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
permit persons to whom the Software is furnished to do so, subject to the following 
conditions:
The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''



from scipy import stats
import scipy.io
from matplotlib import cm
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import  DateFormatter
import numpy as np
import astropy.constants as const
import time
import pickle
import seaborn as sns
import os
import urllib
import json
import importlib
import heliopy.data.spice as spicedata
import heliopy.spice as spice
import astropy
import pandas as pd
import openpyxl
import heliosat
from sunpy.time import parse_time
import datetime
import seaborn as sns
import copy
from numba import njit



from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging

#where the final data are located
#data_path='/nas/helio/data/insitu_python/'

data_path='data/'


###################################### ICMECAT operations ################################


def load_helcats_icmecat_master_from_excel(file):

    print('load HELCATS ICMECAT from file:', file)
    ic=pd.read_excel(file)

    #convert times to datetime objects
    for i in np.arange(0,ic.shape[0]):    
    
        a=str(ic.icme_start_time[i]).strip() #remove leading and ending blank spaces if any
        ic.icme_start_time.loc[i]=parse_time(a).datetime

        a=str(ic.mo_start_time[i]).strip() #remove leading and ending blank spaces if any
        ic.mo_start_time.loc[i]=parse_time(a).datetime 

        a=str(ic.mo_end_time[i]).strip() #remove leading and ending blank spaces if any
        ic.mo_end_time.loc[i]=parse_time(a).datetime 

        
        a=str(ic.icme_end_time[i]).strip() #remove leading and ending blank spaces if any
        if a!= '9999-99-99T99:99Z':
            ic.icme_end_time.loc[i]=parse_time(a).datetime 
        else: ic.icme_end_time.loc[i]=np.nan

    return ic


def dynamic_pressure(density, speed):
   '''
   make dynamic pressure from density and speed
   assume pdyn is only due to protons
   pdyn=np.zeros(len([density])) #in nano Pascals
   '''
   proton_mass=1.6726219*1e-27  #kg
   pdyn=np.multiply(np.square(speed*1e3),density)*1e6*proton_mass*1e9  #in nanoPascal
   return pdyn



def get_cat_parameters(sc, sci, ic,name,sctime_num):
    '''
    get parameters
    sc - spacecraft data recarray
    sci - indices for this spacecraft in icmecat
    ic - icmecat pandas dataframe
    '''
    fileind='icmecat/ICMECAT_indices_'+name+'.p'

    #extract indices of ICMEs in the respective data (time consuming)
    
    make_indices=0

    if make_indices > 0:
        #### get all ICMECAT times for this spacecraft as datenum
        sc_icme_start=mdates.date2num(ic.icme_start_time[sci])   
        sc_mo_start=mdates.date2num(ic.mo_start_time[sci])
        sc_mo_end=mdates.date2num(ic.mo_end_time[sci])
        if name=='Wind': sc_icme_end=mdates.date2num(ic.icme_end_time[sci])

        ### arrays containing the indices of where the ICMEs are in the data
        icme_start_ind=np.zeros(len(sci),dtype=int) 
        mo_start_ind=np.zeros(len(sci),dtype=int)
        mo_end_ind=np.zeros(len(sci),dtype=int)
        if name=='Wind': icme_end_ind=np.zeros(len(sci),dtype=int)

        #this takes some time 
        #sctime_num=mdates.date2num(sc.time)   
        
        #get indices in data for each ICMECAT
        for i in np.arange(len(sci))-1:
         
            icme_start_ind[i]=np.where(sctime_num >sc_icme_start[i])[0][0]-1 
            print(icme_start_ind[i])        
            mo_start_ind[i]=np.where(sctime_num >sc_mo_start[i])[0][0]-1   
            mo_end_ind[i]=np.where(sctime_num >sc_mo_end[i])[0][0]-1 
            if name=='Wind': icme_end_ind[i]=np.where(sctime_num >sc_icme_end[i])[0][0]-1 

        if name=='Wind': 
            pickle.dump([icme_start_ind, mo_start_ind,mo_end_ind,icme_end_ind], open(fileind, 'wb'))     
        else:
            pickle.dump([icme_start_ind, mo_start_ind,mo_end_ind], open(fileind, 'wb'))     

    if name=='Wind': 
       [icme_start_ind, mo_start_ind,mo_end_ind,icme_end_ind]=pickle.load(open(fileind, 'rb'))           
    else: 
       [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))           

    ###### get parameters

    #ICME B_max
    for i in np.arange(len(sci))-1:
        if name=='Wind': 
            ic.icme_bmax.loc[sci[i]]=np.round(np.nanmax(sc.bt[icme_start_ind[i]:icme_end_ind[i]]),1)
        else:
            ic.icme_bmax.loc[sci[i]]=np.round(np.nanmax(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)

    #ICME B_mean
    for i in np.arange(len(sci))-1:
        if name=='Wind': 
            ic.icme_bmean.loc[sci[i]]=np.round(np.nanmean(sc.bt[icme_start_ind[i]:icme_end_ind[i]]),1)
        else:
            ic.icme_bmean.loc[sci[i]]=np.round(np.nanmean(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)

    #MO B_max
    for i in np.arange(len(sci))-1:
        ic.mo_bmax.loc[sci[i]]=np.round(np.nanmax(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO B_mean
    for i in np.arange(len(sci))-1:
        ic.mo_bmean.loc[sci[i]]=np.round(np.nanmean(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO B_std
    for i in np.arange(len(sci))-1:
        ic.mo_bstd.loc[sci[i]]=np.round(np.nanstd(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)


    #MO Bz_mean
    for i in np.arange(len(sci))-1:
        ic.mo_bzmean.loc[sci[i]]=np.round(np.nanmean(sc.bz[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO Bz_min
    for i in np.arange(len(sci))-1:
        ic.mo_bzmin.loc[sci[i]]=np.round(np.nanmin(sc.bz[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO heliodistance
    for i in np.arange(len(sci))-1:
        ic.mo_sc_heliodistance.loc[sci[i]]=np.round(np.nanmean(sc.r[mo_start_ind[i]:mo_end_ind[i]]),4)

    #MO longitude
    for i in np.arange(len(sci))-1:
        ic.sc_long_heeq.loc[sci[i]]=np.round(np.nanmean(sc.lon[mo_start_ind[i]:mo_end_ind[i]]),2)

    #MO latitude
    for i in np.arange(len(sci))-1:
        ic.sc_lat_heeq.loc[sci[i]]=np.round(np.nanmean(sc.lat[mo_start_ind[i]:mo_end_ind[i]]),2)


    # ICME duration
    sci_istart=mdates.date2num(ic.icme_start_time[sci])   
    if name=='Wind': 
        sci_iend=mdates.date2num(ic.icme_end_time[sci])   
    else:
        sci_iend=mdates.date2num(ic.mo_end_time[sci])   
    ic.icme_duration.loc[sci]=np.round((sci_iend-sci_istart)*24,2)

    # sheath duration
    #sci_istart=mdates.date2num(ic.icme_start_time[sci])   
    #sci_iend=mdates.date2num(ic.mo_start_time[sci])   
    #ic.icme_duration.loc[sci]=np.round((sci_iend-sci_istart)*24,2)


    # MO duration
    sci_istart=mdates.date2num(ic.mo_start_time[sci])   
    sci_iend=mdates.date2num(ic.mo_end_time[sci])   
    ic.mo_duration.loc[sci]=np.round((sci_iend-sci_istart)*24,2)

    
    return ic










    
##########################################################################################
######################################## MAIN PROGRAM ####################################
##########################################################################################

##################################### (1) load new data with HelioSat and heliocats.data


    
load_data=1

if load_data >0:

    print('load new Wind, STEREO-A, MAVEN, and ParkerProbe data')


    #MAVEN
    #filemav='maven_2014_2018.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )

    #filemav='maven_2014_2018_removed.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )
    
    filemav='maven_2014_2018_removed_smoothed.p'
    [mav,hmav]=pickle.load(open(data_path+filemav, 'rb' ) )
    

    #Wind
    filewin="wind_2018_2020.p" 
    #for updating data
    #start=datetime.datetime(2018, 1, 1)
    #end=datetime.datetime.utcnow()
    #hd.save_wind_data(data_path,filewin,start,end)
    [win2,hwin2]=pickle.load(open(data_path+filewin, "rb" ) )  


    #STEREO-A    
    filesta2='sta_2018_2019_beacon.p'
    #start=datetime.datetime(2018, 1, 1)
    #end=datetime.datetime(2019, 12, 31)
    #hd.save_stereoa_beacon_data(data_path,filesta,start,end)
   
    sta2=pickle.load(open(data_path+filesta2, "rb" ) )  

    #Parker Solar Probe
    filepsp='psp_2018_2019.p'
    [psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) )  


    # ADD BepiColombo  
    
    
    # ADD Solar Orbiter



    ##################################### (2) load HELCATS DATACAT

    [vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly]=hd.load_helcats_datacat(data_path+'helcats_all_data_removed.p') 



################################ (3) make ICMECAT 

print('data loaded')


ic=load_helcats_icmecat_master_from_excel('icmecat/HELCATS_ICMECAT_v20_master.xlsx')

#get indices for all spacecraft
wini=np.where(ic.sc_insitu == 'Wind')[:][0] 
vexi=np.where(ic.sc_insitu == 'VEX')[:][0]  
mesi=np.where(ic.sc_insitu == 'MESSENGER')[:][0]   
stai=np.where(ic.sc_insitu == 'STEREO-A')[:][0]    
stbi=np.where(ic.sc_insitu == 'STEREO-B')[:][0]    
mavi=np.where(ic.sc_insitu == 'MAVEN')[:][0]    
ulyi=np.where(ic.sc_insitu == 'ULYSSES')[:][0]    



filetimes='icmecat/ICMECAT_numtimes.p'
'''
#save times as mdates
wintime_num=mdates.date2num(win.time) 
statime_num=mdates.date2num(sta.time) 
stbtime_num=mdates.date2num(stb.time) 
mestime_num=mdates.date2num(mes.time) 
vextime_num=mdates.date2num(vex.time) 
ulytime_num=mdates.date2num(uly.time) 
mavtime_num=mdates.date2num(mav.time) 
pickle.dump([wintime_num,statime_num,stbtime_num,mestime_num,vextime_num,ulytime_num,mavtime_num], open(filetimes, 'wb'))     
print('times as num saved')
'''


[wintime_num,statime_num,stbtime_num,mestime_num,vextime_num,ulytime_num,mavtime_num]=pickle.load(open(filetimes,'rb')) 




#pspi=np.where(ic.sc_insitu == 'ParkerSolarProbe')[:][0]    


#get parameters for all spacecraft one after another
ic=get_cat_parameters(win,wini,ic,'Wind',wintime_num)
ic=get_cat_parameters(sta,stai,ic,'STEREO-A',statime_num)
ic=get_cat_parameters(stb,stbi,ic,'STEREO_B',stbtime_num)
ic=get_cat_parameters(mes,mesi,ic,'MESSENGER',mestime_num)
ic=get_cat_parameters(vex,vexi,ic,'VEX',vextime_num)
ic=get_cat_parameters(uly,ulyi,ic,'ULYSSES',ulytime_num)
ic=get_cat_parameters(mav,mavi,ic,'MAVEN',mavtime_num)



################################ (4) save ICMECAT #################################

ic3=copy.deepcopy(ic)  

#pickle, excel, json, csv, txt (cdf? votable?)

#save as pickle with datetime
file='icmecat/HELCATS_ICMECAT_v20.p'
pickle.dump(ic, open(file, 'wb'))



#use date and time format from master table
ic2=pd.read_excel('icmecat/HELCATS_ICMECAT_v20_master.xlsx')
ic3.icme_start_time=ic2.icme_start_time
ic3.mo_start_time=ic2.mo_start_time
ic3.mo_end_time=ic2.mo_end_time
ic3.icme_end_time=ic2.icme_end_time
del(ic2)

#save as Excel
file='icmecat/HELCATS_ICMECAT_v20.xlsx'
ic3.to_excel(file,sheet_name='ICMECATv2.0')

#save as json
file='icmecat/HELCATS_ICMECAT_v20.json'
ic3.to_json(file)

#save as csv
file='icmecat/HELCATS_ICMECAT_v20.csv'
ic3.to_csv(file)


#save as hdf needs pip install tables
#file='icmecat/HELCATS_ICMECAT_v20.hdf'
#ic.to_hdf(file,key='icmecat')


#save as .mat does not work yet
#ile='icmecat/HELCATS_ICMECAT_v20.mat'
#icdict=ic.to_dict()
#scipy.io.savemat(file,ic.values)


#save as txt
file='icmecat/HELCATS_ICMECAT_v20.txt'
np.savetxt(file, ic3.values.astype(str), fmt='%s' )

print('ICMECAT saved as '+file)



sys.exit()


#icl=pickle.load(open(file, 'rb' ) )

######################################################################################
################################### END MAIN #########################################
######################################################################################




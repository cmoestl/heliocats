#stats.py
#statistics stuff for heliocats
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy
from sunpy.time import parse_time
import copy
import matplotlib.dates as mdates
import matplotlib
import seaborn as sns
import datetime
import urllib
import json
import os
import pdb
import scipy.io
import pickle
import sys
import cdflib
import matplotlib.pyplot as plt
import heliosat

from input import *



###################################### ICMECAT operations ################################




def load_helcats_icmecat_master_from_excel(file):
    ''' convert excel master file to pandas dataframe and convert times
        to datetime objects
    '''

    print('load HELCATS ICMECAT from file:', file)
    ic=pd.read_excel(file)

    #convert all times to datetime objects
    for i in np.arange(0,ic.shape[0]):    
    
        #remove leading and ending blank spaces if any and write datetime object into dataframe
        ic.at[i,'icme_start_time']= parse_time(str(ic.icme_start_time[i]).strip()).datetime 
        ic.at[i,'mo_start_time']=parse_time(str(ic.mo_start_time[i]).strip()).datetime
        ic.at[i,'mo_end_time']=parse_time(str(ic.mo_end_time[i]).strip()).datetime
        
        #icme end time only at Wind
        a=str(ic.icme_end_time[i]).strip() #remove leading and ending blank spaces if any
        if a!= '9999-99-99T99:99Z':
            ic.at[i,'icme_end_time']=parse_time(a).datetime
        else: ic.at[i,'icme_end_time']=np.nan

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



def get_cat_parameters(sc, sci, ic, name):
    '''
    get parameters
    sc - spacecraft data recarray
    sci - indices for this spacecraft in icmecat
    ic - icmecat pandas dataframe
    '''
    fileind='data/indices_icmecat/ICMECAT_indices_'+name+'.p'

    #### extract indices of ICMEs in the respective data (time consuming, so do it once)
    
    if os.path.isfile(fileind) == False:
    
        print('extract indices of ICMEs in '+ name+ ' data')
        #### get all ICMECAT times for this spacecraft as datenum
        sc_icme_start=ic.icme_start_time[sci]
        sc_mo_start=ic.mo_start_time[sci]
        sc_mo_end=ic.mo_end_time[sci]
        if name=='Wind': sc_icme_end=ic.icme_end_time[sci]
                
        ### arrays containing the indices of where the ICMEs are in the data
        icme_start_ind=np.zeros(len(sci),dtype=int) 
        mo_start_ind=np.zeros(len(sci),dtype=int)
        mo_end_ind=np.zeros(len(sci),dtype=int)
        if name=='Wind': icme_end_ind=np.zeros(len(sci),dtype=int)

        #this takes some time, get indices in data for each ICMECAT
        for i in np.arange(sci[0],sci[-1]+1):
        
            print(i-sci[0])
        
            icme_start_ind[i-sci[0]]=np.where(sc.time  > sc_icme_start[i])[0][0]-1 
            #print(icme_start_ind[i])        
            mo_start_ind[i-sci[0]]=np.where(sc.time > sc_mo_start[i])[0][0]-1   
            mo_end_ind[i-sci[0]]=np.where(sc.time   > sc_mo_end[i])[0][0]-1 
            if name=='Wind': icme_end_ind[i-sci[0]]=np.where(sc.time > sc_icme_end[i])[0][0]-1

        if name=='Wind': 
            pickle.dump([icme_start_ind, mo_start_ind,mo_end_ind,icme_end_ind], open(fileind, 'wb'))     
            print(name+' indices done.')
        else:
            pickle.dump([icme_start_ind, mo_start_ind,mo_end_ind], open(fileind, 'wb'))
            
                
    ##### load files               
    if name=='Wind': 
       [icme_start_ind, mo_start_ind,mo_end_ind,icme_end_ind]=pickle.load(open(fileind, 'rb'))           
    else: 
       [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))           

    ###### get parameters
    #ICME B_max
    for i in np.arange(len(sci))-1:
        if name=='Wind': 
            ic.at[sci[i],'icme_bmax']=np.round(np.nanmax(sc.bt[icme_start_ind[i]:icme_end_ind[i]]),1)
        else:
            ic.at[sci[i],'icme_bmax']=np.round(np.nanmax(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)

    #ICME B_mean
    for i in np.arange(len(sci))-1:
        if name=='Wind': 
            ic.at[sci[i],'icme_bmean']=np.round(np.nanmean(sc.bt[icme_start_ind[i]:icme_end_ind[i]]),1)
        else:
            ic.at[sci[i],'icme_bmean']=np.round(np.nanmean(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)

    #MO B_max
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_bmax']=np.round(np.nanmax(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO B_mean
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_bmean']=np.round(np.nanmean(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO B_std
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_bstd']=np.round(np.nanstd(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO Bz_mean
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_bzmean']=np.round(np.nanmean(sc.bz[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO Bz_min
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_bzmin']=np.round(np.nanmin(sc.bz[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO heliodistance
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_sc_heliodistance']=np.round(np.nanmean(sc.r[mo_start_ind[i]:mo_end_ind[i]]),4)

    #MO longitude
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'sc_long_heeq']=np.round(np.nanmean(sc.lon[mo_start_ind[i]:mo_end_ind[i]]),2)

    #MO latitude
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'sc_lat_heeq']=np.round(np.nanmean(sc.lat[mo_start_ind[i]:mo_end_ind[i]]),2)

    # ICME duration
    sci_istart=mdates.date2num(ic.icme_start_time[sci])   
    #if name=='Wind': 
    #    sci_iend=mdates.date2num(ic.icme_end_time[sci])   
    #else:
    sci_iend=mdates.date2num(ic.mo_end_time[sci])   
    ic.at[sci,'icme_duration']=np.round((sci_iend-sci_istart)*24,2)

    # sheath duration
    #sci_istart=mdates.date2num(ic.icme_start_time[sci])   
    #sci_iend=mdates.date2num(ic.mo_start_time[sci])   
    #ic.icme_duration.loc[sci]=np.round((sci_iend-sci_istart)*24,2)

    # MO duration
    sci_istart=mdates.date2num(ic.mo_start_time[sci])   
    sci_iend=mdates.date2num(ic.mo_end_time[sci])   
    ic.at[sci,'mo_duration']=np.round((sci_iend-sci_istart)*24,2)

    
    return ic







'''
def get_cat_parameters_old(sc, sci, ic,name,sctime_num):

    #get parameters
    #sc - spacecraft data recarray
    #sci - indices for this spacecraft in icmecat
    #ic - icmecat pandas dataframe
    fileind='icmecat/ICMECAT_indices_'+name+'.p'

    #extract indices of ICMEs in the respective data (time consuming)
    
    
    #check whether file is there*****************************
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
            print('indices done')
        else:
            pickle.dump([icme_start_ind, mo_start_ind,mo_end_ind], open(fileind, 'wb'))
            
                
    #load files               

    if name=='Wind': 
       [icme_start_ind, mo_start_ind,mo_end_ind,icme_end_ind]=pickle.load(open(fileind, 'rb'))           
    else: 
       [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))           

    ###### get parameters

    #ICME B_max
    for i in np.arange(len(sci))-1:
        if name=='Wind': 
            ic.at[sci[i],'icme_bmax']=np.round(np.nanmax(sc.bt[icme_start_ind[i]:icme_end_ind[i]]),1)
        else:
            ic.at[sci[i],'icme_bmax']=np.round(np.nanmax(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)

    #ICME B_mean
    for i in np.arange(len(sci))-1:
        if name=='Wind': 
            ic.at[sci[i],'icme_bmean']=np.round(np.nanmean(sc.bt[icme_start_ind[i]:icme_end_ind[i]]),1)
        else:
            ic.at[sci[i],'icme_bmean']=np.round(np.nanmean(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)

    #MO B_max
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_bmax']=np.round(np.nanmax(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO B_mean
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_bmean']=np.round(np.nanmean(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO B_std
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_bstd']=np.round(np.nanstd(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO Bz_mean
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_bzmean']=np.round(np.nanmean(sc.bz[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO Bz_min
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_bzmin']=np.round(np.nanmin(sc.bz[mo_start_ind[i]:mo_end_ind[i]]),1)

    #MO heliodistance
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_sc_heliodistance']=np.round(np.nanmean(sc.r[mo_start_ind[i]:mo_end_ind[i]]),4)

    #MO longitude
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'sc_long_heeq']=np.round(np.nanmean(sc.lon[mo_start_ind[i]:mo_end_ind[i]]),2)

    #MO latitude
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'sc_lat_heeq']=np.round(np.nanmean(sc.lat[mo_start_ind[i]:mo_end_ind[i]]),2)

    # ICME duration
    sci_istart=mdates.date2num(ic.icme_start_time[sci])   
    #if name=='Wind': 
    #    sci_iend=mdates.date2num(ic.icme_end_time[sci])   
    #else:
    sci_iend=mdates.date2num(ic.mo_end_time[sci])   
    ic.at[sci,'icme_duration']=np.round((sci_iend-sci_istart)*24,2)

    # sheath duration
    #sci_istart=mdates.date2num(ic.icme_start_time[sci])   
    #sci_iend=mdates.date2num(ic.mo_start_time[sci])   
    #ic.icme_duration.loc[sci]=np.round((sci_iend-sci_istart)*24,2)

    # MO duration
    sci_istart=mdates.date2num(ic.mo_start_time[sci])   
    sci_iend=mdates.date2num(ic.mo_end_time[sci])   
    ic.at[sci,'mo_duration']=np.round((sci_iend-sci_istart)*24,2)

    
    return ic



'''




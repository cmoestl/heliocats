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

from input import data_path



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
       
   
    return ic


def make_icmecat_header(ic):
    ''' todo
    '''
  
    
    header='header'
   
    return header



def pdyn(density, speed):
    '''
    make dynamic pressure from density []# ccm-3] and speed [km/s]
    assume pdyn is only due to protons
    pdyn=np.zeros(len([density])) #in nano Pascals
    '''
    proton_mass=1.6726219*1e-27  #kg
    pdyn=np.multiply(np.square(speed*1e3),density)*1e6*proton_mass*1e9  #in nanoPascal

    return pdyn
    
    
def load_pickle(file):    

    ic=pickle.load( open(file, 'rb'))    
    
    return ic


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

    
        ### arrays containing the indices of where the ICMEs are in the data
        icme_start_ind=np.zeros(len(sci),dtype=int) 
        mo_start_ind=np.zeros(len(sci),dtype=int)
        mo_end_ind=np.zeros(len(sci),dtype=int)
   
        #this takes some time, get indices in data for each ICMECAT
        for i in np.arange(sci[0],sci[-1]+1):
        
            print(i-sci[0])

            icme_start_ind[i-sci[0]]=np.where(sc.time  > sc_icme_start[i])[0][0]-1 
            #print(icme_start_ind[i])        
            mo_start_ind[i-sci[0]]=np.where(sc.time > sc_mo_start[i])[0][0]-1   
            mo_end_ind[i-sci[0]]=np.where(sc.time   > sc_mo_end[i])[0][0]-1 

        pickle.dump([icme_start_ind, mo_start_ind,mo_end_ind], open(fileind, 'wb'))
    ############################################            
                
    [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))           
    
    
    #plasma available?
    if name=='Wind': plasma=True
    if name=='STEREO-A': plasma=True
    if name=='STEREO-B': plasma=True
    if name=='ULYSSES': plasma=True
    if name=='MAVEN': plasma=True
    if name=='PSP': plasma=True
    if name=='VEX': plasma=False
    if name=='MESSENGER': plasma=False

    print('Get parameters for ',name)
    

    ####### position

    #MO heliodistance
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_sc_heliodistance']=np.round(sc.r[mo_start_ind[i]],4)

        #MO longitude
        ic.at[sci[i],'mo_sc_long_heeq']=np.round(sc.lon[mo_start_ind[i]],2)

        #MO latitude
        ic.at[sci[i],'mo_sc_lat_heeq']=np.round(sc.lat[mo_start_ind[i]],2)


    ############ ICME    
    # ICME duration
    sci_istart=mdates.date2num(ic.icme_start_time[sci])   
    sci_iend=mdates.date2num(ic.mo_end_time[sci])   
    ic.at[sci,'icme_duration']=np.round((sci_iend-sci_istart)*24,2)

    for i in np.arange(len(sci))-1:

        #ICME B_max
        ic.at[sci[i],'icme_bmax']=np.round(np.nanmax(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)

        #ICME B_mean
        ic.at[sci[i],'icme_bmean']=np.round(np.nanmean(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)

        #icme_bstd
        ic.at[sci[i],'icme_bstd']=np.round(np.nanstd(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)
        
    if plasma==True:        
        #ICME speed_mean and std
        for i in np.arange(len(sci))-1:
            ic.at[sci[i],'icme_speed_mean']=np.round(np.nanmean(sc.vt[icme_start_ind[i]:mo_end_ind[i]]),1)
            ic.at[sci[i],'icme_speed_std']=np.round(np.nanstd(sc.vt[icme_start_ind[i]:mo_end_ind[i]]),1)
    else: #set nan    
        for i in np.arange(len(sci))-1:
            ic.at[sci[i],'icme_speed_mean']=np.nan
            ic.at[sci[i],'icme_speed_std']=np.nan

        
    ########### MO
    # MO duration
    sci_istart=mdates.date2num(ic.mo_start_time[sci])   
    sci_iend=mdates.date2num(ic.mo_end_time[sci])   
    ic.at[sci,'mo_duration']=np.round((sci_iend-sci_istart)*24,2)        

    for i in np.arange(len(sci))-1:
    
        #MO B_max
        ic.at[sci[i],'mo_bmax']=np.round(np.nanmax(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)
    
        #MO B_mean
        ic.at[sci[i],'mo_bmean']=np.round(np.nanmean(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)
    
        #MO B_std
        ic.at[sci[i],'mo_bstd']=np.round(np.nanstd(sc.bt[mo_start_ind[i]:mo_end_ind[i]]),1)

        #MO Bz_mean
        ic.at[sci[i],'mo_bzmean']=np.round(np.nanmean(sc.bz[mo_start_ind[i]:mo_end_ind[i]]),1)

        #MO Bz_min
        ic.at[sci[i],'mo_bzmin']=np.round(np.nanmin(sc.bz[mo_start_ind[i]:mo_end_ind[i]]),1)

         #MO Bz_std
        ic.at[sci[i],'mo_bzstd']=np.round(np.nanstd(sc.bz[mo_start_ind[i]:mo_end_ind[i]]),1)

        #MO By_mean
        ic.at[sci[i],'mo_bymean']=np.round(np.nanmean(sc.by[mo_start_ind[i]:mo_end_ind[i]]),1)

        #MO By_std
        ic.at[sci[i],'mo_bystd']=np.round(np.nanstd(sc.by[mo_start_ind[i]:mo_end_ind[i]]),1)

    
    if plasma==True:   
         
        for i in np.arange(len(sci))-1:
        
            #mo speed_mean and std
            ic.at[sci[i],'mo_speed_mean']=np.round(np.nanmean(sc.vt[mo_start_ind[i]:mo_end_ind[i]]),1)
            ic.at[sci[i],'mo_speed_std']=np.round(np.nanstd(sc.vt[mo_start_ind[i]:mo_end_ind[i]]),1)
            
            ic.at[sci[i],'mo_expansion_speed']=np.round( (sc.vt[mo_start_ind[i]]-sc.vt[mo_end_ind[i]])/2 ,1 )

            ic.at[sci[i],'mo_density_mean']=np.round(np.nanmean(sc.np[mo_start_ind[i]:mo_end_ind[i]]),1)
            ic.at[sci[i],'mo_density_std']=np.round(np.nanstd(sc.np[mo_start_ind[i]:mo_end_ind[i]]),1)

            ic.at[sci[i],'mo_temperature_mean']=np.round(np.nanmean(sc.tp[mo_start_ind[i]:mo_end_ind[i]]),1)
            ic.at[sci[i],'mo_temperature_std']=np.round(np.nanstd(sc.tp[mo_start_ind[i]:mo_end_ind[i]]),1)

            pdyn_i=pdyn(sc.np[mo_start_ind[i]:mo_end_ind[i]],sc.vt[mo_start_ind[i]:mo_end_ind[i]])
            
            ic.at[sci[i],'mo_pdyn_mean']=np.round(np.nanmean(pdyn_i),1)
            ic.at[sci[i],'mo_pdyn_std']=np.round(np.nanstd(pdyn_i),1)
            
            
            #icme speed_mean and std
            ic.at[sci[i],'sheath_speed_mean']=np.round(np.nanmean(sc.vt[icme_start_ind[i]:mo_start_ind[i]]),1)
            ic.at[sci[i],'sheath_speed_std']=np.round(np.nanstd(sc.vt[icme_start_ind[i]:mo_start_ind[i]]),1)         

            ic.at[sci[i],'sheath_density_mean']=np.round(np.nanmean(sc.np[icme_start_ind[i]:mo_start_ind[i]]),1)
            ic.at[sci[i],'sheath_density_std']=np.round(np.nanstd(sc.np[icme_start_ind[i]:mo_start_ind[i]]),1)        

            pdyn_i=pdyn(sc.np[icme_start_ind[i]:mo_start_ind[i]],sc.vt[icme_start_ind[i]:mo_start_ind[i]])

            ic.at[sci[i],'sheath_pdyn_mean']=np.round(np.nanmean(pdyn_i),1)
            ic.at[sci[i],'sheath_pdyn_std']=np.round(np.nanstd(pdyn_i),1)

            
            
    else: #set nan    
    
        for i in np.arange(len(sci))-1:
            ic.at[sci[i],'mo_speed_mean']=np.nan
            ic.at[sci[i],'mo_speed_std']=np.nan
            
            ic.at[sci[i],'mo_expansion_speed']=np.nan
    
            ic.at[sci[i],'mo_density_mean']=np.nan
            ic.at[sci[i],'mo_density_std']=np.nan
    
            ic.at[sci[i],'mo_temperature_mean']=np.nan
            ic.at[sci[i],'mo_temperature_std']=np.nan
            
            ic.at[sci[i],'mo_pdyn_mean']=np.nan
            ic.at[sci[i],'mo_pdyn_std']=np.nan

            ic.at[sci[i],'sheath_pdyn_mean']=np.nan
            ic.at[sci[i],'sheath_pdyn_std']=np.nan

    

    
    return ic










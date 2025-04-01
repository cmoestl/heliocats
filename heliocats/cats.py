#cats.py

#catalog creation for heliocats

#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy
import copy
import matplotlib.dates as mdates
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import datetime
import urllib
import json
import os
import pdb
import scipy.io
import pickle
import sys
import importlib
import cdflib


import spiceypy

import sunpy

from sunpy.coordinates import frames, get_horizons_coord, HeliographicStonyhurst
from sunpy.time import parse_time

import astropy
from astropy.constants import au
import astropy.units as u
from astropy.io.votable import parse_single_table

import astroquery

#import astrospice

from heliocats import data as hd
importlib.reload(hd) #reload again while debugging

#define AU in km
AU=au.value/1e3







################################ HI arrival catalog ARRCAT operations ##############################


def load_higeocat_vot(file):
    #read HIGEOCAT from https://www.helcats-fp7.eu/catalogues/wp3_cat.html
    #https://docs.astropy.org/en/stable/io/votable/
   
    table = parse_single_table('data/HCME_WP3_V06.vot')
    higeocat = table.array
    #usage e.g.
    #higeocat['Date']=parse_time(higeocat['Date'][10]).datetime
    #access data
    #a=table.array['HM HEEQ Long'][10]
    
    return higeocat




   
def cart2sphere_emma_rad(x,y,z):
    
    r = np.sqrt(x**2+ y**2 + z**2) /1.495978707E8         
    theta = np.arctan2(z,np.sqrt(x**2+ y**2)) * 360 / 2 / np.pi
    phi = np.arctan2(y,x) * 360 / 2 / np.pi    
    
    theta=np.deg2rad(theta)
    phi=np.deg2rad(phi)
    
    return (r, theta, phi)









############# for spiceypy positions


#JUNO
def furnish_juno(kernel_path):

    #need to load both directories
    juno_kernel_path = kernel_path+'juno'

    kernels = os.listdir(juno_kernel_path)
    for kernel in kernels:
        spiceypy.furnsh(os.path.join(juno_kernel_path, kernel))

    
    generic_path = kernel_path+'generic/'
    generic_kernels = os.listdir(generic_path)
    for kernel in generic_kernels:
        spiceypy.furnsh(os.path.join(generic_path, kernel))





#STEREO-A and STEREO-B
def furnish_stereo(kernel_path,aorb):

    #need to load both directories
    stereo_kernel_path = kernel_path+'stereo'+aorb+'_predicted/'

    stereoa_kernels = os.listdir(stereo_kernel_path)
    for kernel in stereoa_kernels:
        spiceypy.furnsh(os.path.join(stereo_kernel_path, kernel))
    
    stereo_kernel_path = kernel_path+'stereo'+aorb+'/'

    stereoa_kernels = os.listdir(stereo_kernel_path)
    for kernel in stereoa_kernels:
        spiceypy.furnsh(os.path.join(stereo_kernel_path, kernel))
    
    
    
    generic_path = kernel_path+'generic/'
    generic_kernels = os.listdir(generic_path)
    for kernel in generic_kernels:
        spiceypy.furnsh(os.path.join(generic_path, kernel))



def furnish(kernel_path,kernel_file):

    spiceypy.furnsh(os.path.join(kernel_path, kernel_file))
    print(kernel_path, kernel_file)

    generic_path = kernel_path+'generic/'
    generic_kernels = os.listdir(generic_path)
    for kernel in generic_kernels:
        spiceypy.furnsh(os.path.join(generic_path, kernel))
        

def get_pos(t,name):    

    pos = spiceypy.spkpos(name, spiceypy.datetime2et(t), "HEEQ", "NONE", "SUN")[0]
    r, lat, lon = cart2sphere_emma_rad(pos[0],pos[1],pos[2])
    position = t, pos[0], pos[1], pos[2], r, lat, lon
    return position


def get_sc_pos(time, name):
    
    position = get_pos(time, name)
    return position






def get_position(time1,name):
        
    #this works if the respective kernels are loaded beforehand
    
    insitu_exist=True
    
    #for missions, cut off times directly that do not exist:
    if name=='PSP': 
        #exclude if time before launch time
        if parse_time(time1).plot_date < parse_time(datetime.datetime(2018, 8, 13)).plot_date:
            insitu_exist=False 

    if name=='SolarOrbiter': 
        if parse_time(time1).plot_date < parse_time(datetime.datetime(2020, 4, 15)).plot_date:            
            insitu_exist=False 
        
    if name=='BepiColombo': 
        if parse_time(time1).plot_date < parse_time(datetime.datetime(2018, 10, 24)).plot_date:
            insitu_exist=False               
            
    if name=='STEREO-A': 
        if parse_time(time1).plot_date < parse_time(datetime.datetime(2007, 1, 1)).plot_date:
            insitu_exist=False  
                 
    if name=='STEREO-B': 
        if parse_time(time1).plot_date > parse_time(datetime.datetime(2014, 9, 27)).plot_date:
            insitu_exist=False  
            
    if name=='MESSENGER': 
        if parse_time(time1).plot_date > parse_time(datetime.datetime(2015, 4, 30)).plot_date:
            insitu_exist=False                          
                               
    if name=='JUICE': 
        if parse_time(time1).plot_date < parse_time(datetime.datetime(2023, 4, 14)).plot_date:
            insitu_exist=False                  
                
    if name=='VEX': 
        if parse_time(time1).plot_date > parse_time(datetime.datetime(2014, 12, 16)).plot_date:
            insitu_exist=False                   
              
    if name=='JUNO': 
        if parse_time(time1).plot_date < parse_time(datetime.datetime(2011,8, 12)).plot_date:
            insitu_exist=False                   
                
                
    if name=='Ulysses': 
        #cut off ulysses when no decent in situ data is available anymore
        if parse_time(time1).plot_date > parse_time(datetime.datetime(2008, 5, 1)).plot_date:
            insitu_exist=False              

    if insitu_exist == True:


        if name=='PSP':            
            #time, spice code
            coords=get_sc_pos(time1,"PARKER SOLAR PROBE")
            
        if name=='SolarOrbiter':
            coords=get_sc_pos(time1,"SOLAR ORBITER")

        if name=='BepiColombo':
            coords=get_sc_pos(time1,"BEPICOLOMBO MPO")            

        if name=='STEREO-A':
            coords=get_sc_pos(time1,"STEREO AHEAD")
            
        if name=='STEREO-B':
            coords=get_sc_pos(time1,"STEREO BEHIND")
            
        if name=='Mercury':
            coords=get_sc_pos(time1,'MERCURY_BARYCENTER')
            
        if name=='Venus':
            coords=get_sc_pos(time1,'VENUS_BARYCENTER')

        if name=='Mars':
            coords=get_sc_pos(time1,'MARS_BARYCENTER')

        if name=='Earth_L1':
            coords=get_sc_pos(time1,'EARTH_BARYCENTER')
            
        if name=='JUICE':
            coords=get_sc_pos(time1,'JUICE')

        if name=='MESSENGER':
            coords=get_sc_pos(time1,'MESSENGER')

        if name=='Ulysses':
            coords=get_sc_pos(time1,'ULYSSES')
            
        if name=='JUNO':
            coords=get_sc_pos(time1,'JUNO')

        if name=='VEX':
            coords=get_sc_pos(time1,'VENUS EXPRESS')

        
        insitu_r = coords[4]
        insitu_lat =np.rad2deg(coords[5])
        insitu_lon =np.rad2deg(coords[6])

        #correct roughly for L1 position
        if name=='Earth_L1': insitu_r=insitu_r-0.01

        insitu_time=time1


        
    else:
        insitu_time=np.nan
        insitu_r=np.nan
        insitu_lat=np.nan
        insitu_lon=np.nan    

    return [insitu_time,insitu_r,insitu_lat,insitu_lon]






def calculate_arrival(vsse,delta,lamda,rdist,t0_num):
   
    #calculate arrival time after Möstl and Davies 2013 but using ta=t0+Ri/Visse equivalent to ta=t0+Risse/Vsse    
   
    visse=vsse * (  np.cos(np.radians(delta))   \
                  + np.sqrt(  np.sin(np.radians(lamda))**2-np.sin(np.radians(delta))**2 ) ) \
                    /(1+np.sin(np.radians(lamda)) )      
                              
    
    #arrival time: convert AU to km  and seconds to days                
    ta=t0_num+(rdist*AU/visse)/(3600*24)        
                              
    return [mdates.num2date(ta),visse]












def make_arrival_catalog_insitu_ssef30(higeocat,arrcat, target_name, column_list,kernels_path):
    
    #get parameters from HIGEOCAT for arrival catalog

    higeocat_time=parse_time(higeocat['Date']).datetime    #first HI observation
    higeocat_t0=parse_time(higeocat['SSE Launch']).datetime   #backprojected launch time
    higeocat_t0_num=parse_time(higeocat_t0).plot_date
    higeocat_vsse=np.array(higeocat['SSE Speed'])
    higeocat_vsse_err=np.array(higeocat['SSE Speed Err'])
    higeocat_sse_lon=np.array(higeocat['SSE HEEQ Long' ])
    higeocat_sse_lat=np.array(higeocat['SSE HEEQ Lat' ])
    higeocat_id=np.array(higeocat['ID'])
    higeocat_sc=np.array(higeocat['SC'])
    higeocat_pan=np.array(higeocat['PA-N'])
    higeocat_pas=np.array(higeocat['PA-S'])
    higeocat_pafit=np.array(higeocat['PA-fit'])
    higeocat_pacenter=abs((higeocat_pan+higeocat_pas)/2)
    
    if target_name=='PSP':
        kernels_file='psp/spp_nom_20180812_20300101_v042_PostV7.bsp'
        print(kernels_file)
        furnish(kernels_path,kernels_file)    
  
    if target_name=='SolarOrbiter':
        kernels_file='solo/solo_ANC_soc-orbit-stp_20200210-20301120_353_V1_00424_V01.bsp'
        print(kernels_file)
        furnish(kernels_path,kernels_file)            
                
    if target_name=='BepiColombo':        
        kernels_file='bepi/bc_mpo_fcp_00199_20181020_20270407_v02.bsp'
        print(kernels_file)
        furnish(kernels_path,kernels_file)    
        
    if target_name=='STEREO-A':
        furnish_stereo(kernels_path,'a')           
        
    if target_name=='STEREO-B':
        furnish_stereo(kernels_path,'b')        
        
    if target_name=='MESSENGER':
        kernels_file='messenger/msgr_040803_150430_150430_od431sc_2.bsp'
        print(kernels_file)
        furnish(kernels_path,kernels_file)    

    if target_name=='Ulysses':
        kernels_file='ulysses/ulysses_1990_2009_2050.bsp'
        print(kernels_file)
        furnish(kernels_path,kernels_file)    
        
    if target_name=='JUICE':
        kernels_file='JUICE/juice_orbc_000080_230414_310721_v01.bsp'
        print(kernels_file)
        furnish(kernels_path,kernels_file)     
        
    if target_name=='VEX':
        kernels_file='vex/ORVM_T19___________00001.BSP'
        print(kernels_file)
        furnish(kernels_path,kernels_file)            
        
    if target_name=='JUNO':
        furnish_juno(kernels_path)  

        
    #half width for SSEF30
    lamda=30.0

    #ARRCAT with iteration
    arrcat_insitu_list = []

    ############
    #go through all HIGEOCAT CME events and check for hit at insitu, with 4 iterations in total
    for i in np.arange(len(higeocat_time)):

        #get insitu position for launch time t0            
        [insitu_time,insitu_r,insitu_lat,insitu_lon]=get_position(higeocat_t0[i], target_name)            
        delta=abs(higeocat_sse_lon[i]-insitu_lon)
        #print([insitu_time,insitu_r,insitu_lat,insitu_lon])

        if delta < 30:               

            #calculate arrival time 
            #print(delta,lamda,insitu_r)
            [ta,visse]=calculate_arrival(higeocat_vsse[i],delta, lamda, insitu_r,higeocat_t0_num[i])                
            
  
            [insitu_time2,insitu_r2,insitu_lat2,insitu_lon2]=get_position(ta, target_name) 
            #print(insitu_lon-insitu_lon2)               
            delta2=abs(higeocat_sse_lon[i]-insitu_lon2)
            
            if delta2 <30:

                [ta2,visse2]=calculate_arrival(higeocat_vsse[i],delta2, lamda, insitu_r2,higeocat_t0_num[i])
                #print(int((parse_time(ta2).plot_date-parse_time(ta).plot_date)*24))

                [insitu_time3,insitu_r3,insitu_lat3,insitu_lon3]=get_position(ta2, target_name)       
                delta3=abs(higeocat_sse_lon[i]-insitu_lon3)

                if delta3 <30:
                    
                    [ta3,visse3]=calculate_arrival(higeocat_vsse[i],delta3, lamda, insitu_r3,higeocat_t0_num[i])
                    #print(np.round((parse_time(ta3).plot_date-parse_time(ta2).plot_date)*24,1),int(delta3))

                    [insitu_time4,insitu_r4,insitu_lat4,insitu_lon4]=get_position(ta3, target_name)       
                    delta4=abs(higeocat_sse_lon[i]-insitu_lon4)

                    if delta4 <30:
                        
                        #calculate finally iterated arrival time
                        [ta4,visse4]=calculate_arrival(higeocat_vsse[i],delta4, lamda, insitu_r4,higeocat_t0_num[i])
                        #print(np.round((parse_time(ta4).plot_date-parse_time(ta3).plot_date)*24,1),int(delta4))  
 
                        #print(int(delta4-delta))                                            
                                                
                        #estimate error bar on arrival time adding or subtracting the error in the Vsse speed
                        [ta4_low,visse4_low]=calculate_arrival(higeocat_vsse[i]-higeocat_vsse_err[i],delta4, lamda, insitu_r4,higeocat_t0_num[i])
                        [ta4_high,visse4_high]=calculate_arrival(higeocat_vsse[i]+higeocat_vsse_err[i],delta4, lamda, insitu_r4,higeocat_t0_num[i])
                        
                        #calculate difference in ours high / low to original arrival time and convert to hours
                        ta4_err_low=abs(parse_time(ta4).plot_date-parse_time(ta4_low).plot_date)*24
                        ta4_err_high=abs(parse_time(ta4).plot_date-parse_time(ta4_high).plot_date)*24
                        ta4_err=np.round(np.mean([ta4_err_high,ta4_err_low]),1)
                        #print(ta4_err_low,ta4_err_high,ta4_err)
                 
                  
                        #same for arrival speed error
                        visse4_err_low=abs(visse4_low-visse4)
                        visse4_err_high=abs(visse4_high-visse4)
                        visse4_err=int(np.rint(np.mean([visse4_err_high,visse4_err_low])))
                        #print(visse4_err_low,visse4_err_high,visse4_err,higeocat_vsse_err[i])
                        #print()
                        #
                        #print(ta4)
                        #print(str(ta4[0].isoformat())[0:16])
                        #print(str(higeocat_t0[i].isoformat())[0:16]+'Z')
                        #print(higeocat_t0[i])
                        
                        list1=[higeocat_id[i],higeocat_sc[i],target_name,\
                                str(higeocat_t0[i].isoformat())[0:16]+'Z',str(ta4.isoformat())[0:16]+'Z',ta4_err,\
                                np.round(insitu_r4,3), np.round(insitu_lon4,2), np.round(insitu_lat4,2),np.round(insitu_lon4-higeocat_sse_lon[i],1),\
                                higeocat_sse_lon[i],higeocat_sse_lat[i],higeocat_vsse[i],\
                                higeocat_vsse_err[i], int(np.rint(visse4)),visse4_err,higeocat_pafit[i],higeocat_pan[i],higeocat_pas[i],higeocat_pacenter[i]]
                        #print(list1)
                        arrcat_insitu_list.append(list1)
    ###############
                    

    #arrcat_insitu=np.array(arrcat_insitu_list)    
    #print(arrcat_insitu_list)

    #make dataframe out of list
    ac1 = pd.DataFrame(arrcat_insitu_list, columns = column_list)    
    #arrcat=arrcat.append(ac1)   
    arrcat=pd.concat([arrcat,ac1])
    
    
    print('SSEF30 events: ',len(arrcat_insitu_list)   ) 
    print(target_name,' SSEF30 arrival catalog finished.')
    print()
        
    
    return arrcat
















###################################### SIRCAT operations ################################



def load_helio4cast_sircat_master_from_excel(file):
    ''' convert excel master file to pandas dataframe and convert times
        to datetime objects
    '''
    print('load HELCATS SIRCAT from file:', file)
    sc=pd.read_excel(file)
    sc=sc.drop(columns='Unnamed: 0')
    
    
    #get beginning of tags for STA to identify allen and jian events
    tag_list_sta=[]
    for i in np.arange(0,len(sc)):
        tag_list_sta.append(sc.sircat_id[i][13]) #j

    #get beginning of tags for STA to identify allen and grandin events
    tag_list_wind=[]
    for i in np.arange(0,len(sc)):
        tag_list_wind.append(sc.sircat_id[i][9]) #j

        

    #convert all times to datetime objects
    for i in np.arange(0,sc.shape[0]):            
        
 
        #for STEREO and MAVEN same
        if sc.sc_insitu[i] == 'STEREO-A':
            
            #jian events
            if tag_list_sta[i] =='J':
                #remove leading and ending blank spaces if any and write datetime object into dataframe
                sc.at[i,'sir_start_time']= parse_time(str(sc.sir_start_time[i]).strip()).datetime 
                #sc.at[i,'hss_start_time']= parse_time(str(sc.hss_start_time[i]).strip()).datetime       
                sc.at[i,'sir_end_time']= parse_time(str(sc.sir_end_time[i]).strip()).datetime       

            #allen events
            if tag_list_sta[i] =='A':
                #remove leading and ending blank spaces if any and write datetime object into dataframe
                sc.at[i,'sir_start_time']= parse_time(str(sc.sir_start_time[i]).strip()).datetime       
                sc.at[i,'sir_end_time']= parse_time(str(sc.sir_end_time[i]).strip()).datetime       

                
            
        #for Wind PSP convert different - check PSP wind different sources if needed (Allen and Grandin)
        if sc.sc_insitu[i] == 'Wind': 
            
            
            #Allen events
            if tag_list_wind[i] =='A':
            
                sc.at[i,'sir_start_time']= parse_time(str(sc.sir_start_time[i]).strip()).datetime 
                sc.at[i,'sir_end_time']=parse_time(str(sc.sir_end_time[i]).strip()).datetime     
            
            #Grandin events
            if tag_list_wind[i] =='G':
            
                sc.at[i,'hss_start_time']= parse_time(str(sc.hss_start_time[i]).strip()).datetime 
                sc.at[i,'hss_end_time']=parse_time(str(sc.hss_end_time[i]).strip()).datetime     

            

        if sc.sc_insitu[i] == 'PSP': 
            sc.at[i,'sir_start_time']= parse_time(str(sc.sir_start_time[i]).strip()).datetime 
            sc.at[i,'sir_end_time']=parse_time(str(sc.sir_end_time[i]).strip()).datetime     

            
        if sc.sc_insitu[i] == 'STEREO-B':
            #remove leading and ending blank spaces if any and write datetime object into dataframe
            sc.at[i,'sir_start_time']= parse_time(str(sc.sir_start_time[i]).strip()).datetime 
            #sc.at[i,'hss_start_time']= parse_time(str(sc.hss_start_time[i]).strip()).datetime       
            sc.at[i,'sir_end_time']= parse_time(str(sc.sir_end_time[i]).strip()).datetime       

            
        if sc.sc_insitu[i] == 'MAVEN':
            #remove leading and ending blank spaces if any and write datetime object into dataframe
            sc.at[i,'sir_start_time']= parse_time(str(sc.sir_start_time[i]).strip()).datetime 
            sc.at[i,'hss_start_time']= parse_time(str(sc.hss_start_time[i]).strip()).datetime       
            sc.at[i,'sir_end_time']= parse_time(str(sc.sir_end_time[i]).strip()).datetime       



            
            
    return sc




def get_sircat_parameters(sc, sci, scat, name):
    '''
    get parameters
    sc - spacecraft data recarray
    sci - indscates for this spacecraft in sircat
    scat - scatmecat pandas dataframe
    '''
    fileind='sircat/indices_sircat/SIRCAT_indices_'+name+'.p'    
    

    ################ extract indices of ICMEs in the respective data (time consuming, so do it once and save)
    
    if os.path.isfile(fileind) == False:
    
        print('extract indices of SIRs in '+ name+ ' data')
        #### get all ICMECAT times for this spacecraft as datenum
        sc_sir_start=scat.sir_start_time[sci]
        sc_hss_start=scat.hss_start_time[sci]
        sc_sir_end=scat.sir_end_time[sci]
        sc_hss_end=scat.hss_end_time[sci]


        ### arrays containing the indices of where the SIRs are in the data
        sir_start_ind=np.zeros(len(sci),dtype=int)
        hss_start_ind=np.zeros(len(sci),dtype=int)
        sir_end_ind=np.zeros(len(sci),dtype=int)
        hss_end_ind=np.zeros(len(sci),dtype=int)

        #check where vt is < or > 450 km/s
        vt_lt_450=np.where(sc.vt < 450)[0]
        vt_gt_450=np.where(sc.vt > 450)[0]

        
        #check where vt is < or > 350 km/s
        vt_lt_350=np.where(sc.vt < 350)[0]
        vt_gt_350=np.where(sc.vt > 350)[0]



        
        #this takes some time, get indices in data for each SIRCAT time
        for i in np.arange(sci[0],sci[-1]+1):
        
            print(i-sci[0])

            if (name== 'STEREO-A'):
        
                tag=scat.sircat_id[i][13]
             
                if tag=='J': #Jian events
                    print('J', sc_sir_start[i] )
                    sir_start_ind[i-sci[0]]=np.where(sc.time > sc_sir_start[i])[0][0]-1   
                    #hss_start_ind[i-sci[0]]=np.where(sc.time > sc_hss_start[i])[0][0]-1   
                    sir_end_ind[i-sci[0]]=np.where(sc.time   > sc_sir_end[i])[0][0]-1  
                    
                if tag=='A': #Allen events
                    print('A', sc_sir_start[i])
                    sir_start_ind[i-sci[0]]=np.where(sc.time > sc_sir_start[i])[0][0]-1   
                    sir_end_ind[i-sci[0]]=np.where(sc.time   > sc_sir_end[i])[0][0]-1     
                



            if (name== 'STEREO-B'):

                sir_start_ind[i-sci[0]]=np.where(sc.time > sc_sir_start[i])[0][0]-1   
                hss_start_ind[i-sci[0]]=np.where(sc.time > sc_hss_start[i])[0][0]-1   
                sir_end_ind[i-sci[0]]=np.where(sc.time   > sc_sir_end[i])[0][0]-1 
                
                #here the hss_end_time needs to be extracted - criteria similar to Grandin et al. 2018 
                #where stream goes back to (< 450 km/s) after hss start time
                #check the indices in the 450 array that are greater than the hss_start index +0.5 days           
                #24*60 data points
                #and take the first one           
                
                
                #take next data point > 450 km/s after hss_start + 6 hours (for getting rid of rapid variations)
                #next450=np.where(vt_gt_450 > hss_start_ind[i-sci[0]])[0][0]+6*60
                #print(hss_start_ind[i-sci[0]],vt_gt_450[next450])

                #then take next data point below 450 after this 
                #hss_end_ind[i-sci[0]]=vt_lt_450[ np.where(vt_lt_450 > vt_gt_450[next450])[0][0]   ]              
                                
                #print('hss duration in hours ',(hss_end_ind[i-sci[0]]-hss_start_ind[i-sci[0]])/60)           
                #print(hss_start_ind[i-sci[0]],hss_end_ind[i-sci[0]])           
                                
            if name== 'MAVEN':
                
                sir_start_ind[i-sci[0]]=np.where(sc.time > sc_sir_start[i])[0][0]-1   
                hss_start_ind[i-sci[0]]=np.where(sc.time > sc_hss_start[i])[0][0]-1   
                sir_end_ind[i-sci[0]]=np.where(sc.time   > sc_sir_end[i])[0][0]-1 
                
                
                #hss_end_ind[i-sci[0]]=vt_lt_450[np.where(vt_lt_450 > sir_end_ind[i-sci[0]])[0][0]  ]    
                
                      
                #take next data point > 450 km/s after hss_start + 2 orbits (for getting rid of rapid variations)
                #next350=np.where(vt_gt_350 > hss_start_ind[i-sci[0]])[0][0]+2
                #print(hss_start_ind[i-sci[0]],vt_gt_450[next450])

                #then take next data point below 450 after this 
                #hss_end_ind[i-sci[0]]=vt_lt_350[ np.where(vt_lt_350 > vt_gt_350[next350])[0][0]   ]              
                                
                #print('hss duration in hours ',(hss_end_ind[i-sci[0]]-hss_start_ind[i-sci[0]])*4.5)           
                #print(hss_start_ind[i-sci[0]],hss_end_ind[i-sci[0]])                           
                
        
        
            if name=='Wind':      
                
                  
                tag=scat.sircat_id[i][9]
                
                #Allen events
                if tag[i] =='A':

                    sir_start_ind[i-sci[0]]=np.where(sc.time > sc_sir_start[i])[0][0]-1   
                    sir_end_ind[i-sci[0]]=np.where(sc.time   > sc_sir_end[i])[0][0]-1       

                #Grandin events
                if tag[i] =='G':

                    hss_start_ind[i-sci[0]]=np.where(sc.time > sc_hss_start[i])[0][0]-1   
                    hss_end_ind[i-sci[0]]=np.where(sc.time   > sc_hss_end[i])[0][0]-1     

                
                
            if name=='PSP':      
                
                #here only sir start and sir end exist
                sir_start_ind[i-sci[0]]=np.where(sc.time > sc_sir_start[i])[0][0]-1   
                sir_end_ind[i-sci[0]]=np.where(sc.time   > sc_sir_end[i])[0][0]-1   
                
            

        pickle.dump([sir_start_ind,hss_start_ind,sir_end_ind,hss_end_ind], open(fileind, 'wb'))
    ############################################            
                
    [sir_start_ind, hss_start_ind,sir_end_ind,hss_end_ind]=pickle.load(open(fileind, 'rb'))           
    

    print('Get parameters for ',name)
    
   
    print('position')

    #SIR heliodistance
    #    for i in np.arange(len(sci))-1:
    #        scat.at[sci[i],'sc_heliodistance']=np.round(sc.r[hss_start_ind[i]],4)
    #        #SIR longitude
    #        scat.at[sci[i],'sc_long_heeq']=np.round(sc.lon[hss_start_ind[i]],2)
    #        ##SIR latitude
    #        scat.at[sci[i],'sc_lat_heeq']=np.round(sc.lat[hss_start_ind[i]],2)

    print('hss')
        
    if (name=='PSP'):

        #position
        for i in np.arange(0,len(sci)):
            
            scat.at[sci[i],'sc_heliodistance']=np.round(sc.r[sir_start_ind[i]],4)
            #SIR longitude
            scat.at[sci[i],'sc_long_heeq']=np.round(sc.lon[sir_start_ind[i]],2)
            ##SIR latitude
            scat.at[sci[i],'sc_lat_heeq']=np.round(sc.lat[sir_start_ind[i]],2)
            
            sci_istart=mdates.date2num(scat.sir_start_time[sci[i]])       
            sci_iend=mdates.date2num(scat.sir_end_time[sci[i]])   
            scat.at[sci[i],'sir_duration']=np.round((sci_iend-sci_istart)*24,2)

            vmax=np.round(np.nanmax(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)
 
            
            scat.at[sci[i],'sir_vtmax']=vmax
            scat.at[sci[i],'sir_vtmax_time']=sc.time[np.nanargmax(sc.vt[sir_start_ind[i]:sir_end_ind[i]])+sir_start_ind[i]]                     
       
            scat.at[sci[i],'sir_vtmean']=np.round(np.nanmean(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)
            scat.at[sci[i],'sir_vtstd']=np.round(np.nanstd(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)

            scat.at[sci[i],'sir_btmax']=np.round(np.nanmax(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
            scat.at[sci[i],'sir_btmean']=np.round(np.nanmean(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
            scat.at[sci[i],'sir_btstd']=np.round(np.nanstd(sc.bt[hss_start_ind[i]:sir_end_ind[i]]),1)
            scat.at[sci[i],'sir_bzmin']=np.round(np.nanmin(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
            scat.at[sci[i],'sir_bzmean']=np.round(np.nanmean(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
            scat.at[sci[i],'sir_bzstd']=np.round(np.nanstd(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)



        
        
    if (name== 'Wind'): 

        
        for i in np.arange(0,len(sci)):        
            
            
            ############ HSS duration
            sci_istart=mdates.date2num(scat.hss_start_time[sci[i]])       
            sci_hss_iend=mdates.date2num(scat.hss_end_time[sci[i]])   
            scat.at[sci[i],'hss_duration']=np.round((sci_hss_iend-sci_istart)*24,2)


            #print(i)
            #print('hss duration in hours ',(hss_end_ind[i]-hss_start_ind[i])/60)
            tag=scat.sircat_id[i][13]

            #v_max
            scat.at[sci[i],'hss_vtmax']=np.round(np.nanmax(sc.vt[hss_start_ind[i]:hss_end_ind[i]]),1)

            #vtmaxtime - search for index in sliced array and at beginning of array to see the index in the whole dataset
            #scat.at[sci[i],'hss_vtmax_time']=sc.time[np.nanargmax(sc.vt[hss_start_ind[i]:hss_end_ind[i]])]        
            # v_mean
            scat.at[sci[i],'hss_vtmean']=np.round(np.nanmean(sc.vt[hss_start_ind[i]:hss_end_ind[i]]),1)
            #v_bstd
            scat.at[sci[i],'hss_vtstd']=np.round(np.nanstd(sc.vt[hss_start_ind[i]:hss_end_ind[i]]),1)

            #B_max
            scat.at[sci[i],'hss_btmax']=np.round(np.nanmax(sc.bt[hss_start_ind[i]:hss_end_ind[i]]),1)
            # B_mean
            scat.at[sci[i],'hss_btmean']=np.round(np.nanmean(sc.bt[hss_start_ind[i]:hss_end_ind[i]]),1)
            #bstd
            scat.at[sci[i],'hss_btstd']=np.round(np.nanstd(sc.bt[hss_start_ind[i]:hss_end_ind[i]]),1)
            #bz
            scat.at[sci[i],'hss_bzmin']=np.round(np.nanmin(sc.bz[hss_start_ind[i]:hss_end_ind[i]]),1)
            scat.at[sci[i],'hss_bzmean']=np.round(np.nanmean(sc.bz[hss_start_ind[i]:hss_end_ind[i]]),1)
            scat.at[sci[i],'hss_bzstd']=np.round(np.nanstd(sc.bz[hss_start_ind[i]:hss_end_ind[i]]),1)
 
        
    print('sir')
    ###SIR parameters only for STEREO and MAVEN
    
    ############ SIR duration
    
    if (name== 'MAVEN'):


        ########## SIR general parameters

        for i in np.arange(0,len(sci)):
            
            
            sci_istart=mdates.date2num(scat.sir_start_time[sci[i]]) 
            sci_iend=mdates.date2num(scat.sir_end_time[sci[i]])   
            scat.at[sci[i],'sir_duration']=np.round((sci_iend-sci_istart)*24,2)
            
            #vtmaxtime
            #scat.at[sci[i],'sir_vtmax_time']=sc.time[np.nanargmax(sc.vt[sir_start_ind[i]:sir_end_ind[i]])+sir_start_ind[i]]     

            #v_max
            scat.at[sci[i],'sir_vtmax']=np.round(np.nanmax(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)
            # v_mean
            scat.at[sci[i],'sir_vtmean']=np.round(np.nanmean(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)
            #v_bstd
            scat.at[sci[i],'sir_vtstd']=np.round(np.nanstd(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)

            #B_max
            scat.at[sci[i],'sir_btmax']=np.round(np.nanmax(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
            # B_mean
            scat.at[sci[i],'sir_btmean']=np.round(np.nanmean(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
            #bstd
            scat.at[sci[i],'sir_btstd']=np.round(np.nanstd(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
            #bz
            scat.at[sci[i],'sir_bzmin']=np.round(np.nanmin(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
            scat.at[sci[i],'sir_bzmean']=np.round(np.nanmean(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
            scat.at[sci[i],'sir_bzstd']=np.round(np.nanstd(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
    
    
    if (name== 'STEREO-B'):
    
        ########## SIR general parameters

        for i in np.arange(0,len(sci)):
            
            
            sci_istart=mdates.date2num(scat.sir_start_time[sci[i]]) 
            sci_iend=mdates.date2num(scat.sir_end_time[sci[i]])   
            scat.at[sci[i],'sir_duration']=np.round((sci_iend-sci_istart)*24,2)
            
            #vtmaxtime
            scat.at[sci[i],'sir_vtmax_time']=sc.time[np.nanargmax(sc.vt[sir_start_ind[i]:sir_end_ind[i]])+sir_start_ind[i]]     

            #v_max
            scat.at[sci[i],'sir_vtmax']=np.round(np.nanmax(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)
            # v_mean
            scat.at[sci[i],'sir_vtmean']=np.round(np.nanmean(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)
            #v_bstd
            scat.at[sci[i],'sir_vtstd']=np.round(np.nanstd(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)

            #B_max
            scat.at[sci[i],'sir_btmax']=np.round(np.nanmax(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
            # B_mean
            scat.at[sci[i],'sir_btmean']=np.round(np.nanmean(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
            #bstd
            scat.at[sci[i],'sir_btstd']=np.round(np.nanstd(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
            #bz
            scat.at[sci[i],'sir_bzmin']=np.round(np.nanmin(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
            scat.at[sci[i],'sir_bzmean']=np.round(np.nanmean(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
            scat.at[sci[i],'sir_bzstd']=np.round(np.nanstd(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
    
    
    

    if (name== 'STEREO-A'):

        
        for i in np.arange(0,len(sci)):        
            
            #check which catalog
            
            tag=scat.sircat_id[sci[i]][13]
        
            if tag=='J': #Jian events
            
            
                sci_istart=mdates.date2num(scat.sir_start_time[sci[i]])   
                sci_iend=mdates.date2num(scat.sir_end_time[sci[i]])   
                scat.at[sci[i],'sir_duration']=np.round((sci_iend-sci_istart)*24,2) 

                #v_max
                scat.at[sci[i],'sir_vtmax']=np.round(np.nanmax(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)
                # v_mean
                scat.at[sci[i],'sir_vtmean']=np.round(np.nanmean(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)
                #v_bstd
                scat.at[sci[i],'sir_vtstd']=np.round(np.nanstd(sc.vt[sir_start_ind[i]:sir_end_ind[i]]),1)

                #B_max
                scat.at[sci[i],'sir_btmax']=np.round(np.nanmax(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
                # B_mean
                scat.at[sci[i],'sir_btmean']=np.round(np.nanmean(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
                #bstd
                scat.at[sci[i],'sir_btstd']=np.round(np.nanstd(sc.bt[sir_start_ind[i]:sir_end_ind[i]]),1)
                #bz
                scat.at[sci[i],'sir_bzmin']=np.round(np.nanmin(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
                scat.at[sci[i],'sir_bzmean']=np.round(np.nanmean(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
                scat.at[sci[i],'sir_bzstd']=np.round(np.nanstd(sc.bz[sir_start_ind[i]:sir_end_ind[i]]),1)
    



            if tag=='A': #Allen events
            
                ############ HSS duration
                sci_istart=mdates.date2num(scat.hss_start_time[sci[i]])       
                sci_hss_iend=mdates.date2num(scat.hss_end_time[sci[i]])   
                scat.at[sci[i],'hss_duration']=np.round((sci_hss_iend-sci_istart)*24,2)


                #v_max
                scat.at[sci[i],'hss_vtmax']=np.round(np.nanmax(sc.vt[hss_start_ind[i]:hss_end_ind[i]]),1)

                #vtmaxtime - search for index in sliced array and at beginning of array to see the index in the whole dataset
                scat.at[sci[i],'hss_vtmax_time']=sc.time[np.nanargmax(sc.vt[hss_start_ind[i]:hss_end_ind[i]])+hss_start_ind[i]]        
                # v_mean
                scat.at[sci[i],'hss_vtmean']=np.round(np.nanmean(sc.vt[hss_start_ind[i]:hss_end_ind[i]]),1)
                #v_bstd
                scat.at[sci[i],'hss_vtstd']=np.round(np.nanstd(sc.vt[hss_start_ind[i]:hss_end_ind[i]]),1)

                #B_max
                scat.at[sci[i],'hss_btmax']=np.round(np.nanmax(sc.bt[hss_start_ind[i]:hss_end_ind[i]]),1)
                # B_mean
                scat.at[sci[i],'hss_btmean']=np.round(np.nanmean(sc.bt[hss_start_ind[i]:hss_end_ind[i]]),1)
                #bstd
                scat.at[sci[i],'hss_btstd']=np.round(np.nanstd(sc.bt[hss_start_ind[i]:hss_end_ind[i]]),1)
                #bz
                scat.at[sci[i],'hss_bzmin']=np.round(np.nanmin(sc.bz[hss_start_ind[i]:hss_end_ind[i]]),1)
                scat.at[sci[i],'hss_bzmean']=np.round(np.nanmean(sc.bz[hss_start_ind[i]:hss_end_ind[i]]),1)
                scat.at[sci[i],'hss_bzstd']=np.round(np.nanstd(sc.bz[hss_start_ind[i]:hss_end_ind[i]]),1)  


    
    
    return scat





















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








def create_icme_indices(sc,sci,ic,name):

    # extract indices of ICMEs in the respective data and save
    
    #this is the bottleneck
    sc_numtime=mdates.date2num(sc.time)

    sc_icme_start=mdates.date2num(ic.icme_start_time[sci])
    sc_mo_start=mdates.date2num(ic.mo_start_time[sci])
    sc_mo_end=mdates.date2num(ic.mo_end_time[sci])
  
    #search where the minimum is in the array between the given time and the time in the data array
    #this is fast
    
    #for i in np.arange(0,len(sc_icme_start)):
    
    #     print(np.argmin(abs(sc_icme_start[i]-sc_numtime)))
    #     print(np.argmin(abs(sc_mo_start[i]-sc_numtime)))
    #     print(np.argmin(abs(sc_mo_end[i]-sc_numtime)))
    #     print()
    
    icme_start_ind=[np.argmin(abs(sc_icme_start[i]-sc_numtime)) for i in np.arange(0,len(sc_icme_start))]
    mo_start_ind=[np.argmin(abs(sc_mo_start[i]-sc_numtime)) for i in np.arange(0,len(sc_mo_end))]
    mo_end_ind=[np.argmin(abs(sc_mo_end[i]-sc_numtime)) for i in np.arange(0,len(sc_mo_end))]

    
    
    fileind='icmecat/indices_icmecat/ICMECAT_indices_'+name+'.p'
    pickle.dump([icme_start_ind, mo_start_ind,mo_end_ind], open(fileind, 'wb'))
    
    print(name,'indices done')
    









def get_cat_parameters(sc, sci, ic, name):
    '''
    get parameters
    sc - spacecraft data recarray
    sci - indices for this spacecraft in icmecat
    ic - icmecat pandas dataframe
    '''
    fileind='icmecat/indices_icmecat/ICMECAT_indices_'+name+'.p'
                
    [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))           
    
    
    #plasma available?
    if name=='Wind': plasma=True
    if name=='STEREO-A': plasma=True
    if name=='STEREO-B': plasma=True
    if name=='ULYSSES': plasma=True
    if name=='MAVEN': plasma=True
    if name=='PSP': plasma=True
    if name=='SolarOrbiter': plasma=True    
    
    if name=='VEX': plasma=False
    if name=='MESSENGER': plasma=False
    if name=='BepiColombo': plasma=False
    if name=='Juno': plasma=False

    print('Get parameters for ',name)
    

    ####### position

    #MO heliodistance
    for i in np.arange(len(sci))-1:
        ic.at[sci[i],'mo_sc_heliodistance']=np.round(sc.r[mo_start_ind[i]],4)

        #MO longitude
        ic.at[sci[i],'mo_sc_long_heeq']=np.round(sc.lon[mo_start_ind[i]],2)

        #MO latitude
        ic.at[sci[i],'mo_sc_lat_heeq']=np.round(sc.lat[mo_start_ind[i]],2)

    #######durations 
    
    for i in np.arange(len(sci))-1:

        sci_istart=mdates.date2num(ic.icme_start_time[sci[i]])   
        sci_iend=mdates.date2num(ic.mo_end_time[sci[i]])   
        ic.at[sci[i],'icme_duration']=np.round((sci_iend-sci_istart)*24,2)
        
        # MO duration
        sci_istart=mdates.date2num(ic.mo_start_time[sci[i]])   
        sci_iend=mdates.date2num(ic.mo_end_time[sci[i]])   
        ic.at[sci[i],'mo_duration']=np.round((sci_iend-sci_istart)*24,2)      
    
    
    
    #########icme values

    for i in np.arange(len(sci))-1:

        # ########## use for debugging to find event with wrong times
        #print(sci[i])
        #print(ic.mo_start_time[sci[i]] )   
        #print(ic.mo_end_time[sci[i]] )  
        #ICME B_max
        ic.at[sci[i],'icme_bmax']=np.round(np.nanmax(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)

        #ICME B_mean
        ic.at[sci[i],'icme_bmean']=np.round(np.nanmean(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)

        #icme_bstd
        ic.at[sci[i],'icme_bstd']=np.round(np.nanstd(sc.bt[icme_start_ind[i]:mo_end_ind[i]]),1)
            
        #ic.at[sci[i],'icme_bmax']=np.nan
        #ic.at[sci[i],'icme_bmean']=np.nan
        #ic.at[sci[i],'icme_bstd']=np.nan
        
    if plasma==True:        
        #ICME speed_mean and std
        
        ###### check whether all nan values in this interval - set to nan and skip next        
        
        for i in np.arange(len(sci))-1:
            ic.at[sci[i],'icme_speed_mean']=np.round(np.nanmean(sc.vt[icme_start_ind[i]:mo_end_ind[i]]),1)
            ic.at[sci[i],'icme_speed_std']=np.round(np.nanstd(sc.vt[icme_start_ind[i]:mo_end_ind[i]]),1)            
    else: #set nan    
        for i in np.arange(len(sci))-1:
            ic.at[sci[i],'icme_speed_mean']=np.nan
            ic.at[sci[i],'icme_speed_std']=np.nan

       
  

    
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



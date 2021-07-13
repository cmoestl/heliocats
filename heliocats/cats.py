#stats.py
#catalog creation for heliocats
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
import astropy
from astropy.constants import au
import importlib
import cdflib
import matplotlib.pyplot as plt
import heliosat

import heliopy.data.spice as spicedata
import heliopy.spice as spice

from astropy.io.votable import parse_single_table

from config import data_path


from heliocats import data as hd
importlib.reload(hd) #reload again while debugging



#define AU in km
AU=au.value/1e3



######################## general position functions


# def get_mars_position_array():
    
#     ############### Mars position

#     planet_kernel=spicedata.get_kernel('planet_trajectories')
#     starttime = datetime.datetime(2007, 1, 1)
#     endtime = datetime.datetime(2020, 12, 31)
#     res_in_hours=1
#     mars_time = []
#     while starttime < endtime:
#         mars_time.append(starttime)
#         starttime += datetime.timedelta(hours=res_in_hours)
#     mars=spice.Trajectory('4')  
#     frame='HEEQ'
#     mars.generate_positions(mars_time,'Sun',frame)  
#     mars.change_units(astropy.units.AU)  
#     [mars_r, mars_lat, mars_lon]=hd.cart2sphere(mars.x,mars.y,mars.z)
#     print('mars position done') 
    
#     mars_time=np.array(mars_time)
#     mars_r=np.array(mars_r)
#     mars_lat=np.array(mars_lat)
#     mars_lon=np.array(mars_lon)

#     return [mars_time,mars_r,np.degrees(mars_lat),np.degrees(mars_lon)]






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



def get_insitu_position_time(time1,insitu_location_string,insitu_str,insitu_kernel):
    
    
    insitu_exist=True
    
    if insitu_location_string=='PSP': 
        #exclude if time before launch time
        if parse_time(time1).plot_date < parse_time(datetime.datetime(2018, 8, 13)).plot_date:
            insitu_exist=False 

    if insitu_location_string=='Solo': 
        if parse_time(time1).plot_date < parse_time(datetime.datetime(2020, 3, 1)).plot_date:
            insitu_exist=False 
        
    if insitu_location_string=='Bepi': 
        if parse_time(time1).plot_date < parse_time(datetime.datetime(2018, 10, 24)).plot_date:
            insitu_exist=False
                    
    if insitu_location_string=='STB': 
        if parse_time(time1).plot_date > parse_time(datetime.datetime(2014, 9, 27)).plot_date:
            insitu_exist=False  
            
                               
    if insitu_location_string=='Ulysses': 
        #cut off ulysses when no decent in situ data is available anymore
        if parse_time(time1).plot_date > parse_time(datetime.datetime(2008, 5, 1)).plot_date:
            insitu_exist=False              

            
    if insitu_exist == True:
        #insitu_kernel=spicedata.get_kernel('insitu_trajectories')

        #this needs to be an array, so make two similar times and take the first entry later
        insitu_time=[parse_time(time1).datetime,parse_time(time1).datetime]
        insitu=spice.Trajectory(insitu_str)  
        frame='HEEQ'
        insitu.generate_positions(insitu_time,'Sun',frame)  
        insitu.change_units(astropy.units.AU)  
        [insitu_r, insitu_lat, insitu_lon]=hd.cart2sphere(insitu.x,insitu.y,insitu.z)
        
        #Earth position to Earth L1
        if insitu_str=='3': insitu_r[0]=insitu_r[0]-1.5*1e6/AU
       

        insitu_time=np.array(insitu_time)[0]
        insitu_r=np.array(insitu_r)[0]
        insitu_lat=np.array(insitu_lat)[0]
        insitu_lon=np.array(insitu_lon)[0]
        
    else:
        insitu_time=np.nan
        insitu_r=np.nan
        insitu_lat=np.nan
        insitu_lon=np.nan    

    return [insitu_time,insitu_r,np.degrees(insitu_lat),np.degrees(insitu_lon)]




def calculate_arrival(vsse,delta,lamda,rdist,t0_num):
   
    #calculate arrival time after Möstl and Davies 2013 but using ta=t0+Ri/Visse equivalent to ta=t0+Risse/Vsse    
   
    visse=vsse * (  np.cos(np.radians(delta))   \
                  + np.sqrt(  np.sin(np.radians(lamda))**2-np.sin(np.radians(delta))**2 ) ) \
                    /(1+np.sin(np.radians(lamda)) )      
                              
    
    #arrival time: convert AU to km  and seconds to days                
    ta=t0_num+(rdist*AU/visse)/(3600*24)        
                              
    return [mdates.num2date(ta),visse]





def make_arrival_catalog_insitu_ssef30(higeocat,arrcat,ac_old, insitu_location_string, column_list):
    
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
  
    
    
    #load spice here once for each spacecraft
        
    if insitu_location_string=='STB': 
        insitu_str='-235'
        insitu_kernel=spicedata.get_kernel('stereo_b')    
        target_name='STEREO-B'
    
    if insitu_location_string=='STA': 
        insitu_str='-234'
        insitu_kernel=spicedata.get_kernel('stereo_a_pred')
        insitu_kernel2=spicedata.get_kernel('stereo_a')
        spice.furnish(insitu_kernel2)
        target_name='STEREO-A'
       
    if insitu_location_string=='Mercury': 
        insitu_str='1'
        insitu_kernel=spicedata.get_kernel('planet_trajectories')
        target_name='Mercury'
        
    if insitu_location_string=='Venus': 
        insitu_str='2'
        insitu_kernel=spicedata.get_kernel('planet_trajectories')
        target_name='Venus'
       
    if insitu_location_string=='Earth': 
        insitu_str='3'
        insitu_kernel=spicedata.get_kernel('planet_trajectories')
        target_name='Earth_L1'
        
    if insitu_location_string=='Mars': 
        insitu_str='4'
        insitu_kernel=spicedata.get_kernel('planet_trajectories')        
        target_name='Mars'
        
    if insitu_location_string=='PSP': 
        insitu_str='-96'
        insitu_kernel=spicedata.get_kernel('psp_pred')
        target_name='PSP'

    if insitu_location_string=='Solo': 
        insitu_str='Solar Orbiter'
        insitu_kernel=spicedata.get_kernel('solo_2020')   
        target_name='SolarOrbiter'        
        
    if insitu_location_string=='Bepi': 
        insitu_str='BEPICOLOMBO MPO'
        insitu_kernel=spicedata.get_kernel('bepi_pred')
        target_name='BepiColombo'

    if insitu_location_string=='Ulysses': 
        insitu_str='ulysses'
        insitu_kernel=spicedata.get_kernel('ulysses')
        target_name='Ulysses'


        
    spice.furnish(insitu_kernel)
 
           

    #half width for SSEF30
    lamda=30.0

    #new version of ARRCAT with iteration
    arrcat_insitu_list = []
    #old version without iteration
    arrcat_insitu_list_old = []



    #go through all HIGEOCAT CME events and check for hit at insitu, with 4 iterations in total
    for i in np.arange(len(higeocat_time)):

        #get insitu position for launch time t0    
        [insitu_time,insitu_r,insitu_lat,insitu_lon]=get_insitu_position_time(higeocat_t0[i], insitu_location_string,insitu_str, insitu_kernel)            
        delta=abs(higeocat_sse_lon[i]-insitu_lon)
        #print([insitu_time,insitu_r,insitu_lat,insitu_lon])

        if delta < 30:               

            #calculate arrival time 
            #print(delta,lamda,insitu_r)
            [ta,visse]=calculate_arrival(higeocat_vsse[i],delta, lamda, insitu_r,higeocat_t0_num[i])                
            
            #make old version of ARRCAT without iteration and errors
            list_old=[higeocat_id[i].decode(),higeocat_sc[i].decode(),target_name,\
                   parse_time(higeocat_t0[i]).iso[:-7],parse_time(ta).iso[:-7],0,\
                   np.round(insitu_r,3), np.round(insitu_lon,2), np.round(insitu_lat,2),np.round(insitu_lon-higeocat_sse_lon[i],1),\
                   higeocat_sse_lon[i],higeocat_sse_lat[i],higeocat_vsse[i],\
                   higeocat_vsse_err[i], int(np.rint(visse)),0,higeocat_pafit[i],higeocat_pan[i],higeocat_pas[i],higeocat_pacenter[i]]
                   #print(list1)
            arrcat_insitu_list_old.append(list_old)
            
        

            [insitu_time2,insitu_r2,insitu_lat2,insitu_lon2]=get_insitu_position_time(ta, insitu_location_string,insitu_str, insitu_kernel)       
            #print(insitu_lon-insitu_lon2)               
            delta2=abs(higeocat_sse_lon[i]-insitu_lon2)
            if delta2 <30:

                [ta2,visse2]=calculate_arrival(higeocat_vsse[i],delta2, lamda, insitu_r2,higeocat_t0_num[i])
                #print(int((parse_time(ta2).plot_date-parse_time(ta).plot_date)*24))

                [insitu_time3,insitu_r3,insitu_lat3,insitu_lon3]=get_insitu_position_time(ta2, insitu_location_string,insitu_str, insitu_kernel)       
                delta3=abs(higeocat_sse_lon[i]-insitu_lon3)

                if delta3 <30:
                    [ta3,visse3]=calculate_arrival(higeocat_vsse[i],delta3, lamda, insitu_r3,higeocat_t0_num[i])
                    #print(np.round((parse_time(ta3).plot_date-parse_time(ta2).plot_date)*24,1),int(delta3))

                    [insitu_time4,insitu_r4,insitu_lat4,insitu_lon4]=get_insitu_position_time(ta3, insitu_location_string,insitu_str, insitu_kernel)       
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

                        
                        list1=[higeocat_id[i].decode(),higeocat_sc[i].decode(),target_name,\
                                parse_time(higeocat_t0[i]).iso[:-7],parse_time(ta4).iso[:-7],ta4_err,\
                                np.round(insitu_r4,3), np.round(insitu_lon4,2), np.round(insitu_lat4,2),np.round(insitu_lon4-higeocat_sse_lon[i],1),\
                                higeocat_sse_lon[i],higeocat_sse_lat[i],higeocat_vsse[i],\
                                higeocat_vsse_err[i], int(np.rint(visse4)),visse4_err,higeocat_pafit[i],higeocat_pan[i],higeocat_pas[i],higeocat_pacenter[i]]
                        #print(list1)
                        arrcat_insitu_list.append(list1)



                    

    #arrcat_insitu=np.array(arrcat_insitu_list)    
    #print(arrcat_insitu_list)

    
    #make dataframe out of list
    ac_old1 = pd.DataFrame(arrcat_insitu_list_old, columns = column_list)    
    ac_old=ac_old.append(ac_old1)   

    
    #make dataframe out of list
    ac1 = pd.DataFrame(arrcat_insitu_list, columns = column_list)    
    arrcat=arrcat.append(ac1)   
    
    
    print('SSEF30 events: ',len(arrcat_insitu_list)   ) 
    print(insitu_location_string,' SSEF30 arrival catalog finished.')
    print()
        
    
    return [arrcat,ac_old]






###################################### SIRCAT operations ################################



def load_helio4cast_sircat_master_from_excel(file):
    ''' convert excel master file to pandas dataframe and convert times
        to datetime objects
    '''
    print('load HELCATS SIRCAT from file:', file)
    sc=pd.read_excel(file)
    sc=sc.drop(columns='Unnamed: 0')
    
    
    #get beginning of tags for STA to identify allen and jian events
    tag_list=[]
    for i in np.arange(0,len(sc)):
        tag_list.append(sc.sircat_id[i][13]) #j

    #convert all times to datetime objects
    for i in np.arange(0,sc.shape[0]):            
        
 
        #for STEREO and MAVEN same
        if sc.sc_insitu[i] == 'STEREO-A':
            
            #jian events
            if tag_list[i] =='J':
                #remove leading and ending blank spaces if any and write datetime object into dataframe
                sc.at[i,'sir_start_time']= parse_time(str(sc.sir_start_time[i]).strip()).datetime 
                sc.at[i,'hss_start_time']= parse_time(str(sc.hss_start_time[i]).strip()).datetime       
                sc.at[i,'sir_end_time']= parse_time(str(sc.sir_end_time[i]).strip()).datetime       

            #allen events
            if tag_list[i] =='A':
                #remove leading and ending blank spaces if any and write datetime object into dataframe
                sc.at[i,'hss_start_time']= parse_time(str(sc.hss_start_time[i]).strip()).datetime       
                sc.at[i,'hss_end_time']= parse_time(str(sc.hss_end_time[i]).strip()).datetime       

                
            
        #for Wind PSP convert different - check PSP wind different sources if needed (Allen and Grandin)
        if sc.sc_insitu[i] == 'Wind': 
            sc.at[i,'hss_start_time']= parse_time(str(sc.hss_start_time[i]).strip()).datetime 
            sc.at[i,'hss_end_time']=parse_time(str(sc.hss_end_time[i]).strip()).datetime     
            
            
            

        if sc.sc_insitu[i] == 'PSP': 
            sc.at[i,'hss_start_time']= parse_time(str(sc.hss_start_time[i]).strip()).datetime 
            sc.at[i,'hss_end_time']=parse_time(str(sc.hss_end_time[i]).strip()).datetime     

            
        if sc.sc_insitu[i] == 'STEREO-B':
            #remove leading and ending blank spaces if any and write datetime object into dataframe
            sc.at[i,'sir_start_time']= parse_time(str(sc.sir_start_time[i]).strip()).datetime 
            sc.at[i,'hss_start_time']= parse_time(str(sc.hss_start_time[i]).strip()).datetime       
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
                    hss_start_ind[i-sci[0]]=np.where(sc.time > sc_hss_start[i])[0][0]-1   
                    sir_end_ind[i-sci[0]]=np.where(sc.time   > sc_sir_end[i])[0][0]-1  
                    
                if tag=='A': #Allen events
                    print('A', sc_sir_start[i])
                    hss_start_ind[i-sci[0]]=np.where(sc.time > sc_hss_start[i])[0][0]-1   
                    hss_end_ind[i-sci[0]]=np.where(sc.time   > sc_hss_end[i])[0][0]-1     
                



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
                
                #here only hss start and hss end exist
                hss_start_ind[i-sci[0]]=np.where(sc.time > sc_hss_start[i])[0][0]-1   
                hss_end_ind[i-sci[0]]=np.where(sc.time   > sc_hss_end[i])[0][0]-1     
                
                #future update: set hss_start as sir_start, and add time for hss_start by pt max after sir_start
                
            if name=='PSP':      
                
                #here only hss start and hss end exist
                hss_start_ind[i-sci[0]]=np.where(sc.time > sc_hss_start[i])[0][0]-1   
                hss_end_ind[i-sci[0]]=np.where(sc.time   > sc_hss_end[i])[0][0]-1     
                
                #future update: set hss_start as sir_start, and add time for hss_start by pt max after sir_start
                
            
            
            

        pickle.dump([sir_start_ind,hss_start_ind,sir_end_ind,hss_end_ind], open(fileind, 'wb'))
    ############################################            
                
    [sir_start_ind, hss_start_ind,sir_end_ind,hss_end_ind]=pickle.load(open(fileind, 'rb'))           
    
        
    #first make hss end time for STEREO-A/B from hss_end_ind index
    #if (name== 'STEREO-A') or (name== 'STEREO-B') or (name== 'MAVEN'):
    #      for i in np.arange(len(sci))-1:
    #         scat.at[sci[i],'hss_end_time']=sc.time[hss_end_ind[i]]


    print('Get parameters for ',name)
    
    ####### position
    
    print('position')

    #SIR heliodistance
    for i in np.arange(len(sci))-1:
        scat.at[sci[i],'sc_heliodistance']=np.round(sc.r[hss_start_ind[i]],4)
        #SIR longitude
        scat.at[sci[i],'sc_long_heeq']=np.round(sc.lon[hss_start_ind[i]],2)
        ##SIR latitude
        scat.at[sci[i],'sc_lat_heeq']=np.round(sc.lat[hss_start_ind[i]],2)

    print('hss')
        
    if (name=='PSP'):

        
        sci_istart=mdates.date2num(scat.hss_start_time[sci])       
        sci_hss_iend=mdates.date2num(scat.hss_end_time[sci])   
        scat.at[sci,'hss_duration']=np.round((sci_hss_iend-sci_istart)*24,2)

        
        for i in np.arange(0,len(sci)):        

            #print(i)
            #print('hss duration in hours ',(hss_end_ind[i]-hss_start_ind[i])/60)


            #v_max
            scat.at[sci[i],'hss_vtmax']=np.nan

            try:
                vmax=np.round(np.nanmax(sc.vt[hss_start_ind[i]:hss_end_ind[i]]),1)
     
                #if vmax ok:
                if np.isnan(vmax)==False:
                        scat.at[sci[i],'hss_vtmax']=vmax
                        #vtmaxtime - search for index in sliced array and at beginning of array to see the index in the whole dataset
                        scat.at[sci[i],'hss_vtmax_time']=sc.time[np.nanargmax(sc.vt[hss_start_ind[i]:hss_end_ind[i]])+hss_start_ind[i]]     
                
            except:     
              print('vmax nan')
       
            # v_mean
            try:
                scat.at[sci[i],'hss_vtmean']=np.round(np.nanmean(sc.vt[hss_start_ind[i]:hss_end_ind[i]]),1)
            except:
                print()
            #v_bstd
            try:
                scat.at[sci[i],'hss_vtstd']=np.round(np.nanstd(sc.vt[hss_start_ind[i]:hss_end_ind[i]]),1)
            except:
                print()


            try:
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
            except:
                print()



        
        
    if (name== 'Wind'): 

        ############ HSS duration
        sci_istart=mdates.date2num(scat.hss_start_time[sci])       
        sci_hss_iend=mdates.date2num(scat.hss_end_time[sci])   
        scat.at[sci,'hss_duration']=np.round((sci_hss_iend-sci_istart)*24,2)


        
        for i in np.arange(0,len(sci)):        

            #print(i)
            #print('hss duration in hours ',(hss_end_ind[i]-hss_start_ind[i])/60)
            tag=scat.sircat_id[i][13]

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
 
        
    print('sir')
    ###SIR parameters only for STEREO and MAVEN
    
    ############ SIR duration
    
    if (name== 'STEREO-B') or (name== 'MAVEN'):

        sci_istart=mdates.date2num(scat.hss_start_time[sci])   ##***Fehler? sir_start?
        sci_iend=mdates.date2num(scat.sir_end_time[sci])   
        scat.at[sci,'sir_duration']=np.round((sci_iend-sci_istart)*24,2)


        ########## SIR general parameters

        for i in np.arange(0,len(sci)):

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


def get_cat_parameters(sc, sci, ic, name):
    '''
    get parameters
    sc - spacecraft data recarray
    sci - indices for this spacecraft in icmecat
    ic - icmecat pandas dataframe
    '''
    fileind='icmecat/indices_icmecat/ICMECAT_indices_'+name+'.p'

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
    if name=='SolarOrbiter': plasma=False
    if name=='BepiColombo': plasma=False

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
    

    for i in np.arange(0,len(sci)):
        

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
    
    
    #print(sci_istart)
    #print(sci_iend)
    
    #print(mo_start_ind[i])
    #print(mo_end_ind[i])


    
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



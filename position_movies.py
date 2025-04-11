#!/usr/bin/env python
# coding: utf-8

# ### Spacecraft and planet trajectory movies
# 
# 
# https://github.com/cmoestl/heliocats
#  
# Authors: C. Möstl, Eva Weiler, Emma Davies, Austrian Space Weather Office, GeoSphere Austria
# bluesky @chrisoutofspace, https://github.com/cmoestl
# 
# last update: April 2025
#  
# needs the helio4 environment (see README.md)
# 
# spacecraft position files are made with positions.ipynb
#  
# ---
# 
# Published in figshare repository: https://doi.org/10.6084/m9.figshare.25301260
# 
# **Animations of spacecraft trajectories as context for Solar Orbiter (2020-2030) and PUNCH (2025-2029)**
#  
#   
# Starting from the Solar Orbiter launch in April 2020, we present here animations of the trajectories of Solar Orbiter, Parker Solar Probe, BepiColombo, STEREO-A and JUICE, with the 4 inner planets, for 2020/4-12/2029, in 6 hour time resolution.
# 
# Further, there are 4 more movies with 6 hour time resolution for 2025/4 - 12/2029, covering a shorter timerange as context for the PUNCH mission launched in March 2025. 
# 
# One includes the Mars orbit, and the other shows only the inner heliosphere (named "zoom").
# 
# All movies are available with 1080p (default) and 4k resolution (indicated with 4k in the file name).
# 
# Last update of trajectories: March 2025
# 
# The spacecraft positions have all been generated with spiceypy from the mission kernels.
# 
# Austrian Space Weather Office / GeoSphere Austria, Möstl/Davies/Weiler
# 
# 
# **ISSUES:**
# 
# - to be done: include side view on the bottom left so inclination of SolO can be seen better
# - add lagrange points, create with positions.ipynb first
#  

# In[9]:


import os
import datetime
from datetime import datetime, timedelta


############################ directories

animdirectory   = 'results/positions/movies_2025'
outputdirectory = 'results/positions/movies_2025/frames'
#outputdirectory = 'results/positions/movies_2024/frames_zoom'

#for the movie
#animdirectory   = 'icmecat'
#outputdirectory = 'results/positions/movies_2024/frames'
#for the frames
#outputdirectory = 'icmecat/anim_frames'


if os.path.isdir(outputdirectory) == False: os.mkdir(outputdirectory)
if os.path.isdir(animdirectory) == False: os.mkdir(animdirectory)

#movie_filename='positions_punch_2025_2030'
movie_filename='positions_punch_2025_2030_zoom'
#movie_filename='positions_2020_2030_test_lagrange'
#movie_filename=''


####################### Time resolution

#res_in_hours=24
res_in_hours=6
#res_in_hours=24*14


print('time resolution in hours',res_in_hours)

############### make time range

##PUNCH 4.5 years
t_start = datetime(2025,4,1)
t_end   = datetime(2029,12,31)

## Solar Orbiter 10 years
#t_start = datetime(2020,4,1)
#t_end   = datetime(2029,12,31)

## from psp launch onwards
#t_start = datetime(2018,10,1)
#t_end   = datetime(2030,11,19)


## ICME visuals
#t_start = datetime(2023,3,10)
#t_end   = datetime(2023,5,10)

#t_start = datetime(2020,4,1)
#t_end   = datetime(2024,4,1)


########## plots

#rmax=1.72
rmax=1.1

dpisave=100 #for 1080p
#dpisave=200 #for 4K

##multiprocessing

#used=100
used=8


# In[10]:


###########################################################

ffmpeg_path=''


from sunpy.time import parse_time

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.signal import medfilt
import numpy as np
import pdb
import pickle
import seaborn as sns
import sys
import heliopy.data.spice as spicedata
import heliopy.spice as spice
import astropy
import importlib    
import time
import numba
from numba import jit
import multiprocessing as mp

#ignore warnings
import warnings
warnings.filterwarnings('ignore')

#for server
#matplotlib.use('Agg')

##### check for system type
#server
if sys.platform == 'linux': 
    print('system is linux')
    matplotlib.use('Agg') 
    from config_server import data_path
    
    
#mac - make sure the dpi is always set similar to plt.savefig
if sys.platform =='darwin':  
    print('system is mac')
    #for testing
    get_ipython().run_line_magic('matplotlib', 'inline')
    from config_local import data_path
    #matplotlib.use('Agg') 


print(data_path)



#Convert this notebook to a script 
os.system('jupyter nbconvert --to script position_movies.ipynb')    


# In[11]:


#old files
#[psp, solo, sta, stb, bepi, l1, juno, juice, uly, earth, mercury, venus, mars, jupiter, saturn, uranus, neptune]=pickle.load(open(data_path+'/positions_psp_solo_sta_bepi_wind_juno_juice_ulysses_planets_HEEQ_1hour_rad.p', "rb" ) )


##use 10 min version, these are pandas dataframes
print('load positions')
[psp, bepi, solo, sta, juice, earth, mercury, venus, mars, jupiter, saturn, uranus, neptune,l4,l5]=pickle.load( open( 'results/positions/positions_2020_all_HEEQ_1h_rad_cm.p', "rb" ) )    
print('merged positions loaded')


print('load icmecat')
#load icmecat as numpy array 
file='icmecat/HELIO4CAST_ICMECAT_v23_numpy.p'
[ic,ic_np,h,p]=pickle.load( open(file, 'rb'))   
ic=ic.to_records()

print('done')


# In[12]:


l4.lon[0]
#l4.time[0]
#earth.time[0]


# In[ ]:





# In[ ]:





# In[13]:


psp.time.values[-1]
np.max(solo.lon[0])


# In[14]:


def make_frame(k):
    '''
    loop each frame in multiprocessing
    '''
     
    
    if not black: 
        fig=plt.figure(1, figsize=(19.2,10.8), dpi=100) #full hd
        #fig=plt.figure(1, figsize=(19.2*2,10.8*2), dpi=100) #4k
        ax = plt.subplot2grid((5,2), (0, 0), rowspan=5, projection='polar')
        backcolor='black'
        psp_color='black'
        bepi_color='blue'
        solo_color='green'

    if black: #4k resolution        
        fig=plt.figure(1, figsize=(19.2,10.8), dpi=dpisave, facecolor='black', edgecolor='black')
        ax = plt.subplot(111,projection='polar',facecolor='black') 
        #ax = plt.subplot2grid((5,2), (0, 0), rowspan=5, projection='polar')
        backcolor='white'
        psp_color='white'
        bepi_color='skyblue'
        solo_color='springgreen'
        sta_color='salmon'
        juice_color='gold'
        l4_color='lavender'
        l5_color='lavender'


    frame_time_str=str(mdates.num2date(frame_time_num+k*res_in_hours/24))
    
    #print(frame_time_num)
    print( 'current frame_time_num', frame_time_str, '     ',k)


    #all same times
    dct=frame_time_num+k*res_in_hours/24-earth.time
    earth_timeind=np.argmin(abs(dct))
    
    
    
    #these have their own times
    
    #all same times
    dct=frame_time_num+k*res_in_hours/24-sta.time
    sta_timeind=np.argmin(abs(dct))
    
    dct=frame_time_num+k*res_in_hours/24-psp.time
    psp_timeind=np.argmin(abs(dct))

    dct=frame_time_num+k*res_in_hours/24-bepi.time
    bepi_timeind=np.argmin(abs(dct))

    dct=frame_time_num+k*res_in_hours/24-solo.time
    solo_timeind=np.argmin(abs(dct))

    dct=frame_time_num+k*res_in_hours/24-juice.time
    juice_timeind=np.argmin(abs(dct))
    
    dct=frame_time_num+k*res_in_hours/24-l4.time
    l4_timeind=np.argmin(abs(dct))

    dct=frame_time_num+k*res_in_hours/24-l5.time
    l5_timeind=np.argmin(abs(dct))

    
    print(psp_timeind)
    
    

    #plot all positions including text R lon lat for some 

    #white background
    if not black:
        ax.scatter(venus.lon[earth_timeind], venus.r[earth_timeind]*np.cos(venus.lat[earth_timeind]), s=symsize_planet, c='orange', alpha=1,lw=0,zorder=3)
        ax.scatter(mercury.lon[earth_timeind], mercury.r[earth_timeind]*np.cos(mercury.lat[earth_timeind]), s=symsize_planet, c='dimgrey', alpha=1,lw=0,zorder=3)
        ax.scatter(earth.lon[earth_timeind], earth.r[earth_timeind]*np.cos(earth.lat[earth_timeind]), s=symsize_planet, c='mediumseagreen', alpha=1,lw=0,zorder=3)
        ax.scatter(sta.lon[earth_timeind], sta.r[earth_timeind]*np.cos(sta.lat[earth_timeind]), s=symsize_spacecraft, c='red', marker='s', alpha=1,lw=0,zorder=3)
        ax.scatter(mars.lon[earth_timeind], mars.r[earth_timeind]*np.cos(mars.lat[earth_timeind]), s=symsize_planet, c='orangered', alpha=1,lw=0,zorder=3)

        ax.scatter(l4.lon[l4_timeind], l4.r[l4_timeind]*np.cos(l4.lat[l4_timeind]), s=symsize_planet-2, c='yellow', alpha=1,lw=0,zorder=3,marker='p')
        ax.scatter(l5.lon[l5_timeind], l5.r[l5_timeind]*np.cos(l5.lat[l5_timeind]), s=symsize_planet-2, c='yellowgreen', alpha=1,lw=0,zorder=3,marker='p')

        plt.figtext(0.95,0.75,'PSP ', color='black', ha='center',fontsize=fsize+3)
        plt.figtext(0.95,0.5,'Wind', color='mediumseagreen', ha='center',fontsize=fsize+3)
        plt.figtext(0.95,0.25,'STEREO-A', color='red', ha='center',fontsize=fsize+3)
        '''
        plt.figtext(0.9,0.9,'Mercury', color='dimgrey', ha='center',fontsize=fsize+5)
        plt.figtext(0.9	,0.8,'Venus', color='orange', ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.7,'Earth', color='mediumseagreen', ha='center',fontsize=fsize+5)
        #plt.figtext(0.9,0.7,'Mars', color='orangered', ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.6,'STEREO-A', color='red', ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.5,'Parker Solar Probe', color='black', ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.4,'Bepi Colombo', color='blue', ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.3,'Solar Orbiter', color='green', ha='center',fontsize=fsize+5)
        '''

    #black background
    if black:
        ax.scatter(venus.lon[earth_timeind], venus.r[earth_timeind]*np.cos(venus.lat[earth_timeind]), s=symsize_planet, c='orange', alpha=1,lw=0,zorder=3)
        ax.scatter(mercury.lon[earth_timeind], mercury.r[earth_timeind]*np.cos(mercury.lat[earth_timeind]), s=symsize_planet, c='grey', alpha=1,lw=0,zorder=3)
        ax.scatter(earth.lon[earth_timeind], earth.r[earth_timeind]*np.cos(earth.lat[earth_timeind]), s=symsize_planet, c='mediumseagreen', alpha=1,lw=0,zorder=3)
        ax.scatter(mars.lon[earth_timeind], mars.r[earth_timeind]*np.cos(mars.lat[earth_timeind]), s=symsize_planet, c='orangered', alpha=1,lw=0,zorder=3)

        ax.scatter(l4.lon[l4_timeind], l4.r[l4_timeind]*np.cos(l4.lat[l4_timeind]), s=symsize_planet/2, c=l4_color, alpha=1,lw=0,zorder=3,marker='p')
        ax.scatter(l5.lon[l5_timeind], l5.r[l5_timeind]*np.cos(l5.lat[l5_timeind]), s=symsize_planet/2, c=l5_color, alpha=1,lw=0,zorder=3,marker='p')
        
        ax.scatter(sta.lon[sta_timeind], sta.r[sta_timeind]*np.cos(sta.lat[sta_timeind]), s=symsize_spacecraft, c=sta_color, marker='s', alpha=1,lw=0,zorder=3)

        plt.figtext(0.9,0.9,'Mercury', color='grey', ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.8,'Venus', color='orange', ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.7,'Earth', color='mediumseagreen', ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.6,'Mars', color='orangered', ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.5,'STEREO-A', color=sta_color, ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.4,'Parker Solar Probe', color=psp_color, ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.3,'Bepi Colombo', color=bepi_color, ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.2,'Solar Orbiter', color=solo_color, ha='center',fontsize=fsize+5)
        plt.figtext(0.9,0.1,'JUICE', color=juice_color, ha='center',fontsize=fsize+5)

        
    #get 
    if add_icmes: 
        #get all indices for icmes starting before ** days of current time
        cur_icme=np.where(np.logical_and((frame_time_num+k*res_in_hours/24-ic.icme_start_time) > 0,(frame_time_num+k*res_in_hours/24-ic.icme_start_time) < time_diff_icme))[0]
        #print(cur_icme)
        


    #positions text
    f10=plt.figtext(0.01,0.93,'              R     lon     lat', fontsize=fsize+2, ha='left',color=backcolor)

    if frame=='HEEQ': earth_text='Earth: '+str(f'{earth.r[earth_timeind]:6.2f}')+str(f'{0.0:8.1f}')+str(f'{np.rad2deg(earth.lat[earth_timeind]):8.1f}')
    else: earth_text='Earth: '+str(f'{earth.r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(earth.lon[earth_timeind]):8.1f}')+str(f'{np.rad2deg(earth.lat[earth_timeind]):8.1f}')

    mars_text='Mars:  '+str(f'{mars.r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(mars.lon[earth_timeind]):8.1f}')+str(f'{np.rad2deg(mars.lat[earth_timeind]):8.1f}')
    sta_text='STA:   '+str(f'{sta.r[sta_timeind]:6.2f}')+str(f'{np.rad2deg(sta.lon[sta_timeind]):8.1f}')+str(f'{np.rad2deg(sta.lat[sta_timeind]):8.1f}')
    
    
    
    #ICMEs for Earth    
    if add_icmes:                        
            cur_icme_wind=np.where(ic.sc_insitu[cur_icme] == 'Wind')[0]            
            for i in np.arange(len(cur_icme_wind)):
                time_dist=frame_time_num+k*res_in_hours/24-ic.icme_start_time[cur_icme][cur_icme_wind]
                fadedist=1-time_dist/time_diff_icme
                ax.scatter(earth.lon[earth_timeind], earth.r[earth_timeind]*np.cos(earth.lat[earth_timeind]), s=symsize_icme, c='mediumseagreen', marker='o', alpha=fadedist,lw=0,zorder=3)

    #ICMEs for STEREO-A
    if add_icmes:                        
            cur_icme_sta=np.where(ic.sc_insitu[cur_icme] == 'STEREO-A')[0]            
            for i in np.arange(len(cur_icme_sta)):
                #get difference of icme_start_time to current frame time
                time_dist=frame_time_num+k*res_in_hours/24-ic.icme_start_time[cur_icme][cur_icme_sta]
                #define fading
                fadedist=1-time_dist/time_diff_icme
                ax.scatter(sta.lon[sta_timeind], sta.r[sta_timeind]*np.cos(sta.lat[sta_timeind]), s=symsize_icme, c=sta_color, marker='o', alpha=fadedist,lw=0,zorder=3)


    #position and text for PSP only before the kernel ends
    if np.logical_and(psp_timeind > 0,(frame_time_num+k*res_in_hours/24) < psp.time.values[-1] ):
                      
        #plot trajectory
        ax.scatter(psp.lon[psp_timeind], psp.r[psp_timeind]*np.cos(psp.lat[psp_timeind]), s=symsize_spacecraft, c=psp_color, marker='s', alpha=1,lw=0,zorder=3)
        #plot position as text
        psp_text='PSP:   '+str(f'{psp.r[psp_timeind]:6.2f}')+str(f'{np.rad2deg(psp.lon[psp_timeind]):8.1f}')+str(f'{np.rad2deg(psp.lat[psp_timeind]):8.1f}')
        f5=plt.figtext(0.01,0.78,psp_text, fontsize=fsize, ha='left',color=psp_color)
        if plot_orbit: 
            fadestart=psp_timeind-fadeind
            if  fadestart < 0: fadestart=0
            ax.plot(psp.lon[fadestart:psp_timeind+fadeind], psp.r[fadestart:psp_timeind+fadeind]*np.cos(psp.lat[fadestart:psp_timeind+fadeind]), c=psp_color, alpha=0.6,lw=1,zorder=3)
            
        if add_icmes:                        
            cur_icme_psp=np.where(ic.sc_insitu[cur_icme] == 'PSP')[0]            
            for i in np.arange(len(cur_icme_psp)):
                time_dist=frame_time_num+k*res_in_hours/24-ic.icme_start_time[cur_icme][cur_icme_psp]
                #define fading
                fadedist=1-time_dist/time_diff_icme
                ax.scatter(psp.lon[psp_timeind], psp.r[psp_timeind]*np.cos(psp.lat[psp_timeind]), s=symsize_icme, c=psp_color, marker='o', alpha=fadedist,lw=0,zorder=3)
            
            

    if np.logical_and(bepi_timeind > 0,(frame_time_num+k*res_in_hours/24) < bepi.time.values[-1] ):
        
        ax.scatter(bepi.lon[bepi_timeind], bepi.r[bepi_timeind]*np.cos(bepi.lat[bepi_timeind]), s=symsize_spacecraft, c=bepi_color, marker='s', alpha=1,lw=0,zorder=3)
        bepi_text='Bepi:   '+str(f'{bepi.r[bepi_timeind]:6.2f}')+str(f'{np.rad2deg(bepi.lon[bepi_timeind]):8.1f}')+str(f'{np.rad2deg(bepi.lat[bepi_timeind]):8.1f}')
        f6=plt.figtext(0.01,0.74,bepi_text, fontsize=fsize, ha='left',color=bepi_color)
        if plot_orbit: 
            fadestart=bepi_timeind-fadeind
            if  fadestart < 0: fadestart=0            
            ax.plot(bepi.lon[fadestart:bepi_timeind+fadeind], bepi.r[fadestart:bepi_timeind+fadeind]*np.cos(bepi.lat[fadestart:bepi_timeind+fadeind]), c=bepi_color, alpha=0.6,lw=1,zorder=3)            
        if add_icmes:                        
            cur_icme_bepi=np.where(ic.sc_insitu[cur_icme] == 'BepiColombo')[0]            
            for i in np.arange(len(cur_icme_bepi)):
                time_dist=frame_time_num+k*res_in_hours/24-ic.icme_start_time[cur_icme][cur_icme_bepi]
                #define fading
                fadedist=1-time_dist/time_diff_icme
                ax.scatter(bepi.lon[bepi_timeind], bepi.r[bepi_timeind]*np.cos(bepi.lat[bepi_timeind]), s=symsize_icme, c=bepi_color, marker='o', alpha=fadedist,lw=0,zorder=3)
            
                        
            
    #after Bepi kernel is done, mercury position if needed
    if (frame_time_num+k*res_in_hours/24) > bepi.time.values[-1]:
        
        ax.scatter(mercury.lon[earth_timeind], mercury.r[earth_timeind]*np.cos(mercury.lat[earth_timeind]), s=symsize_spacecraft, c=bepi_color, marker='s', alpha=1,lw=0,zorder=3)
        bepi_text='Bepi:   '+str(f'{mercury.r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(mercury.lon[earth_timeind]):8.1f}')+str(f'{np.rad2deg(mercury.lat[earth_timeind]):8.1f}')
        f6=plt.figtext(0.01,0.74,bepi_text, fontsize=fsize, ha='left',color=bepi_color)
        if plot_orbit: 
            fadestart=earth_timeind-fadeind
            if  fadestart < 0: fadestart=0            
            ax.plot(mercury.lon[fadestart:earth_timeind+fadeind], mercury.r[fadestart:earth_timeind+fadeind]*np.cos(mercury.lat[fadestart:earth_timeind+fadeind]), c=bepi_color, alpha=0.6,lw=1,zorder=3)
      

    if solo_timeind > 0:
        ax.scatter(solo.lon[solo_timeind], solo.r[solo_timeind]*np.cos(solo.lat[solo_timeind]), s=symsize_spacecraft, c=solo_color, marker='s', alpha=1,lw=0,zorder=3)
        solo_text='SolO:  '+str(f'{solo.r[solo_timeind]:6.2f}')+str(f'{np.rad2deg(solo.lon[solo_timeind]):8.1f}')+str(f'{np.rad2deg(solo.lat[solo_timeind]):8.1f}')
        f7=plt.figtext(0.01,0.7,solo_text, fontsize=fsize, ha='left',color=solo_color)
        if plot_orbit: 
            fadestart=solo_timeind-fadeind
            if  fadestart < 0: fadestart=0            
            ax.plot(solo.lon[fadestart:solo_timeind+fadeind], solo.r[fadestart:solo_timeind+fadeind]*np.cos(solo.lat[fadestart:solo_timeind+fadeind]), c=solo_color, alpha=0.6,lw=1,zorder=3)
         
        if add_icmes:                        
            cur_icme_solo=np.where(ic.sc_insitu[cur_icme] == 'SolarOrbiter')[0]            
            #plot current solo cmes
            for i in np.arange(len(cur_icme_solo)):
                time_dist=frame_time_num+k*res_in_hours/24-ic.icme_start_time[cur_icme][cur_icme_solo]
                #define fading
                fadedist=1-time_dist/time_diff_icme
                ax.scatter(solo.lon[solo_timeind], solo.r[solo_timeind]*np.cos(solo.lat[solo_timeind]), s=symsize_icme, c=solo_color, marker='o', alpha=fadedist,lw=0,zorder=3)

 
            
    if juice_timeind > 0:
        ax.scatter(juice.lon[juice_timeind], juice.r[juice_timeind]*np.cos(juice.lat[juice_timeind]), s=symsize_spacecraft, c=juice_color, marker='s', alpha=1,lw=0,zorder=3)
        juice_text='JUICE:  '+str(f'{juice.r[juice_timeind]:6.2f}')+str(f'{np.rad2deg(juice.lon[juice_timeind]):8.1f}')+str(f'{np.rad2deg(juice.lat[juice_timeind]):8.1f}')
        f7=plt.figtext(0.01,0.66,juice_text, fontsize=fsize, ha='left',color=juice_color)
        if plot_orbit: 
            fadestart=juice_timeind-fadeind
            if  fadestart < 0: fadestart=0            
            ax.plot(juice.lon[fadestart:juice_timeind+fadeind], juice.r[fadestart:juice_timeind+fadeind]*np.cos(juice.lat[fadestart:juice_timeind+fadeind]), c=juice_color, alpha=0.6,lw=1,zorder=3)
         
            
    f10=plt.figtext(0.01,0.9,earth_text, fontsize=fsize, ha='left',color='mediumseagreen')
    f9=plt.figtext(0.01,0.86,mars_text, fontsize=fsize, ha='left',color='orangered')
    f8=plt.figtext(0.01,0.82,sta_text, fontsize=fsize, ha='left',color=sta_color)
    
   
    #parker spiral
    if plot_parker:
        for q in np.arange(0,12):
            omega=2*np.pi/(sun_rot*60*60*24) #solar rotation in seconds
            v=400/AUkm #km/s
            r0=695000/AUkm
            r=v/omega*theta+r0*7
            if not black: 
                ax.plot(-theta+np.deg2rad(0+(360/24.47)*res_in_hours/24*k+360/12*q), r, alpha=0.4, lw=0.5,color='grey',zorder=2)
            if black: 
                ax.plot(-theta+np.deg2rad(0+(360/24.47)*res_in_hours724*k+360/12*q), r, alpha=0.7, lw=0.7,color='grey',zorder=2)

    #set axes and grid
    ax.set_theta_zero_location('E')
    #plt.thetagrids(range(0,360,45),(u'0\u00b0 '+frame+' longitude',u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'- 135\u00b0',u'- 90\u00b0',u'- 45\u00b0'), ha='right', fmt='%d',fontsize=fsize-1,color=backcolor, alpha=0.9)
    plt.thetagrids(range(0,360,45),(u'0\u00b0',u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'- 135\u00b0',u'- 90\u00b0',u'- 45\u00b0'), ha='center', fmt='%d',fontsize=fsize-1,color=backcolor, alpha=0.9,zorder=4)


    #plt.rgrids((0.10,0.39,0.72,1.00,1.52),('0.10','0.39','0.72','1.0','1.52 AU'),angle=125, fontsize=fsize,alpha=0.9, color=backcolor)
    plt.rgrids((0.1,0.3,0.5,0.7,1.0,1.3),('0.10','0.3','0.5','0.7','1.0','1.3 AU'),angle=125, fontsize=fsize-3,alpha=0.5, color=backcolor)

    #ax.set_ylim(0, 1.75) #with Mars
    ax.set_ylim(0, rmax) #Mars at 1.66 max 

    #Sun
    ax.scatter(0,0,s=100,c='yellow',alpha=1, edgecolors='black', linewidth=0.3)

    
    #signature
    plt.figtext(0.01,0.02,'Austrian Space Weather Office, GeoSphere Austria', fontsize=fsize, ha='left',color=backcolor) 
   
    
    logo = plt.imread('logo/GSA_Basislogo_NegativAufMidnightGreen_RGB_XXS.png')
    newax = fig.add_axes([0.91,0.91,0.08,0.08], anchor='NE', zorder=1)
    newax.imshow(logo)
    newax.axis('off')
    
    
    
    ############ date
    #plot text for date extra so it does not move    
    xclock=0.65
    yclock=0.06
    f1=plt.figtext(xclock,yclock,frame_time_str[0:4],  ha='center',color=backcolor,fontsize=fsize+6)
    #month
    f2=plt.figtext(xclock+0.04,yclock,frame_time_str[5:7], ha='center',color=backcolor,fontsize=fsize+6)
    #day
    f3=plt.figtext(xclock+0.08,yclock,frame_time_str[8:10], ha='center',color=backcolor,fontsize=fsize+6)
    #hours
    f4=plt.figtext(xclock+0.12,yclock,frame_time_str[11:13], ha='center',color=backcolor,fontsize=fsize+6)
        
    plt.tight_layout()


    #save figure
    framestr = '%05i' % (k)  
    filename=outputdirectory+'/pos_anim_'+framestr+'.jpg'  
    if k==0: print(filename)
    plt.savefig(filename,dpi=dpisave,facecolor=fig.get_facecolor(), edgecolor='none')
    #plt.clf()
    plt.close('all')


# In[15]:


plt.close('all')

####################### SETTINGS

#Coordinate System
#frame='HCI'
frame='HEEQ'
print(frame)

#sidereal solar rotation rate
if frame=='HCI': sun_rot=24.47
#synodic
if frame=='HEEQ': sun_rot=26.24

AUkm=149597870.7   

#black background on or off
#black=True
black=True

#animation settings
plot_orbit=True
#plot_orbit=False
plot_parker=False
#plot_parker=False

high_res_mode=False


# for adding ICMECAT indicators for ICMEs
add_icmes=False
#this is the time how long an ICME is visible as a faded circle
time_diff_icme=7 
#how big the symbol for the ICME is
symsize_icme=450
        


def datetime_range(start, end, resolution_hours):
    current = start
    while current < end:
        yield current
        current += timedelta(hours=resolution_hours)

time_array = [dt for dt in datetime_range(t_start, t_end, res_in_hours)]

k_all=np.size(time_array)
counter=[i for i in range(k_all)]



print(time_array[0])
print(time_array[-1])
print('number of frames:',k_all)
#print(time_array)






############################ 

#animation start time in matplotlib format
frame_time_num=parse_time(t_start).plot_date

sns.set_context('talk')

if not black: sns.set_style('darkgrid'),#{'grid.linestyle': ':', 'grid.color': '.35'}) 

if black: sns.set_style('white',{'grid.linestyle': ':', 'grid.color': '.35'})   

# animation settings 
fsize=13

#how long the trajectory tracks are
#1 hour resolution of the position file
fadeind=int(24*60)

print('days for trajectory',fadeind/24)

symsize_planet=110
symsize_spacecraft=80



#for parker spiral   
theta=np.arange(0,np.deg2rad(180),0.01)


# ### single processing

#print()
#print('make animation')
#print()

#k_all=1

#for debugging
k_all=1000
make_frame(1)
#for i in np.arange(1,10,1):
#    make_frame(i)
    
#os.system(ffmpeg_path+'ffmpeg -r 25 -i '+str(outputdirectory)+'/pos_anim_%05d.jpg -b 5000k \
#    -r 25 '+str(animdirectory)+'/positions_punch.mp4 -y -loglevel quiet')   


# In[ ]:


# ### Multiprocessing


print()
print('make animation')
print()

#number of processes depends on your machines memory; check with command line "top"
#how much memory is used by all your processes


print('Using multiprocessing, nr of cores',mp.cpu_count(), \
      'with nr of processes used: ',used)


#define pool using fork and number of processes
pool=mp.get_context('fork').Pool(processes=used)
# Map the worker function onto the parameters    
t0 = time.time()
pool.map(make_frame, counter) #or use apply_async?,imap
pool.close()
pool.join()     
t1 = time.time()


print('time in min: ',np.round((t1-t0)/60))
print('plots done, frames saved in ',outputdirectory)


os.system(ffmpeg_path+'ffmpeg -r 25 -i '+str(outputdirectory)+'/pos_anim_%05d.jpg -b 5000k \
    -r 30 '+str(animdirectory)+'/'+movie_filename+'.mp4 -y -loglevel quiet')    
print('movie done, saved in ',animdirectory)



# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





#data.py
#read and save data for heliocats
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy
import scipy.io
import scipy.signal
import copy
import os
import sys
import glob
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.image as mpimg
import time
import datetime
import urllib
from urllib.request import urlopen
import requests
import json
from sunpy.time import parse_time
from sunpy.coordinates import HeliocentricInertial, HeliographicStonyhurst, HeliocentricEarthEcliptic
import astropy
import astropy.units as u
from astropy.constants import au
from astropy.time import Time, TimeDelta
import pickle
import cdflib
import spiceypy
import astrospice
from astrospice.net.reg import RemoteKernel, RemoteKernelsBase
from numba import njit
import math
import h5py
from bs4 import BeautifulSoup 



#MIT LICENSE
#Copyright 2020-2025, Christian Moestl, Rachel L. Bailey, Emma E. Davies, Eva Weiler
#Permission is hereby granted, free of charge, to any person obtaining a copy of this 
#software and associated documentation files (the "Software"), to deal in the Software
#without restriction, including without limitation the rights to use, copy, modify, 
#merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
#permit persons to whom the Software is furnished to do so, subject to the following 
#conditions:
#The above copyright notice and this permission notice shall be included in all copies 
#or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
#CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
#OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.




#################################################################################



def get_noaa_xray(noaa_path,data_path,xraypickle, xraypickle2):

    
    
    ################ get primary GOES satellite data
    
    #xray data
    xray='https://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json'

    ################ same for secondary GOES satellite

    #xray2 data
    xray2='https://services.swpc.noaa.gov/json/goes/secondary/xrays-7-day.json'

    #flare list
    #https://services.swpc.noaa.gov/json/goes/primary/xray-flares-7-day.json

    print('download NOAA Xrays')
    datestr=str(datetime.datetime.utcnow().strftime("%Y-%m-%d"))
    print(datestr+' UTC')

    try: 
        urllib.request.urlretrieve(xray, noaa_path+'xray/xrays-7-day_'+datestr+'.json')
        print(noaa_path+'xray/xray-7-day_'+datestr+'.json')
        success=True
    except urllib.error.URLError as e:
        print(' ', xray,' ',e.reason)
        
    try: 
        urllib.request.urlretrieve(xray2, noaa_path+'xray/xrays2-7-day_'+datestr+'.json')
        print(noaa_path+'xray/xray2-7-day_'+datestr+'.json')
        success2=True
    except urllib.error.URLError as e:
        print(' ', xray2,' ',e.reason)
        
        

    if success: 
        #convert json to pickle

        xrayfile=open(noaa_path+'xray/xrays-7-day_'+datestr+'.json','r')

        xd = json.loads(xrayfile.read())

        #this is in the time data    
        #{'time_tag': '2024-03-05T11:43:00Z',
        #'satellite': 16,
        #'flux': 1.2437139957910404e-06,
        #'observed_flux': 1.259415967069799e-06,
        #'electron_correction': 1.5701953515190326e-08,
        #'electron_contaminaton': False,
        #'energy': '0.1-0.8nm'},

        #energy 1 0.1-0.8nm indices
        xdc1=np.zeros(len(xd),dtype=[('time',object),('flux', float)]) 
        xdc1=xdc1.view(np.recarray)   

        #energy 2 0.05-0.4nm indices
        xdc2=np.zeros(len(xd),dtype=[('time',object),('flux', float)]) 
        xdc2=xdc2.view(np.recarray)   

        #indices for energy 1    
        e1_ind=np.arange(1,len(xd),2)

        #indices for energy 2
        e2_ind=np.arange(0,len(xd),2)

        for i in e1_ind:        
            xdc1[i].time=parse_time(xd[i]['time_tag']).datetime
            xdc1[i].flux=xd[i]['flux'] 

        #get all data points from flux, and use those indices to remove the 0s from time array
        datind=np.where(xdc1.flux>1e-7)[0]
        xdc1=xdc1[datind]

        for i in e2_ind:
            xdc2[i].time=parse_time(xd[i]['time_tag']).datetime
            xdc2[i].flux=xd[i]['flux']         

        datind=np.where(xdc2.flux>1e-9)
        xdc2=xdc2[datind]

        pickle.dump([xdc1,xdc2], open(data_path+xraypickle, "wb"))
        print('saved as',data_path+xraypickle,'pickle')    

        
        
    if success2: 
        #convert json to pickle

        xrayfile2=open(noaa_path+'xray/xrays2-7-day_'+datestr+'.json','r')

        xd2 = json.loads(xrayfile2.read())

        #energy 1 0.1-0.8nm indices
        #xray data 2 channel 1
        xd2c1=np.zeros(len(xd2),dtype=[('time',object),('flux', float)]) 
        xd2c1=xd2c1.view(np.recarray)   

        #energy 2 0.05-0.4nm indices
        #xray data 2 channel 2
        xd2c2=np.zeros(len(xd2),dtype=[('time',object),('flux', float)]) 
        xd2c2=xd2c2.view(np.recarray)   

        #indices for energy 1    
        e21_ind=np.arange(1,len(xd2),2)

        #indices for energy 2
        e22_ind=np.arange(0,len(xd2),2)

        for i in e21_ind:        
            xd2c1[i].time=parse_time(xd2[i]['time_tag']).datetime
            xd2c1[i].flux=xd2[i]['flux'] 

        #get all data points from flux, and use those indices to remove the 0s from time array
        datind2=np.where(xd2c1.flux>1e-7)[0]
        xd2c1=xd2c1[datind2]

        for i in e22_ind:
            xd2c2[i].time=parse_time(xd2[i]['time_tag']).datetime
            xd2c2[i].flux=xd2[i]['flux']         

        datind2=np.where(xd2c2.flux>1e-9)
        xd2c2=xd2c2[datind2]
        
               
        pickle.dump([xd2c1,xd2c2], open(data_path+xraypickle2, "wb"))
        print('saved as',data_path+xraypickle2,'pickle')    
        
        
        
        
        
        
        
        
        
        

 

def get_sdo_realtime_image(data_path_sun):
    """Downloads latest SDO image."""

    get_193_file=False
    get_193_date=False
    

    #-------------------------- AIA

  
    
    
    
    ########### 193 Image with timestamp
    
    sdo_latest='https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_0193.jpg'
    #sdo_latest='https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_0193pfss.jpg'
    
   
    
    #download the file
    try: 
        urllib.request.urlretrieve(sdo_latest,data_path_sun+'latest_1024_0193.jpg')
        print('saved ',data_path_sun+'latest_1024_0193.jpg')       
        get_193_file=True

    except urllib.error.URLError as e:
        print('Failed downloading 193 image', sdo_latest,' ',e)

 
    #get the file creation date
    sdo_date= 'https://sdo.gsfc.nasa.gov/assets/img/latest/times0193.txt'
     
    try: 
        with urllib.request.urlopen(sdo_date) as response:
            file_contents = response.read().decode('utf-8')  # Decode the bytes to a string using UTF-8
        date1=file_contents[6:14]
        time1=file_contents[15:]

        #date1=file_contents[6:10]
        timestamp=date1[0:4]+'-'+date1[4:6]+'-'+date1[6:8]+' '+time1[0:2]+':'+time1[2:4]+ ' UT'
        print(timestamp)
        get_193_date=True
        
    except urllib.error.URLError as e:        
        print('Failed downloading 193 date file', sdo_latest,' ',e)

    #if the date file and the image file could both be downloaded add the timestamp
    if get_193_date and get_193_file:
        image_path = data_path_sun+'latest_1024_0193.jpg'
        image = mpimg.imread(image_path)


        font_size = 24
        text_color = 'white'
        text_bg_color = 'black'

        # Create a figure and axes
        fig=plt.figure(figsize=(1024 / 80,1024 / 80), dpi=80)
        ax = fig.add_subplot(111)

        # Display the image
        ax.imshow(image)

        # Add the timestamp as text
        ax.text(20, 50, timestamp, fontsize=font_size, color=text_color)
        ax.text(770, 50, 'SDO/AIA 193', fontsize=font_size, color=text_color)

        # Remove the axis ticks and labels
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')


        # Save the modified image

        output_path = data_path_sun+'latest_1024_0193.jpg'
        plt.savefig(output_path,bbox_inches='tight', pad_inches=0,dpi=110)

        # Show the modified image
        #plt.show()

    
    
    
    #--------------- HMI
    
    get_HMI_file=False
    get_HMI_date=False
    
    sdo_latest='https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_HMIB.jpg'        
    try: 
        urllib.request.urlretrieve(sdo_latest,data_path_sun+'latest_1024_HMIB.jpg')
        get_HMI_file=True
    except urllib.error.URLError as e:
        print('Failed downloading ', sdo_latest,' ',e)
        
    print('saved ',data_path_sun+'latest_1024_HMIB.jpg')    

    
 
    #get the file creation date
    sdo_date= 'https://sdo.gsfc.nasa.gov/assets/img/latest/timesHMIB.txt'
     
    try: 
        with urllib.request.urlopen(sdo_date) as response:
            file_contents = response.read().decode('utf-8')  # Decode the bytes to a string using UTF-8
        date1=file_contents[6:14]
        time1=file_contents[15:]

        #date1=file_contents[6:10]
        timestamp=date1[0:4]+'-'+date1[4:6]+'-'+date1[6:8]+' '+time1[0:2]+':'+time1[2:4]+ ' UT'
        print(timestamp)
        get_HMI_date=True
        
    except urllib.error.URLError as e:        
        print('Failed downloading 193 date file', sdo_latest,' ',e)

    #if the date file and the image file could both be downloaded add the timestamp
    if get_HMI_date and get_HMI_file:
        image_path = data_path_sun+'latest_1024_HMIB.jpg'
        image = mpimg.imread(image_path)


        font_size = 24
        text_color = 'white'
        text_bg_color = 'black'

        # Create a figure and axes
        fig=plt.figure(figsize=(1024 / 80,1024 / 80), dpi=80)
        ax = fig.add_subplot(111)

        # Display the image
        ax.imshow(image,cmap='gray')

        # Add the timestamp as text
        ax.text(20, 50, timestamp, fontsize=font_size, color=text_color)
        ax.text(770, 50, 'SDO/HMI', fontsize=font_size, color=text_color)

        # Remove the axis ticks and labels
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')


        # Save the modified image

        output_path = data_path_sun+'latest_1024_HMIB.jpg'
        plt.savefig(output_path,bbox_inches='tight', pad_inches=0,dpi=110)

        # Show the modified image
        #plt.show()
        
        
        
        
        #download 2 more HMI files
        
        hmicc='https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_HMIIC.jpg'

        try: 
            urllib.request.urlretrieve(hmicc,data_path_sun+'latest_1024_HMIIC.jpg')
            print('saved ',data_path_sun+'latest_1024_HMIIC.jpg')        

        except urllib.error.URLError as e:
            print('Failed downloading HMICC image', hmic,' ',e)

        hmipfss='https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_HMIBpfss.jpg'
        
        try: 
            urllib.request.urlretrieve(hmipfss,data_path_sun+'latest_1024_HMIBpfss.jpg')
            print('saved ',data_path_sun+'latest_1024_HMIBpfss.jpg')        

        except urllib.error.URLError as e:
            print('Failed downloading HMIBpfss image', sdo_latest,' ',e)


        
    
        
    '''    
    #convert to png
    #check if ffmpeg is available locally in the folder or systemwide
    if os.path.isfile('ffmpeg'):
        os.system('./ffmpeg -i latest_1024_0193.jpg latest_1024_0193.png -loglevel quiet -y')
        ffmpeg_avail=True
        logger.info('downloaded SDO latest_1024_0193.jpg converted to png')
        os.system('rm latest_1024_0193.jpg')
    else:
        os.system('ffmpeg -i latest_1024_0193.jpg latest_1024_0193.png -loglevel quiet -y')
        os.system('rm latest_1024_0193.jpg')
    '''    















####################################### get new data ####################################


#general



@njit
def cart2sphere(x,y,z):
    r = np.sqrt(x**2+ y**2 + z**2)           
    theta = np.arctan2(z,np.sqrt(x**2+ y**2))
    phi = np.arctan2(y,x)                    
    return (r, theta, phi)
    


@njit
def sphere2cart(r,lat,lon):
    x = r * np.sin( lat ) * np.cos( lon )
    y = r * np.sin( lat ) * np.sin( lon )
    z = r * np.cos( lat )
    return (x, y,z)
    




def sphere2cart_emma(r, lat, lon):
    x = r*np.cos(lat*(np.pi/180))*np.cos(lon*(np.pi/180))
    y = r*np.cos(lat*(np.pi/180))*np.sin(lon*(np.pi/180))
    z = r*np.sin(lat*(np.pi/180))
    r_au = r/1.495978707E8
    return x.value, y.value, z.value, r_au.value

def cart2sphere_emma(x,y,z):
    r = np.sqrt(x**2+ y**2 + z**2) /1.495978707E8         
    theta = np.arctan2(z,np.sqrt(x**2+ y**2)) * 360 / 2 / np.pi
    phi = np.arctan2(y,x) * 360 / 2 / np.pi                   
    return (r, theta, phi)




















########### to be done

def convert_GSM_to_RTN(sc_in):

    sc=copy.deepcopy(sc_in)

    print('GSM to RTN for NOAA L1 GSM data to RTN for comparison to STEREO-A needed TBD')    
    
    #first GSM to GSE
    #GSE to HEE sign flip
    #HEE to RTN
    
    return sc



##########################################

def convert_GSE_to_GSM(sc_in):

    sc=copy.deepcopy(sc_in)

    print('conversion GSE to GSM')                                

    jd=np.zeros(len(sc))
    mjd=np.zeros(len(sc))

    for i in np.arange(0,len(sc)):
        #get all dates right
        jd[i]=parse_time(sc.time[i]).jd
        mjd[i]=float(int(jd[i]-2400000.5))  #use modified julian date 
        
        #define T00 and UT
        T00=(mjd[i]-51544.5)/36525.0
        dobj=sc.time[i]
        UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours    
        
        #define position of geomagnetic pole in GEO coordinates
        pgeo=78.8+4.283*((mjd[i]-46066)/365.25)*0.01 #in degrees
        lgeo=289.1-1.413*((mjd[i]-46066)/365.25)*0.01 #in degrees
        
        #GEO vector
        Qg=[np.cos(pgeo*np.pi/180)*np.cos(lgeo*np.pi/180), np.cos(pgeo*np.pi/180)*np.sin(lgeo*np.pi/180), np.sin(pgeo*np.pi/180)]
    
        #CREATE T1, T00, UT is known from above
        zeta=(100.461+36000.770*T00+15.04107*UT)*np.pi/180
        T1=np.matrix([[np.cos(zeta), np.sin(zeta),  0], [-np.sin(zeta) , np.cos(zeta) , 0], [0,  0,  1]]) #angle for transpose
        
        #lambda_sun in Hapgood, equation 5, here in degrees
        LAMBDA=280.460+36000.772*T00+0.04107*UT
        M=357.528+35999.050*T00+0.04107*UT
        lt2=(LAMBDA+(1.915-0.0048*T00)*np.sin(M*np.pi/180)+0.020*np.sin(2*M*np.pi/180))*np.pi/180
        
        #CREATE T2, LAMBDA, M, lt2 known from above
        ##################### lamdbda und Z
        t2z=np.matrix([[np.cos(lt2), np.sin(lt2),  0], [-np.sin(lt2) , np.cos(lt2) , 0], [0,  0,  1]])
        et2=(23.439-0.013*T00)*np.pi/180
        
        ###################### epsilon und x
        t2x=np.matrix([[1,0,0],[0,np.cos(et2), np.sin(et2)], [0, -np.sin(et2), np.cos(et2)]])
        T2=np.dot(t2z,t2x)  #equation 4 in Hapgood 1992
        #matrix multiplications   
        T2T1t=np.dot(T2,np.matrix.transpose(T1))
        ################
        Qe=np.dot(T2T1t,Qg) #Q=T2*T1^-1*Qq
        psigsm=np.arctan(Qe.item(1)/Qe.item(2)) #arctan(ye/ze) in between -pi/2 to +pi/2
        
        T3=np.matrix([[1,0,0],[0,np.cos(-psigsm), np.sin(-psigsm)], [0, -np.sin(-psigsm), np.cos(-psigsm)]])
        b_gse=np.matrix([[sc.bx[i]],[sc.by[i]],[sc.bz[i]]]) 
        b_gsm=np.dot(T3,b_gse)   #equation 6 in Hapgood
        sc.bx[i]=b_gsm[0]
        sc.by[i]=b_gsm[1]
        sc.bz[i]=b_gsm[2]
    
    print('conversion GSE to GSM done')   
    
    return sc




def convert_RTN_to_GSE_sta_l1(sc_in):

    sc=copy.deepcopy(sc_in)

    print('conversion RTN to GSE') 
    
    heeq_bx=np.zeros(len(sc))
    heeq_by=np.zeros(len(sc))
    heeq_bz=np.zeros(len(sc))
    
    jd=np.zeros(len(sc))
    mjd=np.zeros(len(sc))
    
     
    ########## first RTN to HEEQ 
    
    #go through all data points
    for i in np.arange(0,len(sc)):
        #print(sc.time[i])
        #HEEQ vectors
        X_heeq=[1,0,0]
        Y_heeq=[0,1,0]
        Z_heeq=[0,0,1]

        #normalized X RTN vector
        Xrtn=[sc.x[i],sc.y[i],sc.z[i]]/np.linalg.norm([sc.x[i],sc.y[i],sc.z[i]])
        #solar rotation axis at 0, 0, 1 in HEEQ
        Yrtn=np.cross(Z_heeq,Xrtn)/np.linalg.norm(np.cross(Z_heeq,Xrtn))
        Zrtn=np.cross(Xrtn, Yrtn)/np.linalg.norm(np.cross(Xrtn, Yrtn))
        
        #project into new system
        heeq_bx[i]=np.dot(np.dot(sc.bx[i],Xrtn)+np.dot(sc.by[i],Yrtn)+np.dot(sc.bz[i],Zrtn),X_heeq)
        heeq_by[i]=np.dot(np.dot(sc.bx[i],Xrtn)+np.dot(sc.by[i],Yrtn)+np.dot(sc.bz[i],Zrtn),Y_heeq)
        heeq_bz[i]=np.dot(np.dot(sc.bx[i],Xrtn)+np.dot(sc.by[i],Yrtn)+np.dot(sc.bz[i],Zrtn),Z_heeq)
    
        #then HEEQ to GSE
        jd[i]=parse_time(sc.time[i]).jd
        mjd[i]=float(int(jd[i]-2400000.5))  #use modified julian date  
        
        #then lambda_sun
        T00=(mjd[i]-51544.5)/36525.0
        dobj=sc.time[i]
        UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   
        LAMBDA=280.460+36000.772*T00+0.04107*UT
        M=357.528+35999.050*T00+0.04107*UT
        
        #lt2 is lambdasun in Hapgood, equation 5, here in rad
        lt2=(LAMBDA+(1.915-0.0048*T00)*np.sin(M*np.pi/180)+0.020*np.sin(2*M*np.pi/180))*np.pi/180
        
        #note that some of these equations are repeated later for the GSE to GSM conversion
        S1=np.matrix([[np.cos(lt2+np.pi), np.sin(lt2+np.pi),  0], [-np.sin(lt2+np.pi) , np.cos(lt2+np.pi) , 0], [0,  0,  1]])
        #create S2 matrix with angles with reversed sign for transformation HEEQ to HAE
        omega_node=(73.6667+0.013958*((mjd[i]+3242)/365.25))*np.pi/180 #in rad
        S2_omega=np.matrix([[np.cos(-omega_node), np.sin(-omega_node),  0], [-np.sin(-omega_node) , np.cos(-omega_node) , 0], [0,  0,  1]])
        inclination_ecl=7.25*np.pi/180
        S2_incl=np.matrix([[1,0,0],[0,np.cos(-inclination_ecl), np.sin(-inclination_ecl)], [0, -np.sin(-inclination_ecl), np.cos(-inclination_ecl)]])
        #calculate theta
        theta_node=np.arctan(np.cos(inclination_ecl)*np.tan(lt2-omega_node)) 

        #quadrant of theta must be opposite lt2 - omega_node Hapgood 1992 end of section 5   
        #get lambda-omega angle in degree mod 360   
        lambda_omega_deg=np.mod(lt2-omega_node,2*np.pi)*180/np.pi
        x = np.cos(np.deg2rad(lambda_omega_deg))
        y = np.sin(np.deg2rad(lambda_omega_deg))
        
        #get theta_node in deg
        theta_node_deg=theta_node*180/np.pi
        x_theta = np.cos(theta_node)
        y_theta = np.sin(theta_node)
        #if in same quadrant, then theta_node = theta_node +pi   
        #if abs(lambda_omega_deg-theta_node_deg) < 180: theta_node=theta_node+np.pi ------> diese Zeile mit den if-Schleifen drunter ersetzen.
        if (x>=0 and y>=0):
            if (x_theta>=0 and y_theta>=0): theta_node = theta_node - np.pi
            elif (x_theta<=0 and y_theta<=0): theta_node = theta_node
            elif (x_theta>=0 and y_theta<=0): theta_node = theta_node - np.pi/2
            elif (x_theta<=0 and y_theta>=0): theta_node = np.pi+(theta_node-np.pi/2)
            
        elif (x<=0 and y<=0):
            if (x_theta>=0 and y_theta>=0): theta_node = theta_node
            elif (x_theta<=0 and y_theta<=0): theta_node = theta_node + np.pi
            elif (x_theta>=0 and y_theta<=0): theta_node = theta_node + np.pi/2
            elif (x_theta<=0 and y_theta>=0): theta_node = theta_node-np.pi/2
            
        elif (x>=0 and y<=0):
            if (x_theta>=0 and y_theta>=0): theta_node = theta_node + np.pi/2
            elif (x_theta<=0 and y_theta<=0): theta_node = theta_node-np.pi/2
            elif (x_theta>=0 and y_theta<=0): theta_node = theta_node + np.pi
            elif (x_theta<=0 and y_theta>=0): theta_node = theta_node

        elif (x<0 and y>0):
            if (x_theta>=0 and y_theta>=0): theta_node = theta_node - np.pi/2
            elif (x_theta<=0 and y_theta<=0): theta_node = theta_node + np.pi/2
            elif (x_theta>=0 and y_theta<=0): theta_node = theta_node
            elif (x_theta<=0 and y_theta>=0): theta_node = theta_node -np.pi
        
        
        S2_theta=np.matrix([[np.cos(-theta_node), np.sin(-theta_node),  0], [-np.sin(-theta_node) , np.cos(-theta_node) , 0], [0,  0,  1]])

        #make S2 matrix
        S2=np.dot(np.dot(S2_omega,S2_incl),S2_theta)
        #this is the matrix S2^-1 x S1
        HEEQ_to_HEE_matrix=np.dot(S1, S2)
        #convert HEEQ components to HEE
        HEEQ=np.matrix([[heeq_bx[i]],[heeq_by[i]],[heeq_bz[i]]]) 
        HEE=np.dot(HEEQ_to_HEE_matrix,HEEQ)
        #change of sign HEE X / Y to GSE is needed
        sc.bx[i]=-HEE[0]
        sc.by[i]=-HEE[1]
        sc.bz[i]=HEE[2]
    
    print('conversion RTN to GSE done') 
   
    return sc




def convert_HEEQ_to_RTN(sc_in):
    '''
    for all spacecraft
    '''

    print('conversion HEEQ to RTN start for both V and B')                                
    
    sc=copy.deepcopy(sc_in)
    
    sc_len=len(sc.bt)

    #unit vectors of HEEQ basis
    heeq_x=[1,0,0]
    heeq_y=[0,1,0]
    heeq_z=[0,0,1]

    #project into new system RTN
    for i in np.arange(0,sc_len):

        #make unit vectors of RTN in basis of HEEQ
        rtn_r=[sc.x[i],sc.y[i],sc.z[i]]/np.linalg.norm([sc.x[i],sc.y[i],sc.z[i]])
        rtn_t=np.cross(heeq_z,rtn_r)/np.linalg.norm(np.cross(heeq_z,rtn_r))
        rtn_n=np.cross(rtn_r,rtn_t)/np.linalg.norm(np.cross(rtn_r,rtn_t))

        sc.bx[i]=sc_in.bx[i]*np.dot(heeq_x,rtn_r)+sc_in.by[i]*np.dot(heeq_y,rtn_r)+sc_in.bz[i]*np.dot(heeq_z,rtn_r)
        sc.by[i]=sc_in.bx[i]*np.dot(heeq_x,rtn_t)+sc_in.by[i]*np.dot(heeq_y,rtn_t)+sc_in.bz[i]*np.dot(heeq_z,rtn_t)
        sc.bz[i]=sc_in.bx[i]*np.dot(heeq_x,rtn_n)+sc_in.by[i]*np.dot(heeq_y,rtn_n)+sc_in.bz[i]*np.dot(heeq_z,rtn_n)
   
        sc.vx[i]=sc_in.vx[i]*np.dot(heeq_x,rtn_r)+sc_in.vy[i]*np.dot(heeq_y,rtn_r)+sc_in.vz[i]*np.dot(heeq_z,rtn_r)
        sc.vy[i]=sc_in.vx[i]*np.dot(heeq_x,rtn_t)+sc_in.vy[i]*np.dot(heeq_y,rtn_t)+sc_in.vz[i]*np.dot(heeq_z,rtn_t)
        sc.vz[i]=sc_in.vx[i]*np.dot(heeq_x,rtn_n)+sc_in.vy[i]*np.dot(heeq_y,rtn_n)+sc_in.vz[i]*np.dot(heeq_z,rtn_n)

    print('conversion HEEQ to RTN done')                                

    return sc



def convert_HEEQ_to_RTN_mag(sc_in):
    '''
    for all spacecraft 
    '''

    print('conversion HEEQ to RTN start for B only')                                
    
    sc=copy.deepcopy(sc_in)
    
    sc_len=len(sc.bt)

    #unit vectors of HEEQ basis
    heeq_x=[1,0,0]
    heeq_y=[0,1,0]
    heeq_z=[0,0,1]

    #project into new system RTN
    for i in np.arange(0,sc_len):

        #make unit vectors of RTN in basis of HEEQ
        rtn_r=[sc.x[i],sc.y[i],sc.z[i]]/np.linalg.norm([sc.x[i],sc.y[i],sc.z[i]])
        rtn_t=np.cross(heeq_z,rtn_r)/np.linalg.norm(np.cross(heeq_z,rtn_r)) 
        rtn_n=np.cross(rtn_r,rtn_t)/np.linalg.norm(np.cross(rtn_r,rtn_t)) 

        sc.bx[i]=sc_in.bx[i]*np.dot(heeq_x,rtn_r)+sc_in.by[i]*np.dot(heeq_y,rtn_r)+sc_in.bz[i]*np.dot(heeq_z,rtn_r)
        sc.by[i]=sc_in.bx[i]*np.dot(heeq_x,rtn_t)+sc_in.by[i]*np.dot(heeq_y,rtn_t)+sc_in.bz[i]*np.dot(heeq_z,rtn_t)
        sc.bz[i]=sc_in.bx[i]*np.dot(heeq_x,rtn_n)+sc_in.by[i]*np.dot(heeq_y,rtn_n)+sc_in.bz[i]*np.dot(heeq_z,rtn_n)
   
    print('conversion HEEQ to RTN done')                                

    return sc



def convert_HEE_to_HEEQ(sc_in):
    '''
    for Wind positions: convert HEE to HAE to HEEQ
    '''

    print('conversion HEE to HEEQ')                                
    
    sc=copy.deepcopy(sc_in)
    
    jd=np.zeros(len(sc))
    mjd=np.zeros(len(sc))
        

    for i in np.arange(0,len(sc)):

        jd[i]=parse_time(sc.time[i]).jd
        mjd[i]=float(int(jd[i]-2400000.5)) #use modified julian date    
        
       
        w_hee=[sc.x[i],sc.y[i],sc.z[i]]
        
        #HEE to HAE        
        
        #define T00 and UT
        T00=(mjd[i]-51544.5)/36525.0          
        dobj=sc.time[i]
        UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   

        #lambda_sun in Hapgood, equation 5, here in rad
        M=np.radians(357.528+35999.050*T00+0.04107*UT)
        LAMBDA=280.460+36000.772*T00+0.04107*UT        
        lambda_sun=np.radians( (LAMBDA+(1.915-0.0048*T00)*np.sin(M)+0.020*np.sin(2*M)) )
        
        #S-1 Matrix equation 12 hapgood 1992, change sign in lambda angle for inversion HEE to HAE instead of HAE to HEE
        c, s = np.cos(-(lambda_sun+np.radians(180))), np.sin(-(lambda_sun+np.radians(180)))
        Sm1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
        w_hae=np.dot(Sm1,w_hee)

        #HAE to HEEQ
        
        iota=np.radians(7.25)
        omega=np.radians((73.6667+0.013958*((mjd[i]+3242)/365.25)))                      
        theta=np.arctan(np.cos(iota)*np.tan(lambda_sun-omega))                       
                      
    
        #quadrant of theta must be opposite lambda_sun minus omega; Hapgood 1992 end of section 5   
        #get lambda-omega angle in degree mod 360 and theta in degrees
        lambda_omega_deg=np.mod(np.degrees(lambda_sun)-np.degrees(omega),360)
        theta_node_deg=np.degrees(theta)


        ##if the 2 angles are close to similar, so in the same quadrant, then theta_node = theta_node +pi           
        if np.logical_or(abs(lambda_omega_deg-theta_node_deg) < 1, abs(lambda_omega_deg-360-theta_node_deg) < 1): theta=theta+np.pi                                                                                                          
        
        #rotation around Z by theta
        c, s = np.cos(theta), np.sin(theta)
        S2_1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))

        #rotation around X by iota  
        iota=np.radians(7.25)
        c, s = np.cos(iota), np.sin(iota)
        S2_2 = np.array(( (1,0,0), (0,c, s), (0, -s, c)) )
                
        #rotation around Z by Omega  
        c, s = np.cos(omega), np.sin(omega)
        S2_3 = np.array( ((c,s, 0), (-s, c, 0), (0, 0, 1)) )
        
        #matrix multiplication to go from HAE to HEEQ components                
        [x_heeq,y_heeq,z_heeq]=np.dot(  np.dot(   np.dot(S2_1,S2_2),S2_3), w_hae) 
        
        sc.x[i]=x_heeq
        sc.y[i]=y_heeq
        sc.z[i]=z_heeq

    
    print('HEE to HEEQ done')  
    
    return sc



  
















def convert_E2K_to_HEE(sc_in,kernels_path):
    '''
    for Bepi magnetic field components: convert E2K to HEE by projection
    '''

    
    sc=copy.deepcopy(sc_in)
    
    print('conversion E2K to HEE')     
    
    
    #old heliopy code   #get Earth trajectory in ECLIPJ2000
    #earth=spice.Trajectory('399')  #399 for Earth, not barycenter (because of moon)
    #earth.generate_positions(sc_in.time,'Sun','ECLIPJ2000')
    #earth.change_units(astropy.units.AU)  
    #earth_x_ej2=earth.x.value
    #earth_y_ej2=earth.y.value
    #earth_z_ej2=earth.z.value

    
    #get x y z in ECLIPJ2000 in au
    #for sc_in.time
    
    
    
    ###### use spice kernels and spiceypy    

    #spice_files = [kernels_path+'generic/naif0012.tls', 
    #               kernels_path+'generic/de430.bsp',
    #               kernels_path+'generic/pck00010.tpc']    
    #spiceypy.furnsh(spice_files)
    
    
    generic_path = kernels_path+'generic/'  
    generic_kernels = os.listdir(generic_path)
    for kernel in generic_kernels:
        spiceypy.furnsh(os.path.join(generic_path, kernel))
    
    time_spice = [spiceypy.str2et(t.strftime('%Y-%m-%d %H:%M')) for t in sc_in.time]

    positions_earth, lightTimes_earth = spiceypy.spkezr('Earth', time_spice, 'ECLIPJ2000', 'NONE', 'Sun')

    pos_earth = np.array(positions_earth)[:, :3] * u.km

    earth_x_ej2 = pos_earth[:, 0].to(u.au).value
    earth_y_ej2 = pos_earth[:, 1].to(u.au).value
    earth_z_ej2 = pos_earth[:, 2].to(u.au).value
    
   
    
    

    
    sc_len=len(sc_in.bt)

    #unit vectors of EJ2 basis
    ej_x=[1,0,0]
    ej_y=[0,1,0]
    ej_z=[0,0,1]

    #project into new system HEE
    for i in np.arange(0,sc_len):   
        #make unit vectors of HEE in basis of EJ2
        hee_x=[earth_x_ej2[i],earth_y_ej2[i],earth_z_ej2[i]]/np.linalg.norm([earth_x_ej2[i],earth_y_ej2[i],earth_z_ej2[i]])
        hee_z=[0,0,1]
        hee_y=np.cross(hee_z,hee_x)

        sc.bx[i]=sc_in.bx[i]*np.dot(ej_x,hee_x)+sc_in.by[i]*np.dot(ej_y,hee_x)+sc_in.bz[i]*np.dot(ej_z,hee_x)
        sc.by[i]=sc_in.bx[i]*np.dot(ej_x,hee_y)+sc_in.by[i]*np.dot(ej_y,hee_y)+sc_in.bz[i]*np.dot(ej_z,hee_y)
        sc.bz[i]=sc_in.bx[i]*np.dot(ej_x,hee_z)+sc_in.by[i]*np.dot(ej_y,hee_z)+sc_in.bz[i]*np.dot(ej_z,hee_z)
        
    print('conversion E2K to HEE done')  

    return sc




 
def convert_GSE_to_HEEQ(sc_in):
    '''
    for Wind magnetic field components: convert GSE to HEE to HAE to HEEQ
    '''
    
    
    sc=copy.deepcopy(sc_in)

    print('conversion GSE to HEEQ start for both B and V')                                

    jd=np.zeros(len(sc))
    mjd=np.zeros(len(sc))
        

    for i in np.arange(0,len(sc)):

        jd[i]=parse_time(sc.time[i]).jd
        mjd[i]=float(int(jd[i]-2400000.5)) #use modified julian date    
        
        #GSE to HEE
        #Hapgood 1992 rotation by 180 degrees, or simply change sign in bx by    
        #rotangle=np.radians(180)
        #c, s = np.cos(rotangle), np.sin(rotangle)
        #T1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
        #[bx_hee,by_hee,bz_hee]=T1[sc.bx[i],sc.by[i],sc.bz[i]]        
        b_hee=[-sc.bx[i],-sc.by[i],sc.bz[i]]
        v_hee=[-sc.vx[i],-sc.vy[i],sc.vz[i]]        
        
        
        #HEE to HAE        
        
        #define T00 and UT
        T00=(mjd[i]-51544.5)/36525.0          
        dobj=sc.time[i]
        UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   

        #lambda_sun in Hapgood, equation 5, here in rad
        M=np.radians(357.528+35999.050*T00+0.04107*UT)
        LAMBDA=280.460+36000.772*T00+0.04107*UT        
        lambda_sun=np.radians( (LAMBDA+(1.915-0.0048*T00)*np.sin(M)+0.020*np.sin(2*M)) )
        
        #S-1 Matrix equation 12 hapgood 1992, change sign in lambda angle for inversion HEE to HAE instead of HAE to HEE
        c, s = np.cos(-(lambda_sun+np.radians(180))), np.sin(-(lambda_sun+np.radians(180)))
        Sm1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
        b_hae=np.dot(Sm1,b_hee)
        v_hae=np.dot(Sm1,v_hee)

        #HAE to HEEQ
        
        iota=np.radians(7.25)
        omega=np.radians((73.6667+0.013958*((mjd[i]+3242)/365.25)))                      
        theta=np.arctan(np.cos(iota)*np.tan(lambda_sun-omega))                       
                      
    
        #quadrant of theta must be opposite lambda_sun minus omega; Hapgood 1992 end of section 5   
        #get lambda-omega angle in degree mod 360 and theta in degrees
        lambda_omega_deg=np.mod(np.degrees(lambda_sun)-np.degrees(omega),360)
        theta_node_deg=np.degrees(theta)


        ##if the 2 angles are close to similar, so in the same quadrant, then theta_node = theta_node +pi           
        if np.logical_or(abs(lambda_omega_deg-theta_node_deg) < 1, abs(lambda_omega_deg-360-theta_node_deg) < 1): theta=theta+np.pi                                                                                                          
        
        #rotation around Z by theta
        c, s = np.cos(theta), np.sin(theta)
        S2_1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))

        #rotation around X by iota  
        iota=np.radians(7.25)
        c, s = np.cos(iota), np.sin(iota)
        S2_2 = np.array(( (1,0,0), (0,c, s), (0, -s, c)) )
                
        #rotation around Z by Omega  
        c, s = np.cos(omega), np.sin(omega)
        S2_3 = np.array( ((c,s, 0), (-s, c, 0), (0, 0, 1)) )
        
        #matrix multiplication to go from HAE to HEEQ components                
        [bx_heeq,by_heeq,bz_heeq]=np.dot(  np.dot(   np.dot(S2_1,S2_2),S2_3), b_hae) 
        
        sc.bx[i]=bx_heeq
        sc.by[i]=by_heeq
        sc.bz[i]=bz_heeq
        
        [vx_heeq,vy_heeq,vz_heeq]=np.dot(  np.dot(   np.dot(S2_1,S2_2),S2_3), v_hae) 
        
        sc.vx[i]=vx_heeq
        sc.vy[i]=vy_heeq
        sc.vz[i]=vz_heeq

        
   
    print('conversion GSE to HEEQ done')                                
    return sc

    

    
    
    
    
    
    
    
    
    
    
    
    
################################################################# OMNI2


def omni_loader():
   '''
   downloads all omni2 data into the "data" folder
   '''

   #if overwrite>0: 
   #      print('download OMNI2 again')
   #      if os.path.exists('data/omni2_all_years.dat'): os.remove('data/omni2_all_years.dat')
  
   #if not os.path.exists('data/omni2_all_years.dat'):
      #see http://omniweb.gsfc.nasa.gov/html/ow_data.html
   print('load OMNI2 .dat into "data" directory from')
   omni2_url='https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
   print(omni2_url)
   try: urllib.request.urlretrieve(omni2_url, 'data/omni2_all_years.dat')
   except urllib.error.URLError as e:
         print(' ', omni2_url,' ',e.reason)
         sys.exit()


def save_omni_data(path,file):
    '''
    save variables from OMNI2 dataset as pickle

    documentation https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2.text
    omni2_url='https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
    '''        
  
    print('start omni')
    
    omni_loader()
    #check how many rows exist in this file
    f=open('data/omni2_all_years.dat')
    dataset= len(f.readlines())
    
    #make array
    o=np.zeros(dataset,dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bygsm', float),('bzgsm', float),('bt', float),\
                ('vt', float),('np', float),('tp', float),('alpha', float),\
                ('dst', float),('kp', float),('spot', float),\
                ('ae', float),('ap', float),('f107', float),\
                ('pcn', float),('al', float),('au', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])
                
    o=o.view(np.recarray)  
    print(dataset, ' datapoints')   #for reading data from OMNI file
    
    j=0
    with open('data/omni2_all_years.dat') as f:
        for line in f:
            line = line.split() # to deal with blank 
            
            #time - need to convert from year doy hour to datetime object
            o.time[j]=datetime.datetime(int(line[0]), 1, 1) + datetime.timedelta(int(line[1]) - 1) \
                              + datetime.timedelta(hours=int(line[2])) 

            #25 is bulkspeed F6.0, in km/s
            o.vt[j]=line[24]
            if o.vt[j] == 9999: o.vt[j]=np.NaN
            
            #24 in file, index 23 proton density /ccm
            o.np[j]=line[23]
            if o.np[j] == 999.9: o.np[j]=np.NaN
            
            #23 in file, index 22 Proton temperature  /ccm
            o.tp[j]=line[22]
            if o.tp[j] == 9999999.: o.tp[j]=np.NaN

            #28 in file, index 27 alpha to proton ratio
            o.alpha[j]=line[27]
            if o.alpha[j] == 9.999: o.alpha[j]=np.NaN
           
            #9 is total B  F6.1 also fill ist 999.9, in nT
            o.bt[j]=line[9]
            if o.bt[j] == 999.9: o.bt[j]=np.NaN

            #GSE components from 13 to 15, so 12 to 14 index, in nT
            o.bx[j]=line[12]
            if o.bx[j] == 999.9: o.bx[j]=np.NaN
            o.by[j]=line[13]
            if o.by[j] == 999.9: o.by[j]=np.NaN
            o.bz[j]=line[14]
            if o.bz[j] == 999.9: o.bz[j]=np.NaN
          
            #GSM
            o.bygsm[j]=line[15]
            if o.bygsm[j] == 999.9: o.bygsm[j]=np.NaN
          
            o.bzgsm[j]=line[16]
            if o.bzgsm[j] == 999.9: o.bzgsm[j]=np.NaN    
          

            o.kp[j]=line[38]
            if o.kp[j] == 99: o.kp[j]=np.nan
           
            o.spot[j]=line[39]
            if o.spot[j] ==  999: o.spot[j]=np.nan

            o.dst[j]=line[40]
            if o.dst[j] == 99999: o.dst[j]=np.nan

            o.ae[j]=line[41]
            if o.ae[j] ==  9999: o.ae[j]=np.nan
            
            o.ap[j]=line[49]
            if o.ap[j] ==  999: o.ap[j]=np.nan
            
            o.f107[j]=line[50]
            if o.f107[j] ==  999.9  : o.f107[j]=np.nan
            
            o.pcn[j]=line[51]
            if o.pcn[j] ==  999.9  : o.pcn[j]=np.nan

            o.al[j]=line[52]
            if o.al[j] ==   99999  : o.al[j]=np.nan

            o.au[j]=line[53]
            if o.au[j] ==  99999 : o.au[j]=np.nan
          
          
            j=j+1     
            
    print('position start')
    

    frame = HeliographicStonyhurst()

    coords_earth = astrospice.generate_coords('Earth', o.time)
    coords_earth = coords_earth.transform_to(frame)
    o.r = coords_earth.radius.to(u.au).value
    o.lon = coords_earth.lon.value #degrees
    o.lat = coords_earth.lat.value

    [o.x, o.y, o.z] = sphere2cart(o.r, np.deg2rad(-o.lat+90), np.deg2rad(o.lon))
    
    
    header='Near Earth OMNI2 1 hour solar wind and geomagnetic indices data since 1963. ' + \
    'Obtained from https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/  '+ \
    'Timerange: '+o.time[0].strftime("%Y-%b-%d %H:%M")+' to '+o.time[-1].strftime("%Y-%b-%d %H:%M")+'. '+\
    'The data are available in a numpy recarray, fields can be accessed by o.time, o.bx, o.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(o.size)+'. '+\
    'For units and documentation see: https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2.text, the '+\
    'heliospheric position of Earth was added and is given in x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_omni_data .'+\
    'By C. Möstl, A. J. Weiss. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    
    pickle.dump([o,header], open(path+file, 'wb') )
    
    print('done omni')
    print()

    










############# BepiColombo


##Bepi POSITION FUNCTIONS with spiceypy, preferred method
#https://naif.jpl.nasa.gov/pub/naif/BEPICOLOMBO/kernels/spk/
def bepi_furnish(kernels_path):
    """Main"""
    bepi_path = kernels_path+'bepi/'
  
    #put the latest file here manually
    bepi_kernel = 'bc_mpo_fcp_00199_20181020_20270407_v02.bsp'
    spiceypy.furnsh(os.path.join(bepi_path, bepi_kernel))
    print(bepi_kernel)

    #for kernel in solo_kernels:
    #    spiceypy.furnsh(os.path.join(solo_path, kernel))
    #spiceypy.furnsh(solo_kernels)

    
    #bepi_kernels = astrospice.SPKKernel(bepi_path+'bc_mpo_fcp_00199_20181020_20270407_v02.bsp')    
    
    generic_path = kernels_path+'generic/'  
    generic_kernels = os.listdir(generic_path)
    print(generic_kernels)
    for kernel in generic_kernels:
        spiceypy.furnsh(os.path.join(generic_path, kernel))


def get_bepi_pos(t,kernels_path):
    if spiceypy.ktotal('ALL') < 1:
        bepi_furnish(kernels_path)
    pos = spiceypy.spkpos("BEPICOLOMBO MPO", spiceypy.datetime2et(t), "HEEQ", "NONE", "SUN")[0]
    r, lat, lon = cart2sphere_emma(pos[0],pos[1],pos[2])
    
    position = t, pos[0], pos[1], pos[2], r, lat, lon
    return position


def get_bepi_positions(time_series,kernels_path):
    positions = []
    for t in time_series:
        position = get_bepi_pos(t,kernels_path)
        positions.append(position)
    df_positions = pd.DataFrame(positions, columns=['time', 'x', 'y', 'z', 'r', 'lat', 'lon'])
    return df_positions

    




def create_bepi_pickle(start_date,end_date,finalfile,bepi_path, sensor, kernels_path):
    
    
    t_start = start_date
    t_end = end_date
    
    #create an array with 1 minute resolution between t start and end
    time1 = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*24*60))] 
    time1_mat=mdates.date2num(time1) 
    print(time1[0],time1[-1])
    print(len(time1))
    
    
    #make array with length time1
    bepi=np.zeros(len(time1),dtype=[('time',object),('bx', float),('by', float),\
                    ('bz', float),('bt', float),('np', float),('vt', float),('vx', float),\
                          ('vy', float),('vz', float),\
                          ('tp', float),('x', float),('y', float),('z', float),\
                    ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    bepi = bepi.view(np.recarray)  
    
    
    bepi.time=time1
    
    
    if sensor=='outbound':
        print('data from outbound directory')
        tag='ob'
        bepi_path=bepi_path+'cruise_ob/'        
        print(bepi_path)
        
    if sensor=='inbound':
        print('data from inbound directory')
        tag='ib'
        bepi_path=bepi_path+'cruise_ib/'        
        print(bepi_path)
        
        

    #Tab File looks like
    #SC POSITION X COMPONENT, ECLIPJ2000 COORDINATES. VALUE IS GIVEN IN KM. THE SC POSITION IS STATED AS THE DISTANCE TO SUN.
    #MAGNETIC FIELD X COMPONENT, 1-SECOND-AVERAGED CALIBRATED DATA,   ECLIPJ2000 COORDINATES, OB SENSOR. VALUE IS GIVEN IN NANOTESLA.

    ######## create array of filenames needed for the current timerange
    
    #array only for days for time1
    time1_days = [ t_start + datetime.timedelta(days=1*n) for n in range(int ((t_end - t_start).days))] 

    data_filename=[]
    
    for i in np.arange(0,len(time1_days)):
    
        #make zeros for day and month
        month=str(np.char.zfill(str(time1_days[i].month), 2))
        day=str(np.char.zfill(str(time1_days[i].day), 2))
        data_filename.append('mag_der_sc_'+tag+'_a001_e2k_00000_'+str(time1_days[i].year)+month+day+'.tab')
        
    
    
    #make one big array of the raw data and then interpolate
    
    #make array with length for all days with 1 second resolution
    bepi_raw=np.zeros(len(data_filename)*24*60*60,dtype=[('time',float),('bx', float),('by', float),\
                    ('bz', float)])   

    #convert to recarray
    bepi_raw = bepi_raw.view(np.recarray)  
    
    
    counter=0
        
    for filename in data_filename:
        print(bepi_path+filename)


        #check if file exists
        try: 
            braw=np.loadtxt(bepi_path+filename,delimiter=",",dtype=[('time','<U30'),('tsc','<U30'),('x', float),\
                                                    ('y', float),('z', float),('bx', float),('by', float),\
                                                    ('bz', float),('t1', float),('t2', float),('t3', float)])
            
            braw = braw.view(np.recarray) 
            
            #write in bepi_raw
            bepi_raw[counter:counter+len(braw)].time=mdates.date2num(Time(braw.time).datetime)
            bepi_raw.bx[counter:counter+len(braw)]=braw.bx
            bepi_raw.by[counter:counter+len(braw)]=braw.by
            bepi_raw.bz[counter:counter+len(braw)]=braw.bz
            
            counter=counter+len(braw)
    
        except:

            print(filename,'does not exist or something went wrong when reading it')

    bepi_raw=bepi_raw[0:counter]
 
    bepi.bx = np.interp(time1_mat, bepi_raw.time, bepi_raw.bx)
    bepi.by = np.interp(time1_mat, bepi_raw.time, bepi_raw.by)
    bepi.bz = np.interp(time1_mat, bepi_raw.time, bepi_raw.bz)
    bepi.bt=np.sqrt(bepi.bx**2+bepi.by**2+bepi.bz**2)                                                      
        
    
    #get positions
    print('positions start')
    
    #with astrospice
    ###### coords_bepi,times_bepi=get_positions_bepi(t_start,t_end,10)
    #use this
    #coords_bepi,times_bepi=get_positions_bepi(time1)
    
    #with spiceypy  
    bepi_furnish(kernels_path)    
    coords_bepi = get_bepi_positions(time1,kernels_path)
    print('positions end')
    
    print(u.au)
            
    bepi.r=coords_bepi.r
    bepi.lon=coords_bepi.lon
    bepi.lat=coords_bepi.lat
    bepi.x=coords_bepi.x #/(au.value*1e-3) to AU
    bepi.y=coords_bepi.y #/(au.value*1e-3)
    bepi.z=coords_bepi.z #/(au.value*1e-3)

    #bepi.x,bepi.y,bepi.z,r_au=sphere2cart_emma(coords_bepi.r,  bepi.lat, bepi.lon)
    
    bepi.np=np.nan    
    bepi.vt=np.nan
    bepi.vx=np.nan
    bepi.vy=np.nan
    bepi.vz=np.nan
    bepi.tp=np.nan
    
    
        
        
    #indicate inbound or outbound sensor    
    header='Bepi Colombo cruise magnetic field data, from TU Braunschweig, sensor '+sensor+'. ' + \
    'Timerange: '+bepi.time[0].strftime("%Y-%b-%d %H:%M")+' to '+bepi.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(bepi.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by bepi.time, bepi.bx, bepi.bt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(bepi.size)+'. '+\
    'Units are btxyz [nT, E2K], heliospheric position x/y/z [km] r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with heliocats/create_bepi_pickle  '+\
    'By Christian Möstl, Eva Weiler, Emma Davies. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([bepi,header], open(finalfile, "wb"))
    
    print('file written',finalfile)
    
    print('bepi update done')
    print()
    

def bepi_rtn_header(bepi, sensor):
    
    header='Bepi Colombo cruise magnetic field data, from TU Braunschweig, sensor '+sensor+'. ' + \
    'Timerange: '+bepi.time[0].strftime("%Y-%b-%d %H:%M")+' to '+bepi.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(bepi.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by bepi.time, bepi.bx, bepi.bt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(bepi.size)+'. '+\
    'Units are btxyz [nT, RTN], heliospheric position x/y/z [km] r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with heliocats/create_bepi_pickle  '+\
    'By Christian Möstl, Eva Weiler, Emma Davies. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    return header







###################### Wind





def wind_download_ascii(start_year, wind_path):

    #download MFI and SWE in ASCII from SPDF
        
    print('downloading Wind ascii data to ', wind_path)
        
    read_data_end_year=datetime.datetime.utcnow().year
    read_data_end_month=datetime.datetime.utcnow().month
    
    wind_years_strings=[]
    for j in np.arange(start_year,read_data_end_year+1):
        wind_years_strings.append(str(j))

    print('for years', wind_years_strings)
    
  
    ######## MFI

    mfi_url='https://spdf.gsfc.nasa.gov/pub/data/wind/mfi/ascii/1min_ascii/'
    print(mfi_url)
    
    
    #for all years
    for i in np.arange(0,len(wind_years_strings)-1):    

        for k in np.arange(1,13):    

            a=str(k).zfill(2) #add leading zeros            
            filewind=wind_years_strings[i]+a+'_wind_mag_1min.asc'
            print(filewind)
            try: urllib.request.urlretrieve(mfi_url+filewind, wind_path+'mfi_1min_ascii/'+filewind)
            except urllib.error.URLError as e:
                print(' ', mfi_url,' ',e.reason) 

    
    
    #for latest year
    for k in np.arange(1,read_data_end_month):    

        a=str(k).zfill(2) #add leading zeros
        filewind=wind_years_strings[-1]+a+'_wind_mag_1min.asc'
        print(filewind)
        try: urllib.request.urlretrieve(mfi_url+filewind, wind_path+'mfi_1min_ascii/'+filewind)
        except urllib.error.URLError as e:
            print(' ', mfi_url,' ',e.reason) 
                      

    
    ############## SWE
        
    swe_url='https://spdf.gsfc.nasa.gov/pub/data/wind/swe/ascii/swe_kp_unspike/'
    print(swe_url)
        
    for i in np.arange(0,len(wind_years_strings)):    

        filewind='wind_kp_unspike'+wind_years_strings[i]+'.txt'
        print(filewind)
        try: urllib.request.urlretrieve(swe_url+filewind, wind_path+'swe_92sec_ascii/'+filewind)
        except urllib.error.URLError as e:
            print(' ', swe_url,' ',e.reason) 


      
    
    
    






def save_wind_data_ascii(start_date,end_date,path,finalfile,coord):
    
    '''
    description of data sources used in this function
    
    SWE 92 sec
    https://spdf.gsfc.nasa.gov/pub/data/wind/swe/ascii/swe_kp_unspike
    
    MFI 1 min    
    https://spdf.gsfc.nasa.gov/pub/data/wind/mfi/ascii/1min_ascii/
    
    examples:
    
    https://spdf.gsfc.nasa.gov/pub/data/wind/swe/ascii/swe_kp_unspike/wind_kp_unspike1996.txt    
    https://spdf.gsfc.nasa.gov/pub/data/wind/mfi/ascii/1min_ascii/201908_wind_mag_1min.asc
    '''
    
    wind_data_path_mag=path+'mfi_1min_ascii/'
    wind_data_path_pla=path+'swe_92sec_ascii/'

    
    t_start = start_date
    t_end = end_date
    
    #create an array with 2 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=2*n) for n in range(int ((t_end - t_start).days*30*24))]  
    time_mat=mdates.date2num(time) 
    
    
    read_data_end_year=datetime.datetime.utcnow().year
    read_data_end_month=datetime.datetime.utcnow().month
        
    
            
    #############mfi
    
    wind_years_strings=[]
    for j in np.arange(start_date.year,end_date.year+1):
        wind_years_strings.append(str(j))

    print(wind_years_strings)
    #array for 40 years
    win_mag=np.zeros(60*24*365*40,dtype=[('time',object),('bx', float),('by', float),\
                    ('bz', float),('bt', float),('np', float),('vt', float),('vx', float),\
                          ('vy', float),('vz', float),\
                          ('tp', float),('x', float),('y', float),('z', float),\
                    ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    win_mag = win_mag.view(np.recarray)  


    counter=0

    #previous years
    for i in np.arange(0,len(wind_years_strings)-1):   
        

        for k in np.arange(1,13):    

            a=str(k).zfill(2) #add leading zeros

            file=wind_data_path_mag+wind_years_strings[i]+a+'_wind_mag_1min.asc' 
            print(file)
            #get data from file
            mfi_data=np.genfromtxt(file,dtype="i8,i8,i8,i8,i8,i8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,")


            #put data in array
            for p in np.arange(0,len(mfi_data)):
                #time
                win_mag.time[p+counter]=datetime.datetime(mfi_data[p][0],mfi_data[p][1],mfi_data[p][2],mfi_data[p][3],mfi_data[p][4],mfi_data[p][5])
                win_mag.bx[p+counter]=mfi_data[p][6]
                
                if coord=='GSE' or coord=='HEEQ':
                    win_mag.by[p+counter]=mfi_data[p][7]
                    win_mag.bz[p+counter]=mfi_data[p][8]
                    win_mag.bt[p+counter]=mfi_data[p][11]

                if coord=='GSM':
                    win_mag.by[p+counter]=mfi_data[p][9]
                    win_mag.bz[p+counter]=mfi_data[p][10]
                    win_mag.bt[p+counter]=mfi_data[p][11]

    
                    
                #gse position
                win_mag.x[p+counter]=mfi_data[p][17]
                win_mag.y[p+counter]=mfi_data[p][18]
                win_mag.z[p+counter]=mfi_data[p][19]

            counter=counter+len(mfi_data)    
            
            
    #current year
    for k in np.arange(1,read_data_end_month-1):    

        a=str(k).zfill(2) #add leading zeros

        file=wind_data_path_mag+wind_years_strings[-1]+a+'_wind_mag_1min.asc' 
        print(file)
        #get data from file
        mfi_data=np.genfromtxt(file,dtype="i8,i8,i8,i8,i8,i8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,")


        #put data in array
        for p in np.arange(0,len(mfi_data)):
                #time
                win_mag.time[p+counter]=datetime.datetime(mfi_data[p][0],mfi_data[p][1],mfi_data[p][2],mfi_data[p][3],mfi_data[p][4],mfi_data[p][5])
                win_mag.bx[p+counter]=mfi_data[p][6]
                
                if coord=='GSE' or coord=='HEEQ':
                    win_mag.by[p+counter]=mfi_data[p][7]
                    win_mag.bz[p+counter]=mfi_data[p][8]
                    win_mag.bt[p+counter]=mfi_data[p][11]
                    
                    

                if coord=='GSM':
                    win_mag.by[p+counter]=mfi_data[p][9]
                    win_mag.bz[p+counter]=mfi_data[p][10]
                    win_mag.bt[p+counter]=mfi_data[p][11]

    
                    
                #gse position
                win_mag.x[p+counter]=mfi_data[p][17]
                win_mag.y[p+counter]=mfi_data[p][18]
                win_mag.z[p+counter]=mfi_data[p][19]

        counter=counter+len(mfi_data)             
            

    #cutoff        
    win_mag2=win_mag[0:counter]    
    win_time2=mdates.date2num(win_mag2.time)

    #set missing data to nan
    win_mag2.bt[np.where(win_mag2.bt == -1e31)]=np.nan  
    win_mag2.bx[np.where(win_mag2.bx == -1e31)]=np.nan  
    win_mag2.by[np.where(win_mag2.by == -1e31)]=np.nan  
    win_mag2.bz[np.where(win_mag2.bz == -1e31)]=np.nan  
    
    win_mag2.x[np.where(win_mag2.x == -1e31)]=np.nan  
    win_mag2.y[np.where(win_mag2.y == -1e31)]=np.nan  
    win_mag2.z[np.where(win_mag2.z == -1e31)]=np.nan  
    
    
    
    ############################swe

    
    #array for 40 years
    win_swe=np.zeros(60*24*365*40,dtype=[('time',object),('bx', float),('by', float),\
                    ('bz', float),('bt', float),('np', float),('vt', float),('vx', float),\
                          ('vy', float),('vz', float),\
                          ('tp', float),('x', float),('y', float),('z', float),\
                    ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    win_swe = win_swe.view(np.recarray)  

    
    
    counter=0


    for i in np.arange(0,len(wind_years_strings)):    

            file=wind_data_path_pla+'wind_kp_unspike'+wind_years_strings[i]+'.txt'
            print(file)
            swe_data=np.genfromtxt(file,dtype="i8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8")

            firstjan=mdates.date2num(datetime.datetime(int(wind_years_strings[i]),1,1))-1

            for p in np.arange(0,len(swe_data)):

                #datenum
                win_swe.time[p+counter]=firstjan+swe_data[p][1]
                win_swe.vt[p+counter]=swe_data[p][2]
                win_swe.vx[p+counter]=swe_data[p][3]
                win_swe.vy[p+counter]=swe_data[p][4]
                win_swe.vz[p+counter]=swe_data[p][5]                
                win_swe.tp[p+counter]=swe_data[p][6]
                win_swe.np[p+counter]=swe_data[p][7]            

            counter=counter+len(swe_data) 


   
    win_swe=win_swe[0:counter]    
    win_swe_time=np.array(win_swe.time,dtype='float')


    win_swe.vt[np.where(win_swe.vt == 99999.9)]=np.nan
    win_swe.vx[np.where(win_swe.vx == 99999.9)]=np.nan
    win_swe.vy[np.where(win_swe.vy == 99999.9)]=np.nan
    win_swe.vz[np.where(win_swe.vz == 99999.9)]=np.nan
    win_swe.tp[np.where(win_swe.tp ==9999999.)]=np.nan  
    win_swe.np[np.where(win_swe.np ==999.99)]=np.nan  


    
    ###################### interpolate
    
     
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, win_time2, win_mag2.bx )
    by = np.interp(time_mat, win_time2, win_mag2.by )
    bz = np.interp(time_mat, win_time2, win_mag2.bz  )
    bt = np.sqrt(bx**2+by**2+bz**2)
        
    den = np.interp(time_mat, win_swe_time, win_swe.np)
    vt = np.interp(time_mat, win_swe_time,win_swe.vt)
    vx = np.interp(time_mat, win_swe_time,win_swe.vx)
    vy = np.interp(time_mat, win_swe_time,win_swe.vy)
    vz = np.interp(time_mat, win_swe_time,win_swe.vz)
    tp = np.interp(time_mat, win_swe_time,win_swe.tp)
    
    
    #make array
    win=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('np', float),('vt', float),('vx', float),('vy', float),('vz', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    win = win.view(np.recarray)  

        
    #fill with data
    win.time=time
    win.bx=bx
    win.by=by
    win.bz=bz 
    win.bt=bt

    win.np=den    
    win.vt=vt
    win.vx=vx
    win.vy=vy
    win.vz=vz
    
    win.tp=tp
    

    print('position start') 
    
    au_km = au.to('km').value

    #interpolate the GSE position over full data range
    x_gse = np.interp(time_mat,  win_time2, win_mag2.x)*6378.1/au_km
    y_gse = np.interp(time_mat,  win_time2, win_mag2.y)*6378.1/au_km
    z_gse = np.interp(time_mat,  win_time2, win_mag2.z)*6378.1/au_km
    
        
    coords_earth = astrospice.generate_coords('Earth', time)
    #earth_hee_lon, earth_hee_lat, earth_hee_r = coords_earth.transform_to(HeliocentricEarthEcliptic())  #returns degrees, km
    coords_earth_hee=coords_earth.transform_to(HeliocentricEarthEcliptic())  #returns degrees, km
    #change angles to -90/90 latitude, distance to AU
    [earth_hee_x, earth_hee_y, earth_hee_z] = sphere2cart(coords_earth_hee.distance.value/au_km,\
                                                              np.deg2rad(-coords_earth_hee.lat.value+90) ,np.deg2rad(coords_earth_hee.lon.value))

    win.x=earth_hee_x-x_gse  #both values are in au, sign reversed for xy for GSE compared to HEE
    win.y=earth_hee_y-y_gse
    win.z=earth_hee_z+z_gse
    
    #convert position to HEEQ
    win = convert_HEE_to_HEEQ(win)
    
    [win.r, win.lat, win.lon]=cart2sphere(win.x,win.y,win.z)       
    win.lat=np.degrees(win.lat)
    win.lon=np.degrees(win.lon)

    print('position end ')
    
    
    
    
    ############ spike removal
    
    spike_removal=0
    
    if spike_removal > 0:
            
        #plasma    
        win.np[np.where(win.np> 500)]=1000000
        #get rid of all single spikes with scipy signal find peaks
        peaks, properties = scipy.signal.find_peaks(win.np, height=500,width=(1, 250))
        #go through all of them and set to nan according to widths
        for i in np.arange(len(peaks)):
            #get width of current peak
            width=int(np.ceil(properties['widths']/2)[i])
            #remove data
            win.np[peaks[i]-width-2:peaks[i]+width+2]=np.nan

        win.tp[np.where(win.tp> 1e8)]=1e11
        #get rid of all single spikes with scipy signal find peaks
        peaks, properties = scipy.signal.find_peaks(win.tp, height=1e8,width=(1, 250))
        #go through all of them and set to nan according to widths
        for i in np.arange(len(peaks)):
            #get width of current peak
            width=int(np.ceil(properties['widths']/2)[i])
            #remove data
            win.tp[peaks[i]-width-2:peaks[i]+width+2]=np.nan

        win.vt[np.where(win.vt> 3000)]=1e11
        #get rid of all single spikes with scipy signal find peaks
        peaks, properties = scipy.signal.find_peaks(win.vt, height=1e8,width=(1, 250))
        #go through all of them and set to nan according to widths
        for i in np.arange(len(peaks)):
            #get width of current peak
            width=int(np.ceil(properties['widths']/2)[i])
            #remove data
            win.vt[peaks[i]-width-2:peaks[i]+width+2]=np.nan




        #magnetic field    
        peaks, properties = scipy.signal.find_peaks(win.bt, prominence=30,width=(1, 10))
        #go through all of them and set to nan according to widths
        for i in np.arange(len(peaks)):
            #get width of current peak
            width=int(np.ceil(properties['widths'])[i])
            #remove data
            win.bt[peaks[i]-width-5:peaks[i]+width+5]=np.nan    

        peaks, properties = scipy.signal.find_peaks(abs(win.bx), prominence=30,width=(1, 10))
        for i in np.arange(len(peaks)):
            width=int(np.ceil(properties['widths'])[i])
            win.bx[peaks[i]-width-5:peaks[i]+width+5]=np.nan    

        peaks, properties = scipy.signal.find_peaks(abs(win.by), prominence=30,width=(1, 10))
        for i in np.arange(len(peaks)):
            width=int(np.ceil(properties['widths'])[i])
            win.by[peaks[i]-width-5:peaks[i]+width+5]=np.nan    

        peaks, properties = scipy.signal.find_peaks(abs(win.bz), prominence=30,width=(1, 10))
        for i in np.arange(len(peaks)):
            width=int(np.ceil(properties['widths'])[i])
            win.bz[peaks[i]-width-5:peaks[i]+width+5]=np.nan    



        #manual spike removal for magnetic field
        if t_start < datetime.datetime(2018, 7, 19, 16, 25):    
            if t_end > datetime.datetime(2018, 7, 19, 16, 25):         

                remove_start=datetime.datetime(2018, 7, 19, 16, 25)
                remove_end=datetime.datetime(2018, 7, 19, 17, 35)
                remove_start_ind=np.where(remove_start<win.time)[0][0]
                remove_end_ind=np.where(remove_end<win.time)[0][0] 

                win.bt[remove_start_ind:remove_end_ind]=np.nan
                win.bx[remove_start_ind:remove_end_ind]=np.nan
                win.by[remove_start_ind:remove_end_ind]=np.nan
                win.bz[remove_start_ind:remove_end_ind]=np.nan

        if t_start < datetime.datetime(2018, 8, 29, 19, 00):    
            if t_end > datetime.datetime(2018, 8, 29, 19, 00):         

                remove_start=datetime.datetime(2018, 8, 29, 19, 00)
                remove_end=datetime.datetime(2018,8, 30, 5, 00)
                remove_start_ind=np.where(remove_start<win.time)[0][0]
                remove_end_ind=np.where(remove_end<win.time)[0][0] 

                win.bt[remove_start_ind:remove_end_ind]=np.nan
                win.bx[remove_start_ind:remove_end_ind]=np.nan
                win.by[remove_start_ind:remove_end_ind]=np.nan
                win.bz[remove_start_ind:remove_end_ind]=np.nan


        if t_start < datetime.datetime(2019, 8, 8, 22, 45):    
            if t_end > datetime.datetime(2019, 8, 8, 22, 45):         

                remove_start=datetime.datetime(2019, 8, 8, 22, 45)
                remove_end=datetime.datetime(2019,   8, 9, 17, 00)
                remove_start_ind=np.where(remove_start<win.time)[0][0]
                remove_end_ind=np.where(remove_end<win.time)[0][0] 

                win.bt[remove_start_ind:remove_end_ind]=np.nan
                win.bx[remove_start_ind:remove_end_ind]=np.nan
                win.by[remove_start_ind:remove_end_ind]=np.nan
                win.bz[remove_start_ind:remove_end_ind]=np.nan

        if t_start < datetime.datetime(2019, 8, 21, 22, 45):    
            if t_end > datetime.datetime(2019, 8, 21, 22, 45):         

                remove_start=datetime.datetime(2019, 8, 20, 18, 0)
                remove_end=datetime.datetime(2019,   8, 21, 12, 0)
                remove_start_ind=np.where(remove_start<win.time)[0][0]
                remove_end_ind=np.where(remove_end<win.time)[0][0] 

                win.bt[remove_start_ind:remove_end_ind]=np.nan
                win.bx[remove_start_ind:remove_end_ind]=np.nan
                win.by[remove_start_ind:remove_end_ind]=np.nan
                win.bz[remove_start_ind:remove_end_ind]=np.nan            

        if t_start < datetime.datetime(2019, 8, 21, 22, 45):    
            if t_end > datetime.datetime(2019, 8, 21, 22, 45):         

                remove_start=datetime.datetime(2019, 8, 22, 1, 0)
                remove_end=datetime.datetime(2019,   8, 22, 9, 0)
                remove_start_ind=np.where(remove_start<win.time)[0][0]
                remove_end_ind=np.where(remove_end<win.time)[0][0] 

                win.bt[remove_start_ind:remove_end_ind]=np.nan
                win.bx[remove_start_ind:remove_end_ind]=np.nan
                win.by[remove_start_ind:remove_end_ind]=np.nan
                win.bz[remove_start_ind:remove_end_ind]=np.nan            
   
    
     
    #convert magnetic field to SCEQ
    #if coord=='HEEQ':
    #    win=convert_GSE_to_HEEQ(win)
   
    
        
    ###################### 

    
    header='Wind magnetic field (MAG instrument) and plasma data (SWE), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/wind/ ascii files  '+ \
    'Timerange: '+win.time[0].strftime("%Y-%b-%d %H:%M")+' to '+win.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(win.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by win.time, win.bx, win.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(win.size)+'. '+\
    'Units are btxyz [nT, '+coord+'], vtxyz [km s^-1,'+coord+'], np[cm^-3], tp [K], heliospheric position x/y/z [km] r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with heliocats/save_wind_data_ascii.  '+\
    'By Christian Möstl, Emma Davies, Eva Weiler. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    pickle.dump([win,header], open(finalfile, "wb"))
    
    print('file written',finalfile)
    
    
    print('wind update done')
    print()
    

def wind_heeq_header(data):
    
    
    header_heeq='Wind magnetic field (MAG instrument) and plasma data (SWE), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/wind/ ascii files  '+ \
    'Timerange: '+data.time[0].strftime("%Y-%b-%d %H:%M")+' to '+data.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(data.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by win.time, win.bx, win.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(data.size)+'. '+\
    'Units are btxyz [nT, HEEQ], vtxyz [km s^-1, HEEQ], np[cm^-3], tp [K], heliospheric position x/y/z [km] r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with heliocats/save_wind_data_ascii.  '+\
    'By Christian Möstl, Emma Davies, Eva Weiler. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    return header_heeq



def wind_rtn_header(data):
    
    
    header_rtn='Wind magnetic field (MAG instrument) and plasma data (SWE), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/wind/ ascii files  '+ \
    'Timerange: '+data.time[0].strftime("%Y-%b-%d %H:%M")+' to '+data.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(data.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by win.time, win.bx, win.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(data.size)+'. '+\
    'Units are btxyz [nT, RTN], vtxyz [km s^-1, RTN], np[cm^-3], tp [K], heliospheric position x/y/z [km] r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with heliocats/save_wind_data_ascii.  '+\
    'By Christian Möstl, Emma Davies, Eva Weiler. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    return header_rtn



      

      
  

    
    
    
    
    
    
    
    
    















################## PSP


def download_pspmag_1min(start_timestamp,end_timestamp, psp_path):
    
    
    
    ### fix here that the server timeout may happen when there is no output for too long when new files are searched
    
    path=psp_path+'fields/level2/'
    
    start = start_timestamp.date()
    end = end_timestamp.date()
    #end = datetime.utcnow().date() + timedelta(days=1)
    while start < end:
        year = start.year
        date_str = f'{year}{start.month:02}{start.day:02}'
        data_url = f'https://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_1min/{year}/'
        soup = BeautifulSoup(urlopen(data_url), 'html.parser')
        for link in soup.find_all('a'):
            href = link.get('href')
            if href is not None and href.startswith('psp_fld_l2_mag_rtn_1min_'+date_str):
                filename = href
                if os.path.isfile(f"{path}{filename}") == True:
                    print(f'{filename} has already been downloaded.')
                else:
                    urllib.request.urlretrieve(data_url+filename, f"{path}{filename}")
                    print(f'Successfully downloaded {filename}')
        start+= datetime.timedelta(days=1)
        


def download_pspplas(start_timestamp, end_timestamp,psp_path):
    
    path=psp_path+'sweap/spi/level2/'
    start = start_timestamp.date()
    #end = datetime.utcnow().date() + timedelta(days=1)
    end = end_timestamp.date()
    
    while start < end:
        year = start.year
        date_str = f'{year}{start.month:02}{start.day:02}'
        data_url = f'https://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spi/l3/spi_sf00_l3_mom/{year}/'

        soup = BeautifulSoup(urlopen(data_url), 'html.parser')
        for link in soup.find_all('a'):
            href = link.get('href')
            if href is not None and href.startswith('psp_swp_spi_sf00_l3_mom_'+date_str):
                filename = href
                if os.path.isfile(f"{path}{filename}") == True:
                    print(f'{filename} has already been downloaded.')
                else:
                    urllib.request.urlretrieve(data_url+filename, f"{path}{filename}")
                    print(f'Successfully downloaded {filename}')
        start+= datetime.timedelta(days=1)


"""
LOAD IN PSP DATA FUNCTIONS: from datapath, and arranges into large dataframes for timerange
"""


def get_pspmag_1min(fp):
    """raw = rtn"""
    try:
        cdf = cdflib.CDF(fp)
        t1 = cdflib.cdfepoch.to_datetime(cdf.varget('epoch_mag_RTN_1min'))
        df = pd.DataFrame(t1, columns=['time'])
        bx, by, bz = cdf['psp_fld_l2_mag_RTN_1min'][:].T
        df['bx'] = bx
        df['by'] = by
        df['bz'] = bz
        df['bt'] = np.linalg.norm(df[['bx', 'by', 'bz']], axis=1)
    except Exception as e:
        print('ERROR:', e, fp)
        df = None
    return df


def get_pspmag_range_1min(start_timestamp, end_timestamp,psp_path ):

    """Pass two datetime objects and grab .cdf files between dates, from
    directory given."""

    path=psp_path+'fields/level2/'

    df = None
    start = start_timestamp.date()
    end = end_timestamp.date() + datetime.timedelta(days=1)
    while start < end:
        date_str = f'{start.year}{start.month:02}{start.day:02}'
        try: 
            fn = glob.glob(path+f'psp_fld_l2_mag_rtn_1min_{date_str}*')[0]
            _df = get_pspmag_1min(fn)
            if _df is not None:
                if df is None:
                    df = _df.copy(deep=True)
                else:
                    df = pd.concat([df, _df])
        except Exception as e:
            print('ERROR:', e, f'{date_str} does not exist')
        start += datetime.timedelta(days=1)       
    return df


def get_pspspi_mom(fp):
    """raw = rtn"""
    try:
        cdf = cdflib.CDF(fp)
        t1 = cdflib.cdfepoch.to_datetime(cdf.varget('Epoch'))
        df = pd.DataFrame(t1, columns=['time'])
        df['np'] = cdf['DENS']
        df['tp'] = cdf['TEMP']*(1.602176634*1e-19)/(1.38064852*1e-23) # from ev to K 
        vx, vy, vz = cdf['VEL_RTN_SUN'][:].T
        df['vx'] = vx
        df['vy'] = vy
        df['vz'] = vz
        df['vt'] = np.linalg.norm(df[['vx', 'vy', 'vz']], axis=1)
    except Exception as e:
        print('ERROR:', e, fp)
        df = None
    return df


def get_pspspi_range_mom(start_timestamp, end_timestamp,psp_path):
    """Pass two datetime objects and grab .cdf files between dates, from
    directory given."""
    
    path=psp_path+'sweap/spi/level2/'
    df = None
    start = start_timestamp.date()
    end = end_timestamp.date() + datetime.timedelta(days=1)
    while start < end:
        date_str = f'{start.year}{start.month:02}{start.day:02}'
        try: 
            fn = glob.glob(path+f'psp_swp_spi_sf00_l3_mom_{date_str}*')[0]
            _df = get_pspspi_mom(fn)
            if _df is not None:
                if df is None:
                    df = _df.copy(deep=True)
                else:
                    df = pd.concat([df, _df])
        except Exception as e:
            print('ERROR:', e, f'{date_str} does not exist')
        start += datetime.timedelta(days=1)       
    return df


"""
PSP POSITION FUNCTIONS: coord maths, call position for each timestamp using astrospice
"""



def get_psp_positions_old(time_series):
    kernels_psp = astrospice.registry.get_kernels('psp', 'predict')
    frame = HeliographicStonyhurst()
    coords_psp = astrospice.generate_coords('Solar probe plus', time_series)
    coords_psp = coords_psp.transform_to(frame)
    x, y, z, r_au = sphere2cart_emma(coords_psp.radius, coords_psp.lat.value, coords_psp.lon.value)
    lat = coords_psp.lat.value
    lon = coords_psp.lon.value
    t = [element.to_pydatetime() for element in list(time_series)]
    positions = np.array([t, x, y, z, r_au, lat, lon])
    df_positions = pd.DataFrame(positions.T, columns=['time', 'x', 'y', 'z', 'r', 'lat', 'lon'])
    return df_positions



#################### new functions without astrospice ************* 
########### only spiceypy



def psp_furnish(kernel_path):

    psp_kernel_path = kernel_path+'psp/'
    generic_path = kernel_path+'generic/'
    
    #psp_kernels = os.listdir(psp_kernel_path)
    #for kernel in psp_kernels:
    #    spiceypy.furnsh(os.path.join(psp_kernel_path, kernel))

    psp_kernel = 'spp_nom_20180812_20300101_v042_PostV7.bsp'
    spiceypy.furnsh(os.path.join(psp_kernel_path, psp_kernel))
    
    generic_kernels = os.listdir(generic_path)
    for kernel in generic_kernels:
        spiceypy.furnsh(os.path.join(generic_path, kernel))
        

def get_psp_pos(t,kernel_path):    
    #if spiceypy.ktotal('ALL') < 1:
    #    psp_furnish(kernel_path)        

    pos = spiceypy.spkpos("PARKER SOLAR PROBE", spiceypy.datetime2et(t), "HEEQ", "NONE", "SUN")[0]
    r, lat, lon = cart2sphere_emma(pos[0],pos[1],pos[2])
    position = t, pos[0], pos[1], pos[2], r, lat, lon
    return position


def get_psp_positions(time_series,kernel_path):
    positions = []
   
    for t in time_series:
        position = get_psp_pos(t,kernel_path)
        positions.append(position)
    df_positions = pd.DataFrame(positions, columns=['time', 'x', 'y', 'z', 'r', 'lat', 'lon'])
    return df_positions



##################################################



























"""
FINAL FUNCTION TO CREATE PICKLE FILE: uses all above functions to create pickle file of 
data from input timestamp to now. 
Can be read in to DataFrame using:
obj = pd.read_pickle('psp_rtn.p')
df = pd.DataFrame.from_records(obj)
"""


def create_psp_pkl(start_time, end_time,psp_file,psp_path,kernels_path):

    #load in mag data to DataFrame and resample, create empty mag and resampled DataFrame if no data
    # if empty, drop time column ready for concat
    df_mag = get_pspmag_range_1min(start_time,end_time,psp_path)
    if df_mag is None:
        print(f'PSP FIELDS data is empty for this timerange')
        df_mag = pd.DataFrame({'time':[], 'bt':[], 'bx':[], 'by':[], 'bz':[]})
        mag_rdf = df_mag.drop(columns=['time'])
    else:
        mag_rdf = df_mag.set_index('time').resample('1min').mean().reset_index(drop=False)
        mag_rdf.set_index(pd.to_datetime(mag_rdf['time']), inplace=True)

    #load in plasma data to DataFrame and resample, create empty plasma and resampled DataFrame if no data
    #only drop time column if MAG DataFrame is not empty
    df_plas = get_pspspi_range_mom(start_time,end_time,psp_path)
    if df_plas is None:
        print(f'PSP SPI/MOM data is empty for this timerange')
        df_plas = pd.DataFrame({'time':[], 'vt':[], 'vx':[], 'vy':[], 'vz':[], 'np':[], 'tp':[]})
        plas_rdf = df_plas
    else:
        plas_rdf = df_plas.set_index('time').resample('1min').mean().reset_index(drop=False)
        plas_rdf.set_index(pd.to_datetime(plas_rdf['time']), inplace=True)
        if mag_rdf.shape[0] != 0:
            plas_rdf = plas_rdf.drop(columns=['time'])

    print('combine mag and plasma')
    #need to combine mag and plasma dfs to get complete set of timestamps for position calculation
    magplas_rdf = pd.concat([mag_rdf, plas_rdf], axis=1)
    #some timestamps may be NaT so after joining, drop time column and reinstate from combined index col
    magplas_rdf = magplas_rdf.drop(columns=['time'])
    magplas_rdf['time'] = magplas_rdf.index

    #get psp positions for corresponding timestamps
    print('get positions')
   
    psp_furnish(kernels_path)
    psp_pos = get_psp_positions(magplas_rdf['time'],kernels_path)   
    psp_pos.set_index(pd.to_datetime(psp_pos['time']), inplace=True)
    psp_pos = psp_pos.drop(columns=['time'])

    #produce final combined DataFrame with correct ordering of columns 
    comb_df = pd.concat([magplas_rdf, psp_pos], axis=1)

    print('produce recarray')

    #produce recarray with correct datatypes
    time_stamps = comb_df['time']
    dt_lst= [element.to_pydatetime() for element in list(time_stamps)] #extract timestamps in datetime.datetime format

    psp=np.zeros(len(dt_lst),dtype=[('time',object),('bx', float),('by', float),('bz', float),('bt', float),\
                ('vx', float),('vy', float),('vz', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float), ('r', float),('lat', float),('lon', float)])
    psp = psp.view(np.recarray) 

    psp.time=dt_lst
    psp.bx=comb_df['bx']
    psp.by=comb_df['by']
    psp.bz=comb_df['bz']
    psp.bt=comb_df['bt']
    psp.vx=comb_df['vx']
    psp.vy=comb_df['vy']
    psp.vz=comb_df['vz']
    psp.vt=comb_df['vt']
    psp.np=comb_df['np']
    psp.tp=comb_df['tp']
    psp.x=comb_df['x']
    psp.y=comb_df['y']
    psp.z=comb_df['z']
    psp.r=comb_df['r']
    psp.lat=comb_df['lat']
    psp.lon=comb_df['lon']
    
    #dump to pickle file
    print('save pickle')

    header='Science level 2 solar wind magnetic field (FIELDS) and plasma data (SWEAP/SPI/MOM) from Parker Solar Probe, ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_1min and https://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spi/l3/spi_sf00_l3_mom/  '+ \
    'Timerange: '+psp.time[0].strftime("%Y-%b-%d %H:%M")+' to '+psp.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', resampled to a time resolution of 1 min. '+\
    'The data are available in a numpy recarray, fields can be accessed by psp.time, psp.bx, psp.vt etc. '+\
    'Total number of data points: '+str(psp.size)+'. '+\
    'Units are btxyz [nT, RTN], vtxyz  [km s^-1, RTN], np[cm^-3], tp [K], heliospheric position x/y/z [km] or r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with heliocats/data_update_web_science.ipynb, by Emma Davies and Christian Möstl. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    pickle.dump([psp,header], open(psp_file, "wb"))

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





#################################################################### SolO

def download_solomag_1min(start_timestamp,end_timestamp, solo_path):
    
    path=solo_path+'mag/level2/'
    start = start_timestamp.date()
    end = end_timestamp.date() 
    #end=datetime.datetime.utcnow().date() + datetime.timedelta(days=1)
    while start < end:
        date_str = f'{start.year}{start.month:02}{start.day:02}'
        data_item_id = f'solo_L2_mag-rtn-normal-1-minute_{date_str}'
        if os.path.isfile(f"{path}{data_item_id}.cdf") == True:
            print(f'{data_item_id}.cdf has already been downloaded.')
            start += datetime.timedelta(days=1)
        else:
            try:
                data_url = f'http://soar.esac.esa.int/soar-sl-tap/data?retrieval_type=PRODUCT&data_item_id={data_item_id}&product_type=SCIENCE'
                urllib.request.urlretrieve(data_url, f"{path}{data_item_id}.cdf")
                print(f'Successfully downloaded {data_item_id}.cdf')
                start += datetime.timedelta(days=1)
            except Exception as e:
                print('ERROR', e, data_item_id)
                start += datetime.timedelta(days=1)


def download_soloplas(start_timestamp, end_timestamp,solo_path):
    
    path=solo_path+'swa/level2/'
    
    start = start_timestamp.date()
    end = end_timestamp.date() 
    #end = datetime.utcnow().date() + timedelta(days=1)
    while start < end:
        date_str = f'{start.year}{start.month:02}{start.day:02}'
        data_item_id = f'solo_L2_swa-pas-grnd-mom_{date_str}'
        if os.path.isfile(f"{path}{data_item_id}.cdf") == True:
            print(f'{data_item_id}.cdf has already been downloaded.')
            start +=datetime. timedelta(days=1)
        else:
            try:
                data_url = f'http://soar.esac.esa.int/soar-sl-tap/data?retrieval_type=PRODUCT&data_item_id={data_item_id}&product_type=SCIENCE'
                urllib.request.urlretrieve(data_url, f"{path}{data_item_id}.cdf")
                print(f'Successfully downloaded {data_item_id}.cdf')
                start += datetime.timedelta(days=1)
            except Exception as e:
                print('ERROR', e, data_item_id)
                start += datetime.timedelta(days=1)


                

"""
LOAD IN SOLO DATA FUNCTIONS: from datapath, and arranges into large dataframes for timerange
"""


def get_solomag(fp):
    """raw = rtn"""
    try:
        cdf = cdflib.CDF(fp)
        t1 = cdflib.cdfepoch.to_datetime(cdf.varget('EPOCH'))
        df = pd.DataFrame(t1, columns=['time'])
        bx, by, bz = cdf['B_RTN'][:].T
        df['bx'] = bx
        df['by'] = by
        df['bz'] = bz
        df['bt'] = np.linalg.norm(df[['bx', 'by', 'bz']], axis=1)
        cols = ['bx', 'by', 'bz', 'bt']
        for col in cols:
            df[col].mask(df[col] < -1E9 , pd.NA, inplace=True)
    except Exception as e:
        print('ERROR:', e, fp)
        df = None
    return df
                
        
        
def get_soloplas(fp):
    """raw = rtn"""
    try:
        cdf = cdflib.CDF(fp)
        t1 = cdflib.cdfepoch.to_datetime(cdf.varget('EPOCH'))
        df = pd.DataFrame(t1, columns=['time'])
        df['np'] = cdf['N']
        df['tp'] = cdf['T']
        vx, vy, vz = cdf['V_RTN'][:].T
        df['vx'] = vx
        df['vy'] = vy
        df['vz'] = vz
        df['vt'] = np.linalg.norm(df[['vx', 'vy', 'vz']], axis=1)
        cols = ['np', 'tp', 'vx', 'vy', 'vz', 'vt']
        for col in cols:
            df[col].mask(df[col] < -9.999E29 , pd.NA, inplace=True)
    except Exception as e:
        print('ERROR:', e, fp)
        df = None
    return df        
                
    

def get_solomag_range_1min(start_timestamp, end_timestamp, solo_path ):
    """Pass two datetime objects and grab .cdf files between dates, from
    directory given."""
    path=solo_path+'mag/level2/'    
    
    df = None
    start = start_timestamp.date()
    end = end_timestamp.date() + datetime.timedelta(days=1)
    while start < end:
        date_str = f'{start.year}{start.month:02}{start.day:02}'
        fn = f'{path}solo_L2_mag-rtn-normal-1-minute_{date_str}.cdf'
        _df = get_solomag(fn)
        if _df is not None:
            if df is None:
                df = _df.copy(deep=True)
            else:
                df = pd.concat([df, _df])
        start += datetime.timedelta(days=1)
    return df


def get_soloplas_range(start_timestamp, end_timestamp, solo_path):

        
    """Pass two datetime objects and grab .cdf files between dates, from
    directory given."""
    
    
    path=solo_path+'swa/level2/'    

    df = None
    start = start_timestamp.date()
    end = end_timestamp.date() + datetime.timedelta(days=1)
    
    while start < end:
        date_str = f'{start.year}{start.month:02}{start.day:02}'
        fn = f'{path}solo_L2_swa-pas-grnd-mom_{date_str}.cdf'
        _df = get_soloplas(fn)
        if _df is not None:
            if df is None:
                df = _df.copy(deep=True)
            else:
                df = pd.concat([df, _df])
        start += datetime.timedelta(days=1)
    return df            
                
                
                
"""
SOLO POSITION FUNCTIONS: coord maths, furnish kernels, and call position for each timestamp
"""



#http://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels/fk/ for solo_ANC_soc-sci-fk_V08.tf
#http://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels/spk/ for solo orbit .bsp

def solo_furnish(kernels_path):
    """Main"""
    solo_path = kernels_path+'solo/'
    generic_path = kernels_path+'generic/'
    #put the latest file here manually
    solo_kernels = astrospice.SPKKernel(solo_path+'solo_ANC_soc-orbit-stp_20200210-20301120_353_V1_00424_V01.bsp')
    generic_kernels = os.listdir(generic_path)
    print(solo_kernels)
    print(generic_kernels)
    #for kernel in solo_kernels:
    #    spiceypy.furnsh(os.path.join(solo_path, kernel))
    #spiceypy.furnsh(solo_kernels)
    for kernel in generic_kernels:
        spiceypy.furnsh(os.path.join(generic_path, kernel))


def get_solo_pos(t,kernels_path):
    if spiceypy.ktotal('ALL') < 1:
        solo_furnish(kernels_path)
    pos = spiceypy.spkpos("SOLAR ORBITER", spiceypy.datetime2et(t), "HEEQ", "NONE", "SUN")[0]
    r, lat, lon = cart2sphere_emma(pos[0],pos[1],pos[2])
    position = t, pos[0], pos[1], pos[2], r, lat, lon
    return position


def get_solo_positions(time_series,kernels_path):
    positions = []
    for t in time_series:
        position = get_solo_pos(t,kernels_path)
        positions.append(position)
    df_positions = pd.DataFrame(positions, columns=['time', 'x', 'y', 'z', 'r', 'lat', 'lon'])
    return df_positions

           
                
    
"""
FINAL FUNCTION TO CREATE PICKLE FILE: uses all above functions to create pickle file of 
data from input timestamp to now. 
Can be read in to DataFrame using:
obj = pd.read_pickle('solo_rtn.p')
df = pd.DataFrame.from_records(obj)
"""


def create_solo_pkl(start_timestamp,end_timestamp,solo_file,solo_path,kernels_path):

    #load in mag data to DataFrame and resample, create empty mag and resampled DataFrame if no data
    # if empty, drop time column ready for concat
    df_mag = get_solomag_range_1min(start_timestamp,end_timestamp,solo_path)
    if df_mag is None:
        print(f'SolO MAG data is empty for this timerange')
        df_mag = pd.DataFrame({'time':[], 'bt':[], 'bx':[], 'by':[], 'bz':[]})
        mag_rdf = df_mag.drop(columns=['time'])
    else:
        mag_rdf = df_mag.set_index('time').resample('1min').mean().reset_index(drop=False)
        mag_rdf.set_index(pd.to_datetime(mag_rdf['time']), inplace=True)
        
    #load in plasma data to DataFrame and resample, create empty plasma and resampled DataFrame if no data
    #only drop time column if MAG DataFrame is not empty
    df_plas = get_soloplas_range(start_timestamp,end_timestamp,solo_path)
    if df_plas is None:
        print(f'SolO SWA data is empty for this timerange')
        df_plas = pd.DataFrame({'time':[], 'vt':[], 'vx':[], 'vy':[], 'vz':[], 'np':[], 'tp':[]})
        plas_rdf = df_plas
    else:
        plas_rdf = df_plas.set_index('time').resample('1min').mean().reset_index(drop=False)
        plas_rdf.set_index(pd.to_datetime(plas_rdf['time']), inplace=True)
        if mag_rdf.shape[0] != 0:
            plas_rdf = plas_rdf.drop(columns=['time'])

    #need to combine mag and plasma dfs to get complete set of timestamps for position calculation
    magplas_rdf = pd.concat([mag_rdf, plas_rdf], axis=1)
    #some timestamps may be NaT so after joining, drop time column and reinstate from combined index col
    magplas_rdf = magplas_rdf.drop(columns=['time'])
    magplas_rdf['time'] = magplas_rdf.index
     
    solo_furnish(kernels_path)    
    #get solo positions for corresponding timestamps
    solo_pos = get_solo_positions(magplas_rdf['time'],kernels_path)
    solo_pos.set_index(pd.to_datetime(solo_pos['time']), inplace=True)
    solo_pos = solo_pos.drop(columns=['time'])

    #produce final combined DataFrame with correct ordering of columns 
    comb_df = pd.concat([magplas_rdf, solo_pos], axis=1)

    #produce recarray with correct datatypes
    time_stamps = comb_df['time']
    dt_lst= [element.to_pydatetime() for element in list(time_stamps)] #extract timestamps in datetime.datetime format

    solo=np.zeros(len(dt_lst),dtype=[('time',object),('bx', float),('by', float),('bz', float),('bt', float),\
                ('vx', float),('vy', float),('vz', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float), ('r', float),('lat', float),('lon', float)])
    solo = solo.view(np.recarray) 

    solo.time=dt_lst
    solo.bx=comb_df['bx']
    solo.by=comb_df['by']
    solo.bz=comb_df['bz']
    solo.bt=comb_df['bt']
    solo.vx=comb_df['vx']
    solo.vy=comb_df['vy']
    solo.vz=comb_df['vz']
    solo.vt=comb_df['vt']
    solo.np=comb_df['np']
    solo.tp=comb_df['tp']*(1.602176634*1e-19)/(1.38064852*1e-23) # from ev to K 
    solo.x=comb_df['x']
    solo.y=comb_df['y']
    solo.z=comb_df['z']
    solo.r=comb_df['r']
    solo.lat=comb_df['lat']
    solo.lon=comb_df['lon']
    
    
    
    #remove SolO Earth flyby
    
    
    #delete Earth flyby
    nt1=parse_time('2021-11-27T02:00Z').datetime 
    nt2=parse_time('2021-11-27T08:00Z').datetime
    gap=np.where(np.logical_and(solo.time > nt1,solo.time < nt2 ))[0]
    
    solo.bt[gap]=np.nan
    solo.bx[gap]=np.nan
    solo.by[gap]=np.nan
    solo.bz[gap]=np.nan
    
    header='Science level 2 solar wind magnetic field (MAG) and plasma data (SWA) from Solar Orbiter, ' + \
    'obtained from https://soar.esac.esa.int/soar/  '+ \
    'Timerange: '+solo.time[0].strftime("%Y-%b-%d %H:%M")+' to '+solo.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', resampled to a time resolution of 1 min. '+\
    'The data are available in a numpy recarray, fields can be accessed by solo.time, solo.bx, solo.vt etc. '+\
    'Total number of data points: '+str(solo.size)+'. '+\
    'Units are btxyz [nT, RTN], vtxyz [km s^-1, RTN], np[cm^-3], tp [K], heliospheric position x/y/z [km], r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with heliocats/data_update_web_science.ipynb, by Emma Davies and Christian Möstl. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    
    #dump to pickle file 
    pickle.dump([solo,header], open(solo_file, "wb"))
         
                
                
                
                
################################## STEREO-A            


"""
DOWNLOAD STEREOA DATA FUNCTIONS: MERGED 1MIN IMPACT
"""

def download_stereoa_merged(start_timestamp, end_timestamp, stereoa_path):
        
    path=stereoa_path+'impact/merged/level2/'
    
    start = start_timestamp.year
    end = end_timestamp.year +1
    while start < end:
        year = start
        date_str = f'{year}0101'
        try: 
            data_url = f'https://spdf.gsfc.nasa.gov/pub/data/stereo/ahead/l2/impact/magplasma/1min/{year}/'
            soup = BeautifulSoup(urlopen(data_url), 'html.parser')
            for link in soup.find_all('a'):
                href = link.get('href')
                if href is not None and href.startswith('sta_l2_magplasma_1m_'+date_str):
                    filename = href
                    if os.path.isfile(f"{path}{filename}") == True:
                        print(f'{filename} has already been downloaded.')
                    else:
                        urllib.request.urlretrieve(data_url+filename, f"{path}{filename}")
                        print(f'Successfully downloaded {filename}')
        except Exception as e:
            print('ERROR', e, f'.File for {year} does not exist.')
        start+=1


"""
LOAD IN STEREOA DATA FUNCTIONS: from datapath, and arranges into large dataframes for timerange
"""


def get_stereoa_merged(fp):
    """raw = rtn"""
    try:
        cdf = cdflib.CDF(fp)
        t1 = cdflib.cdfepoch.to_datetime(cdf.varget('Epoch'))
        df = pd.DataFrame(t1, columns=['time'])
        bx, by, bz = cdf['BFIELDRTN'][:].T
        df['bx'] = bx
        df['by'] = by
        df['bz'] = bz
        df['bt'] = cdf['BTOTAL']
        df['np'] = cdf['Np']
        df['tp'] = cdf['Tp']
        df['vt'] = cdf['Vp']
        cols = ['bx', 'by', 'bz', 'bt', 'np', 'tp', 'vt']
        for col in cols:
            df[col].mask(df[col] < -9.999E29 , pd.NA, inplace=True)
        df['vx'] = cdf['Vr_Over_V_RTN']*df['vt']
        df['vy'] = cdf['Vt_Over_V_RTN']*df['vt']
        df['vz'] = cdf['Vn_Over_V_RTN']*df['vt']
        v_cols = ['vx', 'vy', 'vz']
        for v_col in v_cols:
            df[v_col].mask(df[v_col] < -9.999E29 , pd.NA, inplace=True)
    except Exception as e:
        print('ERROR:', e, fp)
        df = None
    return df


def get_stereoa_merged_range(start_timestamp, end_timestamp, stereoa_path):
    """Pass two datetime objects and grab .cdf files between dates, from
    directory given."""
    
    path=stereoa_path+'impact/merged/level2/'
    
    df=None
    start = start_timestamp.year
    end = end_timestamp.year+1
    while start < end:
        year = start
        date_str = f'{year}0101'
        try: 
            fn = glob.glob(path+f'sta_l2_magplasma_1m_{date_str}*')[0]
            print(fn)
            _df = get_stereoa_merged(fn)
            if _df is not None:
                if df is None:
                    df = _df.copy(deep=True)
                else:
                    df = pd.concat([df, _df])
        except Exception as e:
            print('ERROR:', e, f'{date_str} does not exist')
        start += 1
    timemask = (df['time']>=start_timestamp) & (df['time']<=end_timestamp)
    df = df[timemask]
    return df


"""
STEREOA POSITION FUNCTIONS: coord maths, call position for each timestamp using astrospice
"""



#def get_stereoa_positions(time_series):
#    kernels_sta = astrospice.registry.get_kernels('stereo-a', 'predict')
#    frame = HeliographicStonyhurst()
#    coords_sta = astrospice.generate_coords('Stereo ahead', time_series)
#    coords_sta = coords_sta.transform_to(frame)
#    x, y, z, r_au = sphere2cart_emma(coords_sta.radius, coords_sta.lat.value, coords_sta.lon.value)
#    lat = coords_sta.lat.value
#    lon = coords_sta.lon.value
#    t = [element.to_pydatetime() for element in list(time_series)]
#    positions = np.array([t, x, y, z, r_au, lat, lon])
#    df_positions = pd.DataFrame(positions.T, columns=['time', 'x', 'y', 'z', 'r', 'lat', 'lon'])
#    return df_positions




def stereoa_furnish(kernel_path):
    """Main"""
    stereo_kernel_path = kernel_path+'stereoa/'
    generic_path = kernel_path+'generic/'
    stereoa_kernels = os.listdir(stereo_kernel_path)
    generic_kernels = os.listdir(generic_path)
    for kernel in stereoa_kernels:
        spiceypy.furnsh(os.path.join(stereo_kernel_path, kernel))
    for kernel in generic_kernels:
        spiceypy.furnsh(os.path.join(generic_path, kernel))


def get_stereoa_pos(t,kernel_path):
    if spiceypy.ktotal('ALL') < 1:
        stereoa_furnish(kernel_path)
    pos = spiceypy.spkpos("STEREO AHEAD", spiceypy.datetime2et(t), "HEEQ", "NONE", "SUN")[0]
    r, lat, lon = cart2sphere_emma(pos[0],pos[1],pos[2])
    position = t, pos[0], pos[1], pos[2], r, lat, lon
    return position


def get_stereoa_positions(time_series,kernel_path):
    positions = []
    for t in time_series:
        position = get_stereoa_pos(t,kernel_path)
        positions.append(position)
    df_positions = pd.DataFrame(positions, columns=['time', 'x', 'y', 'z', 'r', 'lat', 'lon'])
    return df_positions



"""
FINAL FUNCTION TO CREATE PICKLE FILE: uses all above functions to create pickle file of 
data from input timestamp to now. 
Can be read in to DataFrame using:
obj = pd.read_pickle('stereoa_rtn.p')
df = pd.DataFrame.from_records(obj)
"""


def create_stereoa_pkl(start_timestamp,end_timestamp,stereoa_file,stereoa_path,kernel_path):

    # #download stereo-a merged magplasma data up to endtime
    #download_stereoa_merged(start_timestamp,end_timestamp)

    
    
    #load in merged mag and plasma data to DataFrame, create empty DataFrame if no data
    # if empty, drop time column ready for concat
    df_magplas = get_stereoa_merged_range(start_timestamp,end_timestamp,stereoa_path)
    
    if df_magplas is None:
        print(f'STEREO A data is empty for this timerange')
        df_magplas = pd.DataFrame({'time':[], 'bt':[], 'bx':[], 'by':[], 'bz':[], 'vt':[], 'vx':[], 'vy':[], 'vz':[], 'np':[], 'tp':[]})
        df_magplas = df_magplas.drop(columns=['time'])
    else:
        df_magplas.set_index(pd.to_datetime(df_magplas['time']), inplace=True)

    #get stereoa positions for corresponding timestamps
    print(kernel_path)
    stereoa_furnish(kernel_path)  
    sta_pos = get_stereoa_positions(df_magplas['time'],kernel_path)
    sta_pos.set_index(pd.to_datetime(sta_pos['time']), inplace=True)
    sta_pos = sta_pos.drop(columns=['time'])

    #produce final combined DataFrame with correct ordering of columns 
    comb_df = pd.concat([df_magplas, sta_pos], axis=1)

    #produce recarray with correct datatypes
    time_stamps = comb_df['time']
    dt_lst= [element.to_pydatetime() for element in list(time_stamps)] #extract timestamps in datetime.datetime format

    stereoa=np.zeros(len(dt_lst),dtype=[('time',object),('bx', float),('by', float),('bz', float),('bt', float),\
                ('vx', float),('vy', float),('vz', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float), ('r', float),('lat', float),('lon', float)])
    stereoa = stereoa.view(np.recarray) 

    stereoa.time=dt_lst
    stereoa.bx=comb_df['bx']
    stereoa.by=comb_df['by']
    stereoa.bz=comb_df['bz']
    stereoa.bt=comb_df['bt']
    stereoa.vx=comb_df['vx']
    stereoa.vy=comb_df['vy']
    stereoa.vz=comb_df['vz']
    stereoa.vt=comb_df['vt']
    stereoa.np=comb_df['np']
    stereoa.tp=comb_df['tp']
    stereoa.x=comb_df['x']
    stereoa.y=comb_df['y']
    stereoa.z=comb_df['z']
    stereoa.r=comb_df['r']
    stereoa.lat=comb_df['lat']
    stereoa.lon=comb_df['lon']
    
    #dump to pickle file
    
    header='Science level 2 solar wind magnetic field and plasma data (IMPACT) from STEREO A, ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/stereo/ahead/l2/impact/magplasma/1min/ '+ \
    'Timerange: '+stereoa.time[0].strftime("%Y-%b-%d %H:%M")+' to '+stereoa.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', with a time resolution of 1 min. '+\
    'The data are available in a numpy recarray, fields can be accessed by stereoa.time, stereoa.bx, stereoa.vt etc. '+\
    'Total number of data points: '+str(stereoa.size)+'. '+\
    'Units are btxyz [nT, RTN], vtxy  [km s^-1, RTN], np[cm^-3], tp [K], heliospheric position x/y/z [km], r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with heliocats/data_update_web_science.ipynb by Emma Davies and Christian Möstl. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    pickle.dump([stereoa,header], open(stereoa_file, "wb"))
    

























###################################################################




def get_noaa_ephem(ephemfile):
    
    data = json.loads(ephemfile.read())    
    
    #make array for 7 days
    ephem=np.zeros(len(data)-1,dtype=[('time',object),('xgse', float),('ygse', float),('zgse', float)]) 
    ephem=ephem.view(np.recarray)    

    ephem.time = [datetime.datetime.strptime(x[0], "%Y-%m-%d %H:%M:%S.%f")  for x in data[1:]]
    #print(ephem.time)
    ephem.xgse=[x[1] for x in data[1:]]
    ephem.ygse=[x[2] for x in data[1:]]
    ephem.zgse=[x[3] for x in data[1:]]

    return ephem


def get_noaa_dst(dstfile):
    
    dm = json.loads(dstfile.read())
    
    
    #make array for 7 days
    dst=np.zeros(len(dm)-1,dtype=[('time',object),('dst', float)]) 
    dst=dst.view(np.recarray)    

    dst.time = [datetime.datetime.strptime(x[0], "%Y-%m-%d %H:%M:%S")  for x in dm[1:]]
    dst.dst=[x[1] for x in dm[1:]]

    return dst


def get_noaa_json(magfile, plasmafile):
    
    
    
    # Read magnetic field data:
    dm = json.loads(magfile.read())
    
    dmn = [[np.nan if x == None else x for x in d] for d in dm]     # Replace None w NaN
    dtype=[(x, 'float') for x in dmn[0]]
    datesm = [datetime.datetime.strptime(x[0], "%Y-%m-%d %H:%M:%S.%f")  for x in dmn[1:]]
    mdatesm=mdates.date2num(datesm)
    dm_ = [tuple([d]+[float(y) for y in x[1:]]) for d, x in zip(mdatesm, dm[1:])] 
    noaa_m = np.array(dm_, dtype=dtype)
    
     
    #print(noaa_m)
    
    
    # Read plasma data:
    dp = json.loads(plasmafile.read())
    dpn = [[np.nan if x == None else x for x in d] for d in dp]     # Replace None w NaN
    dtype=[(x, 'float') for x in dp[0]]
    datesp = [datetime.datetime.strptime(x[0], "%Y-%m-%d %H:%M:%S.%f")  for x in dpn[1:]]
    #convert datetime to matplotlib times
    mdatesp=mdates.date2num(datesp)
    dp_ = [tuple([d]+[float(y) for y in x[1:]]) for d, x in zip(mdatesp, dpn[1:])] 
    noaa_p = np.array(dp_, dtype=dtype)

    #print(noaa_p)

  
    
    #use mag for times
    t_start=datesm[0]
    t_end=datesm[-1]

    #1 minute resolution between t_start and t_end for given mag file
    itime = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int (((t_end - t_start).days+1)*60*24))]  
    #print(itime[0])
    #print(itime[-1])
    
    #make date numbers
    itime_num=mdates.date2num(itime)
    #interpolate all onto itime_num
    
    rbtot_m = np.interp(itime_num, noaa_m['time_tag'], noaa_m['bt'])
    rbxgsm_m = np.interp(itime_num, noaa_m['time_tag'], noaa_m['bx_gsm'])
    rbygsm_m = np.interp(itime_num, noaa_m['time_tag'], noaa_m['by_gsm'])
    rbzgsm_m = np.interp(itime_num, noaa_m['time_tag'], noaa_m['bz_gsm'])
    rpv_m = np.interp(itime_num, noaa_p['time_tag'], noaa_p['speed'])
    rpn_m = np.interp(itime_num, noaa_p['time_tag'], noaa_p['density'])
    rpt_m = np.interp(itime_num, noaa_p['time_tag'], noaa_p['temperature'])

    #make array
    noaa_data=np.zeros(np.size(rbtot_m),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('vx', float),('vy', float),('vz', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
                        
    #convert to recarray
    noaa_data = noaa_data.view(np.recarray)                             

    noaa_data.time=itime
    noaa_data.bt=rbtot_m
    noaa_data.bx=rbxgsm_m
    noaa_data.by=rbygsm_m
    noaa_data.bz=rbzgsm_m

    noaa_data.vt=rpv_m
    noaa_data.np=rpn_m
    noaa_data.tp=rpt_m

    noaa_data.vx=np.nan
    noaa_data.vy=np.nan
    noaa_data.vz=np.nan
    
    
    #print('NOAA data read completed for file with end time: ',itime[-1])
    
    
      
    return noaa_data


























def save_noaa_rtsw_data(data_path,noaa_path,filenoaa,filedst, cutoff):


    print(' ')
    print('convert NOAA real time solar wind archive to pickle file')
    
    print('directories for the json data')
    
    #get mag data files
    print(noaa_path+'mag/')
    items=os.listdir(noaa_path+'mag/')  
    maglist = [] 
    for names in items: 
       if names.endswith(".json"):
            maglist.append(names)
    maglist=np.sort(maglist)
    #print(maglist)    
    
    #cutoff last N files
    maglist=maglist[-cutoff:]
    print('Sorted file list to be read with cutoff ',cutoff,' files. ')
    print(maglist)
    
    #---------------------------------

    #get plasma data files
    print(noaa_path+'plasma/')
    items=os.listdir(noaa_path+'plasma/')  
    plasmalist = [] 
    for names in items: 
       if names.endswith(".json"):
            plasmalist.append(names)
    plasmalist=np.sort(plasmalist)        

    #cutoff last N files
    plasmalist=plasmalist[-cutoff:]

    print(plasmalist)
    print()
        
    #-------------------------        
    #get ephemerides files
    print(noaa_path+'ephem/')
    items=os.listdir(noaa_path+'ephem/')  
    ephemlist = [] 
    for names in items: 
       if names.endswith(".json"):
            ephemlist.append(names)
    ephemlist=np.sort(ephemlist)        

    #cutoff last N files
    ephemlist=ephemlist[-cutoff:]

    print(ephemlist)
    print()
               
   
    #make array 
    noaa=np.zeros(int(1e7),dtype=[('time',object),('bx', float),('by', float),\
                    ('bz', float),('bt', float),('vt', float),('vx', float),('vy', float),('vz', float),('np', float),('tp', float),\
                    ('x', float),('y', float),('z', float),\
                    ('r', float),('lat', float),('lon', float)])   

 

    #counter    
    k=0
    
    
    #plasma and mag

    for i in np.arange(len(maglist))-1:

        #read in data of corresponding files
        
        m1=open(noaa_path+'mag/'+maglist[i],'r')
        #print(noaa_path+'mag/'+maglist[i])
        p1=open(noaa_path+'plasma/'+plasmalist[i],'r')
        #print(noaa_path+'plasma/'+plasmalist[i])
        
      
        #extract data from files
        try: 
            d1=get_noaa_json(m1, p1)
            #save in large array
            noaa[k:k+np.size(d1)]=d1
            k=k+np.size(d1) 
        except:
            print('one of these json could not be loaded')
            print(m1)
            print(p1)
            

    #cut zeros, sort, convert to recarray, and find unique times and data

    noaa_cut=noaa[0:k]
    noaa_cut.sort()
    
    nu=noaa_cut.view(np.recarray)
    [dum,ind]=np.unique(nu.time,return_index=True)  
    nf=nu[ind]
    
    
    #add positions to the final array
    
    
    print('position start')
    
    ephem=np.zeros(int(1e6),dtype=[('time',object),('xgse', float),('ygse', float),('zgse', float)]) 
    ephem=ephem.view(np.recarray)  

    
    # exact position from ephemerides.json files, downloaded each day        
    k=0
    
    for i in np.arange(len(ephemlist))-1:

        #read in data of corresponding files
        
        e1=open(noaa_path+'ephem/'+ephemlist[i],'r')
      
        #extract data from files
        try: 
            e2=get_noaa_ephem(e1)
            ephem[k:k+np.size(e2)]=e2
            k=k+np.size(e2) 
        except:
            print('one of these json could not be loaded')
            print(e1)
            

    ephem_cut=ephem[0:k]
    ephem_cut.sort()   
             
    eu=ephem_cut.view(np.recarray)
    [dum,ind]=np.unique(eu.time,return_index=True)  
    ef=eu[ind]

    
    au_km = au.to('km').value
    
    #interpolate hourly ephemerides onto 1 minute data 
    #but only for those times where it is available, otherwise 
    #set x_gse to 0.01 AU and thats it
    
    nf_mat=mdates.date2num(nf.time)
    ef_mat=mdates.date2num(ef.time)
    
    
    #first set all xyz_gse to 0.01 and 0 for L1
    x_gse=np.zeros(len(nf_mat))+0.01
    y_gse=np.zeros(len(nf_mat))
    z_gse=np.zeros(len(nf_mat))

    #ef_mat is the earliest available ephemerides time
    nf_eph_ind=np.where(nf_mat > ef_mat[0])[0]
    print(nf_eph_ind)

    #only interpolate the nf_mat times after this date
    x_gse[nf_eph_ind] = np.interp(nf_mat[nf_eph_ind], ef_mat,ef.xgse )/au_km 
    y_gse[nf_eph_ind] = np.interp(nf_mat[nf_eph_ind], ef_mat,ef.ygse )/au_km 
    z_gse[nf_eph_ind] = np.interp(nf_mat[nf_eph_ind], ef_mat,ef.zgse )/au_km 

    #print(x_gse)
    #print(y_gse)
    #print(z_gse)

    coords_earth = astrospice.generate_coords('Earth', nf.time)
    #coords_earth_heeq=coords_earth.transform_to(HeliographicStonyhurst())  #returns degrees, km
    coords_earth_hee=coords_earth.transform_to(HeliocentricEarthEcliptic())  #returns degrees, km

    #change angles to -90/90 latitude, distance to AU
    [earth_hee_x, earth_hee_y, earth_hee_z] = sphere2cart(coords_earth_hee.distance.value/au_km,\
                                                              np.deg2rad(-coords_earth_hee.lat.value+90),\
                                                          np.deg2rad(coords_earth_hee.lon.value))
    
 
    nf.x=earth_hee_x-x_gse  #both values are in au, 0.01 for L1
    nf.y=earth_hee_y-y_gse
    nf.z=earth_hee_z+z_gse
    
    #convert position to HEEQ  for xyz
    nf = convert_HEE_to_HEEQ(nf)
    
    [nf.r, nf.lat, nf.lon]=cart2sphere(nf.x,nf.y,nf.z)       
    nf.lat=np.rad2deg(nf.lat)
    nf.lon=np.rad2deg(nf.lon)  



    
    print('position end ')
       
        
    header='Real time solar wind magnetic field and plasma data from NOAA, ' + \
        'obtained daily from https://services.swpc.noaa.gov/products/solar-wind/  '+ \
        'Timerange: '+nf.time[0].strftime("%Y-%b-%d %H:%M")+' to '+nf.time[-1].strftime("%Y-%b-%d %H:%M")+\
        ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(nf.time)).seconds)+' seconds. '+\
        'The data are available in a numpy recarray, fields can be accessed by nf.time, nf.bx, nf.vt etc. '+\
        'Total number of data points: '+str(nf.size)+'. '+\
        'Units are btxyz [nT, GSM], vt  [km s^-1], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
        'Made with https://github.com/cmoestl/heliocats save_noaa_rtsw_data  '+\
        'By C. Moestl (twitter @chrisoutofspace) and R. Bailey. File creation date: '+\
        datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    
    print('file saved ',data_path+filenoaa)
    pickle.dump([nf,header], open(data_path+filenoaa, "wb"))
    
    
    
    #################### dst file
    
    print()
    print()
    print('-------- make Dst file ----')

    
    #get dst data files
    print(noaa_path+'dst/')
    items=os.listdir(noaa_path+'dst/')  
    dstlist = [] 
    for names in items: 
       if names.endswith(".json"):
            dstlist.append(names)
    dstlist=np.sort(dstlist)
    #cutoff last N files
    dstlist=dstlist[-cutoff:]

    print(dstlist)

        
    #make array 
    #dst=np.zeros(int(1e7),dtype=[('time',object),('dst', float),('nc',float)])   
    dst=np.zeros(int(1e7),dtype=[('time',object),('dst', float)])

    #counter    
    k=0
    
    for i in np.arange(len(dstlist))-1:

        #read in data of corresponding files
        
        dstfile1=open(noaa_path+'dst/'+dstlist[i],'r')
 
        #extract data from files
        try: 
            dst1=get_noaa_dst(dstfile1)
            dst[k:k+np.size(dst1)]=dst1
            k=k+np.size(dst1) 
        except:
            print('json could not be loaded')
            print(dstfile1)
      
           
    #cut zeros, sort, convert to recarray, and find unique times and data

    dst_cut=dst[0:k]
    dst_cut.sort()    
             
    dstu=dst_cut.view(np.recarray)
    [dum,ind]=np.unique(dstu.time,return_index=True)  
    dstf=dstu[ind]
    
    
    #calculate newell coupling from the solar wind
    #TBD
    
    
    
    #put both in a recarray
    
    dst_nc=np.zeros(len(dstf),dtype=[('time',object),('dst', float),('nc',float)])   
    dst_nc=dst_nc.view(np.recarray)
    
    dst_nc.time=dstf.time
    dst_nc.dst=dstf.dst
    dst_nc.nc=np.zeros(len(dstf))
    
    print('file saved ',data_path+filedst)
    pickle.dump(dst_nc, open(data_path+filedst, "wb"))
    print(' ')
    print(' ')
    



    
    
def stereoa_download_beacon(start, end, stereoa_path):

    #download MFI and SWE in ASCII from SPDF
        
    print('download STEREO-A beacon data to ', stereoa_path)
        
        
        
    #define utc now as end date   
    read_data_end_year=end.year
    read_data_end_month=end.month
    read_data_end_day=end.day
    
    
    #sta_years_strings=[]
    #for j in np.arange(start.year,read_data_end_year+1):
    #    sta_years_strings.append(str(j))

    #print('for years', sta_years_strings)
    
    ## BEACON SOURCE location  
    
    impact_url='https://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/impact/'
    plastic_url='https://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/plastic/'
    
    #make a list of all dates from start year,month,day onwards until now
    tstart1=datetime.datetime(start.year,start.month,start.day,1)
    
    time_1=[]
    while tstart1 < end:
        time_1.append(tstart1)  
        tstart1 += datetime.timedelta(days=1)
    
    for i in np.arange(0,len(time_1)):
        
        yearstr=str(time_1[i].year)        
        monthstr=str(time_1[i].month).zfill(2) #add leading zeros 
        daystr=str(time_1[i].day).zfill(2) #add leading zeros 

        filesta='STA_LB_IMPACT_'+yearstr+monthstr+daystr+'_V02.cdf'
        print(filesta)        
        try: urllib.request.urlretrieve(impact_url+yearstr+'/'+monthstr+'/'+filesta, stereoa_path+'beacon/impact/'+filesta)
        except urllib.error.URLError as e:
                 print(' ', impact_url,' ',e.reason)
                
                
        filesta='STA_LB_PLASTIC_'+yearstr+monthstr+daystr+'_V14.cdf'    
        print(filesta) 
        try: urllib.request.urlretrieve(plastic_url+yearstr+'/'+monthstr+'/'+filesta, stereoa_path+'beacon/plastic/'+filesta)
        except urllib.error.URLError as e:
                   print(' ', plastic_url,' ',e.reason)
      
            
            

def save_stereoa_beacon_data(data_path,path,file_rtn,file_gsm,t_start,t_end,coord):
            
        

    #round tstart to nearest minute
    t_start=datetime.datetime(t_start.year, t_start.month, t_start.day,t_start.hour,t_start.minute)
    t_start1=copy.deepcopy(t_start)
    time_1=[]
        
    #make 1 min datetimes
    while t_start1 < t_end:
        time_1.append(t_start1)  
        t_start1 += datetime.timedelta(minutes=1)
    print(time_1[0])    
    print(time_1[-1])

    #make array for 1 min data
    sta=np.zeros(len(time_1),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vx', float),('vy', float),('vz', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    sta = sta.view(np.recarray)  
    sta.time=time_1
   

    #make data file names
    t_start1=copy.deepcopy(t_start)
    days_sta = []
    days_str = []
    i=0
    while t_start < t_end:
        days_sta.append(t_start)  
        days_str.append(str(days_sta[i])[0:4]+str(days_sta[i])[5:7]+str(days_sta[i])[8:10])        
        i=i+1
        t_start +=datetime.timedelta(days=1)

    print(days_str)
    

    bx=np.zeros(int(1e7))
    by=np.zeros(int(1e7))
    bz=np.zeros(int(1e7))
    t2=[]
    i=0
    
    #go through all files    
    for days_date in days_str:
        
        cdf_file = 'STA_LB_IMPACT_{}_V02.cdf'.format(days_date)
        
        if os.path.exists(path+'beacon/impact/'+cdf_file):
            
            print(cdf_file)
            f1 = cdflib.CDF(path+'beacon/impact/'+cdf_file)
            #print(f1)
            #convert from epoch to datetime
            t1 = cdflib.cdfepoch.to_datetime(f1.varget('Epoch_MAG')) 
            
            t2.extend(t1)
            bfield=f1.varget('MAGBField') #RTN
            #bt[i:i+len(bfield[:,3])]=bfield[:,3]
            bx[i:i+len(bfield[:,0])]=bfield[:,0]
            by[i:i+len(bfield[:,1])]=bfield[:,1]
            bz[i:i+len(bfield[:,2])]=bfield[:,2]
            
            i=i+len(bfield[:,0])
            
    print('done') 

    
    #cut array

    bx=bx[0:i]
    by=by[0:i]
    bz=bz[0:i]
    
    tm2=mdates.date2num(t2)
    time_mat=mdates.date2num(time_1)


    #linear interpolation to time_mat times    
    sta.bx = np.interp(time_mat, tm2, bx )
    sta.by = np.interp(time_mat, tm2, by )
    sta.bz = np.interp(time_mat, tm2, bz )
    #sta.bt = np.sqrt(sta.bx**2+sta.by**2+sta.bz**2)
    
    
    #set missing data to NaN
    
    #round first each original time to full minutes   original data at 30sec
    #get rid of datetime64
    tround=copy.deepcopy(t2)

    format_str = '%Y-%m-%d %H:%M'  
    for k in np.arange(np.size(t2)):
         tround[k] = datetime.datetime.strptime(datetime.datetime.strftime(t2[k].astype(datetime.datetime), format_str), format_str) 
    tm2_round=mdates.date2num(tround)

    #which values are not in original data compared to full time range
    isin=np.isin(time_mat,tm2_round)      
    setnan=np.where(isin==False)
    #set to to nan that is not in original data
    sta.bx[setnan]=np.nan
    sta.by[setnan]=np.nan
    sta.bz[setnan]=np.nan
    sta.bt = np.sqrt(sta.bx**2+sta.by**2+sta.bz**2)
    
    
    
    #plt.plot(t2,bx)
    #plt.plot(t2,by)
    #plt.plot(t2,bz)
    #plt.plot(t2,bt,'-k')
    
    
    
    ########## PLASTIC

    print(days_str)
    
    #set variables
    vx=np.zeros(int(1e7))
    vy=np.zeros(int(1e7))
    vz=np.zeros(int(1e7))
    vt=np.zeros(int(1e7))
    
    den=np.zeros(int(1e7))
    
    t2=[]
    i=0
    
    #go through all files    
    for days_date in days_str:
        
        cdf_file = 'STA_LB_PLASTIC_{}_V14.cdf'.format(days_date)
        print(cdf_file)
        
        if os.path.exists(path+'beacon/plastic/'+cdf_file):
            

            f1 = cdflib.CDF(path+'beacon/plastic/'+cdf_file)
            
            #print(f1.cdf_info())
            #print(f1)
            #convert from epoch to datetime
            t1 = cdflib.cdfepoch.to_datetime(f1.varget('Epoch1')) 
            
            #Velocity_RTN
            #Density get attributes for this variabel f1.varattsget('Density')
            
            t2.extend(t1)
            vfield=f1.varget('Velocity_RTN') #RTN
            vx[i:i+len(vfield[:,0])]=vfield[:,0]
            vy[i:i+len(vfield[:,1])]=vfield[:,1]
            vz[i:i+len(vfield[:,2])]=vfield[:,2]
            vt[i:i+len(vfield[:,3])]=vfield[:,3]
            
            den1=f1.varget('Density')             
            den[i:i+len(den1)]=den1            
            
            i=i+len(vfield[:,0])
            

            
    print('PLASTIC done') 

    
    #cut arrays

    vx=vx[0:i]
    vy=vy[0:i]
    vz=vy[0:i]    
    vt=vt[0:i]
    den=den[0:i]
    
    tp2=mdates.date2num(t2)

    #set missing data to nan    
    vx[np.where(vx < -1e30)]=np.nan  
    vy[np.where(vy < -1e30)]=np.nan  
    vz[np.where(vy < -1e30)]=np.nan  
    vt[np.where(vt < -1e30)]=np.nan  
    den[np.where(den < -1e30)]=np.nan  


    #linear interpolation to time_mat times    
    sta.vx = np.interp(time_mat, tp2, vx )
    sta.vy = np.interp(time_mat, tp2, vy )
    sta.vz = np.interp(time_mat, tp2, vz )
    sta.vt = np.interp(time_mat, tp2, vt )
    sta.np = np.interp(time_mat, tp2, den )
    
    
    #POSITION

    print('position start with astrospice stereo-A predict')
    
    frame = HeliographicStonyhurst()

    kernels_stereo_pred = astrospice.registry.get_kernels('stereo-a', 'predict')
    stereo_kernel_pred = kernels_stereo_pred[0]
    coverage_stereo_pred = stereo_kernel_pred.coverage('Stereo ahead')
    times_stereo_pred = sta.time
    coords_stereo_pred = astrospice.generate_coords('Stereo ahead', times_stereo_pred)
    coords_stereo_pred = coords_stereo_pred.transform_to(frame)

    sta.r = coords_stereo_pred.radius.to(u.au).value
    sta.lon = coords_stereo_pred.lon.value #degrees
    sta.lat = coords_stereo_pred.lat.value

    [sta.x, sta.y, sta.z] = sphere2cart(sta.r, np.deg2rad(-sta.lat+90), np.deg2rad(sta.lon))

    print('position end ')

    sta_gse=convert_RTN_to_GSE_sta_l1(sta) #here is the bug in the conversion 
    sta_gsm=convert_GSE_to_GSM(sta_gse)
    
    header='STEREO-A beacon magnetic field (IMPACT instrument) and beacon plasma data (PLASTIC), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/impact/ and   '+ \
    'https://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/plastic/ '+ \
    'Timerange: '+sta.time[0].strftime("%Y-%b-%d %H:%M")+' to '+sta.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', with an average time resolution of '+str(np.mean(np.diff(sta.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by sta.time, sta.bx, sta.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(sta.size)+'. '+\
    'Units are btxyz [nT, RTN], vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats '+\
    'By C. Moestl, E. Weiler, E. E. Davies. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    print('save pickle file')
    pickle.dump([sta,header], open(data_path+file_rtn, "wb"))
    
    
    
    header_gsm='STEREO-A beacon magnetic field (IMPACT instrument) and beacon plasma data (PLASTIC), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/impact/ and   '+ \
    'https://stereo-ssc.nascom.nasa.gov/data/beacon/ahead/plastic/ '+ \
    'Timerange: '+sta.time[0].strftime("%Y-%b-%d %H:%M")+' to '+sta.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', with an average time resolution of '+str(np.mean(np.diff(sta.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by sta.time, sta.bx, sta.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(sta.size)+'. '+\
    'Units are btxyz [nT, GSM], vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats '+\
    'By C. Moestl, E. Weiler, E. E. Davies. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    print('save pickle file')
    pickle.dump([sta_gsm,header_gsm], open(data_path+file_gsm, "wb"))
    
    
    
    print('done sta')
    
            
            
            
            
            

def download_stereoa_science_merged(start_timestamp, end_timestamp, path="/Users/emmadavies/Documents/Data-test"):
    
    ######## TO DO adjust from Emmas code
    
    
    from datetime import datetime, timedelta
    import urllib.request
    import os.path
    
    start = start_timestamp.year
    end = end_timestamp.year + 1
    while start < end:
        year = start
        date_str = f'{year}0101'
        for v in range(6, 0, -1):
            data_item_id = f'sta_l2_magplasma_1m_{date_str}_v{v:02}'
            if os.path.isfile(f"{path}/{data_item_id}.cdf") == True:
                print(f'{data_item_id}.cdf has already been downloaded.')
            else:
                try:
                    data_url = f'https://spdf.gsfc.nasa.gov/pub/data/stereo/ahead/l2/impact/magplasma/1min/{year}/{data_item_id}.cdf'
                    urllib.request.urlretrieve(data_url, f"{path}/{data_item_id}.cdf")
                    print(f'Successfully downloaded {data_item_id}.cdf')
                except Exception as e:
                    print('ERROR', e, data_item_id)
        start += 1  
            

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
#---------------------------------------------------------------            
            
            

    

def wind_download_ascii_old_wget_curl():

    #download MFI and SWE in ASCII from SPDF
        
    print('downloading Wind ascii data')
        
   
    read_data_end_year=datetime.datetime.utcnow().year
    read_data_end_month=datetime.datetime.utcnow().month
    
    wind_years_strings=[]
    for j in np.arange(2022,read_data_end_year+1):
        wind_years_strings.append(str(j))

    print(wind_years_strings)
        
    wind_data_path='/perm/aswo/data/wind/mfi_1min_ascii'
    if mac   == 1: wind_data_path='/Users/chris/python/data/wind/mfi_1min_ascii'

    os.chdir(wind_data_path)

    mfi_url='https://spdf.gsfc.nasa.gov/pub/data/wind/mfi/ascii/1min_ascii/'
    print(mfi_url)
    
    #for all years
    for i in np.arange(0,len(wind_years_strings)-1):    

        for k in np.arange(1,13):    

            a=str(k).zfill(2) #add leading zeros
            print(wind_years_strings[i],' ',a)

            
#                omni2_url='https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
#    print(omni2_url)
#    try: urllib.request.urlretrieve(omni2_url, 'data/omni2_all_years.dat')
#    except urllib.error.URLError as e:
#        print(' ', omni2_url,' ',e.reason) 

            
            
            if linux == 1: os.system('wget -nc '+mfi_url+wind_years_strings[i]+a+'_wind_mag_1min.asc') #199504_wind_mag_1min.asc	
            if mac   == 1: os.system('curl -s -O -C - '+mfi_url+wind_years_strings[i]+a+'_wind_mag_1min.asc') 

    #for latest year
    for k in np.arange(1,read_data_end_month):    

        a=str(k).zfill(2) #add leading zeros
        print(wind_years_strings[-1],' ',a)
                      
            
        if linux == 1: os.system('wget -nc '+mfi_url+wind_years_strings[i]+a+'_wind_mag_1min.asc') 
        if mac   == 1: os.system('curl -s -O -C - '+mfi_url+wind_years_strings[-1]+a+'_wind_mag_1min.asc') 


    
    #download plasma ascii data

    #linux
    #wind_data_path='/perm/aswo/data/wind/swe_92sec_ascii'
        
    #mac
    wind_data_path='/Users/chris/python/data/wind/swe_92sec_ascii'
    
    os.chdir(wind_data_path)

    swe_url='https://spdf.gsfc.nasa.gov/pub/data/wind/swe/ascii/swe_kp_unspike/'
    print(swe_url)

    for i in np.arange(0,len(wind_years_strings)):    
        
            print(wind_years_strings[i])

            if linux == 1: os.system('wget -nc '+swe_url+'wind_kp_unspike'+wind_years_strings[i]+'.txt')
            if mac   == 1: os.system('curl -s -O -C - '+swe_url+'wind_kp_unspike'+wind_years_strings[i]+'.txt')

    if linux == 1: os.chdir('/export/home/aswo/heliofc/code/heliocats')
    if mac   == 1: os.chdir('/Users/chris/python/heliocats')
        
        
    
    






def remove_wind_spikes_gaps(data):
    
    #nan intervals

    nt1=parse_time('2020-04-20 17:06').datetime 
    nt2=parse_time('2020-04-20 17:14').datetime
    gapind1=np.where(np.logical_and(data.time >= nt1,data.time <= nt2 ))[0]

    nt1=parse_time('2020-04-21 01:20').datetime 
    nt2=parse_time('2020-04-21 01:22').datetime
    gapind2=np.where(np.logical_and(data.time >= nt1,data.time <= nt2 ))[0]

    nt1=parse_time('2020-11-09T16:04Z').datetime 
    nt2=parse_time('2020-11-09T17:08Z').datetime
    gapind3=np.where(np.logical_and(data.time >= nt1,data.time <= nt2 ))[0]

    nt1=parse_time('2020-08-31T16:58Z').datetime 
    nt2=parse_time('2020-08-31T18:32Z').datetime
    gapind4=np.where(np.logical_and(data.time >= nt1,data.time <= nt2 ))[0]
    
    nt1=parse_time('2021-02-01T12:32Z').datetime 
    nt2=parse_time('2021-02-01T14:04Z').datetime
    gapind5=np.where(np.logical_and(data.time >= nt1,data.time <= nt2 ))[0]
    
    nt1=parse_time('2022-02-04T07:40Z').datetime 
    nt2=parse_time('2022-02-04T08:14Z').datetime
    gapind6=np.where(np.logical_and(data.time >= nt1,data.time <= nt2 ))[0]
    
    nt1=parse_time('2021-12-10T17:18Z').datetime 
    nt2=parse_time('2021-12-11T14:52Z').datetime
    gapind7=np.where(np.logical_and(data.time >= nt1,data.time <= nt2 ))[0]    
    
    nt1=parse_time('2021-12-14T17:02Z').datetime 
    nt2=parse_time('2021-12-15T13:12Z').datetime
    gapind8=np.where(np.logical_and(data.time >= nt1,data.time <= nt2 ))[0]

    nt1=parse_time('2021-12-29T16:56Z').datetime 
    nt2=parse_time('2021-12-30T02:54Z').datetime
    gapind9=np.where(np.logical_and(data.time >= nt1,data.time <= nt2 ))[0]

    
    

    data.bt[np.hstack([gapind1,gapind2,gapind3,gapind4,gapind5,gapind6,gapind7,gapind8,gapind9])]=np.nan
    data.bx[np.hstack([gapind1,gapind2,gapind3,gapind4,gapind5,gapind6,gapind7,gapind8,gapind9])]=np.nan
    data.by[np.hstack([gapind1,gapind2,gapind3,gapind4,gapind5,gapind6,gapind7,gapind8,gapind9])]=np.nan
    data.bz[np.hstack([gapind1,gapind2,gapind3,gapind4,gapind5,gapind6,gapind7,gapind8,gapind9])]=np.nan
    
    
    #data.bt[np.hstack([gapind1,gapind2,gapind3,gapind4,gapind5,gapind6,gapind7,gapind8,gapind9])]=np.nan
    #data.bx[np.hstack([gapind1,gapind2,gapind3,gapind4,gapind5,gapind6,gapind7,gapind8,gapind9])]=np.nan
    #data.by[np.hstack([gapind1,gapind2,gapind3,gapind4,gapind5,gapind6,gapind7,gapind8,gapind9])]=np.nan
    #data.bz[np.hstack([gapind1,gapind2,gapind3,gapind4,gapind5,gapind6,gapind7,gapind8,gapind9])]=np.nan
    
    
    return data



    
    
    
    



    
    
    
    
    
    
    
    
    
    
    
    
    
def read_wind_icme_catalog():    
    
    #read Wind ICME catalog copied into text file: 


    #load master file for ICMECAT    
    icm=hc.load_helcats_icmecat_master_from_excel('icmecat/HELIO4CAST_ICMECAT_v21_master.xlsx')
    lenicm=len(icm)

    wind_mc_cat='data/wind_mc_cat/wind_mc_cat.txt'

    #extract 3 times for icmecat
    cat=np.genfromtxt(wind_mc_cat)

    file1 = open(wind_mc_cat, 'r')
    lines = file1.readlines()

    for i in np.arange(0,len(lines)):
        t1=lines[i][0:16]
        t11=t1[0:4]+'-'+t1[5:7]+'-'+t1[8:17]
        t1d=parse_time(t11).datetime

        t2=lines[i][19:35]
        t21=t2[0:4]+'-'+t2[5:7]+'-'+t2[8:17]
        t2d=parse_time(t21).datetime

        if lines[i][36]=='E':
            t3=lines[i][38:54]
        else:
            t3=lines[i][39:55]

        t31=t3[0:4]+'-'+t3[5:7]+'-'+t3[8:17]
        t3d=parse_time(t31).datetime

        #print('ICME_WIND_NASA_'+t11[0:4]+t11[5:7]+t11[8:10]+'_01')
        #append to datafrme
        icm.loc[lenicm+i,'icmecat_id']='ICME_Wind_NASA_'+t11[0:4]+t11[5:7]+t11[8:10]+'_01'
        icm.loc[lenicm+i,'sc_insitu']='Wind'
        icm.loc[lenicm+i,'icme_start_time']=t1d
        icm.loc[lenicm+i,'mo_start_time']=t2d
        icm.loc[lenicm+i,'mo_end_time']=t3d




        #icm.loc[len(icm),1] = [' ']#],'Wind',t11,t21,t31]    
        #icm.loc[len(icm),2] = ['Wind']



        #print(parse_time(t11).datetime)
        #print(parse_time(t21).datetime)
        #print(parse_time(t31).datetime)
        #print('---')


    #check if there is a double id


    ################ save to different formats

    ic=icm
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


    #save as Excel with all wind events appended
    file='icmecat/HELIO4CAST_ICMECAT_v21_master_allwind_from_cat.xlsx'
    ic_copy.to_excel(file,sheet_name='ICMECATv2.1')
    print('ICMECAT saved as '+file)




    
    
    
    
    
    
    





    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


























def save_stereoa_science_data_merge_rtn(data_path,file):
    
    print('STEREO-A science data merging')
    filesta="stereoa_2007_2019_rtn.p"
    [sta0,hsta0]=pickle.load(open(data_path+filesta, "rb" ) )  

    filesta="stereoa_2020_april_rtn.p" 
    [sta1,hsta1]=pickle.load(open(data_path+filesta, "rb" ) )  

    filesta="stereoa_2020_may_july_rtn.p" 
    [sta2,hsta2]=pickle.load(open(data_path+filesta, "rb" ) )  

    #beacon data
    #filesta='stereoa_2019_now_sceq_beacon.p'
    #[sta3,hsta3]=pickle.load(open(data_path+filesta2, "rb" ) )  
    #sta2=sta2[np.where(sta2.time >= parse_time('2020-Aug-01 00:00').datetime)[0]]


    #make array
    sta=np.zeros(np.size(sta0.time)+np.size(sta1.time)+np.size(sta2.time),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    sta = sta.view(np.recarray)  
    sta.time=np.hstack((sta0.time,sta1.time,sta2.time))
    sta.bx=np.hstack((sta0.bx,sta1.bx,sta2.bx))
    sta.by=np.hstack((sta0.by,sta1.by,sta2.by))
    sta.bz=np.hstack((sta0.bz,sta1.bz,sta2.bz))
    sta.bt=np.hstack((sta0.bt,sta1.bt,sta2.bt))
    sta.vt=np.hstack((sta0.vt,sta1.vt,sta2.vt))
    sta.np=np.hstack((sta0.np,sta1.np,sta2.np))
    sta.tp=np.hstack((sta0.tp,sta1.tp,sta2.tp))
    sta.x=np.hstack((sta0.x,sta1.x,sta2.x))
    sta.y=np.hstack((sta0.y,sta1.y,sta2.y))
    sta.z=np.hstack((sta0.z,sta1.z,sta2.z))
    sta.r=np.hstack((sta0.r,sta1.r,sta2.r))
    sta.lon=np.hstack((sta0.lon,sta1.lon,sta2.lon))
    sta.lat=np.hstack((sta0.lat,sta1.lat,sta2.lat))


    pickle.dump(sta, open(data_path+file, "wb"))
    print('STEREO-A merging done')
   
    
    
    return 0


def save_stereoa_science_data_merge_sceq(data_path,file):
    
    print('STEREO-A science data merging')
    filesta="stereoa_2007_2019_sceq.p"
    [sta0,hsta0]=pickle.load(open(data_path+filesta, "rb" ) )  

    filesta="stereoa_2020_april_sceq.p" 
    [sta1,hsta1]=pickle.load(open(data_path+filesta, "rb" ) )  

    filesta="stereoa_2020_may_july_sceq.p" 
    [sta2,hsta2]=pickle.load(open(data_path+filesta, "rb" ) )  

    #beacon data
    #filesta='stereoa_2019_now_sceq_beacon.p'
    #[sta3,hsta3]=pickle.load(open(data_path+filesta2, "rb" ) )  
    #sta2=sta2[np.where(sta2.time >= parse_time('2020-Aug-01 00:00').datetime)[0]]


    #make array
    sta=np.zeros(np.size(sta0.time)+np.size(sta1.time)+np.size(sta2.time),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    sta = sta.view(np.recarray)  
    sta.time=np.hstack((sta0.time,sta1.time,sta2.time))
    sta.bx=np.hstack((sta0.bx,sta1.bx,sta2.bx))
    sta.by=np.hstack((sta0.by,sta1.by,sta2.by))
    sta.bz=np.hstack((sta0.bz,sta1.bz,sta2.bz))
    sta.bt=np.hstack((sta0.bt,sta1.bt,sta2.bt))
    sta.vt=np.hstack((sta0.vt,sta1.vt,sta2.vt))
    sta.np=np.hstack((sta0.np,sta1.np,sta2.np))
    sta.tp=np.hstack((sta0.tp,sta1.tp,sta2.tp))
    sta.x=np.hstack((sta0.x,sta1.x,sta2.x))
    sta.y=np.hstack((sta0.y,sta1.y,sta2.y))
    sta.z=np.hstack((sta0.z,sta1.z,sta2.z))
    sta.r=np.hstack((sta0.r,sta1.r,sta2.r))
    sta.lon=np.hstack((sta0.lon,sta1.lon,sta2.lon))
    sta.lat=np.hstack((sta0.lat,sta1.lat,sta2.lat))




    pickle.dump(sta, open(data_path+file, "wb"))
    print('STEREO-A merging done')
   
    
    



    
def save_stereoa_science_data(path,file,t_start, t_end,sceq):


    #impact https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/ahead/ 
    #download with heliosat
    #-------------------    
    #print('start STA')
    #sta_sat = heliosat.STA()
     
    #create an array with 1 minute resolution between t start and end
    #time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    #time_mat=mdates.date2num(time) 
    
    #tm, mag = sta_sat.get_data_raw(t_start, t_end, "sta_impact_l1")
    #print('download complete')
    #---------------------------
    
    
    #2020 PLASTIC download manually
    #https://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/Derived_from_1D_Maxwellian/ASCII/1min/A/2020/
            
    sta_impact_path='/nas/helio/data/heliosat/data/sta_impact_l1/'
    sta_plastic_path='/nas/helio/data/heliosat/data/sta_plastic_l2_ascii/'


    t_start1=copy.deepcopy(t_start)
    time_1=[]
    #make 1 min datetimes
    while t_start1 < t_end:
        time_1.append(t_start1)  
        t_start1 += datetime.timedelta(minutes=1)


    #make array for 1 min data
    sta=np.zeros(len(time_1),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    sta = sta.view(np.recarray)  
    sta.time=time_1

    #make data file names
    t_start1=copy.deepcopy(t_start)
    days_sta = []
    days_str = []
    i=0
    while t_start < t_end:
        days_sta.append(t_start)  
        days_str.append(str(days_sta[i])[0:4]+str(days_sta[i])[5:7]+str(days_sta[i])[8:10])        
        i=i+1
        t_start +=datetime.timedelta(days=1)

    #go through all files

    bt=np.zeros(int(1e9))
    bx=np.zeros(int(1e9))
    by=np.zeros(int(1e9))
    bz=np.zeros(int(1e9))
    t2=[]
    i=0
    for days_date in days_str:
        cdf_file = 'STA_L1_MAG_RTN_{}_V06.cdf'.format(days_date)
        
        if os.path.exists(sta_impact_path+cdf_file):
            print(cdf_file)
            f1 = cdflib.CDF(sta_impact_path+cdf_file)
            t1=parse_time(f1.varget('Epoch'),format='cdf_epoch').datetime
            t2.extend(t1)
            bfield=f1.varget('BFIELD')
            bt[i:i+len(bfield[:,3])]=bfield[:,3]
            bx[i:i+len(bfield[:,0])]=bfield[:,0]
            by[i:i+len(bfield[:,1])]=bfield[:,1]
            bz[i:i+len(bfield[:,2])]=bfield[:,2]
            i=i+len(bfield[:,3])

    #cut array
    bt=bt[0:i]
    bx=bx[0:i]
    by=by[0:i]
    bz=bz[0:i]

    tm2=mdates.date2num(t2)
    time_mat=mdates.date2num(time_1)


    #linear interpolation to time_mat times    
    sta.bx = np.interp(time_mat, tm2, bx )
    sta.by = np.interp(time_mat, tm2, by )
    sta.bz = np.interp(time_mat, tm2, bz )
    #sta.bt = np.sqrt(sta.bx**2+sta.by**2+sta.bz**2)
    
    
    
    #round first each original time to full minutes   original data at 30sec
    tround=copy.deepcopy(t2)
    format_str = '%Y-%m-%d %H:%M'  
    for k in np.arange(np.size(t2)):
         tround[k] = datetime.datetime.strptime(datetime.datetime.strftime(t2[k], format_str), format_str) 
    tm2_round=parse_time(tround).plot_date

    #which values are not in original data compared to full time range
    isin=np.isin(time_mat,tm2_round)      
    setnan=np.where(isin==False)
    #set to to nan that is not in original data
    sta.bx[setnan]=np.nan
    sta.by[setnan]=np.nan
    sta.bz[setnan]=np.nan
    sta.bt = np.sqrt(sta.bx**2+sta.by**2+sta.bz**2)
    


    ########### get PLASTIC new prel data
    #PLASTIC
    #2019 monthly if needed
    #https://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/Derived_from_1D_Maxwellian/ASCII/1min/A/2019/

    #2020 manually all
    #https://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/Derived_from_1D_Maxwellian/ASCII/1min/A/2020/

    #STA_L2_PLA_1DMax_1min_202004_092_PRELIM_v01.txt
    #STA_L2_PLA_1DMax_1min_202005_122_PRELIM_v01.txt
    #STA_L2_PLA_1DMax_1min_202006_153_PRELIM_v01.txt
    #STA_L2_PLA_1DMax_1min_202007_183_PRELIM_v01.txt





    ########
    pvt=np.zeros(int(1e8))
    pnp=np.zeros(int(1e8))
    ptp=np.zeros(int(1e8))
    pt2=[]




    pfiles=['STA_L2_PLA_1DMax_1min_202004_092_PRELIM_v01.txt',     
         'STA_L2_PLA_1DMax_1min_202005_122_PRELIM_v01.txt',
         'STA_L2_PLA_1DMax_1min_202006_153_PRELIM_v01.txt',
         'STA_L2_PLA_1DMax_1min_202007_183_PRELIM_v01.txt']

    j=0
    for name in pfiles:

        p1=np.genfromtxt(sta_plastic_path+name,skip_header=2)
        print(name)

        vt1=p1[:,8]
        np1=p1[:,9]
        tp1=p1[:,10]

        #YEAR	DOY	hour	min	sec
        year1=p1[:,0]
        doy1=p1[:,1]
        hour1=p1[:,2]
        min1=p1[:,3]
        sec1=p1[:,4]



        p1t=[]
        #make datetime array from year and doy
        for i in np.arange(len(doy1)):
            p1t.append(parse_time(str(int(year1[i]))+'-01-01 00:00').datetime+datetime.timedelta(days=doy1[i]-1)+\
                       +datetime.timedelta(hours=hour1[i]) + datetime.timedelta(minutes=min1[i])  )


        pvt[j:j+len(vt1)]=vt1
        pnp[j:j+len(np1)]=np1
        ptp[j:j+len(tp1)]=tp1

        pt2.extend(p1t)

        j=j+len(vt1)

    #cut array
    pvt=pvt[0:j]
    pnp=pnp[0:j]
    ptp=ptp[0:j]
    pt2=pt2[0:j]



    pt2m=mdates.date2num(pt2)


    #linear interpolation to time_mat times    
    sta.vt = np.interp(time_mat, pt2m, pvt )
    sta.np = np.interp(time_mat, pt2m, pnp )
    sta.tp = np.interp(time_mat, pt2m, ptp )



    #which values are not in original data compared to full time range
    isin=np.isin(time_mat,pt2m)      
    setnan=np.where(isin==False)
    #set to to nan that is not in original data
    sta.vt[setnan]=np.nan
    sta.np[setnan]=np.nan
    sta.tp[setnan]=np.nan

    #add position
    
    print('position start')
    frame='HEEQ'
    kernels = spicedata.get_kernel('stereo_a')
    kernels += spicedata.get_kernel('stereo_a_pred')
    spice.furnish(kernels)
    statra=spice.Trajectory('-234') #STEREO-A SPICE NAIF code
    statra.generate_positions(sta.time,'Sun',frame)
    statra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(statra.x,statra.y,statra.z)
    
    sta.x=statra.x
    sta.y=statra.y
    sta.z=statra.z
    
    sta.r=r
    sta.lat=np.degrees(lat)
    sta.lon=np.degrees(lon)

    print('position end ')

    
        
    coord='RTN'
    #convert magnetic field to SCEQ
    if sceq==True:
        print('convert RTN to SCEQ ')
        coord='SCEQ'
        sta=convert_RTN_to_SCEQ(sta,'STEREO-A')
    
    header='STEREO-A magnetic field (IMPACT instrument, science data) and plasma data (PLASTIC, preliminary science data), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/ahead/ and   '+ \
    'https://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/Derived_from_1D_Maxwellian/ASCII/1min/A/2020/ '+ \
    'Timerange: '+sta.time[0].strftime("%Y-%b-%d %H:%M")+' to '+sta.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', with an average time resolution of '+str(np.mean(np.diff(sta.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by sta.time, sta.bx, sta.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(sta.size)+'. '+\
    'Units are btxyz [nT, '+coord+', vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats '+\
    'and https://github.com/heliopython/heliopy. '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, R. L. Bailey and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    print('save pickle file')
    pickle.dump([sta,header], open(path+file, "wb"))
    
    print('done sta')
    print()

   
    

    return 0












































def save_wsa_hux(filein):
    #load wsa hux 

    windraw = np.loadtxt('data/wsa_hux_mars_aug2014_jan2018.txt', dtype=[('time','<U30'),('time2','<U30'),('time_mat', float),('vt', float)] )
    windraw = windraw.view(np.recarray)  

    wind=np.zeros(len(windraw),dtype=[('time',object),('vt', float)])   
    wind=wind.view(np.recarray) 


    for i in np.arange(len(windraw)):
        wind_time_str=windraw.time[i][8:12]+'-'+windraw.time[i][4:7]+'-'+windraw.time[i][1:3]+' '+windraw.time2[i][0:8]
        wind.time[i]=(parse_time(wind_time_str).datetime)

    wind.vt=windraw.vt
    
    fileout='wsa_hux_mars_aug2014_jan2018.p'      
    pickle.dump(wind, open(data_path+fileout, "wb"))
    
    return 0



def load_mars_wsa_hux():
    
    data_path='data/'
    file='wsa_hux_mars_aug2014_jan2018.p'  
    rad=pickle.load(open(data_path+file, "rb"))
    
    return rad





def load_maven_sir_huang():
    
    
    #Huang et al. 2019 APJ convert PDF to excel with https://pdftoxls.com
    
    
   
    mavensir='sircat/sources/Huang_2019_SIR_MAVEN_table_1.xlsx'
    
    print('load MAVEN Huang SIR catalog from ', mavensir)
    ms=pd.read_excel(mavensir)
    ms=ms.drop(index=[0,1,2])
    ms_num=np.array(ms['No.'])
    ms_start=np.array(ms['Start'])
    ms_end=np.array(ms['End'])
    ms_si=np.array(ms['SI'])
    
    
    
    ms=np.zeros(len(ms_num),dtype=[('start',object),('end',object),('si',object)])   
    ms=ms.view(np.recarray) 
        
        
    #make correct years for start time 
    ms_num[np.where(ms_num< 7)[0]]=2014
    ms_num[np.where(ms_num< 27)[0]]=2015
    ms_num[np.where(ms_num< 64)[0]]=2016
    ms_num[np.where(ms_num< 83)[0]]=2017
    ms_num[np.where(ms_num< 127)[0]]=2018
    
    #make correct years for end and si time
    ms_num2=copy.deepcopy(ms_num)
    ms_num2[3]=2015
    ms_num2[62]=2017

    #transform date of start time
    for t in np.arange(0,len(ms_start)):
        #check for nans in between time strings
        if pd.isna(ms_start[t])==False:                    
            
            
            ####################### start time
            #year
            year=str(ms_num[t])
            
            #month
            datetimestr=ms_start[t]    
            datestr=datetimestr[0:2]
            
            monthfloat=float(datestr)
            month=str(int(np.floor(monthfloat)))
            
            #day            
            if int(month) < 10: day=datetimestr[2:4]
            if int(month) > 9: day=datetimestr[3:5]           
            
            #time
            timestr=datetimestr[-5:]                        
            
            #construct year month day
            datetimestrfin=str(ms_num[t])+'-'+month+'-'+day                    
            #remove white spaces at the end and add time
            finaldatetime=datetimestrfin.strip()+' '+timestr
            #print(ms_start[t])
            #print(finaldatetime)            
            ms.start[t]=parse_time(finaldatetime).datetime
            
            
            
            
            ################### end time            

            #year
            year=str(ms_num2[t])

            #month
            datetimestr=ms_end[t]    
            datestr=datetimestr[0:2]
            
            monthfloat=float(datestr)
            month=str(int(np.floor(monthfloat)))
            
            #day            
            if int(month) < 10: day=datetimestr[2:4]
            if int(month) > 9: day=datetimestr[3:5]           
            
            #time
            timestr=datetimestr[-5:]                        
            
            #construct year month day
            datetimestrfin=str(ms_num2[t])+'-'+month+'-'+day                    
            #remove white spaces at the end and add time
            finaldatetime=datetimestrfin.strip()+' '+timestr
            #print(ms_end[t])
            #print(finaldatetime)            
            ms.end[t]=parse_time(finaldatetime).datetime
            
            
                        
            ############# stream interface time            

            #year
            year=str(ms_num2[t])

            #month
            datetimestr=ms_si[t]    
            datestr=datetimestr[0:2]
            
            monthfloat=float(datestr)
            month=str(int(np.floor(monthfloat)))
            
            #day            
            if int(month) < 10: day=datetimestr[2:4]
            if int(month) > 9: day=datetimestr[3:5]           
            
            #time
            timestr=datetimestr[-5:]                        
            
            #construct year month day
            datetimestrfin=str(ms_num2[t])+'-'+month+'-'+day                    
            #remove white spaces at the end and add time
            finaldatetime=datetimestrfin.strip()+' '+timestr
            #print(ms_si[t])
            #print(finaldatetime)            
            ms.si[t]=parse_time(finaldatetime).datetime
            #print()


    #get rid of zeros where the years where stated in the original data        
    ms2 = ms[np.argwhere(ms)]

    return ms2





def save_msl_rad():

    #convert file 
    # year, doy, sol, doseE hourly [uGy/day], doseE sol-filtered [uGy/day]
    raw=np.loadtxt('data/doseE_sol_filter_2019.dat')
    
    rad=np.zeros(len(raw),dtype=[('time',object),('sol', float),('dose_hour', float),('dose_sol', float)])   
    rad = rad.view(np.recarray) 
    
    rad.sol=raw[:,2]
    rad.dose_hour=raw[:,3]
    rad.dose_sol=raw[:,4]


    #make datetime array from year and doy
    for i in np.arange(len(rad)):
        rad[i].time=parse_time(str(int(raw[i,0]))+'-01-01 00:00').datetime+datetime.timedelta(days=raw[i,1]-1)
        print(rad[i].time)
        
    file='msl_2012_2019_rad.p'  
    pickle.dump(rad, open(data_path+file, "wb"))
    
    return 0



def load_msl_rad(data_path):

    file='msl_2012_2019_rad.p'  
    rad=pickle.load(open(data_path+file, "rb"))
    
    return rad



def save_psp_data_mag_only(path, file, sceq):
        
    print('save PSP data mag only')

    t_start = datetime.datetime(2022, 1, 1)
    t_end = datetime.datetime(2022, 7, 31)
    psp=get_psp_data_mag_only(t_start,t_end)
    
    
    #convert magnetic field to SCEQ
    coord='RTN'
    if sceq==True:
        coord='SCEQ'
        psp=convert_RTN_to_SCEQ(psp,'PSP')
    
    
    header='PSP magnetic field (FIELDS instrument), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/psp/  '+ \
    'Timerange: '+psp.time[0].strftime("%Y-%b-%d %H:%M")+' to '+psp.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(psp.time)).seconds)+' seconds. '+\
    'The data are put in a numpy recarray, fields can be accessed by psp.time, psp.bx, psp.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(psp.size)+'. '+\
    'Units are btxyz [nT,'+coord+'], vtxyz [km/s, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([psp,header], open(path+file, "wb"))



def get_psp_data_mag_only(t_start,t_end):
     
    print('start PSP')
     
    psp_sat = heliosat.PSP()
    
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    tm, mag = psp_sat.get_data(time, "psp_fields_l2",return_datetimes=True)
    #tp, pro = psp_sat.get_data(time, "psp_spc_l3",return_datetimes=True)
    
    
    #tm=parse_time(tm,format='unix').datetime 
    #tp=parse_time(tp,format='unix').datetime 
    
    print('download complete')
    
    print('start nan or interpolate')
    
    print('field')
    #round first each original time to full minutes   original data at 30sec
    tround=copy.deepcopy(tm)
    format_str = '%Y-%m-%d %H:%M'  
    for k in np.arange(np.size(tm)):
         tround[k] = datetime.datetime.strptime(datetime.datetime.strftime(tm[k], format_str), format_str) 
    tm_mat=parse_time(tround).plot_date
    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    
   
    #which values are not in original data compared to full time range
    isin=np.isin(time_mat,tm_mat)      
    setnan=np.where(isin==False)
    #set to to nan that is not in original data
    bx[setnan]=np.nan
    by[setnan]=np.nan
    bz[setnan]=np.nan

    bt = np.sqrt(bx**2+by**2+bz**2)
    

    
    #    print('plasma')
        
    #for plasma round first each original time to full minutes
    #tround=copy.deepcopy(tp)
    #format_str = '%Y-%m-%d %H:%M'  
    #for k in np.arange(np.size(tp)):
    #     tround[k] = datetime.datetime.strptime(datetime.datetime.strftime(tp[k], format_str), format_str) 
    #tp_mat=mdates.date2num(tround) 
    
    #isin=np.isin(time_mat,tp_mat)      
    #setnan=np.where(isin==False)
        
    #den = np.interp(time_mat, tp_mat, pro[:,0])
    #vx = np.interp(time_mat, tp_mat, pro[:,1])
    #vy = np.interp(time_mat, tp_mat, pro[:,2])
    #vz = np.interp(time_mat, tp_mat, pro[:,3])
    #temp = np.interp(time_mat, tp_mat, pro[:,4])
    
    den = np.zeros(np.size(time_mat))
    vx = np.zeros(np.size(time_mat))
    vy = np.zeros(np.size(time_mat))
    vz = np.zeros(np.size(time_mat))
    temp = np.zeros(np.size(time_mat))
    
    den[:]=np.nan
    temp[:]=np.nan
    vx[:]=np.nan
    vy[:]=np.nan
    vz[:]=np.nan
  
    vt=np.sqrt(vx**2+vy**2+vz**2)

    print('end nan or interpolate')

        
    print('position start')
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('psp_pred'))
    psptra=spice.Trajectory('SPP')
    psptra.generate_positions(time,'Sun',frame)
    psptra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(psptra.x,psptra.y,psptra.z)
    print('PSP pos')    
    print('position end')

        
    #make array
    psp=np.zeros(np.size(bx),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('np', float),('tp', float),('x', float),('y', float),\
                ('z', float),('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    psp = psp.view(np.recarray)  

    #fill with data
    psp.time=time
    psp.bx=bx
    psp.by=by
    psp.bz=bz 
    psp.bt=bt
    
    psp.x=psptra.x
    psp.y=psptra.y
    psp.z=psptra.z
    
    psp.r=r
    psp.lat=np.degrees(lat)
    psp.lon=np.degrees(lon)
    
    psp.vt=vt
    psp.vx=vx    
    psp.vy=vy  
    psp.vz=vz
    psp.np=den
    #https://en.wikipedia.org/wiki/Thermal_velocity convert from km/s to K
    #from astropy.constants import m_p,k_B
    #psp.tp=np.pi*m_p*((temp*1e3)**2)/(8*k_B) 
    
    #remove spikes
    #psp.tp[np.where(psp.tp > 1e10)]=np.nan
    
    print('done get psp')
    print()
    
    return psp

    
    
    





def save_psp_data(path, file, sceq):
    
    
    print('save PSP data')

    t_start = datetime.datetime(2018, 10, 6)
    t_end = datetime.datetime(2019, 4, 24) #  UNTIL ERROR on Apr 25
    psp1=get_psp_data(t_start,t_end)

    t_start = datetime.datetime(2019, 4, 26)    
    #t_end = datetime.datetime(2019, 4, 30)    
    #t_end = datetime.datetime(2019, 10, 15)
    
    t_end = datetime.datetime(2022, 7, 31)    
    psp2=get_psp_data(t_start,t_end)

    #add both
    
    psp=np.zeros(np.size(psp1.time)+np.size(psp2.time),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('vx', float),('vy', float),('vz', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    psp = psp.view(np.recarray)  
    psp.time=np.hstack((psp1.time,psp2.time))
    psp.bx=np.hstack((psp1.bx,psp2.bx))
    psp.by=np.hstack((psp1.by,psp2.by))
    psp.bz=np.hstack((psp1.bz,psp2.bz))
    psp.bt=np.hstack((psp1.bt,psp2.bt))
    psp.vt=np.hstack((psp1.vt,psp2.vt))
    psp.vx=np.hstack((psp1.vx,psp2.vx))
    psp.vy=np.hstack((psp1.vy,psp2.vy))
    psp.vz=np.hstack((psp1.vz,psp2.vz))
    psp.np=np.hstack((psp1.np,psp2.np))
    psp.tp=np.hstack((psp1.tp,psp2.tp))
    psp.x=np.hstack((psp1.x,psp2.x))
    psp.y=np.hstack((psp1.y,psp2.y))
    psp.z=np.hstack((psp1.z,psp2.z))
    psp.r=np.hstack((psp1.r,psp2.r))
    psp.lon=np.hstack((psp1.lon,psp2.lon))
    psp.lat=np.hstack((psp1.lat,psp2.lat))

    print('Merging done')
    
    
    #convert magnetic field to SCEQ
    coord='RTN'
    if sceq==True:
        coord='SCEQ'
        psp=convert_RTN_to_SCEQ(psp,'PSP')
    
    
    header='PSP magnetic field (FIELDS instrument) and plasma data (SWEAP), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/psp/  '+ \
    'Timerange: '+psp.time[0].strftime("%Y-%b-%d %H:%M")+' to '+psp.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(psp.time)).seconds)+' seconds. '+\
    'The data are put in a numpy recarray, fields can be accessed by psp.time, psp.bx, psp.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(psp.size)+'. '+\
    'Units are btxyz [nT,'+coord+'], vtxyz [km/s, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([psp,header], open(path+file, "wb"))




def get_psp_data(t_start,t_end):
     
    print('start PSP')
     
    psp_sat = heliosat.PSP()
    
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    tm, mag = psp_sat.get_data(time, "psp_fields_l2",return_datetimes=True)
    tp, pro = psp_sat.get_data(time, "psp_spc_l3",return_datetimes=True)
    
    
    #tm=parse_time(tm,format='unix').datetime 
    #tp=parse_time(tp,format='unix').datetime 
    
    print('download complete')
    
    print('start nan or interpolate')
    
    print('field')
    #round first each original time to full minutes   original data at 30sec
    tround=copy.deepcopy(tm)
    format_str = '%Y-%m-%d %H:%M'  
    for k in np.arange(np.size(tm)):
         tround[k] = datetime.datetime.strptime(datetime.datetime.strftime(tm[k], format_str), format_str) 
    tm_mat=parse_time(tround).plot_date
    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    
   
    #which values are not in original data compared to full time range
    isin=np.isin(time_mat,tm_mat)      
    setnan=np.where(isin==False)
    #set to to nan that is not in original data
    bx[setnan]=np.nan
    by[setnan]=np.nan
    bz[setnan]=np.nan

    bt = np.sqrt(bx**2+by**2+bz**2)
    

    
    print('plasma')
        
    #for plasma round first each original time to full minutes
    tround=copy.deepcopy(tp)
    format_str = '%Y-%m-%d %H:%M'  
    for k in np.arange(np.size(tp)):
         tround[k] = datetime.datetime.strptime(datetime.datetime.strftime(tp[k], format_str), format_str) 
    tp_mat=mdates.date2num(tround) 
    
    isin=np.isin(time_mat,tp_mat)      
    setnan=np.where(isin==False)
        
    den = np.interp(time_mat, tp_mat, pro[:,0])
    vx = np.interp(time_mat, tp_mat, pro[:,1])
    vy = np.interp(time_mat, tp_mat, pro[:,2])
    vz = np.interp(time_mat, tp_mat, pro[:,3])
    temp = np.interp(time_mat, tp_mat, pro[:,4])
    
    den[setnan]=np.nan
    temp[setnan]=np.nan
    vx[setnan]=np.nan
    vy[setnan]=np.nan
    vz[setnan]=np.nan
  
    vt=np.sqrt(vx**2+vy**2+vz**2)

    print('end nan or interpolate')

        
    print('position start')
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('psp_pred'))
    psptra=spice.Trajectory('SPP')
    psptra.generate_positions(time,'Sun',frame)
    psptra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(psptra.x,psptra.y,psptra.z)
    print('PSP pos')    
    print('position end')

        
    #make array
    psp=np.zeros(np.size(bx),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('np', float),('tp', float),('x', float),('y', float),\
                ('z', float),('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    psp = psp.view(np.recarray)  

    #fill with data
    psp.time=time
    psp.bx=bx
    psp.by=by
    psp.bz=bz 
    psp.bt=bt
    
    psp.x=psptra.x
    psp.y=psptra.y
    psp.z=psptra.z
    
    psp.r=r
    psp.lat=np.degrees(lat)
    psp.lon=np.degrees(lon)
    
    psp.vt=vt
    psp.vx=vx    
    psp.vy=vy  
    psp.vz=vz
    psp.np=den
    #https://en.wikipedia.org/wiki/Thermal_velocity convert from km/s to K
    from astropy.constants import m_p,k_B
    psp.tp=np.pi*m_p*((temp*1e3)**2)/(8*k_B) 
    
    #remove spikes
    psp.tp[np.where(psp.tp > 1e10)]=np.nan
    
    print('done get psp')
    print()
    
    return psp

    
    
    
      
    
    
    

    
def save_psp_data_non_merged(path, file):
    '''
    *** TO DO
    save PSP data as pickle file with 3 separate arrays for orbit, magnetic field and plasma data    
    '''
     
    print('start PSP')
     
    psp_sat = heliosat.PSP()
    t_start = datetime.datetime(2018, 10, 14,14,14, 30)
    #t_end = datetime.datetime(2018, 12, 12,23,59,30)
    t_end = datetime.datetime(2019, 4, 23,23,59,30)

    #t_end = datetime.datetime(2019, 5, 31,23,59,30)
    #t_end = datetime.datetime(2019, 5, 1,23,59,30)

    
    timeb, mag = psp_sat.get_data(t_start, t_end, "psp_fields_l2")
    timep, pro = psp_sat.get_data(t_start, t_end, "psp_spc_l3")
    print('download complete')
    
    #create an array with 1 minute resolution between t start and end for position
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    
    
    #make separate arrays for orbit plasma and mag
    psp_orbit=np.zeros(np.size(time),dtype=[('time',object),('x', float),('y', float),\
                ('z', float),('r', float),('lat', float),('lon', float)])   
    #convert to recarray
    psp_orbit = psp_orbit.view(np.recarray)  

    psp_mag=np.zeros(np.size(timeb),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float)])   
    psp_mag = psp_mag.view(np.recarray)  

    psp_plasma=np.zeros(np.size(timep),dtype=[('time',object),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('np', float),('tp', float)])
    psp_plasma = psp_plasma.view(np.recarray)  


    psp_orbit.time=time
    psp_mag.time=parse_time(timeb,format='unix').datetime 
    psp_plasma.time=parse_time(timep,format='unix').datetime 
    print('time convert done')


    print('position start')
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('psp_pred'))
    psptra=spice.Trajectory('SPP')
    psptra.generate_positions(time,'Sun',frame)
    psptra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(psptra.x,psptra.y,psptra.z)
    psp_orbit.x=psptra.x
    psp_orbit.y=psptra.y
    psp_orbit.z=psptra.z
    psp_orbit.r=r
    psp_orbit.lat=np.rad2deg(lat)
    psp_orbit.lon=np.rad2deg(lon)
    print('position end')
   
    #fields
    psp_mag.bx = mag[:,0] 
    psp_mag.by = mag[:,1] 
    psp_mag.bz = mag[:,2] 
    psp_mag.bt = np.sqrt(psp_mag.bx**2+psp_mag.by**2+psp_mag.bz**2)
     
    #sweap
    from astropy.constants import m_p,k_B
    psp_plasma.np = pro[:,0]
    #https://en.wikipedia.org/wiki/Thermal_velocity convert from km/s to K

    psp_plasma.tp = np.pi*m_p*((pro[:,4]*1e3)**2)/(8*k_B) 

    psp_plasma.vx = pro[:,1]
    psp_plasma.vy = pro[:,2]
    psp_plasma.vz = pro[:,3]
    psp_plasma.vt=np.sqrt(psp_plasma.vx**2+psp_plasma.vy**2+psp_plasma.vz**2)

    
    header='PSP magnetic field (FIELDS instrument) and plasma data (SWEAP), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/psp/  '+ \
    'Timerange: '+psp_orbit.time[0].strftime("%Y-%b-%d %H:%M")+' to '+psp_orbit.time[-1].strftime("%Y-%b-%d %H:%M")+\
    '. The data are put in 3 numpy recarrays, fields can be accessed by psp_plasma.timep (for plasma), psp.vt etc.; psp_mag.timeb (for magnetic field),psp.bt, etc.; psp_orbit.time (for position) psp.r, psp.lon, ... '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(psp_plasma.size)+' (plasma), '+str(psp_mag.size)+' (mag). '+\
    'Units are btxyz [nT, RTN], vtxyz [km/s, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_psp_data_non_merged (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([psp_orbit,psp_mag,psp_plasma,header], open(path+file, "wb"))

    print('done psp')
    print()









def save_wind_data(path,file,start_date,end_date,heeq):
    
    '''
    description of data sources used in heliosat:
    https://cdaweb.sci.gsfc.nasa.gov/misc/NotesW.html
    '''
  
    print('start wind update')
    wind_sat = heliosat.WIND()
    t_start = start_date
    t_end = end_date
    
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=2*n) for n in range(int ((t_end - t_start).days*30*24))]  
    time_mat=mdates.date2num(time) 
    
    #print(parse_time(time[0:10]).iso)
    
    tm, mag = wind_sat.get_data_raw(t_start, t_end, "wind_mfi_k0", extra_columns=["PGSE"])#, return_datetimes=True)
    #tm, mag = wind_sat.get_data_raw(t_start, t_end, "wind_mfi_h0") k0
    tp, pro = wind_sat.get_data_raw(t_start, t_end, "wind_swe_h1")#,return_datetimes=True)
    
    tm=parse_time(tm,format='unix').datetime 
    tp=parse_time(tp,format='unix').datetime 
   
    
    #print(parse_time(tm[0:10]).iso)
    #print(parse_time(tp[0:10]).iso)
    
    tm=parse_time(tm).plot_date
    tp=parse_time(tp).plot_date
    
    print('download complete')
     
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, tm, mag[:,0] )
    by = np.interp(time_mat, tm, mag[:,1] )
    bz = np.interp(time_mat, tm, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
        
    den = np.interp(time_mat, tp, pro[:,0])
    vt = np.interp(time_mat, tp, pro[:,1])
    tp = np.interp(time_mat, tp, pro[:,2])
        
    #interpolate the GSE position over full data range
    x_gse = np.interp(time_mat, tm, mag[:,3])*6378.1/149597870.7*astropy.units.AU #earth radii to km to AU
    y_gse = np.interp(time_mat, tm, mag[:,4])*6378.1/149597870.7*astropy.units.AU
    z_gse = np.interp(time_mat, tm, mag[:,5])*6378.1/149597870.7*astropy.units.AU
    
    print('position start')    
    frame='HEEQ'
    planet_kernel=spicedata.get_kernel('planet_trajectories')
    earth=spice.Trajectory('399')  #399 for Earth
    earth.generate_positions(time,'Sun',frame)
    earth.change_units(astropy.units.AU)       #from km to AU    
    #add gse position to Earth position
    x=earth.x-x_gse  #earth radii to km
    y=earth.y-y_gse
    z=earth.z+z_gse
    [r, lat, lon]=cart2sphere(x,y,z)    
    
    
    #wind_pos=heliosat.WIND().trajectory(time, frame="HEEQ")
    #x=wind._pos[:,0]
    #y=wind_pos[:,1]
    #z=wind_pos[:,2]
    #[r, lat, lon]=hd.cart2sphere(wind_pos[:,0],wind_pos[:,1],wind_pos[:,2])
    #lon=np.rad2deg(lon) #convert to degree
    #lat=np.rad2deg(lat)


    
    
    print('position end ')
    
    #make array
    win=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('np', float),('vt', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    win = win.view(np.recarray)  

    #fill with data
    win.time=time
    win.bx=bx
    win.by=by
    win.bz=bz 
    win.bt=bt

    win.x=x
    win.y=y
    win.z=z
    
    win.r=r
    win.lat=np.degrees(lat)
    win.lon=np.degrees(lon)
    
    win.np=den
    win.vt=vt
    #https://en.wikipedia.org/wiki/Thermal_velocity convert from km/s to K
    from astropy.constants import m_p,k_B
    win.tp=np.pi*m_p*((tp*1e3)**2)/(8*k_B)
    
    
    ############ spike removal
            
    #plasma    
    win.np[np.where(win.np> 500)]=1000000
    #get rid of all single spikes with scipy signal find peaks
    peaks, properties = scipy.signal.find_peaks(win.np, height=500,width=(1, 250))
    #go through all of them and set to nan according to widths
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        win.np[peaks[i]-width-2:peaks[i]+width+2]=np.nan

    win.tp[np.where(win.tp> 1e8)]=1e11
    #get rid of all single spikes with scipy signal find peaks
    peaks, properties = scipy.signal.find_peaks(win.tp, height=1e8,width=(1, 250))
    #go through all of them and set to nan according to widths
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        win.tp[peaks[i]-width-2:peaks[i]+width+2]=np.nan

    win.vt[np.where(win.vt> 3000)]=1e11
    #get rid of all single spikes with scipy signal find peaks
    peaks, properties = scipy.signal.find_peaks(win.vt, height=1e8,width=(1, 250))
    #go through all of them and set to nan according to widths
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        win.vt[peaks[i]-width-2:peaks[i]+width+2]=np.nan

        
        
        
    #magnetic field    
    peaks, properties = scipy.signal.find_peaks(win.bt, prominence=30,width=(1, 10))
    #go through all of them and set to nan according to widths
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths'])[i])
        #remove data
        win.bt[peaks[i]-width-5:peaks[i]+width+5]=np.nan    

    peaks, properties = scipy.signal.find_peaks(abs(win.bx), prominence=30,width=(1, 10))
    for i in np.arange(len(peaks)):
        width=int(np.ceil(properties['widths'])[i])
        win.bx[peaks[i]-width-5:peaks[i]+width+5]=np.nan    

    peaks, properties = scipy.signal.find_peaks(abs(win.by), prominence=30,width=(1, 10))
    for i in np.arange(len(peaks)):
        width=int(np.ceil(properties['widths'])[i])
        win.by[peaks[i]-width-5:peaks[i]+width+5]=np.nan    

    peaks, properties = scipy.signal.find_peaks(abs(win.bz), prominence=30,width=(1, 10))
    for i in np.arange(len(peaks)):
        width=int(np.ceil(properties['widths'])[i])
        win.bz[peaks[i]-width-5:peaks[i]+width+5]=np.nan    



    #manual spike removal for magnetic field
    if t_start < datetime.datetime(2018, 7, 19, 16, 25):    
        if t_end > datetime.datetime(2018, 7, 19, 16, 25):         

            remove_start=datetime.datetime(2018, 7, 19, 16, 25)
            remove_end=datetime.datetime(2018, 7, 19, 17, 35)
            remove_start_ind=np.where(remove_start<win.time)[0][0]
            remove_end_ind=np.where(remove_end<win.time)[0][0] 

            win.bt[remove_start_ind:remove_end_ind]=np.nan
            win.bx[remove_start_ind:remove_end_ind]=np.nan
            win.by[remove_start_ind:remove_end_ind]=np.nan
            win.bz[remove_start_ind:remove_end_ind]=np.nan

    if t_start < datetime.datetime(2018, 8, 29, 19, 00):    
        if t_end > datetime.datetime(2018, 8, 29, 19, 00):         

            remove_start=datetime.datetime(2018, 8, 29, 19, 00)
            remove_end=datetime.datetime(2018,8, 30, 5, 00)
            remove_start_ind=np.where(remove_start<win.time)[0][0]
            remove_end_ind=np.where(remove_end<win.time)[0][0] 

            win.bt[remove_start_ind:remove_end_ind]=np.nan
            win.bx[remove_start_ind:remove_end_ind]=np.nan
            win.by[remove_start_ind:remove_end_ind]=np.nan
            win.bz[remove_start_ind:remove_end_ind]=np.nan

            
    if t_start < datetime.datetime(2019, 8, 8, 22, 45):    
        if t_end > datetime.datetime(2019, 8, 8, 22, 45):         

            remove_start=datetime.datetime(2019, 8, 8, 22, 45)
            remove_end=datetime.datetime(2019,   8, 9, 17, 00)
            remove_start_ind=np.where(remove_start<win.time)[0][0]
            remove_end_ind=np.where(remove_end<win.time)[0][0] 

            win.bt[remove_start_ind:remove_end_ind]=np.nan
            win.bx[remove_start_ind:remove_end_ind]=np.nan
            win.by[remove_start_ind:remove_end_ind]=np.nan
            win.bz[remove_start_ind:remove_end_ind]=np.nan

    if t_start < datetime.datetime(2019, 8, 21, 22, 45):    
        if t_end > datetime.datetime(2019, 8, 21, 22, 45):         

            remove_start=datetime.datetime(2019, 8, 20, 18, 0)
            remove_end=datetime.datetime(2019,   8, 21, 12, 0)
            remove_start_ind=np.where(remove_start<win.time)[0][0]
            remove_end_ind=np.where(remove_end<win.time)[0][0] 

            win.bt[remove_start_ind:remove_end_ind]=np.nan
            win.bx[remove_start_ind:remove_end_ind]=np.nan
            win.by[remove_start_ind:remove_end_ind]=np.nan
            win.bz[remove_start_ind:remove_end_ind]=np.nan            

    if t_start < datetime.datetime(2019, 8, 21, 22, 45):    
        if t_end > datetime.datetime(2019, 8, 21, 22, 45):         

            remove_start=datetime.datetime(2019, 8, 22, 1, 0)
            remove_end=datetime.datetime(2019,   8, 22, 9, 0)
            remove_start_ind=np.where(remove_start<win.time)[0][0]
            remove_end_ind=np.where(remove_end<win.time)[0][0] 

            win.bt[remove_start_ind:remove_end_ind]=np.nan
            win.bx[remove_start_ind:remove_end_ind]=np.nan
            win.by[remove_start_ind:remove_end_ind]=np.nan
            win.bz[remove_start_ind:remove_end_ind]=np.nan            


    
     
    coord='GSE'
    #convert magnetic field to SCEQ
    if heeq==True:
        win=convert_GSE_to_HEEQ(win)
        coord='HEEQ'
   
    
        
    ###################### 

    
    header='Wind magnetic field (MAG instrument) and plasma data (SWE), ' + \
    'obtained from https://spdf.gsfc.nasa.gov/pub/data/wind/  '+ \
    'Timerange: '+win.time[0].strftime("%Y-%b-%d %H:%M")+' to '+win.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(win.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by win.time, win.bx, win.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(win.size)+'. '+\
    'Units are btxyz [nT, '+coord+'], vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_wind_data (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    pickle.dump([win,header], open(path+file, "wb"))
    
    
    print('wind update done')
    print()
    



    
    
    
    
    
    


def load_stereoa_science_1min():

    varnames = ['Epoch', 'Vp', 'Vr_Over_V_RTN', 'Np', 'Tp', 'BFIELDRTN']
    alldata = {k: [] for k in varnames}
    if not os.path.exists(heliosat_data_path+'sta_magplasma_outside_heliosat'):
        os.mkdir(heliosat_data_path+'sta_magplasma_outside_heliosat')
    for year in ['2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014','2015','2016','2017','2018','2019']:
        print('get STEREO-A yearly 1min data file for ',year)
        cdf_write = heliosat_data_path+'sta_magplasma_outside_heliosat/STA_L2_MAGPLASMA_1m_{}_V01.cdf'.format(year)
        if not os.path.exists(cdf_write):
            cdf_url = ("https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/ahead/magplasma/STA_L2_MAGPLASMA_1m_{}_V01.cdf".format(year))
            cdf_file = requests.get(cdf_url)
            open(cdf_write, 'wb').write(cdf_file.content)
        cdf = cdflib.CDF(cdf_write)
        #cdf.cdf_info() shows all variable names and attributes
        for var in varnames:
            data = cdf[var][...]
            data[np.where(data < cdf.varattsget(var)['VALIDMIN'][0])] = np.NaN
            data[np.where(data > cdf.varattsget(var)['VALIDMAX'][0])] = np.NaN
            alldata[var].append(data)
    arrays = {}
    for var in varnames:
        arrays[var] = np.concatenate(alldata[var])
        
    return arrays



def save_all_stereoa_science_data(path,file,sceq):
    ''' 
    saves all STEREO-Ahead science data btxyz 
    vt np tp x y z r lat lon 1 min resolution as pickle
    sceq=True -> convert RTN to SCEQ coordinates for magnetic field components
    
    filesta_all='stereoa_2007_2019_sceq.p'
    hd.save_all_stereoa_science_data(data_path, filesta_all,sceq=True)

    filesta_all='stereoa_2007_2019_rtn.p'
    hd.save_all_stereoa_science_data(data_path, filesta_all,sceq=False)    
    
    
    [sta_t,hsta_t]=pickle.load(open(data_path+filesta_all, "rb" ) )
    '''
    
    #load all data with function
    s1=load_stereoa_science_1min()    
    print('download complete')

    #make array
    sta=np.zeros(len(s1['Epoch']),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
    
    #convert to recarray
    sta = sta.view(np.recarray)  
    
    #parse Epoch time to datetime objects
    sta.time=parse_time(s1['Epoch'],format='cdf_epoch').datetime  
    print('time conversion complete')

    sta.bx=s1['BFIELDRTN'][:,0]
    sta.by=s1['BFIELDRTN'][:,1]
    sta.bz=s1['BFIELDRTN'][:,2]
    sta.bt=np.sqrt(sta.bx**2+sta.by**2+sta.bz**2)
    
    sta.vt=s1['Vp']    
    sta.np=s1['Np']
    sta.tp=s1['Tp']    
    
    print('parameters into array complete')
    
    print('position start')
    frame='HEEQ'
    kernels = spicedata.get_kernel('stereo_a')
    kernels += spicedata.get_kernel('stereo_a_pred')
    spice.furnish(kernels)
    statra=spice.Trajectory('-234') #STEREO-A SPICE NAIF code
    statra.generate_positions(sta.time,'Sun',frame)
    statra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(statra.x,statra.y,statra.z)
    
    sta.x=statra.x
    sta.y=statra.y
    sta.z=statra.z
    
    sta.r=r
    sta.lat=np.degrees(lat)
    sta.lon=np.degrees(lon)

    print('position end ')
    
    
    
    #remove spike in magnetic field in 2015
    spike_ind=np.where(sta.bt >300)[0]
    if len(spike_ind) > 0:
        sta.bt[spike_ind[0]-77:spike_ind[-1]+5]=np.nan
        sta.bx[spike_ind[0]-77:spike_ind[-1]+5]=np.nan
        sta.by[spike_ind[0]-77:spike_ind[-1]+5]=np.nan
        sta.bz[spike_ind[0]-77:spike_ind[-1]+5]=np.nan



    
    coord='RTN'
    #convert magnetic field to SCEQ
    if sceq==True:
        print('convert RTN to SCEQ ')
        coord='SCEQ'
        sta=convert_RTN_to_SCEQ(sta,'STEREO-A')
    
    header='STEREO-A magnetic field (IMPACT instrument) and plasma data (PLASTIC, science), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/ahead/magplasma   '+ \
    'Timerange: '+sta.time[0].strftime("%Y-%b-%d %H:%M")+' to '+sta.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', with an average time resolution of '+str(np.mean(np.diff(sta.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by sta.time, sta.bx, sta.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(sta.size)+'. '+\
    'Units are btxyz [nT, '+coord+', vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats '+\
    'and https://github.com/heliopython/heliopy. '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, R. L. Bailey and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    print('save pickle file')
    pickle.dump([sta,header], open(path+file, "wb"))
    
    print('done sta')
    print()

   








def load_stereob_science_1min():

    varnames = ['Epoch', 'Vp', 'Vr_Over_V_RTN', 'Np', 'Tp', 'BFIELDRTN']
    alldata = {k: [] for k in varnames}
    if not os.path.exists(heliosat_data_path+'stb_magplasma_outside_heliosat'):
        os.mkdir(heliosat_data_path+'stb_magplasma_outside_heliosat')
    for year in ['2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014']:
    #for year in ['2007']:
        print('get STEREO-B yearly 1min data file for ',year)
        cdf_write = heliosat_data_path+'stb_magplasma_outside_heliosat/STB_L2_MAGPLASMA_1m_{}_V01.cdf'.format(year)
        if not os.path.exists(cdf_write):
            cdf_url = ("https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/behind/magplasma/STB_L2_MAGPLASMA_1m_{}_V01.cdf".format(year))
            cdf_file = requests.get(cdf_url)
            open(cdf_write, 'wb').write(cdf_file.content)
        cdf = cdflib.CDF(cdf_write)
        #cdf.cdf_info() shows all variable names and attributes
        for var in varnames:
            data = cdf[var][...]
            #fillval = cdf[var].attrs['FILLVAL']
            #fillval=cdf.varattsget(var)['FILLVAL'][0]
            data[np.where(data < cdf.varattsget(var)['VALIDMIN'][0])] = np.NaN
            data[np.where(data > cdf.varattsget(var)['VALIDMAX'][0])] = np.NaN
            alldata[var].append(data)
    arrays = {}
    for var in varnames:
        arrays[var] = np.concatenate(alldata[var])
        
    return arrays


def save_all_stereob_science_data(path,file,sceq):
    ''' 
    saves all STEREO-Behind science data btxyz 
    vt np tp x y z r lat lon 1 min resolution as pickle
    sceq=True -> convert RTN to SCEQ coordinates for magnetic field components
    
    use as:
    filestb_all='stereob_2007_2014_sceq.p'
    hd.save_all_stereob_science_data(data_path, filestb_all,sceq=True)
    
    filestb_all='stereob_2007_2014.p'
    hd.save_all_stereob_science_data(data_path, filestb_all,sceq=False)  
    
    [stb_t,hstb_t]=pickle.load(open(data_path+filestb_all, "rb" ) )    
    '''
    
    #load all data with function
    s1=load_stereob_science_1min()    
    print('download complete')

    #make array
    stb=np.zeros(len(s1['Epoch']),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
    
    #convert to recarray
    stb = stb.view(np.recarray)  
   
    
    #parse Epoch time to datetime objects
    stb.time=parse_time(s1['Epoch'],format='cdf_epoch').datetime  
    
    print('time conversion complete')

    stb.bx=s1['BFIELDRTN'][:,0]
    stb.by=s1['BFIELDRTN'][:,1]
    stb.bz=s1['BFIELDRTN'][:,2]
    stb.bt=np.sqrt(stb.bx**2+stb.by**2+stb.bz**2)

    
    stb.vt=s1['Vp']    
    stb.np=s1['Np']
    stb.tp=s1['Tp']    
    
    print('parameters into array complete')
    
    print('position start')
    frame='HEEQ'
    kernels = spicedata.get_kernel('stereo_b')
    spice.furnish(kernels)
    stbtra=spice.Trajectory('-235') #STEREO-A SPICE NAIF code
    stbtra.generate_positions(stb.time,'Sun',frame)
    stbtra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(stbtra.x,stbtra.y,stbtra.z)
    
    stb.x=stbtra.x
    stb.y=stbtra.y
    stb.z=stbtra.z
    
    stb.r=r
    stb.lat=np.rad2deg(lat)
    stb.lon=np.rad2deg(lon)

    print('position end ')

    
    coord='RTN'
    #convert magnetic field to SCEQ
    if sceq==True:
        print('convert RTN to SCEQ ')
        coord='SCEQ'
        stb=convert_RTN_to_SCEQ(stb,'STEREO-B')
    
    
    header='STEREO-B magnetic field (IMPACT instrument) and plasma data (PLASTIC, science), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/behind/magplasma   '+ \
    'Timerange: '+stb.time[0].strftime("%Y-%b-%d %H:%M")+' to '+stb.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', with a an average time resolution of '+str(np.mean(np.diff(stb.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by stb.time, stb.bx, stb.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(stb.size)+'. '+\
    'Units are btxyz [nT, '+coord+'], vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats '+\
    'and https://github.com/heliopython/heliopy. '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, R. L. Bailey and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
  
    print('save pickle file')
    pickle.dump([stb,header], open(path+file, "wb"))
    
    print('done stb')
    print()


 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


def save_noaa_rtsw_data_predstorm(data_path,noaa_path,filenoaa):


    print(' ')
    print('convert NOAA real time solar wind from predstorm h5 file to pickle file')

    hf = h5py.File(noaa_path, 'r')

    #make array
    noaa=np.zeros(len(np.array(hf.get('time'))),dtype=[('time',object),('bx', float),('by', float),\
                    ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                    ('x', float),('y', float),('z', float),\
                    ('r', float),('lat', float),('lon', float)])   

    #convert to recarray
    noaa = noaa.view(np.recarray)  
    noaa.time =  mdates.num2date( np.array(hf.get('time')))
    noaa.bt=np.array(hf.get('bt'))
    noaa.bx=np.array(hf.get('bx_gsm'))
    noaa.by=np.array(hf.get('by_gsm'))
    noaa.bz=np.array(hf.get('bz_gsm'))
    noaa.vt=np.array(hf.get('speed'))
    noaa.np=np.array(hf.get('density'))
    noaa.tp=np.array(hf.get('temperature'))

    header='Real time solar wind magnetic field and plasma data from NOAA, ' + \
            'obtained daily from https://services.swpc.noaa.gov/products/solar-wind/  '+ \
            'Timerange: '+noaa.time[0].strftime("%Y-%b-%d %H:%M")+' to '+noaa.time[-1].strftime("%Y-%b-%d %H:%M")+\
            ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(noaa.time)).seconds)+' seconds. '+\
            'The data are available in a numpy recarray, fields can be accessed by nf.time, nf.bx, nf.vt etc. '+\
            'Total number of data points: '+str(noaa.size)+'. '+\
            'Units are btxyz [nT, RTN], vt  [km s^-1], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
            'Made with https://github.com/cmoestl/heliocats save_noaa_rtsw_data_predstorm  '+\
            'By C. Moestl (twitter @chrisoutofspace) and R. L. Bailey. File creation date: '+\
            datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'

    pickle.dump([noaa,header], open(data_path+filenoaa, "wb"))

    print('NOAA from predstorm done')        
    
    
    
    
    
    
    
    
    
    
    



def save_stereob_beacon_data(path,file,start_time,end_time):

    print('start STB')
    stb_sat = heliosat.STB()
    t_start = start_time
    t_end = end_time
   
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
  
    tm, mag = stb_sat.get_data_raw(t_start, t_end, "mag_beacon")
    tp, pro = stb_sat.get_data_raw(t_start, t_end, "proton_beacon")

    print('download complete')
   
    tm=parse_time(tm,format='unix').datetime 
    tp=parse_time(tp,format='unix').datetime 

    #convert to matplotlib time for linear interpolation
    tm_mat=mdates.date2num(tm) 
    tp_mat=mdates.date2num(tp) 
    
    print('time convert done')
    
    print('position start')
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('stereo_b'))
    stbtra=spice.Trajectory('-235') #STEREO-A SPICE NAIF code
    stbtra.generate_positions(time,'Sun',frame)
    stbtra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(stbtra.x,stbtra.y,stbtra.z)
    print('position end ')
    
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
      
    den = np.interp(time_mat, tp_mat, pro[:,0])
    vt = np.interp(time_mat, tp_mat, pro[:,1])
    tp = np.interp(time_mat, tp_mat, pro[:,2])
    
    #make array
    stb=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    stb = stb.view(np.recarray)  

    #fill with data
    stb.time=time
    stb.bx=bx
    stb.by=by
    stb.bz=bz 
    stb.bt=bt

    stb.x=stbtra.x
    stb.y=stbtra.y
    stb.z=stbtra.z
    
    stb.r=r
    stb.lat=np.rad2deg(lat)
    stb.lon=np.rad2deg(lon)
    
    stb.np=den
    stb.tp=tp    
    stb.vt=vt  
        
    #remove spikes from plasma data
    #median filter
    stb.vt=scipy.signal.medfilt(stb.vt,9)
    #set nans to a high number
    stb.vt[np.where(np.isfinite(stb.vt) == False)]=1e5
    #get rid of all single spikes with scipy signal find peaks (cannot use nan)
    peaks,properties = scipy.signal.find_peaks(stb.vt, prominence=200,width=(1,200))
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        stb.vt[peaks[i]-width-2:peaks[i]+width+2]=np.nan
    #set nan again
    stb.vt[np.where(stb.vt == 1e5)]=np.nan     
    stb.tp[np.where(np.isfinite(stb.vt) == False)]=np.nan
    stb.np[np.where(np.isfinite(stb.vt) == False)]=np.nan   
    
    #remove spikes from magnetic field data
    #median filter
    #set nans to a high number
    stb.bt[np.where(np.isfinite(stb.bt) == False)]=1e5
    #get rid of all single spikes with scipy signal find peaks (cannot use nan)
    peaks,properties = scipy.signal.find_peaks(stb.bt, height=40,width=(1,20))
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths'])[i])
        #remove data
        stb.bt[peaks[i]-width-2:peaks[i]+width+2]=np.nan
        stb.bx[peaks[i]-width-2:peaks[i]+width+2]=np.nan
        stb.by[peaks[i]-width-2:peaks[i]+width+2]=np.nan
        stb.bz[peaks[i]-width-2:peaks[i]+width+2]=np.nan

    #set nan again
    stb.bt[np.where(stb.bt == 1e5)]=np.nan     

    #manual spike removal for speed
    remove_start=datetime.datetime(2007, 7, 18,22, 00)
    remove_end=datetime.datetime(2007, 7, 19, 16, 00)
    remove_start_ind=np.where(remove_start==stb.time)[0][0]
    remove_end_ind=np.where(remove_end==stb.time)[0][0] 

    stb.vt[remove_start_ind:remove_end_ind]=np.nan
    
    
    header='BEACON STEREO-B magnetic field (IMPACT instrument) and plasma data (PLASTIC), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/beacon/behind/  '+ \
    'Timerange: '+stb.time[0].strftime("%Y-%b-%d %H:%M")+' to '+stb.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(stb.time)).seconds)+' seconds. '+\
    'A median filter has been applied (plasma data only) and then spikes were removed with scipy.signal.find_peaks (plasma and field). '+\
    'The data are available in a numpy recarray, fields can be accessed by stb.time, stb.bx, stb.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(stb.size)+'. '+\
    'Units are btxyz [nT, RTN], vt  [km s^-1], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats heliocats.data.save_stereob_beacon_data (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([stb,header], open(path+file, "wb"))
    
    print('done stb')
    print()




    
    
    
    

def save_stereoa_science_data_old(path,file,t_start, t_end,sceq):
    
    '''** TO DO 
    '''

    print('start STA')
    sta_sat = heliosat.STA()
     
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  
    time_mat=mdates.date2num(time) 
    
    #tp, pro = sta_sat.get_data_raw(t_start, t_end, "sta_plastic_l2")
    #tm, mag = sta_sat.get_data_raw(t_start, t_end, "sta_impact_beacon")
    tm, mag = sta_sat.get_data_raw(t_start, t_end, "sta_impact_l1")

    print('download complete')
   
    tm=parse_time(tm,format='unix').datetime 
    tp=parse_time(tp,format='unix').datetime 

    #convert to matplotlib time for linear interpolation
    tm_mat=mdates.date2num(tm) 
    tp_mat=mdates.date2num(tp) 
    
    print('time convert done')
    
    
    print('position start')
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('stereo_a_pred'))
    statra=spice.Trajectory('-234') #STEREO-A SPICE NAIF code
    statra.generate_positions(time,'Sun',frame)
    statra.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(statra.x,statra.y,statra.z)
    print('position end ')
    
    
    #linear interpolation to time_mat times    
    bx = np.interp(time_mat, tm_mat, mag[:,0] )
    by = np.interp(time_mat, tm_mat, mag[:,1] )
    bz = np.interp(time_mat, tm_mat, mag[:,2] )
    bt = np.sqrt(bx**2+by**2+bz**2)
      
      
    den = np.interp(time_mat, tp_mat, pro[:,0])
    vt = np.interp(time_mat, tp_mat, pro[:,1])
    tp = np.interp(time_mat, tp_mat, pro[:,2])
    
    #make array
    sta=np.zeros(np.size(bx),dtype=[('time',object),('bx', float),('by', float),\
                ('bz', float),('bt', float),('vt', float),('np', float),('tp', float),\
                ('x', float),('y', float),('z', float),\
                ('r', float),('lat', float),('lon', float)])   
       
    #convert to recarray
    sta = sta.view(np.recarray)  

    #fill with data
    sta.time=time
    sta.bx=bx
    sta.by=by
    sta.bz=bz 
    sta.bt=bt

    sta.x=statra.x
    sta.y=statra.y
    sta.z=statra.z
    
    sta.r=r
    sta.lat=np.rad2deg(lat)
    sta.lon=np.rad2deg(lon)
    
    sta.vt=vt    
    sta.np=den
    sta.tp=tp    
    
     
     
    #remove spikes from plasma data
    #median filter
    sta.vt=scipy.signal.medfilt(sta.vt,9)
    #set nans to a high number
    sta.vt[np.where(np.isfinite(sta.vt) == False)]=1e5
    #get rid of all single spikes with scipy signal find peaks (cannot use nan)
    peaks,properties = scipy.signal.find_peaks(sta.vt, prominence=200,width=(1,200))
    for i in np.arange(len(peaks)):
        #get width of current peak
        width=int(np.ceil(properties['widths']/2)[i])
        #remove data
        sta.vt[peaks[i]-width-2:peaks[i]+width+2]=np.nan
    #set nan again
    sta.vt[np.where(sta.vt == 1e5)]=np.nan     
    sta.tp[np.where(np.isfinite(sta.vt) == False)]=np.nan
    sta.np[np.where(np.isfinite(sta.vt) == False)]=np.nan   
    
    #manual spike removal for magnetic field
    remove_start=datetime.datetime(2018, 9, 23, 11, 00)
    remove_end=datetime.datetime(2018, 9, 25, 00, 00)
    remove_start_ind=np.where(remove_start==sta.time)[0][0]
    remove_end_ind=np.where(remove_end==sta.time)[0][0] 

    sta.bt[remove_start_ind:remove_end_ind]=np.nan
    sta.bx[remove_start_ind:remove_end_ind]=np.nan
    sta.by[remove_start_ind:remove_end_ind]=np.nan
    sta.bz[remove_start_ind:remove_end_ind]=np.nan

    
    
    #convert magnetic field to SCEQ
    if sceq==True:
        sta=convert_RTN_to_SCEQ(sta,'STEREO-A')
    
    
    
    
    
    header='STEREO-A magnetic field (IMPACT instrument, beacon) and plasma data (PLASTIC, science), ' + \
    'obtained from https://stereo-ssc.nascom.nasa.gov/data/ins_data/   '+ \
    'Timerange: '+sta.time[0].strftime("%Y-%b-%d %H:%M")+' to '+sta.time[-1].strftime("%Y-%b-%d %H:%M")+\
    ', linearly interpolated to a time resolution of '+str(np.mean(np.diff(sta.time)).seconds)+' seconds. '+\
    'The data are available in a numpy recarray, fields can be accessed by sta.time, sta.bx, sta.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(sta.size)+'. '+\
    'Units are btxyz [nT, RTN], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. '+\
    'Made with https://github.com/cmoestl/heliocats (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
  
    pickle.dump([sta,header], open(path+file, "wb"))
    
    print('done sta')
    print()










def convert_MAVEN_mat_original(file_input,filename):

    print('load MAVEN from MAT')
    
    file=file_input
    mavraw = scipy.io.loadmat(file)
    
    #make array 
    mav=np.zeros(np.size(mavraw['BT']),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('tp', float),('np', float),('r', float),('lat', float),\
                ('lon', float),('x', float),('y', float),('z', float),\
                ('ro', float), ('lato', float), ('lono', float),\
                ('xo', float), ('yo', float), ('zo', float)])   
    #convert to recarray
    mav = mav.view(np.recarray)  
    #convert time from matlab to python
    t=mavraw['timeD'][:,0]    
    
    for p in np.arange(np.size(t)):
        mav.time[p]= datetime.datetime.fromordinal(t[p].astype(int) ) + \
        datetime.timedelta(days=t[p]%1) - datetime.timedelta(days = 366) 
        

    mav.bx=mavraw['Bx'][:,0]       
    mav.by=mavraw['By'][:,0] 
    mav.bz=mavraw['Bz'][:,0]      
    mav.bt=mavraw['BT'][:,0]      

    mav.vx=mavraw['Vx'][:,0]      
    mav.vy=mavraw['Vy'][:,0]      
    mav.vz=mavraw['Vz'][:,0]      
    mav.vt=mavraw['VT'][:,0] 

    mav.tp=mavraw['Tp'][:,0]*(1.602176634*1e-19)/(1.38064852*1e-23)        #from ev to K     
    mav.np=mavraw['np'][:,0]      
    
    #add position with respect to Mars center in km in MSO
    print('orbit position start')
    
    insertion=datetime.datetime(2014,9,22,2,24,0)
    #these are the indices of the times for the cruise phase     
    tc=np.where(mdates.date2num(mav.time) < mdates.date2num(insertion))       

    mars_radius=3389.5
    mav.xo=mavraw['Xsc'][:,0]*mars_radius
    mav.yo=mavraw['Ysc'][:,0]*mars_radius
    mav.zo=mavraw['Zsc'][:,0]*mars_radius  
    
    #set to nan for cruise phase
    mav.xo[tc]=np.nan
    mav.yo[tc]=np.nan
    mav.zo[tc]=np.nan
    
    [mav.ro,mav.lato,mav.lono]=cart2sphere(mav.xo,mav.yo,mav.zo)
    mav.lono=np.rad2deg(mav.lono)   
    mav.lato=np.rad2deg(mav.lato)

    print('HEEQ position start')
    frame='HEEQ'
    
    #add position in HEEQ for cruise phase and orbit
    
    #cruise phase
    #use heliopy to load own bsp spice file from MAVEN 
    #obtained through https://naif.jpl.nasa.gov/pub/naif/pds/pds4/maven/maven_spice/spice_kernels/spk/
    spice.furnish(data_path+'input/maven_cru_rec_131118_140923_v1.bsp') 
    cruise=spice.Trajectory('MAVEN') #or NAIF CODE -202 
    cruise.generate_positions(mav.time[tc],'Sun',frame)     
    cruise.change_units(astropy.units.AU)  
    mav.x[tc]=cruise.x
    mav.y[tc]=cruise.y
    mav.z[tc]=cruise.z
    [mav.r[tc], mav.lat[tc], mav.lon[tc]]=cart2sphere(mav.x[tc],mav.y[tc],mav.z[tc])

    #times in orbit
    to=np.where(mdates.date2num(mav.time) > mdates.date2num(insertion))       

    planet_kernel=spicedata.get_kernel('planet_trajectories')
    mars=spice.Trajectory('MARS BARYCENTER')
    mars.generate_positions(mav.time[to],'Sun',frame)
    mars.change_units(astropy.units.AU)  
    mav.x[to]=mars.x
    mav.y[to]=mars.y
    mav.z[to]=mars.z
    [mav.r[to], mav.lat[to], mav.lon[to]]=cart2sphere(mav.x[to],mav.y[to],mav.z[to])

    #convert to degree
    mav.lon=np.rad2deg(mav.lon)   
    mav.lat=np.rad2deg(mav.lat)
    print('position end ')
  
        
    print('save MAVEN as pickle')

    header='MAVEN merged magnetic field and plasma data, obtained from Toulouse and C. Simon Wedlund. '+\
    'Timerange: '+mav.time[0].strftime("%Y-%b-%d %H:%M")+' to '+mav.time[-1].strftime("%Y-%b-%d %H:%M")+'.'+\
    'Mean time resolution: '+str(np.mean(np.diff(mav.time)).seconds)+' seconds. '+\
    'The data are put in a numpy recarray, fields can be accessed by mav.time, mav.bx, mav.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(mav.size)+'. '+\
    'Units are btxyz [nT, MSO], vtxyz [km/s, MSO], np[cm^-3], tp [K], orbital position: '+ \
    'xo/yo/zo/ro/lono/lato [km, degree, MSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    'Made with https://github.com/cmoestl/heliocats heliocats.data.convert_MAVEN_mat_original (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    

    pickle.dump([mav,header], open(filename, "wb"))

   
    



def convert_MAVEN_mat_removed(file_input,filename):

    print('load MAVEN from MAT')

    file=file_input
    mavraw = scipy.io.loadmat(file)
    
    
    #load time data extra
    
    file_input=data_path+'input/Data-MAVEN-MAG_SolarWind_102014-012021_time.mat'
    mavraw_time = scipy.io.loadmat(file_input)
    #mavraw_time['time']

   
    
    #make array 
    mav=np.zeros(np.size(mavraw['BT']),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('tp', float),('np', float),('r', float),('lat', float),\
                ('lon', float),('x', float),('y', float),('z', float),\
                ('ro', float), ('lato', float), ('lono', float),\
                ('xo', float), ('yo', float), ('zo', float)])   
    #convert to recarray
    mav = mav.view(np.recarray)  
    #convert time from matlab to python
    #t=mavraw['timeD'][:,0]
    t=mavraw_time['time']

    for p in np.arange(np.size(t)):
        mav.time[p]= datetime.datetime.fromordinal(t[p][0].astype(int) ) + \
        datetime.timedelta(days=t[p][0]%1) - datetime.timedelta(days = 366) 
        

    mav.bx=mavraw['Bx'][:,0]       
    mav.by=mavraw['By'][:,0] 
    mav.bz=mavraw['Bz'][:,0]      
    mav.bt=mavraw['BT'][:,0]      

    
    
    ### TO DO *** add plasma data
    #file_input=data_path+'input/MAVEN_2014to2018_cyril.mat'
    #mavraw2 = scipy.io.loadmat(file_input)
   
    
    #mav.vx[0:len(mavraw2)]=mavraw2['Vx'][:,0]      
    #mav.vy[0:len(mavraw2)]=mavraw2['Vy'][:,0]      
    #mav.vz[0:len(mavraw2)]=mavraw2['Vz'][:,0]      
    #mav.vt[0:len(mavraw2)]=mavraw2['VT'][:,0] 

    #mav.tp[0:len(mavraw2)]=mavraw2['Tp'][:,0]*(1.602176634*1e-19)/(1.38064852*1e-23)        #from ev to K     
    #mav.np[0:len(mavraw2)]=mavraw2['np'][:,0]      
  
    
    #add position with respect to Mars center in km in MSO
 
    print('orbit position start')
    
    insertion=datetime.datetime(2014,9,22,2,24,0)
    #these are the indices of the times for the cruise phase     
    #tc=np.where(mdates.date2num(mav.time) < mdates.date2num(insertion))       

    mars_radius=3389.5
    mav.xo=mavraw['Xsc'][:,0]*mars_radius
    mav.yo=mavraw['Ysc'][:,0]*mars_radius
    mav.zo=mavraw['Zsc'][:,0]*mars_radius  
    
    #set to nan for cruise phase
    #mav.xo[tc]=np.nan
    #mav.yo[tc]=np.nan
    #mav.zo[tc]=np.nan
    
    [mav.ro,mav.lato,mav.lono]=cart2sphere(mav.xo,mav.yo,mav.zo)
    mav.lono=np.rad2deg(mav.lono)   
    mav.lato=np.rad2deg(mav.lato)


    print('HEEQ position start')
    frame='HEEQ'
    
    #add position in HEEQ for cruise phase and orbit
    
    #cruise phase for plasma file only
    #use heliopy to load own bsp spice file from MAVEN 
    #obtained through https://naif.jpl.nasa.gov/pub/naif/pds/pds4/maven/maven_spice/spice_kernels/spk/
    #spice.furnish(data_path+'input/maven_cru_rec_131118_140923_v1.bsp') 
    #cruise=spice.Trajectory('MAVEN') #or NAIF CODE -202 
    #cruise.generate_positions(mav.time[tc],'Sun',frame)     
    #cruise.change_units(astropy.units.AU)  
    #mav.x[tc]=cruise.x
    #mav.y[tc]=cruise.y
    #mav.z[tc]=cruise.z
    #[mav.r[tc], mav.lat[tc], mav.lon[tc]]=cart2sphere(mav.x[tc],mav.y[tc],mav.z[tc])

    #times in orbit
    to=np.where(mdates.date2num(mav.time) > mdates.date2num(insertion))       

    planet_kernel=spicedata.get_kernel('planet_trajectories')
    mars=spice.Trajectory('MARS BARYCENTER')
    mars.generate_positions(mav.time[to],'Sun',frame)
    mars.change_units(astropy.units.AU)  
    mav.x[to]=mars.x
    mav.y[to]=mars.y
    mav.z[to]=mars.z
    [mav.r[to], mav.lat[to], mav.lon[to]]=cart2sphere(mav.x[to],mav.y[to],mav.z[to])

    #convert to degree
    mav.lon=np.rad2deg(mav.lon)   
    mav.lat=np.rad2deg(mav.lat)
    print('position end ')
    
    
    header='MAVEN merged magnetic field and plasma data, obtained from Toulouse. ' + \
    'The magnetosphere is removed with the Gruesbeck et al. 3D model (by C. Simon Wedlund). '+ \
    'Timerange: '+mav.time[0].strftime("%Y-%b-%d %H:%M")+' to '+mav.time[-1].strftime("%Y-%b-%d %H:%M")+'. '+\
    'Mean time resolution: '+str(np.mean(np.diff(mav.time)).seconds)+' seconds. '+\
    'The data are put in a numpy recarray, fields can be accessed by mav.time, mav.bx, mav.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(mav.size)+'. '+\
    'Units are btxyz [nT, MSO], vtxyz [km/s, MSO], np[cm^-3], tp [K], orbital position: '+ \
    'xo/yo/zo/ro/lono/lato [km, degree, MSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    'Made with https://github.com/cmoestl/heliocats heliocats.data.convert_MAVEN_mat_removed (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
     
    pickle.dump([mav,header], open(filename, "wb"))







def MAVEN_smooth_orbit(filemav,filename):


    [mav,hmav]=pickle.load(open(filemav, 'rb' ) )
    print('loaded ',filemav)
    
    
    ############# smooth over each orbit to extract solar wind signal
    #determine apogees
    #get rid of nans
    mav.ro[np.where(mav.ro)==np.nan]=-1e3 
    peaks,properties = scipy.signal.find_peaks(mav.ro,height=5000,width=(100))
    #bring nans back
    mav.ro[np.where(mav.ro)==-1e3]=np.nan
    
    print('Nr. of orbits in dataset: ',len(peaks))


    #make array 
    mavs=np.zeros(np.size(peaks),dtype=[('time',object),('bt', float),('bx', float),\
                ('by', float),('bz', float),('vt', float),('vx', float),('vy', float),\
                ('vz', float),('tp', float),('np', float),('r', float),('lat', float),\
                ('lon', float),('x', float),('y', float),('z', float),\
                ('ro', float), ('lato', float), ('lono', float),\
                ('xo', float), ('yo', float), ('zo', float)])   
     
    #convert to recarray
    mavs = mavs.view(np.recarray)  

                
    #2h on each side
    window=121
    for i in np.arange(len(peaks)):
          mavs.bt[i]=np.nanmedian(mav.bt[peaks[i]-window:peaks[i]+window])
          mavs.bx[i]=np.nanmedian(mav.bx[peaks[i]-window:peaks[i]+window])
          mavs.by[i]=np.nanmedian(mav.by[peaks[i]-window:peaks[i]+window])
          mavs.bz[i]=np.nanmedian(mav.bz[peaks[i]-window:peaks[i]+window])

          mavs.vt[i]=np.nanmedian(mav.vt[peaks[i]-window:peaks[i]+window])
          mavs.vx[i]=np.nanmedian(mav.vx[peaks[i]-window:peaks[i]+window])
          mavs.vy[i]=np.nanmedian(mav.vy[peaks[i]-window:peaks[i]+window])
          mavs.vz[i]=np.nanmedian(mav.vz[peaks[i]-window:peaks[i]+window])

          mavs.np[i]=np.nanmedian(mav.np[peaks[i]-window:peaks[i]+window])
          mavs.tp[i]=np.nanmedian(mav.tp[peaks[i]-window:peaks[i]+window])
         
          mavs.time[i]=mav.time[peaks[i]] 
          
          mavs.r[i]=mav.r[peaks[i]]
          mavs.lat[i]=mav.lat[peaks[i]]
          mavs.lon[i]=mav.lon[peaks[i]]

          mavs.x[i]=mav.x[peaks[i]]
          mavs.y[i]=mav.y[peaks[i]]
          mavs.z[i]=mav.z[peaks[i]]

          mavs.ro[i]=mav.ro[peaks[i]]
          mavs.lato[i]=mav.lato[peaks[i]]
          mavs.lono[i]=mav.lono[peaks[i]]

          mavs.xo[i]=mav.xo[peaks[i]]
          mavs.yo[i]=mav.yo[peaks[i]]
          mavs.zo[i]=mav.zo[peaks[i]]

    
    '''
    for testing:
    plt.figure(1)
    ax1 = plt.subplot(121)
    #ax1.plot_date(mav.time,mav.bt,'bo') 

    #ax1.plot_date(mav.time,mav.bt,'-r') 
    ax1.plot_date(mav.time[peaks],bt1,'-r') 
    ax1.plot_date(mav.time[peaks],vt1,'-k') 

    #ax1.plot_date(mav.time[peaks],bt2,'bo') 


    #ax1.plot_date(mav.time,g,'-b') 

    #ax1.set_xlim(timeset-days_window*10,timeset+days_window*10)

    ax2 = plt.subplot(122)
    ax2.plot_date(mav.time,mav.ro,'-k') 
    ax2.plot_date(mav.time[peaks],mav.ro[peaks],'bo') 
    #ax2.set_xlim(timeset-days_window,timeset+days_window)

    ax2.set_ylim(7000,10000)
    plt.show()


    plt.figure(2)
    plt.plot_date(mav.time,mav.vt,'-r') 

    plt.plot_date(mav.time[peaks],vt1,'ob') 
    '''
    #pickle.dump(mavs, open(data_path+filename, "wb"))

    
    
    header='MAVEN solar wind merged magnetic field and plasma data. ' + \
    'The magnetosphere was removed with the Gruesbeck et al. 3D model (by C. Simon Wedlund), '+\
    'and a +/-2h median filter around the apogee is used for 1 data point per orbit. '+ \
    'Timerange: '+mavs.time[0].strftime("%Y-%b-%d %H:%M")+' to '+mavs.time[-1].strftime("%Y-%b-%d %H:%M")+'. '+\
    'Mean time resolution: '+str(np.mean(np.diff(mavs.time)).seconds)+' seconds. '+\
    'The data are put in a numpy recarray, fields can be accessed by mav.time, mav.bx, mav.vt etc. '+\
    'Missing data has been set to "np.nan". Total number of data points: '+str(mavs.size)+'. '+\
    'Units are btxyz [nT, MSO], vtxyz [km/s, MSO], np[cm^-3], tp [K], orbital position: '+ \
    'xo/yo/zo/ro/lono/lato [km, degree, MSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    'Made with https://github.com/cmoestl/heliocats heliocats.data.convert_MAVEN_mat_removed_orbit (uses https://github.com/ajefweiss/HelioSat '+\
    'and https://github.com/heliopython/heliopy). '+\
    'By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, and D. Stansby. File creation date: '+\
    datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M")+' UTC'
    
    pickle.dump([mavs,header], open(filename, "wb"))

    


    
    
    
    



########################################## load HISTORIC DATA ############################


def save_helios_data(file):
    
    '''
    **TO DO
    '''

    print('start Helios')
    t_start = datetime.datetime(1975, 1, 1)
    t_end = datetime.datetime(1980, 12, 31)
     
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  


    #h1=heliosdata.corefit(1,t_start,t_end)
    #h2=heliosdata.corefit(2,t_start,t_end)
    h1=heliosdata.merged(1,t_start,t_end)  

    print('end Helios')




def save_cassini_data(file):

    '''
    **TO DO
    '''

    print('start Cassini')
    t_start = datetime.datetime(1999, 8, 16)
    t_end = datetime.datetime(2016, 12, 31)
     
    #create an array with 1 minute resolution between t start and end
    time = [ t_start + datetime.timedelta(minutes=1*n) for n in range(int ((t_end - t_start).days*60*24))]  


    coords='RTN'

    #Cassini Orbiter Magnetometer Calibrated MAG data in 1 minute averages available 
    #covering the period 1999-08-16 (DOY 228) to 2016-12-31 (DOY 366). 
    #The data are provided in RTN coordinates throughout the mission, with Earth, Jupiter, 
    #and Saturn centered coordinates for the respective flybys of those planets.
    cas=cassinidata.mag_hires(t_start,t_end, coords)




def save_ulysses_data(data_path):

   
    print('read Ulysses data from cdf and convert to pickle')   

    datacat_path='/nas/helio/data/DATACAT/'
    #load cdf
    ulycdf = cdflib.CDF(datacat_path+'ulysses_merged_1990_2009_CDAWEB.cdf') 
    #check variables
    #ulycdf.cdf_info()

    #time conversion to datetime      
    time=ulycdf.varget('Epoch')
    t=parse_time(time,format='cdf_epoch').datetime  
    
    
    #cutoff time and later data so that it starts with available position on Oct 6 1990
    t=t[6696:-1]

    
    print('Ulysses position start')
    #position starts on Oct 6 1990
    frame='HEEQ'
    spice.furnish(spicedata.get_kernel('ulysses'))
    upos=spice.Trajectory('-55')
    upos.generate_positions(t,'Sun',frame)
    upos.change_units(astropy.units.AU)  
    [r, lat, lon]=cart2sphere(upos.x,upos.y,upos.z)
    print('position end ')
    
    #make custom array
    uly=np.zeros(len(t),dtype=[('time',object),('bx', float),('by', float), \
     ('bz', float), ('bt', float),('vt', float),('np', float),('tp', float), \
      ('x', float),('y', float), ('z', float),('r', float),('lat', float), \
      ('lon', float)])   
    #convert to recarray
    uly = uly.view(np.recarray)  
    
    uly.time=t
    
    uly.bx=ulycdf.varget('BR')[6696:-1]
    uly.by=ulycdf.varget('BT')[6696:-1]
    uly.bz=ulycdf.varget('BN')[6696:-1]
    uly.bt=ulycdf.varget('ABS_B')[6696:-1]
    uly.vt=ulycdf.varget('plasmaFlowSpeed')[6696:-1]
    uly.np=ulycdf.varget('protonDensity')[6696:-1]
    uly.tp=ulycdf.varget('protonTempLarge')[6696:-1]
    
    
    uly.x=upos.x
    uly.y=upos.y
    uly.z=upos.z
    
    uly.r=r
    uly.lat=np.rad2deg(lat)
    uly.lon=np.rad2deg(lon)
        
    badmag=np.where(uly.bt < -10000)
    uly.bt[badmag]=np.nan  
    uly.bx[badmag]=np.nan  
    uly.by[badmag]=np.nan  
    uly.bz[badmag]=np.nan  
    
    badv=np.where(uly.vt < -100000)
    uly.vt[badv]=np.nan  
    
    badn=np.where(uly.np < -100000)
    uly.np[badn]=np.nan  
    
    badt=np.where(uly.tp < -100000)
    uly.tp[badt]=np.nan  
    
 
    header='Ulysses merged magnetic field and plasma data, obtained from CDAWEB. '+ \
    'Timerange: '+uly.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+uly.time[-1].strftime("%d-%b-%Y %H:%M:%S") +\
    '. Units are btxyz [nT, RTN], vt [km/s], np [cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ].\
Made with https://github.com/cmoestl/heliocats and https://github.com/heliopython/heliopy. By C. Moestl (twitter @chrisoutofspace) and D. Stansby.'
            

    file=data_path+'ulysses_1990_2009_rtn.p'
    pickle.dump([uly,header], open(file, "wb"))
    
    
    print('Ulysses done')




      






############################# HELCATS DATA into single file ###############################



#for reading .sav files
def getcat(filename):

    print('reading CAT '+filename)
    cat=scipy.io.readsav(filename, verbose='true')  
    print('done reading CAT')
    return cat  
  
def convert_sav_to_p():
    start=time.time()
    wind=getcat('/nas/helio/data/DATACAT/WIND_2007to2016_HEEQ.sav')
    end=time.time()
    print('load wind end. time in minutes:', (end-start)/60)
    #save time and components as pickle
    pickle.dump(wind.wind, open( "/nas/helio/data/insitu_python/WIND_2007to2016_HEEQ.p", "wb" ) )


def despike_helcats_speed_wind(vt,vx,vy,vz):


    #set all nan to 0 in the v gradient array
    v1=copy.deepcopy( abs( np.gradient(vt)   ))
    vnan_ind=np.where(np.isnan(v1)==True)[0]
    v1[vnan_ind]=0
 
    peaks, properties = scipy.signal.find_peaks(v1, prominence=40, width=(1, 10))
    vt[np.where(v1>100)]=np.nan
    vx[np.where(v1>100)]=np.nan
    vy[np.where(v1>100)]=np.nan
    vz[np.where(v1>100)]=np.nan


    for i in np.arange(len(peaks)):
         width=int(np.ceil(properties['widths'])[i])
         #print(width)   
         vt[peaks[i]-width:peaks[i]+width]=np.nan
         vx[peaks[i]-width:peaks[i]+width]=np.nan
         vy[peaks[i]-width:peaks[i]+width]=np.nan
         vz[peaks[i]-width:peaks[i]+width]=np.nan


    #print(properties['widths'])
    #plt.plot(win.vt[0:80000],'-b')
    #plt.plot(win.vt[0:80000],'-k',linewidth=1)
    #plt.plot(v,'-r',linewidth=1)
    #plt.plot(v1,'-b',linewidth=1)
    return vt,vx,vy,vz
 
def despike_helcats_density_wind(den):


    #set all nan to 0 in the v gradient array
    den1=copy.deepcopy( abs( np.gradient(den)   ))
    den1nan_ind=np.where(np.isnan(den1)==True)[0]
    den1[den1nan_ind]=0
 
    peaks, properties = scipy.signal.find_peaks(den1, prominence=10, width=(1, 10))
    den[np.where(den1>10)]=np.nan
    den1[np.where(den1>10)]=np.nan

    for i in np.arange(len(peaks)):
         width=int(np.ceil(properties['widths'])[i])
         #print(width)   
         den[peaks[i]-width:peaks[i]+width]=np.nan
            
            
    #print(properties['widths'])
    #plt.plot(win.vt[0:80000],'-b')
    #plt.plot(win.np[1200000:1500000]+50,'-g',linewidth=5)
    #plt.plot(den+1,'-k',linewidth=1)
    #plt.plot(den1,'-b',linewidth=1)

    return den
 

def despike_helcats_temperature_wind(den):

    den=den/1e6    
    #set all nan to 0 in the v gradient array
    den1=copy.deepcopy( abs( np.gradient(den)   ))
    den1nan_ind=np.where(np.isnan(den1)==True)[0]
    den1[den1nan_ind]=0

    peaks, properties = scipy.signal.find_peaks(den1, prominence=0.2, width=(1, 10))
    #den[np.where(den>100)[0]]=np.nan
    den[np.where(den1>0.2)]=np.nan
    den1[np.where(den1>0.2)]=np.nan

    for i in np.arange(len(peaks)):
         width=int(np.ceil(properties['widths'])[i])
         #print(width)   
         den[peaks[i]-width:peaks[i]+width]=np.nan

    return den*1e6



    
    
    
def save_helcats_datacat(data_path,removed):  
    ''' to save all of helcats DATACAT into a single file
    use: hd.save_helcats_datacat(data_path,removed=True)
    '''
      
    print('save all helcats DATACAT into single file')
    datacat_path='/nas/helio/data/DATACAT/'
    print('all data in', datacat_path)
      
      
      
      
      
    print( 'read Wind')
    winin= pickle.load( open(datacat_path+ "WIND_2007to2018_HEEQ.p", "rb" ) )
    winin_time=parse_time(winin.time,format='utime').datetime
    winin=winin.astype([('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('vt', 'float64'), ('vx', 'float64'), ('vy', 'float64'), \
                    ('vz', 'float64'), ('tp', 'float64'), ('np', 'float64'), \
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64')])  
    #make new array with xyz                
    win=np.zeros(np.size(winin),[('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('vt', 'float64'), ('vx', 'float64'), ('vy', 'float64'), \
                    ('vz', 'float64'), ('tp', 'float64'), ('np', 'float64'),\
                    ('x', 'float64'),('y', 'float64'), ('z', 'float64'),\
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64')])      
                    
    win = win.view(np.recarray)                                                       
    
    win.time=winin_time
    win.bx=winin.bx
    win.by=winin.by
    win.bz=winin.bz
    win.bt=winin.bt

    win.vt=winin.vt
    win.vx=winin.vx
    win.vy=winin.vy
    win.vz=winin.vz

    win.np=winin.np
    win.tp=winin.tp
 
    [win.vt,win.vx,win.vy,win.vz]=despike_helcats_speed_wind(win.vt,win.vx,win.vy,win.vz)
    win.np=despike_helcats_density_wind(win.np)
    win.tp=despike_helcats_temperature_wind(win.tp)
    






    win.r=winin.r/(astropy.constants.au.value/1e3)
    win.lat=winin.lat
    win.lon=winin.lon
     
        
    [win.x, win.y, win.z]=sphere2cart(win.r,np.abs(win.lat-np.radians(90)),win.lon)
    win.lon=np.rad2deg(win.lon)   
    win.lat=np.rad2deg(win.lat)

    del(winin)
    
    hwin='Wind merged magnetic field and plasma data, obtained from HELCATS (A. Isavnin). '+ \
    'Timerange: '+win.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+win.time[-1].strftime("%d-%b-%Y %H:%M:%S") +\
    'Units are btxyz [nT, SCEQ], vtxyz [km/s, SCEQ], np [#/cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    
    pickle.dump([win,hwin], open(data_path+ "helcats/wind_2007_2018_helcats.p", "wb" ) )
    print( 'convert Wind done.')   
     
    
    
    
    

    print( 'read STEREO-A')
    stain= pickle.load( open(datacat_path+ "STA_2007to2015_SCEQ.p", "rb" ) )
    stain_time=parse_time(stain.time,format='utime').datetime
    stain=stain.astype([('time', 'object'), ('bt', 'float'),\
                    ('bx', 'float'), ('by', 'float'), ('bz', 'float'), \
                    ('vt', 'float'), ('vx', 'float'), ('vy', 'float'), \
                    ('vz', 'float'), ('tp', 'float'), ('np', 'float'), \
                    ('r', 'float'),('lat', 'float'), ('lon', 'float')])     
    sta=np.zeros(np.size(stain),[('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('vt', 'float64'), ('vx', 'float64'), ('vy', 'float64'), \
                    ('vz', 'float64'), ('tp', 'float64'), ('np', 'float64'),\
                    ('x', 'float64'),('y', 'float64'), ('z', 'float64'),\
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64')])      
    sta = sta.view(np.recarray)                                                       
    
    sta.time=stain_time
    sta.bx=stain.bx
    sta.by=stain.by
    sta.bz=stain.bz
    sta.bt=stain.bt

    sta.vt=stain.vt
    sta.vx=stain.vx
    sta.vy=stain.vy
    sta.vz=stain.vz

    sta.np=stain.np
    sta.tp=stain.tp
 
    sta.r=stain.r/(astropy.constants.au.value/1e3)
    sta.lat=stain.lat
    sta.lon=stain.lon
     
        

    [sta.x, sta.y, sta.z]=sphere2cart(sta.r,np.abs(sta.lat-np.radians(90)),sta.lon)
    sta.lon=np.rad2deg(sta.lon)   
    sta.lat=np.rad2deg(sta.lat)

    del(stain)
  
    hsta='STEREO-A merged magnetic field and plasma data, obtained from HELCATS (A. Isavnin). '+ \
    'Timerange: '+sta.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+sta.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
    'Units are btxyz [nT, SCEQ], vtxyz [km/s, SCEQ], np [#/cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    pickle.dump([sta,hsta], open(data_path+ "helcats/stereoa_2007_2015_helcats.p", "wb" ) )
    print( 'read STA done.')
    
  
  
    
    
    print( 'read STEREO-B')
    stbin= pickle.load( open(datacat_path+ "STB_2007to2014_SCEQ.p", "rb" ) )
    stbin_time=parse_time(stbin.time,format='utime').datetime
    stbin=stbin.astype([('time', 'object'), ('bt', 'float'),\
                    ('bx', 'float'), ('by', 'float'), ('bz', 'float'), \
                    ('vt', 'float'), ('vx', 'float'), ('vy', 'float'), \
                    ('vz', 'float'), ('tp', 'float'), ('np', 'float'), \
                    ('r', 'float'),('lat', 'float'), ('lon', 'float')])     
    stb=np.zeros(np.size(stbin),[('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('vt', 'float64'), ('vx', 'float64'), ('vy', 'float64'), \
                    ('vz', 'float64'), ('tp', 'float64'), ('np', 'float64'),\
                    ('x', 'float64'),('y', 'float64'), ('z', 'float64'),\
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64')])      
    stb = stb.view(np.recarray)                                                       
    
    stb.time=stbin_time
    stb.bx=stbin.bx
    stb.by=stbin.by
    stb.bz=stbin.bz
    stb.bt=stbin.bt

    stb.vt=stbin.vt
    stb.vx=stbin.vx
    stb.vy=stbin.vy
    stb.vz=stbin.vz

    stb.np=stbin.np
    stb.tp=stbin.tp
 
    stb.r=stbin.r/(astropy.constants.au.value/1e3)
    stb.lat=stbin.lat
    stb.lon=stbin.lon
    
    
    
    
     
    [stb.x, stb.y, stb.z]=sphere2cart(stb.r,np.abs(stb.lat-np.radians(90)),stb.lon)
    stb.lon=np.rad2deg(stb.lon)   
    stb.lat=np.rad2deg(stb.lat)
    
    
    #replace missing 2014 plasma data for STEREO-B    
    filestb2='stereob_2013_2014.p'
    [stb2,hstb2]=pickle.load(open(data_path+filestb2, "rb" ) )   
    
    stb_time_mat=parse_time(stb.time).plot_date
    stb2_time_mat=parse_time(stb2.time).plot_date
    
    #interpolate times onto stb.time
    dumvt=np.interp(stb_time_mat, stb2_time_mat,stb2.vt)
    dumtp=np.interp(stb_time_mat, stb2_time_mat,stb2.tp)
    dumnp=np.interp(stb_time_mat, stb2_time_mat,stb2.np)
    #get indices of 1-1-2014 to end
    begin=np.where(stb_time_mat > parse_time('2014-1-1').plot_date )[0][0]
    end=np.size(stb.vt)
    stb.vt[begin:end]=dumvt[begin:end]   
    stb.tp[begin:end]=dumtp[begin:end]   
    stb.np[begin:end]=dumnp[begin:end]   


    

    del(stbin)
    del(dumvt)
    del(dumtp)
    del(dumnp)



    hstb='STEREO-B merged magnetic field and plasma data, obtained from HELCATS (A. Isavnin). '+ \
    'Timerange: '+stb.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+stb.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
    'Units are btxyz [nT, SCEQ], vtxyz [km/s, SCEQ], np [#/cm-3], tp[K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
    pickle.dump([stb,hstb], open(data_path+ "helcats/stereob_2007_2014_helcats.p", "wb" ) )
    print( 'read STB done.')
    
     
      
      
      

    print( 'read MESSENGER')
    #get insitu data from helcats, converted from IDL .sav to pickle
    if removed == True:   
       mesin= pickle.load( open( datacat_path+"MES_2007to2015_SCEQ_removed.p", "rb" ) )
    #non removed dataset
    if removed == False:  
       mesin= pickle.load( open( datacat_path+"MES_2007to2015_SCEQ_non_removed.p", "rb" ) )
    
    
    #time conversion
    mesin_time=parse_time(mesin.time,format='utime').datetime
    #replace mes.time with datetime object
    #new variable names                  
  
    mes=np.zeros(np.size(mesin),[('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('x', 'float64'),('y', 'float64'), ('z', 'float64'),\
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64')])
                    
    #convert to recarray
    mes = mes.view(np.recarray)  
    #set time new  
    mes.time=mesin_time
    mes.bx=mesin.bx
    mes.by=mesin.by
    mes.bz=mesin.bz
    mes.bt=mesin.btot
    
    
    #convert distance from Sun from km to AU, astropy constant is given in m
    mes.r=mesin.mes_radius_in_km_heeq/(astropy.constants.au.value/1e3)
    mes.lon=np.degrees(mesin.mes_longitude_in_radians_heeq.astype('float64'))   
    mes.lat=np.degrees(mesin.mes_latitude_in_radians_heeq.astype('float64'))
    [mes.x, mes.y, mes.z]=sphere2cart(mes.r,np.radians(np.abs(mes.lat-90)),np.radians(mes.lon))
    
    
    if removed == True:   
      hmes='MESSENGER magnetic field data, obtained from NASA PDS. '+ \
      'Timerange: '+mes.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+mes.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
      '. The magnetosphere is removed with a manual magnetopause crossings list (Lydia Philpott, Reka Winslow, Brian Anderson). '+ \
      'Units are btxyz [nT, SCEQ], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
      pickle.dump([mes,hmes], open(data_path+ "helcats/messenger_2007_2015_helcats_removed.p", "wb" ) )
     
    
    if removed == False:  
       hmes='MESSENGER magnetic field data, obtained from NASA PDS. '+ \
       'Timerange: '+mes.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+mes.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
       '. The magnetosphere is removed with a manual magnetopause crossings list (Lydia Philpott, Reka Winslow, Brian Anderson). '+ \
       'Units are btxyz [nT, SCEQ], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
       pickle.dump([mes,hmes], open(data_path+ "helcats/messenger_2007_2015_helcats.p", "wb" ) )
    
    print('convert MESSENGER done.')    
    
        


    print ('read VEX')
    if removed == True:   
         vexin= pickle.load( open(datacat_path+ "VEX_2007to2014_SCEQ_removed.p", "rb" ) )
    if removed == False:  
         vexin= pickle.load( open(datacat_path+ "VEX_2007to2014_SCEQ.p", "rb" ) )
    #time conversion
    vexin_time=parse_time(vexin.time,format='utime').datetime
    vex=np.zeros(np.size(vexin),[('time', 'object'), ('bt', 'float64'),\
                    ('bx', 'float64'), ('by', 'float64'), ('bz', 'float64'), \
                    ('x', 'float64'),('y', 'float64'), ('z', 'float64'),\
                    ('r', 'float64'),('lat', 'float64'), ('lon', 'float64'),\
                    ('xo', 'float64'),('yo', 'float64'), ('zo', 'float64'),\
                    ('ro', 'float64'),('lato', 'float64'), ('lono', 'float64')])      
    vex = vex.view(np.recarray)                                  
 
 
    vex.r=vexin.vex_radius_in_km_heeq/(astropy.constants.au.value/1e3)
    vex.lon=np.rad2deg(vexin.vex_longitude_in_radians_heeq)   
    vex.lat=np.rad2deg(vexin.vex_latitude_in_radians_heeq)

    
    [vex.x, vex.y, vex.z]=sphere2cart(vex.r,np.radians(np.abs(vex.lat-90)),np.radians(vex.lon))
    #convert to degree

    vex.time=vexin_time
    vex.bx=vexin.bx
    vex.by=vexin.by
    vex.bz=vexin.bz
    vex.bt=vexin.btot

       
    #add orbit position 
    #https://www.cosmos.esa.int/web/spice/spice-for-vex
    #or from VEX_2007to2014_VSO.p
    vex2in=pickle.load( open( datacat_path+"VEX_2007to2014_VSO.p", "rb" ) )
  
    vex.xo=vex2in.x
    vex.yo=vex2in.y
    vex.zo=vex2in.z
    
    [vex.ro, vex.lato, vex.lono]=cart2sphere(vex.xo,vex.yo,vex.zo)
    vex.lono=np.rad2deg(vex.lono)   
    vex.lato=np.rad2deg(vex.lato)
    
    del(vexin)
    del(vex2in)

  
    
    if removed == True:   
     hvex='VEX magnetic field data, obtained from the VEX magnetometer PI T. Zhang IWF Graz, Austria. '+ \
     'Timerange: '+vex.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+vex.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
     '. The magnetosphere was removed with the model from Zhang et al. (2008), see Moestl et al. (2017, doi: 10.1002/2017SW001614) for details. '+ \
     'Units are btxyz [nT, SCEQ], orbital position: '+ \
     'xo/yo/zo/ro/lono/lato [km, degree, VSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
     pickle.dump([vex,hvex], open(data_path+ "helcats/vex_2007_2014_helcats_removed.p", "wb" ) )
     
    
    if removed == False:  
     hvex='VEX magnetic field data, obtained from the VEX magnetometer PI T. Zhang IWF Graz, Austria. '+ \
     'Timerange: '+vex.time[0].strftime("%d-%b-%Y %H:%M:%S")+' to '+vex.time[-1].strftime("%d-%b-%Y %H:%M:%S")+\
     '. Units are btxyz [nT, SCEQ], orbital position: '+ \
     'xo/yo/zo/ro/lono/lato [km, degree, VSO], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]'
     pickle.dump([vex,hvex], open(data_path+ "helcats/vex_2007_2014_helcats.p", "wb" ) )
    
    print( 'convert VEX done.')
    
    
    
    

    
    #the Ulysses file has been generated by selecting the merged Ulysses data in CDAWEB
    #and then saved as one cdf 2.7 file
    print('read Ulysses from CDAWEB cdf')    
    save_ulysses_data(data_path)
    #fileuly=data_path+'ulysses_1990_2009_helcats.p'
    #[uly,huly]=pickle.load(open(fileuly, 'rb' ) )
  

    #if removed==True: 
    #    pickle.dump([vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly], open(data_path+ "helcats_all_data_removed.p", "wb" ) )
    #    print('saved as ' +data_path+ 'helcats_all_data_removed.p')

    #if removed==False: 
    #    pickle.dump([vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly], open(data_path+ "helcats_all_data_non_removed.p", "wb" ) )
    #    print('saved as ' +data_path+ 'helcats_all_data_non_removed.p')




  
def load_helcats_datacat(file):  
    ''' to load all of helcats DATACAT from a single file'''
    
    print('load all helcats DATACAT from single file: ', file)
    [vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly]=pickle.load( open(file, "rb" ) )
    print('Use vex,win,sta,stb,mes,uly to access data and position, hvex,hwin, hmes, hsta, hstb, huly for headers.')
    return [vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly]
    
    
    
    
def recarray_to_numpy_array(rec):    
    
    '''convert data recarray to numpy structured array with matplotlib time '''

    #recarr.time=parse_time(recarr.time).plot_date
    #numarr = pd.DataFrame(recarr).to_numpy()
    
    num=copy.deepcopy(np.array(rec))
    num['time']=parse_time(num['time']).plot_date
    
    return num
    
    
  
    
    
    
#################################### MATH ################################################


   
    
  
   
def convert_HEE_to_HEEQ_old(sc_in):
    '''
    for Wind, Bepi magnetic field components: convert HEE to HAE to HEEQ
    '''

    print('conversion HEE to HEEQ')                                
    
    sc=copy.deepcopy(sc_in)
    
    jd=np.zeros(len(sc))
    mjd=np.zeros(len(sc))
        

    for i in np.arange(0,len(sc)):

        jd[i]=parse_time(sc.time[i]).jd
        mjd[i]=float(int(jd[i]-2400000.5)) #use modified julian date    
        
       
        b_hee=[sc.bx[i],sc.by[i],sc.bz[i]]
        
        #HEE to HAE        
        
        #define T00 and UT
        T00=(mjd[i]-51544.5)/36525.0          
        dobj=sc.time[i]
        UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   

        #lambda_sun in Hapgood, equation 5, here in rad
        M=np.radians(357.528+35999.050*T00+0.04107*UT)
        LAMBDA=280.460+36000.772*T00+0.04107*UT        
        lambda_sun=np.radians( (LAMBDA+(1.915-0.0048*T00)*np.sin(M)+0.020*np.sin(2*M)) )
        
        #S-1 Matrix equation 12 hapgood 1992, change sign in lambda angle for inversion HEE to HAE instead of HAE to HEE
        c, s = np.cos(-(lambda_sun+np.radians(180))), np.sin(-(lambda_sun+np.radians(180)))
        Sm1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
        b_hae=np.dot(Sm1,b_hee)

        #HAE to HEEQ
        
        iota=np.radians(7.25)
        omega=np.radians((73.6667+0.013958*((mjd[i]+3242)/365.25)))                      
        theta=np.arctan(np.cos(iota)*np.tan(lambda_sun-omega))                       
                      
    
        #quadrant of theta must be opposite lambda_sun minus omega; Hapgood 1992 end of section 5   
        #get lambda-omega angle in degree mod 360 and theta in degrees
        lambda_omega_deg=np.mod(np.degrees(lambda_sun)-np.degrees(omega),360)
        theta_node_deg=np.degrees(theta)


        ##if the 2 angles are close to similar, so in the same quadrant, then theta_node = theta_node +pi           
        if np.logical_or(abs(lambda_omega_deg-theta_node_deg) < 1, abs(lambda_omega_deg-360-theta_node_deg) < 1): theta=theta+np.pi                                                                                                          
        
        #rotation around Z by theta
        c, s = np.cos(theta), np.sin(theta)
        S2_1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))

        #rotation around X by iota  
        iota=np.radians(7.25)
        c, s = np.cos(iota), np.sin(iota)
        S2_2 = np.array(( (1,0,0), (0,c, s), (0, -s, c)) )
                
        #rotation around Z by Omega  
        c, s = np.cos(omega), np.sin(omega)
        S2_3 = np.array( ((c,s, 0), (-s, c, 0), (0, 0, 1)) )
        
        #matrix multiplication to go from HAE to HEEQ components                
        [bx_heeq,by_heeq,bz_heeq]=np.dot(  np.dot(   np.dot(S2_1,S2_2),S2_3), b_hae) 
        
        sc.bx[i]=bx_heeq
        sc.by[i]=by_heeq
        sc.bz[i]=bz_heeq

    
    print('conversion HEE to HEEQ done')  
    
    return sc





    
    

    
    

   
def convert_HEEQ_to_GSE(sc_in):
    '''
    for solar orbiter magnetic field components: convert HEEQ to HAE to HEE to GSE 
    
    hapgood 1992 https://ui.adsabs.harvard.edu/abs/1992P%26SS...40..711H/abstract
    '''
    
    sc=copy.deepcopy(sc_in)

    print('conversion HEEQ to GSE')                                

    jd=np.zeros(len(sc))
    mjd=np.zeros(len(sc))
        

    for i in np.arange(0,len(sc)):
        
        
        jd[i]=parse_time(sc.time[i]).jd
        mjd[i]=float(int(jd[i]-2400000.5)) #use modified julian date    
        
        
        #general parameters
        
        #define T00 and UT
        T00=(mjd[i]-51544.5)/36525.0          
        dobj=sc.time[i]
        UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   

        #lambda_sun in Hapgood, equation 5, here in rad
        M=np.radians(357.528+35999.050*T00+0.04107*UT)
        LAMBDA=280.460+36000.772*T00+0.04107*UT        
        lambda_sun=np.radians( (LAMBDA+(1.915-0.0048*T00)*np.sin(M)+0.020*np.sin(2*M)) )
        
        

        b_heeq=[sc.bx[i],sc.by[i],sc.bz[i]]

        #HEEQ to HAE - Hapgood equation 13 inverted, angles are set to -angle later, revert order of matrices too
        
        iota=np.radians(7.25)
        omega=np.radians((73.6667+0.013958*((mjd[i]+3242)/365.25)))                      
        theta=np.arctan(np.cos(iota)*np.tan(lambda_sun-omega))                                             
    
        #quadrant of theta must be opposite lambda_sun minus omega; Hapgood 1992 end of section 5   
        #get lambda-omega angle in degree mod 360 and theta in degrees
        lambda_omega_deg=np.mod(np.degrees(lambda_sun)-np.degrees(omega),360)
        theta_node_deg=np.degrees(theta)

        ##if the 2 angles are close to similar, so in the same quadrant, then theta_node = theta_node +pi           
        if np.logical_or(abs(lambda_omega_deg-theta_node_deg) < 1, abs(lambda_omega_deg-360-theta_node_deg) < 1): theta=theta+np.pi                                                                                                          
        
        #rotation around Z by -theta
        c, s = np.cos(-theta), np.sin(-theta)
        St = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))

        #rotation around X by -iota  
        iota=np.radians(7.25)
        c, s = np.cos(-iota), np.sin(-iota)
        Si = np.array(( (1,0,0), (0,c, s), (0, -s, c)) )
                
        #rotation around Z by -Omega  
        c, s = np.cos(-omega), np.sin(-omega)
        So = np.array( ((c,s, 0), (-s, c, 0), (0, 0, 1)) )
        
        #matrix multiplication to go from HAE to HEEQ components - inverted matrix to Hapgood equation 13, 
        #so use inverted order here: omega - iota - theta

        b_hae=np.dot(  np.dot(   np.dot(So,Si),St), b_heeq) 
        
        #HAE to HEE        
        
        #S1 Matrix equation 12 hapgood 1992
        c, s = np.cos((lambda_sun+np.radians(180))), np.sin((lambda_sun+np.radians(180)))
        S1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
        b_hee=np.dot(S1,b_hae)

        
        #HEE to GSE
        #Hapgood 1992 rotation by 180 degrees, or simply change sign in bx by    
        #rotangle=np.radians(180)
        #c, s = np.cos(rotangle), np.sin(rotangle)
        #T1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
        #[bx_hee,by_hee,bz_hee]=T1[sc.bx[i],sc.by[i],sc.bz[i]]        
        [bx_gse,by_gse,bz_gse]=[-b_hee[0],-b_hee[1],b_hee[2]]

         
        sc.bx[i]=bx_gse
        sc.by[i]=by_gse
        sc.bz[i]=bz_gse

    
    print('conversion HEEQ to GSE done')                                
    
    return sc

    
 


   

def convert_HEEQ_to_SCEQ(sc_in):
    '''
    sc is the input recarray for the data, e.g. psp, sta    
    '''
        
    sc=copy.deepcopy(sc_in)
    print('HEEQ to SCEQ')
    
    sc_len=len(sc)
   
    #HEEQ - make to SCEQ for each timestep,   solar rotation axis at 0, 0, 1 in HEEQ
    X_heeq=[1,0,0]
    Y_heeq=[0,1,0]
    Z_heeq=[0,0,1]
        
    # go through all data points    
    for i in np.arange(0,sc_len):        
        
        #rotation for HEEQ to SCEQ vectors, depending on current longitude (correct c, s)
        rotangle=np.radians(sc.lon[i])
        c, s = np.cos(rotangle), np.sin(rotangle)
        #rotation matrix around Z
        R = np.array(((c,-s, 0), (s, c, 0), (0, 0, 1)))
        
        #rotate X and Y  
        X_sceq=np.dot(R,X_heeq)
        Y_sceq=np.dot(R,Y_heeq)         
       
        #project into new system
        sc.bx[i]=np.dot([sc_in.bx[i],sc_in.by[i],sc_in.bz[i]],X_sceq)
        sc.by[i]=np.dot([sc_in.bx[i],sc_in.by[i],sc_in.bz[i]],Y_sceq)
    
    return sc


    

   


  
    
    
   
def convert_RTN_to_SCEQ(sc_in,name):
    '''
    for STEREO-A, B, Solar Orbiter and Parker Solar Probe   
    sc is the input recarray for the data, e.g. psp, sta    
    '''
    
    
    sc=copy.deepcopy(sc_in)
    print('RTN to SCEQ')
    
    sc_len=len(sc)
   
    #HEEQ - make to SCEQ for each timestep,   solar rotation axis at 0, 0, 1 in HEEQ
    X_heeq=[1,0,0]
    Y_heeq=[0,1,0]
    Z_heeq=[0,0,1]
        
    # go through all data points    
    for i in np.arange(0,sc_len):        
        
        #rotation for HEEQ to SCEQ vectors, depending on current longitude (correct c, s)
        rotangle=np.radians(sc.lon[i])
        c, s = np.cos(rotangle), np.sin(rotangle)
        #rotation matrix around Z
        R = np.array(((c,-s, 0), (s, c, 0), (0, 0, 1)))
        
        #rotate X and Y  
        X_sceq=np.dot(R,X_heeq)
        Y_sceq=np.dot(R,Y_heeq)
         
        #make normalized RTN vectors
        Xrtn=[sc.x[i], sc.y[i],sc.z[i]]/np.linalg.norm([sc.x[i], sc.y[i],sc.z[i]])
        Yrtn=np.cross(Z_heeq,Xrtn)/np.linalg.norm(np.cross(Z_heeq,Xrtn))
        Zrtn=np.cross(Xrtn, Yrtn)/np.linalg.norm(np.cross(Xrtn, Yrtn))
        
        #project into new system
        sc.bx[i]=np.dot(np.dot(sc_in.bx[i],Xrtn)+np.dot(sc_in.by[i],Yrtn)+np.dot(sc_in.bz[i],Zrtn),X_sceq)
        sc.by[i]=np.dot(np.dot(sc_in.bx[i],Xrtn)+np.dot(sc_in.by[i],Yrtn)+np.dot(sc_in.bz[i],Zrtn),Y_sceq)
        sc.bz[i]=np.dot(np.dot(sc_in.bx[i],Xrtn)+np.dot(sc_in.by[i],Yrtn)+np.dot(sc_in.bz[i],Zrtn),Z_heeq)
        
        #project into new system
        if name=='PSP':
            sc.vx[i]=np.dot(np.dot(sc_in.vx[i],Xrtn)+np.dot(sc_in.vy[i],Yrtn)+np.dot(sc_in.vz[i],Zrtn),X_sceq)
            sc.vy[i]=np.dot(np.dot(sc_in.vx[i],Xrtn)+np.dot(sc_in.vy[i],Yrtn)+np.dot(sc_in.vz[i],Zrtn),Y_sceq)
            sc.vz[i]=np.dot(np.dot(sc_in.vx[i],Xrtn)+np.dot(sc_in.vy[i],Yrtn)+np.dot(sc_in.vz[i],Zrtn),Z_heeq)


    
    return sc




       
    

def convert_RTN_to_HEEQ(sc_in,name):
    '''
    for STEREO-A, B, Solar Orbiter and Parker Solar Probe   
    sc is the input recarray for the data, e.g. psp, sta    
    '''
    
    sc=copy.deepcopy(sc_in)
    
    sc_len=len(sc)
    
    print('conversion RTN to HEEQ')
    #HEEQ unit vectors (same as spacecraft xyz position)
    X_heeq=[1,0,0]
    Y_heeq=[0,1,0]
    Z_heeq=[0,0,1]
        
    # go through all data points    
    for i in np.arange(0,sc_len):                
       
        #make normalized RTN unit vectors from spacecraft position in HEEQ basis
        Xrtn=[sc.x[i], sc.y[i],sc.z[i]]/np.linalg.norm([sc.x[i], sc.y[i],sc.z[i]])
        Yrtn=np.cross(Z_heeq,Xrtn)/np.linalg.norm(np.cross(Z_heeq,Xrtn))
        Zrtn=np.cross(Xrtn, Yrtn)/np.linalg.norm(np.cross(Xrtn, Yrtn))
        
        #project into new system (HEEQ)
        sc.bx[i]=np.dot(np.dot(sc_in.bx[i],Xrtn)+np.dot(sc_in.by[i],Yrtn)+np.dot(sc_in.bz[i],Zrtn),X_heeq)
        sc.by[i]=np.dot(np.dot(sc_in.bx[i],Xrtn)+np.dot(sc_in.by[i],Yrtn)+np.dot(sc_in.bz[i],Zrtn),Y_heeq)
        sc.bz[i]=np.dot(np.dot(sc_in.bx[i],Xrtn)+np.dot(sc_in.by[i],Yrtn)+np.dot(sc_in.bz[i],Zrtn),Z_heeq)
        
        #project into new system (HEEQ) for speed 
        if name=='PSP':
            sc.vx[i]=np.dot(np.dot(sc_in.vx[i],Xrtn)+np.dot(sc_in.vy[i],Yrtn)+np.dot(sc_in.vz[i],Zrtn),X_heeq)
            sc.vy[i]=np.dot(np.dot(sc_in.vx[i],Xrtn)+np.dot(sc_in.vy[i],Yrtn)+np.dot(sc_in.vz[i],Zrtn),Y_heeq)
            sc.vz[i]=np.dot(np.dot(sc_in.vx[i],Xrtn)+np.dot(sc_in.vy[i],Yrtn)+np.dot(sc_in.vz[i],Zrtn),Z_heeq)

    print('conversion RTN to HEEQ done')     
    
    return sc






       

##########################################################################################################










'''



   
def convert_GSE_to_HEEQ_testing(sc):

    print('conversion GSE to HEEQ start')                                

    jd=np.zeros(len(sc))
    mjd=np.zeros(len(sc))


    angles1=[]

    for i in np.arange(0,len(sc)):

        jd[i]=parse_time(sc.time[i]).jd
        mjd[i]=float(int(jd[i]-2400000.5)) #use modified julian date    

        #GSE to HEE
        #Hapgood 1992 rotation by 180 degrees, or simply change sign in bx by    
        #rotangle=np.radians(180)
        #c, s = np.cos(rotangle), np.sin(rotangle)
        #T1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
        #[bx_hee,by_hee,bz_hee]=T1[sc.bx[i],sc.by[i],sc.bz[i]]        
        b_hee=[-sc.bx[i],-sc.by[i],sc.bz[i]]

        #HEE to HAE     - one rotation   

        #define T00 and UT
        T00=(mjd[i]-51544.5)/36525.0          
        dobj=sc.time[i]
        UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   

        #lambda_sun in Hapgood, equation 5, here in rad
        M=np.radians(357.528+35999.050*T00+0.04107*UT)
        LAMBDA=280.460+36000.772*T00+0.04107*UT        
        lambda_sun=np.radians( (LAMBDA+(1.915-0.0048*T00)*np.sin(M)+0.020*np.sin(2*M)) )

        #S-1 Matrix equation 12 hapgood 1992, change sign in lambda angle
        c, s = np.cos(-(lambda_sun+np.radians(180))), np.sin(-(lambda_sun+np.radians(180)))
        Sm1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
        b_hae=np.dot(Sm1,b_hee)


        #HAE to HEEQ

        iota=np.radians(7.25)
        omega=np.radians((73.6667+0.013958*((mjd[i]+3242)/365.25)))                      
        theta=np.arctan(np.cos(iota)*np.tan(lambda_sun-omega))  


        #quadrant of theta must be opposite lambda_sun minus omega; Hapgood 1992 end of section 5   
        #get lambda-omega angle in degree mod 360 and theta in degrees
        lambda_omega_deg=np.mod(np.degrees(lambda_sun)-np.degrees(omega),360)
        theta_node_deg=np.degrees(theta)

        ##if the 2 angles are close to similar, so in the same quadrant, then theta_node = theta_node +pi           
        if np.logical_or(abs(lambda_omega_deg-theta_node_deg) < 1, abs(lambda_omega_deg-360-theta_node_deg) < 1): theta=theta+np.pi                                                            


        #convert again for array to check    
        theta_node_deg=np.degrees(theta)

        angles1.append([i,np.round(theta_node_deg,1),np.round(lambda_omega_deg,1)])


        #rotation around Z by theta
        c, s = np.cos(theta), np.sin(theta)
        S2_1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))

        #rotation around X by iota  
        iota=np.radians(7.25)
        c, s = np.cos(iota), np.sin(iota)
        S2_2 = np.array(( (1,0,0), (0,c, s), (0, -s, c)) )

        #rotation around Z by Omega  
        c, s = np.cos(omega), np.sin(omega)
        S2_3 = np.array( ((c,s, 0), (-s, c, 0), (0, 0, 1)) )

        #matrix multiplication to go from HAE to HEEQ components                
        [bx_heeq,by_heeq,bz_heeq]=np.dot(  np.dot(   np.dot(S2_1,S2_2),S2_3), b_hae) 

        sc.bx[i]=bx_heeq
        sc.by[i]=by_heeq
        sc.bz[i]=bz_heeq
        
        

    ang=np.array(angles1)


    plt.plot(wing.time,wing.bx,'-k',label='gse')
    plt.plot(sc.time,sc.bx,'-g',label='heeq_new')
    plt.plot(sc.time,ang[:,1],'-r',label='theta')
    plt.plot(sc.time,ang[:,2],'-b',label='lambda - omega mod 360')

    plt.legend()





    print('conversion GSE to HEEQ done')    
########### test Wind HEEQ conversion

print('conversion GSE to HEEQ start')                                

jd=np.zeros(len(sc))
mjd=np.zeros(len(sc))


angles1=[]

for i in np.arange(0,len(sc)):

    jd[i]=parse_time(sc.time[i]).jd
    mjd[i]=float(int(jd[i]-2400000.5)) #use modified julian date    

    #GSE to HEE
    #Hapgood 1992 rotation by 180 degrees, or simply change sign in bx by    
    #rotangle=np.radians(180)
    #c, s = np.cos(rotangle), np.sin(rotangle)
    #T1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
    #[bx_hee,by_hee,bz_hee]=T1[sc.bx[i],sc.by[i],sc.bz[i]]        
    b_hee=[-sc.bx[i],-sc.by[i],sc.bz[i]]

    #HEE to HAE     - one rotation   

    #define T00 and UT
    T00=(mjd[i]-51544.5)/36525.0          
    dobj=sc.time[i]
    UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   

    #lambda_sun in Hapgood, equation 5, here in rad
    M=np.radians(357.528+35999.050*T00+0.04107*UT)
    LAMBDA=280.460+36000.772*T00+0.04107*UT        
    lambda_sun=np.radians( (LAMBDA+(1.915-0.0048*T00)*np.sin(M)+0.020*np.sin(2*M)) )

    #S-1 Matrix equation 12 hapgood 1992, change sign in lambda angle
    c, s = np.cos(-(lambda_sun+np.radians(180))), np.sin(-(lambda_sun+np.radians(180)))
    Sm1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))
    b_hae=np.dot(Sm1,b_hee)


    #HAE to HEEQ

    iota=np.radians(7.25)
    omega=np.radians((73.6667+0.013958*((mjd[i]+3242)/365.25)))                      
    theta=np.arctan(np.cos(iota)*np.tan(lambda_sun-omega))  
    
    

    #quadrant of theta must be opposite lambda_sun minus omega; Hapgood 1992 end of section 5   
    #get lambda-omega angle in degree mod 360 and theta in degrees
    lambda_omega_deg=np.mod(np.degrees(lambda_sun)-np.degrees(omega),360)
    theta_node_deg=np.degrees(theta)
    
       
    ##if the 2 angles are close to similar, so in the same quadrant, then theta_node = theta_node +pi           
    if np.logical_or(abs(lambda_omega_deg-theta_node_deg) < 1, abs(lambda_omega_deg-360-theta_node_deg) < 1): theta=theta+np.pi                                                            
        

    #convert again for array to check    
    theta_node_deg=np.degrees(theta)
    
    angles1.append([i,np.round(theta_node_deg,1),np.round(lambda_omega_deg,1)])

    
    #rotation around Z by theta
    c, s = np.cos(theta), np.sin(theta)
    S2_1 = np.array(((c,s, 0), (-s, c, 0), (0, 0, 1)))

    #rotation around X by iota  
    iota=np.radians(7.25)
    c, s = np.cos(iota), np.sin(iota)
    S2_2 = np.array(( (1,0,0), (0,c, s), (0, -s, c)) )

    #rotation around Z by Omega  
    c, s = np.cos(omega), np.sin(omega)
    S2_3 = np.array( ((c,s, 0), (-s, c, 0), (0, 0, 1)) )

    #matrix multiplication to go from HAE to HEEQ components                
    [bx_heeq,by_heeq,bz_heeq]=np.dot(  np.dot(   np.dot(S2_1,S2_2),S2_3), b_hae) 

    sc.bx[i]=bx_heeq
    sc.by[i]=by_heeq
    sc.bz[i]=bz_heeq


print('conversion GSE to HEEQ done')                                



ang=np.array(angles1)

    
plt.plot(wing.time,wing.bx,'-k',label='gse')
plt.plot(sc.time,sc.bx,'-g',label='heeq_new')
plt.plot(sc.time,ang[:,1],'-r',label='theta')
plt.plot(sc.time,ang[:,2],'-b',label='lambda - omega mod 360')

plt.legend()

'''
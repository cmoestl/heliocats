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
- despike sta stb wind all, new B and V for STA, Wind and PSP converted to SCEQ components, plasma correct for new PSP, wind, sta ...


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

import numpy as np
import scipy.io
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import  DateFormatter
import seaborn as sns
import datetime
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


from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging

from heliocats import data as hd
importlib.reload(hd) #reload again while debugging

from heliocats import cats as hc
importlib.reload(hc) #reload again while debugging

#where the 6 in situ data files are located is read from input.py
#as data_path=....
from input import *


########### make directories first time if not there

resdir='results'
if os.path.isdir(resdir) == False: os.mkdir(resdir)

datadir='data'
if os.path.isdir(datadir) == False: os.mkdir(datadir)

indexdir='data/indices_icmecat' 
if os.path.isdir(indexdir) == False: os.mkdir(indexdir) 

catdir='icmecat'
if os.path.isdir(datadir) == False: os.mkdir(catdir)

    
##########################################################################################
######################################## MAIN PROGRAM ####################################
##########################################################################################

############################## (1) load new data with HelioSat and heliocats.data ########

   
load_data=1

if load_data > 0:

    print('load HELCATS data until 2018 and new Wind, STEREO-A, MAVEN, Parker Solar Probe data')

    # MAVEN
    #filemav='maven_2014_2018.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )

    #filemav='maven_2014_2018_removed.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )
    
    filemav='maven_2014_2018_removed_smoothed.p'
    [mav,hmav]=pickle.load(open(data_path+filemav, 'rb' ) )

    # Wind
    filewin="wind_2018_2020.p" 
    #for updating data
    #start=datetime.datetime(2018, 1, 1)
    #end=datetime.datetime.utcnow()
    #hd.save_wind_data(data_path,filewin,start,end)
    [win2,hwin2]=pickle.load(open(data_path+filewin, "rb" ) )  

    # STEREO-A    
    filesta2='sta_2018_2019_beacon.p'
    #start=datetime.datetime(2018, 1, 1)
    #end=datetime.datetime(2019, 12, 31)
    #hd.save_stereoa_beacon_data(data_path,filesta,start,end)
   
    sta2=pickle.load(open(data_path+filesta2, "rb" ) )  

    # Parker Solar Probe
    filepsp='psp_2018_2019.p'
    [psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) )  


    # ADD BepiColombo  
    
    
    # ADD Solar Orbiter
    
    
    # Ulysses is currently taken from the full 
    # helcats data below, but a file is available on figshare


    # get data file from helcats with headers
    [vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly]=hd.load_helcats_datacat(data_path+'helcats_all_data_removed.p') 



################################ (2) measure new events ##################################








################################ (3) make ICMECAT  #######################################

print('data loaded')
ic=hc.load_helcats_icmecat_master_from_excel('icmecat/HELCATS_ICMECAT_v20_master.xlsx')

####### 3a get indices for all spacecraft
wini=np.where(ic.sc_insitu == 'Wind')[:][0] 
vexi=np.where(ic.sc_insitu == 'VEX')[:][0]  
mesi=np.where(ic.sc_insitu == 'MESSENGER')[:][0]   
stai=np.where(ic.sc_insitu == 'STEREO-A')[:][0]    
stbi=np.where(ic.sc_insitu == 'STEREO-B')[:][0]    
mavi=np.where(ic.sc_insitu == 'MAVEN')[:][0]    
ulyi=np.where(ic.sc_insitu == 'ULYSSES')[:][0]    
#pspi=np.where(ic.sc_insitu == 'ParkerSolarProbe')[:][0]    



####### 3b get parameters for all spacecraft one after another

ic=hc.get_cat_parameters(win,wini,ic,'Wind')
ic=hc.get_cat_parameters(sta,stai,ic,'STEREO-A')
ic=hc.get_cat_parameters(stb,stbi,ic,'STEREO_B')
ic=hc.get_cat_parameters(mes,mesi,ic,'MESSENGER')
ic=hc.get_cat_parameters(vex,vexi,ic,'VEX')
ic=hc.get_cat_parameters(uly,ulyi,ic,'ULYSSES')
ic=hc.get_cat_parameters(mav,mavi,ic,'MAVEN')


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



#icl=pickle.load(open(file, 'rb' ) )

######################################################################################
################################### END MAIN #########################################
######################################################################################




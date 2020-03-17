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
if os.path.isdir(catdir) == False: os.mkdir(catdir)

icplotsdir='results/plots_icmecat/' 
if os.path.isdir(icplotsdir) == False: os.mkdir(icplotsdir) 

    
##########################################################################################
######################################## MAIN PROGRAM ####################################
##########################################################################################

############################## (1) load new data with HelioSat and heliocats.data ########

   
load_data=1

if load_data > 0:

    print('load new Wind, STEREO-A, MAVEN, Parker Solar Probe data')

    # MAVEN
    #filemav='maven_2014_2018.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )

    #filemav='maven_2014_2018_removed.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )
    
    filemav='maven_2014_2018_removed_smoothed.p'
    [mav,hmav]=pickle.load(open(data_path+filemav, 'rb' ) )

    # Wind
    filewin="wind_2018_2019.p" 
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
   
    [sta2,hsta2]=pickle.load(open(data_path+filesta2, "rb" ) )  

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


#for measuring new events use this function from heliocats.plot
#hp.plot_insitu_measure(psp, '2018-Nov-10','2018-Nov-15', 'PSP', 'results/plots_icmecat/')


#for plotting single events
#hp.plot_insitu(psp, ic.icme,'2018-Nov-15', 'PSP', icplotsdir)





################################ (3) make ICMECAT  #######################################

print('data loaded')
ic=hc.load_helcats_icmecat_master_from_excel('icmecat/HELCATS_ICMECAT_v20_master.xlsx')

####### 3a get indices for all spacecraft

wini=np.where(ic.sc_insitu == 'Wind')[:][0] 
stai=np.where(ic.sc_insitu == 'STEREO-A')[:][0]    
stbi=np.where(ic.sc_insitu == 'STEREO-B')[:][0]    
vexi=np.where(ic.sc_insitu == 'VEX')[:][0]  
mesi=np.where(ic.sc_insitu == 'MESSENGER')[:][0]   
ulyi=np.where(ic.sc_insitu == 'ULYSSES')[:][0]    
mavi=np.where(ic.sc_insitu == 'MAVEN')[:][0]    
pspi=np.where(ic.sc_insitu == 'PSP')[:][0]    



####### 3b get parameters for all spacecraft one after another

ic=hc.get_cat_parameters(win,wini,ic,'Wind')
ic=hc.get_cat_parameters(sta,stai,ic,'STEREO-A')
ic=hc.get_cat_parameters(stb,stbi,ic,'STEREO_B')
ic=hc.get_cat_parameters(vex,vexi,ic,'VEX')
ic=hc.get_cat_parameters(mes,mesi,ic,'MESSENGER')
ic=hc.get_cat_parameters(uly,ulyi,ic,'ULYSSES')
ic=hc.get_cat_parameters(mav,mavi,ic,'MAVEN')
ic=hc.get_cat_parameters(psp,pspi,ic,'PSP')


####### 3c make all plots if wanted

'''
matplotlib.use('Agg')
hp.plot_icmecat_events(sta,stai,ic,'STEREO-A',icplotsdir)
hp.plot_icmecat_events(stb,stbi,ic,'STEREO-B',icplotsdir)
hp.plot_icmecat_events(vex,vexi,ic,'VEX',icplotsdir)
hp.plot_icmecat_events(mes,mesi,ic,'MESSENGER',icplotsdir)
hp.plot_icmecat_events(uly,ulyi,ic,'ULYSSES',icplotsdir)
hp.plot_icmecat_events(mav,mavi,ic,'MAVEN',icplotsdir)
hp.plot_icmecat_events(psp,pspi,ic,'PSP',icplotsdir)
hp.plot_icmecat_events(win,wini,ic,'Wind',icplotsdir)
'''


################################ (4) save ICMECAT #################################


#make header
header=hc.make_icmecat_header(ic)
file='icmecat/HELCATS_ICMECAT_v20_header.txt'
with open(file, "w") as text_file:
    text_file.write(header)
print(header)


##### save ICMECAT as pickle with times as datetime objects
file='icmecat/HELCATS_ICMECAT_v20.p'
pickle.dump(ic, open(file, 'wb'))


################ save to different formats

#copy pandas dataframe first to change time format

ic_copy=copy.deepcopy(ic)  
ic_copy.icme_start_time=parse_time(ic.icme_start_time).iso
ic_copy.mo_start_time=parse_time(ic.mo_start_time).iso
ic_copy.mo_end_time=parse_time(ic.mo_end_time).iso

#change time format
for i in np.arange(len(ic)):

    dum=ic_copy.icme_start_time[i] 
    ic_copy.at[i,'icme_start_time']=dum[0:16]
     
    dum=ic_copy.mo_start_time[i] 
    ic_copy.at[i,'mo_start_time']=dum[0:16]
     
    dum=ic_copy.mo_end_time[i] 
    ic_copy.at[i,'mo_end_time']=dum[0:16]


#save as Excel
file='icmecat/HELCATS_ICMECAT_v20.xlsx'
ic_copy.to_excel(file,sheet_name='ICMECATv2.0')

#save as json
file='icmecat/HELCATS_ICMECAT_v20.json'
ic_copy.to_json(file)

#save as csv
file='icmecat/HELCATS_ICMECAT_v20.csv'
ic_copy.to_csv(file)

#save as html
file='icmecat/HELCATS_ICMECAT_v20_simple.html'
ic_copy.to_html(file)

#save as hdf needs pip install tables
#file='icmecat/HELCATS_ICMECAT_v20.hdf'
#ic.to_hdf(file,key='icmecat')

#save as .mat does not work yet
#ile='icmecat/HELCATS_ICMECAT_v20.mat'
#icdict=ic.to_dict()
#scipy.io.savemat(file,ic.values)


#save as txt
file='icmecat/HELCATS_ICMECAT_v20.txt'
np.savetxt(file, ic_copy.values.astype(str), fmt='%s' )

print('ICMECAT saved as '+file)

sys.exit()

############ save as html file

file='icmecat/HELCATS_ICMECAT_v20.p'
ic=pickle.load( open(file, 'rb'))

#save as html
file='icmecat/HELCATS_ICMECAT_v20.html'
#ic.to_html(file,justify='center')

ichtml='{% extends "_base.html" %} \n \n {% block content %} \n \n \n <p> ICMECAT version 2.0 </p>'
ichtml += ic.to_html()
ichtml +='\n \n {% endblock %}'

'''
{% extends "_base.html" %}
{% block content %}
<div class="row">
    <div class="col-12 col-md-12 col-xl-12">
        <p> </br>
        Recent real-time solar and solar wind data from SDO in Earth orbit, 
        from the DSCOVR or ACE spacecraft at the Sun-Earth L1 point provided by NOAA and from STEREO-A.
        Updated twice daily at 05:00 and 17:00 UTC.</br></br>
        </p>
    </div>
</div>
'''



with open(file,'w') as f:
    f.write(ichtml)
    f.close()


######################################################################################
################################### END MAIN #########################################
######################################################################################

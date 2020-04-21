# https://github.com/cmoestl/heliocats  data_update.py

# for updating data every day 

# MIT LICENSE
# Copyright 2020, Christian Moestl 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this 
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify, 
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
# permit persons to whom the Software is furnished to do so, subject to the following 
# conditions:
# The above copyright notice and this permission notice shall be included in all copies 
# or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import pickle
import importlib
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
import sys
import numpy as np
import datetime
import scipy.signal
import urllib
import json
import os   

#import 
from heliocats import data as hd
importlib.reload(hd) #reload again while debugging

from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging

#for server
matplotlib.use('Agg')
#matplotlib.use('qt5agg')

data_path='/nas/helio/data/insitu_python/'
plot_path='/nas/helio/data/insitu_python/plots/'
noaa_path='/nas/helio/data/noaa_rtsw/'
position_path='/nas/helio/data/insitu_python/plots_positions/'


#########################################################################################


#for easier debugging - do not download and process data but do everything else
get_new_data=True

#~/miniconda/envs/helio/bin/python /home/cmoestl/pycode/heliocats/sc_positions_insitu_orbit1.py
#~/miniconda/envs/helio/bin/python /home/cmoestl/pycode/heliocats/sc_positions_insitu_orbit2.py


################################# PSP update

# load PSP data
# go to heliosat directory psp_fields_l2
# wget "ftps://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_1min/2018/*.cdf"
# wget -nc "ftps://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_1min/2019/*.cdf"
# psp_spc_l3
# wget "ftps://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spc/l3/l3i/2018/*.cdf"
# wget "ftps://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spc/l3/l3i/2019/*.cdf"

# print('load PSP data') #from heliosat, converted to SCEQ similar to STEREO-A/B

#change date in hd.save_psp_data
#filepsp='psp_2018_2019_rtn.p'
#hd.save_psp_data(data_path,filepsp, sceq=False)   

#filepsp='psp_2018_2019_sceq.p'
#hd.save_psp_data(data_path,filepsp, sceq=True)   



#print('download PSP  files without overwriting existing files into ',data_path)
#psp_data_path='/nas/helio/data/heliosat/data/psp_spc_l3'
#os.system('wget -nc --directory-prefix='+psp_data_path+' "ftps://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spc/l3/l3i/2019/*.cdf" ')
#wget -nc "ftps://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spc/l3/l3i/2019/*.cdf")




############################################# standard data update each day


######################################### spacecraft positions image
hp.plot_positions(datetime.datetime.utcnow(),position_path, 'HEEQ',now=True)



########################################### SDO images now
hd.get_sdo_realtime_image()



########################################### NOAA real time

print('download NOAA real time solar wind plasma and mag')
datestr=str(datetime.datetime.utcnow().strftime("%Y_%b_%d_%H_%M"))
print(datestr+' UTC')

plasma='http://services.swpc.noaa.gov/products/solar-wind/plasma-7-day.json'
mag='http://services.swpc.noaa.gov/products/solar-wind/mag-7-day.json'

try: urllib.request.urlretrieve(plasma, noaa_path+'plasma-7-day_'+datestr+'.json')
except urllib.error.URLError as e:
  print(' ', plasma,' ',e.reason)

try: urllib.request.urlretrieve(mag, noaa_path+'mag-7-day_'+datestr+'.json')
except urllib.error.URLError as e:
  print(' ', mag,' ',e.reason)
  
print()

filenoaa='noaa_rtsw_jan_2020_now.p'
if get_new_data: hd.save_noaa_rtsw_data(data_path,noaa_path,filenoaa)
[noaa,hnoaa]=pickle.load(open(data_path+filenoaa, "rb" ) ) 

start=noaa.time[-1]-datetime.timedelta(days=14)
end=datetime.datetime.utcnow() #noaa.time[-1]     
hp.plot_insitu_update(noaa, start, end,'NOAA_RTSW',plot_path,now=True)

start=noaa.time[-1]-datetime.timedelta(days=55)
end=datetime.datetime.utcnow() #noaa.time[-1]     
hp.plot_insitu_update(noaa, start, end,'NOAA_RTSW',plot_path,now2=True)

print()



#make quick new plot for ICMEs
# hd.save_noaa_rtsw_data(data_path,noaa_path,filenoaa)
# start=noaa.time[-1]-datetime.timedelta(days=3)
# end=datetime.datetime.utcnow() #noaa.time[-1]     
# hp.plot_insitu(noaa, start, end,'NOAA_RTSW','/home/cmoestl/pycode/heliocats')



########################################### STEREO-A
filesta="stereoa_2019_now_sceq_beacon.p" 
start=datetime.datetime(2019, 1, 1)
end=datetime.datetime.utcnow()
if get_new_data: hd.save_stereoa_beacon_data(data_path,filesta,start,end,sceq=True)
[sta,hsta]=pickle.load(open(data_path+filesta, "rb" ) ) 

start=sta.time[-1]-datetime.timedelta(days=14)
end=datetime.datetime.utcnow()     
hp.plot_insitu_update(sta, start, end,'STEREO-A_beacon',plot_path,now=True)




########################################### Wind

print('download Wind 2020 files without overwriting existing files into ',data_path)
wind_data_path='/nas/helio/data/heliosat/data/wind_mfi_k0'
os.system('wget -nc --directory-prefix='+wind_data_path+' "ftps://spdf.gsfc.nasa.gov/pub/data/wind/mfi/mfi_k0/2020/*.cdf"')
wind_data_path='/nas/helio/data/heliosat/data/wind_swe_h1'
os.system('wget -nc --directory-prefix='+wind_data_path+' "ftps://spdf.gsfc.nasa.gov/pub/data/wind/swe/swe_h1/2020/*.cdf"')


filewin="wind_2018_now_gse.p" 
start=datetime.datetime(2018, 1, 1)
end=datetime.datetime.utcnow()
if get_new_data: hd.save_wind_data(data_path,filewin,start,end,heeq=False)
[win,hwin]=pickle.load(open(data_path+filewin, "rb" ) )  

filewin="wind_2018_now_heeq.p" 
start=datetime.datetime(2018, 1, 1)
end=datetime.datetime.utcnow()
if get_new_data: hd.save_wind_data(data_path,filewin,start,end,heeq=True)

start=win.time[-1]-datetime.timedelta(days=100)
end=datetime.datetime.utcnow()         
hp.plot_insitu_update(win, start, end,'Wind',plot_path,now=True)




####################################### OMNI2
fileomni="omni_1963_now.p"
overwrite=1
if get_new_data: hd.save_omni_data(data_path,fileomni,overwrite)
[o,ho]=pickle.load(open(data_path+fileomni, "rb" ) )  

start=datetime.datetime.utcnow() -datetime.timedelta(days=365)
end=datetime.datetime.utcnow() 
hp.plot_insitu_update(o, start, end,'OMNI2',plot_path,now=True)


#for making quick plots
# plot_path='/home/cmoestl/pycode/heliocats/'
# fileomni="omni_1963_now.p"
# [o,ho]=pickle.load(open(data_path+fileomni, "rb" ) )  
# start=datetime.datetime.utcnow() -datetime.timedelta(days=365)
# end=datetime.datetime.utcnow() 
# hp.plot_insitu_update(o, start, end,'OMNI2',plot_path,now=True)


############### write header file for daily updates
text = open('/nas/helio/data/insitu_python/data_update_headers.txt', 'w')
text.write('Contains headers for the data files which are updated in real time.'+'\n \n')
text.write('File creation date:  '+datetime.datetime.utcnow().strftime("%Y-%b-%d %H:%M") +' \n \n')


text.write('NOAA real time solar wind: '+filenoaa+'\n \n'+ hnoaa+' \n \n')
text.write('load with: >> [noaa,hnoaa]=pickle.load(open("'+data_path+filenoaa+'", "rb"))') 
text.write(' \n \n \n \n')

text.write('STEREO-A beacon: '+filesta+'\n \n'+ hsta+' \n \n')
text.write('load with: >> [sta,hsta]=pickle.load(open("'+data_path+filesta+'", "rb"))') 
text.write(' \n \n \n \n')

text.write('Wind: '+filewin+'\n \n'+ hwin+' \n \n')
text.write('load with: >> [win,hwin]=pickle.load(open("'+data_path+filewin+'", "rb" ))') 
text.write(' \n \n \n \n')


text.write('OMNI2: '+fileomni+'\n \n'+ ho+' \n \n')
text.write('load with: >> [o,ho]=pickle.load(open("'+data_path+fileomni+'", "rb" ))') 
text.write(' \n \n \n \n')

text.close()
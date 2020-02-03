'''
https://github.com/cmoestl/heliocats  data_update.py

for updating data every day 



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



#for easier debugging - do not download and process data but do everything else
get_new_data=True

#~/miniconda/envs/helio/bin/python /home/cmoestl/pycode/heliocats/sc_positions_insitu_orbit1.py
#~/miniconda/envs/helio/bin/python /home/cmoestl/pycode/heliocats/sc_positions_insitu_orbit2.py


'''
filesta2="stereoa_2007_2009_beacon.p"
start=datetime.datetime(2007, 3, 20)
end=datetime.datetime(2009, 12, 31)
hd.save_stereoa_beacon_data(data_path, filesta2,start, end)
#[sta,hsta]=pickle.load(open(data_path+filesta2, "rb" ) ) 

sys.exit()


filesta2="stereoa_2007_2009_beacon.p"
start=datetime.datetime(2007, 3, 20)
end=datetime.datetime(2009, 12, 31)
hd.save_stereoa_beacon_data(data_path, filesta2,start, end)
#[sta,hsta]=pickle.load(open(data_path+filesta2, "rb" ) ) 

sys.exit()
'''

#add real time image of the solar system in data_update


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
print()




##################### standard data update each day


#NOAA
filenoaa='noaa_rtsw_jan_2020_now.p'
if get_new_data: hd.save_noaa_rtsw_data(data_path,noaa_path,filenoaa)
[noaa,hnoaa]=pickle.load(open(data_path+filenoaa, "rb" ) ) 


start=noaa.time[-1]-datetime.timedelta(days=14)
end=datetime.datetime.utcnow() #noaa.time[-1]     
hp.plot_insitu(noaa, start, end,'NOAA_RTSW',plot_path,now=True)

start=noaa.time[-1]-datetime.timedelta(days=32)
end=noaa.time[-1]     
hp.plot_insitu(noaa, start, end,'NOAA_RTSW',plot_path,now2=True)

#STEREO-A
filesta="sta_2018_now_beacon.p" 
start=datetime.datetime(2018, 1, 1)
end=datetime.datetime.utcnow()
if get_new_data: hd.save_stereoa_beacon_data(data_path,filesta,start,end)
[sta,hsta]=pickle.load(open(data_path+filesta, "rb" ) ) 

start=sta.time[-1]-datetime.timedelta(days=14)
end=sta.time[-1]     
hp.plot_insitu(sta, start, end,'STEREO-A_beacon',plot_path,now=True)


#Wind
filewin="wind_2018_now.p" 
start=datetime.datetime(2018, 1, 1)
end=datetime.datetime.utcnow()
if get_new_data: hd.save_wind_data(data_path,filewin,start,end)
[win,hwin]=pickle.load(open(data_path+filewin, "rb" ) )  

start=win.time[-1]-datetime.timedelta(days=100)
end=win.time[-1]     
hp.plot_insitu(win, start, end,'Wind',plot_path,now=True)



#OMNI2
'''
fileomni="omni_1963_now.p"
overwrite=0
if get_new_data: hd.save_omni_data(data_path,fileomni,overwrite)
'''
[o,ho]=pickle.load(open(data_path+fileomni, "rb" ) )  

'''
start=o.time[-1]-datetime.timedelta(days=26*7)
end=o.time[-1]     
hp.plot_insitu(o, start, end,'OMNI2',plot_path,now=True)
'''



#################### write header file for daily updates
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

sys.exit()









####################################

## long jobs


data_path='/nas/helio/data/insitu_python/'

'''
#################### PSP -> ok!
filepsp="psp_2018_2019_merged.p"
hd.save_psp_data(data_path, filepsp)
[psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) ) 

filepsp="psp_2018_2019_non_merged.p"
hd.save_psp_data_non_merged(data_path, filepsp)
[psp_orbit,psp_mag,psp_plasma,header_psp]=pickle.load(open(data_path+filepsp, "rb" ) )  

#####################
















'''


'''
filestb="stereob_2007_2014_beacon.p"
start=datetime.datetime(2007, 3, 20)
end=datetime.datetime(2014, 9, 27)
hd.save_stereob_beacon_data(data_path, filestb,start, end)
#[stb,hstb]=pickle.load(open(data_path+filestb, "rb" ) ) 

filesta2="stereoa_2007_2019_beacon.p"
start=datetime.datetime(2007, 3, 20)
end=datetime.datetime(2019, 12, 31)
hd.save_stereoa_beacon_data(data_path, filesta2,start, end)
#[sta,hsta]=pickle.load(open(data_path+filesta2, "rb" ) ) 


filewin="wind_2007_2019.p" 
start=datetime.datetime(2007, 1, 1)
end=datetime.datetime(2019, 12, 31)
hd.save_wind_data(data_path,filewin, start, end)
#[win,hwin]=pickle.load(open(data_path+filewin, "rb" ) )  



filewin2="wind_2007_2019.p" 
start=datetime.datetime(2007, 1, 1)
end=datetime.datetime(2019, 31, 12)
hd.save_wind_data(data_path,filewin2, start, end)
[win,hwin]=pickle.load(open(data_path+filewin2, "rb" ) )  


filesta2="stereoa_2007_2019_beacon.p"
start=datetime.datetime(2007, 1, 1)
end=datetime.datetime(2019, 31, 12)
hd.save_stereoa_beacon_data(data_path, filesta2,start, end)
[sta,hsta]=pickle.load(open(data_path+filesta2, "rb" ) ) 


sys.exit()



filesta2="stereoa_2018_2019.p"
hd.save_stereoa_science_data(data_path, filesta2)
[sta2,header]=pickle.load(open(data_path+filesta2, "rb" ) ) 



filemav='maven_2014_2018_removed_smoothed.p'
hd.MAVEN_smooth_orbit(data_path,filemav) 

sys.exit()


#filepsp="psp_2018_2019.p"
#hd.save_psp_data(data_path, filepsp)
#sys.exit()


sys.exit()

sys.exit()


    # save files from raw data if necessary for updates
    #hd.save_psp_data(filepsp)
    #hd.save_wind_data(filewin2)
    #hd.save_stereoa_data(filesta2)
    #hd.convert_MAVEN_mat_to_pickle() data from C. S. Wedlund
    # ADD BepiColombo  
    # ADD Solar Orbiter
    #sys.exit()


    #make a single helcats data file if necessary
    #hd.save_helcats_datacat()



filemav='maven_2014_2018.p'
hd.convert_MAVEN_mat_original(data_path,filemav) 
#[mav,hmav]=pickle.load(open(filemav, 'rb' ) )

filemav='maven_2014_2018_removed.p'
hd.convert_MAVEN_mat_removed(data_path,filemav) 
#[mav,hmav]=pickle.load(open(filemav, 'rb' ) )

filemav='maven_2014_2018_removed_smoothed.p'
hd.convert_MAVEN_mat_removed_orbit(data_path,filemav) 
#[mav,hmav]=pickle.load(open(filemav, 'rb' ) )




hd.save_psp_data(data_path, filepsp)

filepsp="psp_2018_2019.p"
data_path='/nas/helio/data/insitu_python/'
[psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) )  

sys.exit()

hd.save_wind_data(data_path,filewin2)
hd.save_stereoa_data(data_path, filesta2)



#filemav=data_path+'maven_2014_2018_removed.p'
#[mav,hmav]=pickle.load(open(filemav, 'rb' ) )

sys.exit()


hd.save_helcats_datacat(data_path,removed=True)
hd.save_helcats_datacat(data_path,removed=False)



file='data/helios.p'
hd.save_helios_data(file)
sys.exit()

filecas='data/cassini_1999_2000.p'
hd.save_cassini_data(filecas)
sys.exit()
'''


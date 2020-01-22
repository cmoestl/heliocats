'''
for updating data every day for Wind and STEREO-A
https://github.com/cmoestl/heliocats


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
import sys
import numpy as np
import datetime
import scipy.signal

from heliocats import data as hd
importlib.reload(hd) #reload again while debugging

from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging



data_path='/nas/helio/data/insitu_python/'


#real time image of the solar system in data_update



##################### standard data update

#STEREO-A
filesta="sta_2018_now_beacon.p" 
start=datetime.datetime(2018, 1, 1)
end=datetime.datetime.utcnow()
hd.save_stereoa_beacon_data(data_path,filesta,start,end)
[sta,hsta]=pickle.load(open(data_path+filesta, "rb" ) ) 


#add omni


'''
#wind - not correct MAG
filewin="wind_2018_now.p" 

start=datetime.datetime(2018, 1, 1)
end=datetime.datetime.utcnow()
hd.save_wind_data(data_path,filewin,start,end)
[win,hwin]=pickle.load(open(data_path+filewin, "rb" ) )  


'''
###########################################


## long jobs


data_path='/nas/helio/data/insitu_python/'


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


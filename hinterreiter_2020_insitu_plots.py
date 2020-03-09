'''
https://github.com/cmoestl/heliocats  hint20.py

plots for hinterreiter et al. 2020


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
from sunpy.time import parse_time
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


#for server use this so no plotting window pops up:
matplotlib.use('Agg')
 

data_path='/nas/helio/data/insitu_python/'
plot_path='/nas/helio/data/insitu_python/hinterreiter2020/'
noaa_path='/nas/helio/data/noaa_rtsw/'


file="wind_2007_2018_helcats.p" 
[win,hwin]=pickle.load(open(data_path+file, "rb" ) ) 


#final mine
mine=[
'2010-02-07 18:04',
'2010-03-23 22:29',
'2010-04-05 07:55',
'2010-04-11 12:20',
'2010-05-28 01:52',
'2010-10-30 09:15',
'2011-02-04 01:50',
'2011-02-18 00:48',
'2011-09-09 11:46',
'2012-01-24 14:36',
'2012-06-16 19:34',
'2012-07-14 17:38']



#helcats and moestl 2014 for 2 events
hm=[
'2010-02-07 18:04',
'2010-03-23 22:33',
'2010-04-05 07:55',
'2010-04-11 12:28',
'2010-05-28 02:23',
'2010-10-31 02:09',
'2011-02-04 01:55',
'2011-02-18 00:43',
'2011-09-09 11:46',
'2012-01-24 14:36',
'2012-06-16 09:07',
'2012-07-14 17:59']




#richardson cane 
rc=[
'2010-02-07 17:00',
'2010-04-05 08:26',
'2010-04-11 13:04',
'2010-05-28 02:58',
'2010-10-30 10:15',
'2011-02-04 13:00',
'2011-02-18 01:30',
'2011-09-09 12:42',
'2012-06-16 20:19',
'2012-07-14 18:09']

wind=[
'2010-02-07 18:04',
'2010-03-23 22:29',
'2010-04-05 07:55',
'2010-04-11 12:20',
'2010-05-28 01:55',
'2010-10-31 02:09',
'2011-02-18 00:49',
'2011-02-18 19:50',
'2012-06-16 09:03',
'2012-06-16 19:34',
'2012-07-14 17:39']


mine=parse_time(mine).datetime          
hm=parse_time(hm).datetime          
rc=parse_time(rc).datetime          
wind=parse_time(wind).datetime     


'''

from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging

#for measuring specific events
window=48
i=5
start=mine[i]-datetime.timedelta(hours=window)
end=mine[i]+datetime.timedelta(hours=window)
hp.plot_insitu_hint20(win, start, end,'Wind',plot_path,mine,rc,wind,lowres=True)
'''


for i in np.arange(len(mine)):

    window=48
    start=mine[i]-datetime.timedelta(hours=window)
    end=mine[i]+datetime.timedelta(hours=window)
    hp.plot_insitu_hint20(win, start, end,'Wind',plot_path,mine,rc,wind,lowres=True)

    window=24
    start=mine[i]-datetime.timedelta(hours=window)
    end=mine[i]+datetime.timedelta(hours=window)
    hp.plot_insitu_hint20(win, start, end,'Wind',plot_path,mine,rc,wind)

    window=3
    start=mine[i]-datetime.timedelta(hours=window)
    end=mine[i]+datetime.timedelta(hours=window)
    hp.plot_insitu_hint20(win, start, end,'Wind',plot_path,mine,rc,wind,highres=True)












































#script for full HIA movie

from scipy import stats
import scipy.io
from matplotlib import cm
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import sunpy.time
import time
import pickle
import seaborn as sns
import math
import urllib





### convert all files to mp4

for q in np.arange(2007,2015):

 year = str(q)
 

 for i in np.arange(1,13):

  month = '%02i' % i
  befehl='/Users/chris/movie/ffmpeg -i HI1B_full/'+year+month+'_hi1b_1_25_m5_stardiff.mpg   HI1B_full/'+year+month+'_hi1b_1_25_m5_stardiff.mp4 -y'
  os.system('/Users/chris/movie/ffmpeg -i HI1B_full/'+year+month+'_hi1b_1_25_m5_stardiff.mpg HI1B_full/'+year+month+'_hi1b_1_25_m5_stardiff.mp4 -y')



sys.exit()


#download all files



names=['January', 'February', 'March', 'April', 'May','June','July','August','September','October','November','December']

for q in np.arange(2008,2015):

 year = str(q)

 for i in np.arange(1,13):
    month = '%02i' % i
    
    url='https://www.ukssdc.ac.uk/solar/stereo/movies/MOVIES/'+year+'/'+month+'_'+names[i-1]+'/'+year+month+'_hi1b_1_25_m5_stardiff.mpg'
    urllib.request.urlretrieve(url, 'HI1B_full/'+year+month+'_hi1b_1_25_m5_stardiff.mpg')







sys.exit()



sys.exit()

for q in np.arange(2015,2018):

 year = str(q)
 

 for i in np.arange(1,13):

  month = '%02i' % i
  befehl='/Users/chris/movie/ffmpeg -i HI1B_full/'+year+month+'_hi1b_1_25_m5_stardiff.mpg   HI1B_full/'+year+month+'_hi1b_1_25_m5_stardiff.mp4 -y'
  os.system('/Users/chris/movie/ffmpeg -i HI1B_full/'+year+month+'_hi1b_1_25_m5_stardiff.mpg HI1B_full/'+year+month+'_hi1b_1_25_m5_stardiff.mp4 -y')





#/Users/chris/movie/ffmpeg -i 201204_hi1a_1_25_m5_stardiff.mpg 201204_hi1a_1_25_m5_stardiff.mp4


#/Users/chris/movie/ffmpeg -i 201205_hi1a_1_25_m5_stardiff.mpg 201205_hi1a_1_25_m5_stardiff.mp4



#/Users/chris/movie/ffmpeg -i https://www.ukssdc.ac.uk/solar/stereo/movies/MOVIES/2012/06_June/201206_hi1a_1_25_m5_stardiff.mpg  201206_hi1a_1_25_m5_stardiff.mp4
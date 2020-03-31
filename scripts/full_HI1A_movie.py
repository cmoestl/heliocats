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



#download all files



for q in np.arange(2006,2015):

 year = str(q)

 for i in np.arange(1,13):
    month = '%02i' % i


testfile = urllib.URLopener()
testfile.retrieve("http://randomsite.com/file.gz", "file.gz")







sys.exit()


### convert all files to mp4

for q in np.arange(2006,2015):

 year = str(q)
 

 for i in np.arange(1,13):

  month = '%02i' % i
  befehl='/Users/chris/movie/ffmpeg -i HIA1_full/'+year+month+'_hi1a_1_25_m5_stardiff.mpg   HIA1_full/'+year+month+'_hi1a_1_25_m5_stardiff.mp4 -y'
  os.system('/Users/chris/movie/ffmpeg -i HIA1_full/'+year+month+'_hi1a_1_25_m5_stardiff.mpg HIA1_full/'+year+month+'_hi1a_1_25_m5_stardiff.mp4 -y')





for q in np.arange(2015,2018):

 year = str(q)
 

 for i in np.arange(1,13):

  month = '%02i' % i
  befehl='/Users/chris/movie/ffmpeg -i HIA1_full/'+year+month+'_hi1a_1_25_m5_stardiff.mpg   HIA1_full/'+year+month+'_hi1a_1_25_m5_stardiff.mp4 -y'
  os.system('/Users/chris/movie/ffmpeg -i HIA1_full/'+year+month+'_hi1a_1_25_m5_stardiff.mpg HIA1_full/'+year+month+'_hi1a_1_25_m5_stardiff.mp4 -y')





#/Users/chris/movie/ffmpeg -i 201204_hi1a_1_25_m5_stardiff.mpg 201204_hi1a_1_25_m5_stardiff.mp4


#/Users/chris/movie/ffmpeg -i 201205_hi1a_1_25_m5_stardiff.mpg 201205_hi1a_1_25_m5_stardiff.mp4



#/Users/chris/movie/ffmpeg -i https://www.ukssdc.ac.uk/solar/stereo/movies/MOVIES/2012/06_June/201206_hi1a_1_25_m5_stardiff.mpg  201206_hi1a_1_25_m5_stardiff.mp4
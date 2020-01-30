#plot.py
#plot stuff for heliocats
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy
import copy
import matplotlib.dates as mdates
import matplotlib
import seaborn as sns
import datetime
import urllib
import json
import os
import pdb
import scipy.io
import pickle
import sys
import cdflib
import matplotlib.pyplot as plt
import heliosat
from numba import njit
from astropy.time import Time
import heliopy.spice as spice
import astropy

data_path='/nas/helio/data/insitu_python/'




####################################### 




def plot_insitu(sc, start, end, sc_label, path):
     '''
     sc = data
    
     '''
     sns.set_style("darkgrid")
     sns.set_context('paper')

     fig=plt.figure(1, figsize=(9,6), dpi=150)
     
     #plt.title(sc_label+' data, start time '+start.strftime("%Y-%b-%d %H:%M"))

     ax1 = plt.subplot(411) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
     
     plt.ylabel('B [nT]')
     plt.legend(loc=3)
     ax1.set_xlim(start,end)
     ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d') )
     plt.ylim((-20, 20))
     ax1.set_xticklabels([])

     plt.title(sc_label+' data, start time '+start.strftime("%Y-%b-%d %H:%M"))


     ax2 = plt.subplot(412) 

     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=0.7)

     plt.ylabel('V [km/s]')
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((250, 800))
     ax2.set_xticklabels([])


     ax3 = plt.subplot(413) 

     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)

     plt.ylabel('N [ccm-3]')
     ax3.set_xlim(start,end)
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((0, 50))
     ax3.set_xticklabels([])

     ax4 = plt.subplot(414) 

     ax4.plot_date(sc.time,sc.tp/1e6,'-k',label='Tp',linewidth=0.7)

     plt.ylabel('T [MK]')
     ax4.set_xlim(start,end)
     ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((0, 0.5))
 
     plt.tight_layout()
     #plt.show()

     plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     plt.savefig(plotfile)
     print('saved as ',plotfile)
         
     
     
     
     
     
     
     
     
##################################
 
''' 
plt.close() 
sns.set_style("darkgrid")
sns.set_context('paper')

 
data_path='/nas/helio/data/insitu_python/'  
plot_path='/nas/helio/data/insitu_python/plots/'
filewin="wind_2018_now.p" 


filewin='wind_2007_2018_helcats.p'
[win,hwin]=pickle.load(open(data_path+filewin, "rb" ) )  

start=datetime.datetime(2010, 10, 30)
end=datetime.datetime(2010, 11, 1)     

plot_insitu(win, start, end,'Wind',plot_path)
'''











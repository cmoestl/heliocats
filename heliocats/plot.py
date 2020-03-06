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




def plot_insitu(sc, start, end, sc_label, path, **kwargs):
     '''
     sc = data
    
     '''
     sns.set_style('darkgrid')
     sns.set_context('paper')

     fig=plt.figure(figsize=(9,6), dpi=150)
     
     #sharex means that zooming in works with all subplots
     ax1 = plt.subplot(411) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
     
     plt.ylabel('B [nT]')
     plt.legend(loc=3,ncol=4,fontsize=8)
     ax1.set_xlim(start,end)
     ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d') )
     plt.ylim((-20, 20))
     #ax1.set_xticklabels([]) does not work with sharex
     #plt.setp(ax1.get_xticklabels(), fontsize=6)
     plt.setp(ax1.get_xticklabels(), visible=False)

     plt.title(sc_label+' data, start: '+start.strftime("%Y-%b-%d %H:%M")+'  end: '+end.strftime("%Y-%b-%d %H:%M"))


     ax2 = plt.subplot(412,sharex=ax1) 
     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=0.7)

     plt.ylabel('V [km/s]')
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((250, 800))
     #ax2.set_xticklabels([])
     plt.setp(ax2.get_xticklabels(), visible=False)



     ax3 = plt.subplot(413,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)

     plt.ylabel('N [ccm-3]')
     ax3.set_xlim(start,end)
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((0, 50))
     #ax3.set_xticklabels([])
     plt.setp(ax3.get_xticklabels(), visible=False)


     ax4 = plt.subplot(414,sharex=ax1) 
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
   

     #if now exists as keyword, save as the file with just now in filename:     
     if 'now' in kwargs:
        plotfile=path+sc_label+'_now.png'
        plt.savefig(plotfile)
        print('saved as ',plotfile)

     #if now2 exists as keyword, save as the file with just now in filename:     
     if 'now2' in kwargs:
        plotfile=path+sc_label+'_now2.png'
        plt.savefig(plotfile)
        print('saved as ',plotfile)
     
     
     
     
     


def plot_insitu_hint20(sc, start, end, sc_label, path, e1,e2,e3,**kwargs):
     '''
     '''
     
     sns.set_style('darkgrid')
     sns.set_context('paper')

     fig=plt.figure(figsize=(9,6), dpi=150)
     
     
     #cut out data
     startindex=np.where(start>sc.time)[0][-1]
     endindex=np.where(end>sc.time)[0][-1]   
          
     sc=sc[startindex:endindex]
         
     
     #sharex means that zooming in works with all subplots
     ax1 = plt.subplot(311) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)

     ax1.plot_date([e1,e1],[-100,100],'-k',linewidth=1)     
     ax1.plot_date([e2,e2],[-100,100],'--k',linewidth=1)
     ax1.plot_date([e3,e3],[-100,100],'-b',linewidth=1)

     
     plt.ylabel('B [nT]')
     plt.legend(loc=3,ncol=4,fontsize=8)
     ax1.set_xlim(start,end)
     ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d') )
     plt.ylim((-np.nanmax(sc.bt)-5,np.nanmax(sc.bt)+5))
     #ax1.set_xticklabels([]) does not work with sharex
     #plt.setp(ax1.get_xticklabels(), fontsize=6)
     plt.setp(ax1.get_xticklabels(), visible=False)



     plt.title(sc_label+' data, start: '+start.strftime("%Y-%b-%d %H:%M")+'  end: '+end.strftime("%Y-%b-%d %H:%M"))

     ax2 = plt.subplot(312,sharex=ax1) 
     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=0.7)

     ax2.plot_date([e1,e1],[0,5000],'-k',linewidth=1)     
     ax2.plot_date([e2,e2],[0,5000],'--k',linewidth=1)
     ax2.plot_date([e3,e3],[0,5000],'-b',linewidth=1)
     

     plt.ylabel('V [km/s]')
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((np.nanmin(sc.vt)-50,np.nanmax(sc.vt)+100))
     #ax2.set_xticklabels([])
     plt.setp(ax2.get_xticklabels(), visible=False)



     ax3 = plt.subplot(313,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)

     ax3.plot_date([e1,e1],[0,500],'-k',linewidth=1)     
     ax3.plot_date([e2,e2],[0,500],'--k',linewidth=1)
     ax3.plot_date([e3,e3],[0,500],'-b',linewidth=1)


     plt.ylabel('N [ccm-3]')
     ax3.set_xlim(start,end)
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((0,np.nanmax(sc.np)+5))
     #ax3.set_xticklabels([])
     #plt.setp(ax3.get_xticklabels(), visible=False)
     
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     
     plt.tight_layout()
     #plt.show()


     '''
     ax4 = plt.subplot(414,sharex=ax1) 
     ax4.plot_date(sc.time,sc.tp/1e6,'-k',linewidth=0.7)

     ax4.plot_date([e1,e1],[0,10],'-k',linewidth=1)     
     ax4.plot_date([e2,e2],[0,10],'--k',linewidth=1)


     plt.ylabel('T [MK]')
     ax4.set_xlim(start,end)
     plt.ylim((0,np.nanmax(sc.tp)+0.2))

     '''



    

     #if highres exists as keyword, save as the file with just now in filename:     
     if 'highres' in kwargs:
        plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'_highres.png'
        plt.savefig(plotfile)
        print('saved as ',plotfile)
        
     #if highres exists as keyword, save as the file with just now in filename:     
     elif 'lowres' in kwargs:
        plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'_lowres.png'
        plt.savefig(plotfile)
        print('saved as ',plotfile)   
        
     else:   
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











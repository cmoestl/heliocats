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
import seaborn as sns

#import 

from heliocats import data as hd
importlib.reload(hd) #reload again while debugging

from heliocats import plot as hp
importlib.reload(hp) #reload again while debugging


#for server use this so no plotting window pops up:
#matplotlib.use('Agg')
 
 
 
 



def plot_insitu(sc, start, end, sc_label, path, lines, **kwargs):
     '''
     sc = data
    
     '''
     
     start=parse_time(start).datetime
     end=parse_time(end).datetime
     #print(start)
     #print(end)
     
     sns.set_style('darkgrid')
     sns.set_context('paper')

     fig=plt.figure(figsize=(9,7), dpi=150)
     
     #sharex means that zooming in works with all subplots
     ax1 = plt.subplot(411) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=1)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=1)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=1)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=1)
     
     plt.ylabel('B [nT]')
     plt.legend(loc=3,ncol=4,fontsize=8)
     ax1.set_xlim(start,end)
     ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )
     plt.ylim((-110, 110))
     #ax1.set_xticklabels([]) does not work with sharex
     #plt.setp(ax1.get_xticklabels(), fontsize=6)
     plt.setp(ax1.get_xticklabels(), visible=False)
     
     ax1.plot_date([lines[0],lines[0]],[-100,100],'-k',linewidth=1) 
     ax1.plot_date([lines[1],lines[1]],[-100,100],'-k',linewidth=1) 
     ax1.plot_date([lines[2],lines[2]],[-100,100],'-k',linewidth=1) 


     plt.title(sc_label+' data, start: '+start.strftime("%Y-%b-%d %H:%M")+'  end: '+end.strftime("%Y-%b-%d %H:%M"))
     


     ax2 = plt.subplot(412,sharex=ax1) 
     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=0.7)

     plt.ylabel('V [km/s]')
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((300, 500))
     #ax2.set_xticklabels([])
     plt.setp(ax2.get_xticklabels(), visible=False)

     ax2.plot_date([lines[0],lines[0]],[0,5000],'-k',linewidth=1) 
     ax2.plot_date([lines[1],lines[1]],[0,5000],'-k',linewidth=1) 
     ax2.plot_date([lines[2],lines[2]],[0,5000],'-k',linewidth=1) 


     ax3 = plt.subplot(413,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)

     plt.ylabel('N [ccm-3]')
     ax3.set_xlim(start,end)
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((0, 330))
     #ax3.set_xticklabels([])
     plt.setp(ax3.get_xticklabels(), visible=False)

     ax3.plot_date([lines[0],lines[0]],[0,5000],'-k',linewidth=1) 
     ax3.plot_date([lines[1],lines[1]],[0,5000],'-k',linewidth=1) 
     ax3.plot_date([lines[2],lines[2]],[0,5000],'-k',linewidth=1) 


     ax4 = plt.subplot(414,sharex=ax1) 
     ax4.plot_date(sc.time,sc.tp/1e6,'-k',label='Tp',linewidth=0.7)

     plt.ylabel('T [MK]')
     ax4.set_xlim(start,end)
     ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((0, 0.3))
     
     
     ax4.plot_date([lines[0],lines[0]],[0,5000],'-k',linewidth=1) 
     ax4.plot_date([lines[1],lines[1]],[0,5000],'-k',linewidth=1) 
     ax4.plot_date([lines[2],lines[2]],[0,5000],'-k',linewidth=1) 

     
     plt.tight_layout()
     #plt.show()

     plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     plt.savefig(plotfile)
     print('saved as ',plotfile)
     
     plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.eps'
     plt.savefig(plotfile)
     print('saved as ',plotfile)
     
     
     
     
     
file='icmecat/HELCATS_ICMECAT_v20.p'
ic=pickle.load( open(file, 'rb'))
   
event=704   
s1=parse_time(ic.icme_start_time.loc[event]).datetime     

s1=parse_time('2018-11-11T22:26Z').datetime     
s2=parse_time(ic.mo_start_time.loc[event]).datetime     
s3=parse_time(ic.mo_end_time.loc[event]).datetime     
     
     
data_path='/nas/helio/data/insitu_python/'

file="psp_2018_2019.p" 
[psp,hpsp]=pickle.load(open(data_path+file, "rb" ) ) 


#final mine


start=parse_time('2018-11-11 16:00').datetime
end=parse_time('2018-11-12 12:00').datetime

plot_path='results/'

plot_insitu(psp, start, end,'Parker Solar Probe',plot_path,[s1,s2,s3])












































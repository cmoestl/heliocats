#plots.py
#plotting routines for heliocats
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy
import copy
import matplotlib.dates as mdates
import matplotlib.image as mpimg
from sunpy.time import parse_time
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
from numba import njit
from astropy.time import Time
import astropy
import importlib

#import heliopy.data.spice as spicedata
#import heliopy.spice as spice


#Plotly imports
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.io as pio
from plotly.offline import iplot, init_notebook_mode
import plotly.express as px



#import heliosat   #not compatible with astrospice, problems with spiceypy in astrospice.generate
#from config import data_path

from heliocats import data as hd
importlib.reload(hd) #reload again while debugging


from heliocats import cats as hc
importlib.reload(hc) #reload again while debugging


'''
MIT LICENSE
Copyright 2020-2023, Christian Moestl, Rachel L. Bailey 
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


    


####################################### 





def plot_noaa_xray(xrayfile,xrayfile2,plot_path): 

    
    #use try except if not available
    
    [xdc1,xdc2]=pickle.load(open(xrayfile, 'rb'))
    
    [xd2c1,xd2c2]=pickle.load(open(xrayfile2, 'rb'))
    
    #remove data points below threshold due to eclipses of the GOES satellites  

    #upper plot
    
    threshold=1*1e-7
    
    rem=np.where(xdc1.flux < threshold)[0]    
    for i in rem:
        xdc1[i-2:i+2]=np.nan
    
    rem=np.where(xd2c1.flux < threshold)[0]    
    for i in rem:
        xd2c1[i-2:i+2]=np.nan

        
    #lower plot     
        
    #rem=np.where(xdc2.flux < 5*1e-9)[0]    
    #for i in rem:
    #    xdc2[i-2:i+2]=np.nan
    
    #rem=np.where(xd2c2.flux < 4*1e-9)[0]    
    #for i in rem:
    #    xd2c2[i-2:i+2]=np.nan
        
    
    sns.set_style('darkgrid')
    sns.set_context('paper')

    fig=plt.figure(3,figsize=(12,6),dpi=100)
    ax=plt.subplot(111)
    ax.set_yscale('log')    
    #ax.set_ylim(1e-9,1e-2)
    ax.set_ylim(1e-8,1e-2)
     
    
    print(xdc1.time)
    print(xdc2.time)
    
    
    ax.plot(xdc1.time,xdc1.flux,'-r',label='GOES-16 0.1-0.8nm')      
    #ax.plot(xdc2.time,xdc2.flux,'-b') 
    
    
    ax.plot(xd2c1.time,xd2c1.flux,color='darkorange', label='GOES-18 0.1-0.8nm')      
    #ax.plot(xd2c2.time,xd2c2.flux,color='purple')
    
    threshold1=1e-4   #for X1 Watts per m^2
    threshold2=5*1e-4 #for X5
    threshold3=1e-3   #for X10
    
    ax.axhline(y=5*1e-5, color='grey', linestyle='--', linewidth=0.5)
    ax.axhline(y=threshold1, color='yellowgreen', linestyle='--',label='X1', linewidth=1.5)
    ax.axhline(y=threshold2, color='orange', linestyle='--',label='X5',linewidth=1.5)
    ax.axhline(y=threshold3, color='red', linestyle='--',label='X10',linewidth=1.5)

    ax.axhline(y=threshold, color='grey', linestyle='--',alpha=0.5)

    
    ax.annotate('NOAA SWPC M5',xy=(datetime.datetime.utcnow()-datetime.timedelta(days=6),3*1e-5),xycoords='data',fontsize=12,ha='left',alpha=0.5)
    
    ax.axvline(x=datetime.datetime.utcnow(), color='k', linestyle='--',alpha=0.5, linewidth=1.0)



    plt.title('GOES X-Ray flux from NOAA',fontsize=16)

    fsize=12
    plt.legend(loc=3,fontsize=13, ncol=5)
    plt.figtext(0.02,0.01,'Austrian Space Weather Office   GeoSphere Austria', color='black', ha='left',fontsize=11, style='italic')
    plt.figtext(0.98,0.01,'helioforecast.space', color='black', ha='right',fontsize=11, style='italic')
    plt.figtext(0.09,0.95,'last update: '+str(datetime.datetime.utcnow())[0:16]+ ' UTC', ha='left', fontsize=10)


    ax.xaxis_date()
    ax.xaxis.set_major_locator(mdates.DayLocator())
    myformat = mdates.DateFormatter('%b %d')
    ax.xaxis.set_major_formatter(myformat)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    ax.set_xlabel(str(datetime.datetime.utcnow().year)+' UTC',fontsize=15)
    ax.set_xlim(datetime.datetime.utcnow()-datetime.timedelta(days=6),datetime.datetime.utcnow()+datetime.timedelta(days=0.5))

    ax.set_ylabel('Watts m$^{-2}$',fontsize=15)


    logo = plt.imread('logo/GSA_Basislogo_Positiv_RGB_XXS.png')
    newax = fig.add_axes([0.86,0.91,0.08,0.08], anchor='NE', zorder=1)
    newax.imshow(logo)
    newax.axis('off')

    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.tick_params(which='both', bottom=True, color='gray')
    

    plt.tight_layout()

    ax.annotate('X10',xy=(datetime.datetime.utcnow()+datetime.timedelta(days=0.15),1.5*1e-3),xycoords='data',fontsize=15,ha='left')
    ax.annotate('X',xy=(datetime.datetime.utcnow()+datetime.timedelta(days=0.15),1.5*1e-4),xycoords='data',fontsize=15,ha='left')
    ax.annotate('M',xy=(datetime.datetime.utcnow()+datetime.timedelta(days=0.15),1.5*1e-5),xycoords='data',fontsize=15,ha='left')
    ax.annotate('C',xy=(datetime.datetime.utcnow()+datetime.timedelta(days=0.15),1.5*1e-6),xycoords='data',fontsize=15,ha='left')
    ax.annotate('B',xy=(datetime.datetime.utcnow()+datetime.timedelta(days=0.15),1.5*1e-7),xycoords='data',fontsize=15,ha='left')
    ax.annotate('A',xy=(datetime.datetime.utcnow()+datetime.timedelta(days=0.15),1.5*1e-8),xycoords='data',fontsize=15,ha='left')

    
    plotfile=plot_path+'latest_xray.jpg'
    plt.savefig(plotfile, dpi=200)
    print('saved as ',plotfile)
        
    ## make plotly html plot
    
    nrows=1
    fig = make_subplots(rows=nrows, cols=1, shared_xaxes=True)

    fig.add_trace(go.Scatter(x=xdc1.time, y=xdc1.flux, name='GOES 16 long', line_color='red' ) )
    
    #fig.add_trace(go.Scatter(x=xdc2.time, y=xdc2.flux, name='GOEs 16 short',line_color='blue') )
    
    fig.add_trace(go.Scatter(x=xd2c1.time, y=xd2c1.flux, name='GOES 18 long', line_color='orange' ) )
    #fig.add_trace(go.Scatter(x=xd2c2.time, y=xd2c2.flux, name='GOEs 18 short',line_color='purple') )


    fig.update_layout(title='GOES Xrays', font=dict(size=20))
    fig.update_layout(xaxis=dict(range=[datetime.datetime.utcnow()-datetime.timedelta(days=6),datetime.datetime.utcnow()+datetime.timedelta(days=0.5)]) )


    fig.update_layout(
        xaxis=dict(
            title=dict(
                text="time",
                font=dict(size=20)  # Adjust the font size as needed
            )
        ),
        yaxis=dict(
            title=dict(
                text="Watts m-2",
                font=dict(size=20)  # Adjust the font size as needed
            )
        )
    )
    # Set the x-axis and y-axis type to log scale
    fig.update_layout(yaxis_type='log')
    fig.update_layout(yaxis_tickvals=[1e-9, 1e-8, 1e-7,1e-6,1e-5,1e-4,1e-3,1e-2], yaxis_ticktext=["10-^9", "10-^8", "10-^7", "10-^6", "10-^5", "10-^4", "10-^3", "10-^2"])
    #fig.update_layout(yaxis=dict(range=[1e-9,1e-2]) )



    plotfile=plot_path+'latest_xray.html'
    fig.write_html(plotfile)
    print('saved as ',plotfile)

    












def data_overview_plot(data,filename):


    sns.set_style('darkgrid')

    plt.figure(figsize=(15,15),dpi=100)


    ax1 = plt.subplot(321) 
    ax1.set_title('Btxyz')
    ax1.set_ylabel('nT')
    ax1.plot(data.time,data.bt,'-k',linewidth=0.5)
    ax1.plot(data.time,data.bx,'-r',linewidth=0.4)
    ax1.plot(data.time,data.by,'-g',linewidth=0.4)
    ax1.plot(data.time,data.bz,'-b',linewidth=0.4)

    ax2 = plt.subplot(322) 
    ax2.set_title('Vtxyz')
    ax2.set_ylabel('km/s')
    ax2.plot(data.time,data.vt,'-k',linewidth=0.5)
    ax2.plot(data.time,data.vx,'-r',linewidth=0.4)
    ax2.plot(data.time,data.vy,'-g',linewidth=0.4)
    ax2.plot(data.time,data.vz,'-b',linewidth=0.4)

    ax3 = plt.subplot(323) 
    ax3.set_title('Tp')
    ax3.set_ylabel('K')
    ax3.plot(data.time,data.tp,'-k',linewidth=0.5)

    ax4 = plt.subplot(324) 
    ax4.set_title('Np')
    ax4.set_ylabel('ccm-3')
    ax4.plot(data.time,data.np,'-k',linewidth=0.5)


    ax5 = plt.subplot(325) 
    ax5.set_title('XYZ HEEQ')
    ax5.set_ylabel('km')
    ax5.plot(data.time,data.x,'-r',linewidth=1,label='x')
    ax5.plot(data.time,data.y,'-g',linewidth=1,label='y')
    ax5.plot(data.time,data.z,'-b',linewidth=1,label='z')
    plt.legend(fontsize=15,loc=2)


    ax6 = plt.subplot(326) 
    ax6.set_title('R lon lat HEEQ')
    ax6.set_ylabel('AU,degree')
    ax6.plot(data.time,data.r*10,'-r',linewidth=1,label='R x 10')
    ax6.plot(data.time,data.lon,'-g',linewidth=1,label='lon')
    ax6.plot(data.time,data.lat,'-b',linewidth=1,label='lat')
    plt.legend(fontsize=15,loc=2)

    plt.tight_layout()

    plt.savefig(filename+'.png',dpi=150)
    #plt.savefig(filename+'.pdf',dpi=150)













def plot_insitu_update(sc, start, end, sc_label, path, **kwargs):
     '''
     sc = data
    
     '''
     sns.set_style('darkgrid')
     sns.set_context('paper')
     
     fsize=10

     fig=plt.figure(figsize=(9,6), dpi=150)
     
     #sharex means that zooming in works with all subplots
     ax1 = plt.subplot(411) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
     
     plt.ylabel('B [nT]',fontsize=fsize)
     plt.legend(loc=3,ncol=4,fontsize=fsize-2)
     ax1.set_xlim(start,end)
     ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d') )
     plt.ylim((-40, 40))
     #ax1.set_xticklabels([]) does not work with sharex
     #plt.setp(ax1.get_xticklabels(), fontsize=6)
     plt.setp(ax1.get_xticklabels(), visible=False)

     plt.title(sc_label+' data       start: '+start.strftime("%Y-%b-%d")+'       end: '+end.strftime("%Y-%b-%d"),fontsize=fsize)
     #plt.title(sc_label+' data       start: '+start.strftime("%Y-%b-%d %H:%M")+'       end: '+end.strftime("%Y-%b-%d %H:%M"),fontsize=fsize)


     ax2 = plt.subplot(412,sharex=ax1) 
     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=0.7)
     plt.ylabel('V [km/s]',fontsize=fsize)
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((250, 900))
     #ax2.set_xticklabels([])
     plt.setp(ax2.get_xticklabels(), visible=False)


     ax3 = plt.subplot(413,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)
     plt.ylabel('N [ccm-3]',fontsize=fsize)
     ax3.set_xlim(start,end)
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((0, 50))
     #ax3.set_xticklabels([])
     plt.setp(ax3.get_xticklabels(), visible=False)


     ax4 = plt.subplot(414,sharex=ax1) 
     ax4.plot_date(sc.time,sc.tp/1e6,'-k',label='Tp',linewidth=0.7)
     plt.ylabel('T [MK]',fontsize=fsize)
     ax4.set_xlim(start,end)
        
     ax4.xaxis.set_major_locator(mdates.MonthLocator())
        
     ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d') )
    
    
    
    
     plt.ylim((0, 1.0))
     
     plt.figtext(0.01,0.01,'Austrian Space Weather Office   GeoSphere Austria', color='black', ha='left',fontsize=fsize-4, style='italic')

     plt.figtext(0.99,0.01,'helioforecast.space', color='black', ha='right',fontsize=fsize-4, style='italic')

     
     plt.tight_layout()
     #plt.show()

     plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     plt.savefig(plotfile)
     print('saved as ',plotfile)

     plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.pdf'
     plt.savefig(plotfile)
     print('saved as ',plotfile)

    

     #if now exists as keyword, save as the file with just now in filename:     
     if 'now' in kwargs:
        plotfile=path+sc_label+'_now.png'
        plt.savefig(plotfile)
        print('saved as ',plotfile)

        plotfile=path+sc_label+'_now.pdf'
        plt.savefig(plotfile)
        print('saved as ',plotfile)

        
     #if now2 exists as keyword, save as the file with just now2 in filename:     
     if 'now2' in kwargs:
        plotfile=path+sc_label+'_now2.png'
        plt.savefig(plotfile)
        print('saved as ',plotfile)
        
        
        
        
def plot_insitu_update_stereoa_noaa(sc1in, sc2in, start, end, sc_label, path, **kwargs):
    '''
    sc = data

    '''

    #cut out data starting with start time that will be plotted, for better scaling
    startind=np.where(sc1in.time > start)[0][0]  
    sc1=sc1in[startind:]    

    startind=np.where(sc2in.time > start)[0][0]  
    sc2=sc2in[startind:]    


    sns.set_style('darkgrid')
    sns.set_context('paper')

    fsize=10

    #take maximum from both arrays
    bscale=np.nanmax(np.hstack((sc1.bt, sc2.bt))) +10
    vscale=np.nanmax(np.hstack((sc1.vt, sc2.vt)))+50


    fig=plt.figure(figsize=(9,6), dpi=150)

    #plt.suptitle(sc_label+' data       start: '+start.strftime("%Y-%b-%d %H:%M UT")+'       end: '+end.strftime("%Y-%b-%d %H:%M UT"),fontsize=fsize+2)

    plt.suptitle('NOAA RTSW and STEREO-A beacon data', fontsize=fsize+2)

    ax1 = plt.subplot(311) 

    ax1.plot_date(sc1.time,sc1.bx,'-r',label='Bx',linewidth=0.8)
    ax1.plot_date(sc1.time,sc1.by,'-g',label='By',linewidth=0.8)
    ax1.plot_date(sc1.time,sc1.bz,'-b',label='Bz',linewidth=0.8)
    ax1.plot_date(sc1.time,sc1.bt,'-k',label='Btotal',lw=0.8)
    plt.ylabel('B [nT] GSM',fontsize=fsize)
    plt.legend(loc=3,ncol=4,fontsize=fsize-2)
    ax1.set_xlim(start,end+datetime.timedelta(days=0.5))
    ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %Hh') )
    ax1.tick_params(which='both', bottom=True, color='gray')
    ax1.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
    ax1.set_ylim((-bscale,bscale))
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.text(sc2.time[40], bscale, 'NOAA L1 real time solar wind', fontsize=fsize, color='black')
    pos_string='R '+str(np.round(sc1.r[-1],3))+' AU   '+ 'lon '+str(np.round(sc1.lon[-1],2))+'  ' + 'lat '+str(np.round(sc1.lat[-1],2))+'' 
    ax1.text(sc1.time[-1350], bscale, pos_string, fontsize=fsize, color='black')

    #ax1.text(sc1.time[-600], bscale-5, 'east   north', fontsize=fsize-2, color='black')
    #ax1.text(sc1.time[-600], -bscale+3, 'west   south', fontsize=fsize-2, color='black')


    ax1.text(0.88,0.89,'east', fontsize=fsize-2, color='green',ha='center', va='center', transform=fig.transFigure)
    ax1.text(0.88,0.68,'west', fontsize=fsize-2, color='green',ha='center', va='center', transform=fig.transFigure)

    ax1.text(0.92,0.89,'north', fontsize=fsize-2, color='blue',ha='center', va='center', transform=fig.transFigure)
    ax1.text(0.92,0.68,'south', fontsize=fsize-2, color='blue',ha='center', va='center', transform=fig.transFigure)



    ax2 = plt.subplot(312,sharex=ax1)  
    ax2.plot_date(sc2.time,sc2.bx,'-r',label='Bx',linewidth=0.8)
    ax2.plot_date(sc2.time,sc2.by,'-g',label='By',linewidth=0.8)
    ax2.plot_date(sc2.time,sc2.bz,'-b',label='Bz',linewidth=0.8)
    ax2.plot_date(sc2.time,sc2.bt,'-k',label='Btotal',lw=0.8)
    plt.ylabel('B [nT] GSM',fontsize=fsize)
    plt.legend(loc=3,ncol=4,fontsize=fsize-2)
    ax2.set_xlim(start,end+datetime.timedelta(days=0.5))    
    ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %Hh') )
    ax2.tick_params(which='both', bottom=True, color='gray')
    ax2.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax2.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))    
    ax2.set_ylim((-bscale,bscale))
    ax2.text(sc2.time[40], bscale, 'STEREO-A beacon data', fontsize=fsize, color='black')
    pos_string='R '+str(np.round(sc2.r[-1],3))+' AU   '+ 'lon '+str(np.round(sc2.lon[-1],2))+'  ' + 'lat '+str(np.round(sc2.lat[-1],2))+''

    ax2.text(sc2.time[-1350], bscale, pos_string, fontsize=fsize, color='black')

    #ax2.text(sc1.time[-600], bscale-5, 'east   north', fontsize=fsize-2, color='black')
    #ax2.text(sc1.time[-600], -bscale+3, 'west   south', fontsize=fsize-2, color='black')


    plt.setp(ax2.get_xticklabels(), visible=False)



    ax3 = plt.subplot(313,sharex=ax1)  
    ax3.plot_date(sc1.time,sc1.vt,'-k',label='NOAA RTSW',linewidth=1)
    ax3.plot_date(sc2.time,sc2.vt,'-r',label='STEREO-A beacon',linewidth=1)
    plt.ylabel('V [km s$^{-1}$]',fontsize=fsize)    
    ax3.set_xlim(start,end+datetime.timedelta(days=0.5))
    ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %Hh') )    
    ax3.tick_params(which='both', bottom=True, color='gray')
    ax3.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax3.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(100))
    #set after maximum in both   
    ax3.set_ylim(250, vscale)
    ax3.text(sc1.time[40], vscale, 'NOAA L1', fontsize=fsize, color='black')
    ax3.text(sc1.time[650], vscale, 'STEREO-A', fontsize=fsize, color='red')



    ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %Hh') )

    ax1.axvline(x=datetime.datetime.utcnow(), color='k', linestyle='--',alpha=0.2, linewidth=1.0)
    ax2.axvline(x=datetime.datetime.utcnow(), color='k', linestyle='--',alpha=0.2, linewidth=1.0)
    ax3.axvline(x=datetime.datetime.utcnow(), color='k', linestyle='--',alpha=0.2, linewidth=1.0)


    plt.figtext(0.03,0.01,'Austrian Space Weather Office   GeoSphere Austria', color='black', ha='left',fontsize=fsize-4, style='italic')
    plt.figtext(0.97,0.01,'helioforecast.space', color='black', ha='right',fontsize=fsize-4, style='italic')
    plt.figtext(0.085,0.94,'last update: '+str(datetime.datetime.utcnow())[0:16]+ ' UTC', ha='left', fontsize=8)


    logo = plt.imread('logo/GSA_Basislogo_Positiv_RGB_XXS.png')
    newax = fig.add_axes([0.85,0.90,0.08,0.08], anchor='NE', zorder=1)
    newax.imshow(logo)
    newax.axis('off')
    
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


    #if now2 exists as keyword, save as the file with just now2 in filename:     
    if 'now2' in kwargs:
        plotfile=path+sc_label+'_now2.png'
        plt.savefig(plotfile)
        print('saved as ',plotfile)

        plotfile=path+sc_label+'_now2.pdf'
        plt.savefig(plotfile)
        print('saved as ',plotfile)





    ## make plotly html plot
    
    
    #cut out data starting with start time that will be plotted, for better scaling
    
    days_going_back=60
    startind=np.where(sc1in.time > end-datetime.timedelta(days=days_going_back))[0][0]  
    sc1=sc1in[startind:]    

    startind=np.where(sc2in.time > end-datetime.timedelta(days=days_going_back))[0][0]  
    sc2=sc2in[startind:]    
    

    nrows=3
    fig = make_subplots(rows=nrows, cols=1, shared_xaxes=True,row_heights=[0.35,0.35, 0.3])

    fig.add_trace(go.Scatter(x=sc1.time, y=sc1.bx, name='Bx',line_color='red'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sc1.time, y=sc1.by, name='By',line_color='green'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sc1.time, y=sc1.bz, name='Bz',line_color='blue'), row=1, col=1)
    fig.add_trace(go.Scatter(x=sc1.time, y=sc1.bt, name='Bt',line_color='black'), row=1, col=1)
    fig.update_yaxes(title_text="L1 B [nT] GSM", row=1, col=1)

    fig.add_trace(go.Scatter(x=sc2.time, y=sc2.bx, name='Bx',line_color='red'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sc2.time, y=sc2.by, name='By',line_color='green'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sc2.time, y=sc2.bz, name='Bz',line_color='blue'), row=2, col=1)
    fig.add_trace(go.Scatter(x=sc2.time, y=sc2.bt, name='Bt',line_color='black'), row=2, col=1)
    fig.update_yaxes(title_text="STEREO-A B [nT] GSM", row=2, col=1)


    fig.add_trace(go.Scatter(x=sc1.time, y=sc1.vt, name='Vt L1',line_color='black'), row=3, col=1)
    fig.add_trace(go.Scatter(x=sc2.time, y=sc2.vt, name='Vt STEREO-A',line_color='red'), row=3, col=1)
    fig.update_yaxes(title_text="V [km/s]", row=3, col=1)


    fig.update_layout(title='NOAA L1 and STEREO-A solar wind', font=dict(size=20))
    fig.update_layout(xaxis=dict(range=[datetime.datetime.utcnow()-datetime.timedelta(days=10),datetime.datetime.utcnow()+datetime.timedelta(days=0.5)]) )


    fig.update_layout(
        xaxis=dict(
            title=dict(
                font=dict(size=20)  # Adjust the font size as needed
            )
        ),
        yaxis=dict(
            title=dict(                
                font=dict(size=20)  # Adjust the font size as needed
            )
        )
    )

    fig.add_annotation(x=0.91, y=1.05, text="East", xref="paper", yref="paper", showarrow=False, font=dict(color='green')  )
    fig.add_annotation(x=0.91, y=0.72, text="West", xref="paper", yref="paper", showarrow=False,  font=dict(color='green')  )
    
    fig.add_annotation(x=0.99, y=1.05, text="North", xref="paper", yref="paper", showarrow=False, font=dict(color='blue')  )
    fig.add_annotation(x=0.99, y=0.72, text="South", xref="paper", yref="paper", showarrow=False,  font=dict(color='blue')  )
    
    
    fig.add_annotation(x=0.82, y=1.05, text="Toward", xref="paper", yref="paper", showarrow=False, font=dict(color='red')  )
    fig.add_annotation(x=0.82, y=0.72, text="Away", xref="paper", yref="paper", showarrow=False,  font=dict(color='red')  )
    
    fig.add_annotation(x=1.1, y=-0.05, text="toward to or away from the Sun", xref="paper", yref="paper", showarrow=False, font=dict(color='red',size=15))

    plotfile=path+sc_label+'_now.html'
    fig.write_html(plotfile)
    print('saved as ',plotfile)





        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

def plot_insitu_update_noaa_rtsw(sc, start, end, sc_label, path, **kwargs):
     '''
     sc = data
    
     '''
        
     #cut out data starting with start time that will be plotted, for better scaling
     startind=np.where(sc.time > start)[0][0]  
     sc=sc[startind:]    
        
     sns.set_style('darkgrid')
     sns.set_context('paper')
     
     fsize=10

     fig=plt.figure(figsize=(9,6), dpi=150)
     
     #sharex means that zooming in works with all subplots
     ax1 = plt.subplot(411) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
     
     plt.ylabel('B [nT] GSM',fontsize=fsize)
     plt.legend(loc=3,ncol=4,fontsize=fsize-2)
     ax1.set_xlim(start,end)
     ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d') )
     ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
     plt.ylim((-np.nanmax(sc.bt)-3,np.nanmax(sc.bt)+3))
     #ax1.set_xticklabels([]) does not work with sharex
     #plt.setp(ax1.get_xticklabels(), fontsize=6)
     plt.setp(ax1.get_xticklabels(), visible=False)

     #plt.title(sc_label+' data       start: '+start.strftime("%Y-%b-%d")+'       end: '+end.strftime("%Y-%b-%d"),fontsize=fsize)
     plt.title(sc_label+' data       start: '+start.strftime("%Y-%b-%d %H:%M UT")+'       end: '+end.strftime("%Y-%b-%d %H:%M UT"),fontsize=fsize)


     ax2 = plt.subplot(412,sharex=ax1) 
     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=0.7)
     plt.ylabel('V [km s$^{-1}$]',fontsize=fsize)
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     ax2.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(100))
     ax2.set_ylim((250, np.nanmax(sc.vt)+50))
     plt.setp(ax2.get_xticklabels(), visible=False)


     ax3 = plt.subplot(413,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)
     plt.ylabel('N [ccm$^{-3}$]',fontsize=fsize)
     ax3.set_xlim(start,end)
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     #ax3.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
     ax3.set_ylim((0, np.nanmax(sc.np)+10))
    
     plt.setp(ax3.get_xticklabels(), visible=False)


     ax4 = plt.subplot(414,sharex=ax1) 
     ax4.plot_date(sc.time,sc.tp/1e6,'-k',label='Tp',linewidth=0.7)
     plt.ylabel('T [MK]',fontsize=fsize)
     ax4.set_xlim(start,end)
     ax4.set_ylim((0, 1))
        
             
     ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %Hh') )
    
     
     plt.figtext(0.04,0.01,'Austrian Space Weather Office   GeoSphere Austria', color='black', ha='left',fontsize=fsize-4, style='italic')
     plt.figtext(0.98,0.01,'helioforecast.space', color='black', ha='right',fontsize=fsize-4, style='italic')

     
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

        
     #if now2 exists as keyword, save as the file with just now2 in filename:     
     if 'now2' in kwargs:
        plotfile=path+sc_label+'_now2.png'
        plt.savefig(plotfile)
        print('saved as ',plotfile)
        
        plotfile=path+sc_label+'_now2.pdf'
        plt.savefig(plotfile)
        print('saved as ',plotfile)
        
     
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
   


def plot_insitu_update_stereoa_beacon(sc, start, end, sc_label, path, coord, **kwargs):
     '''
     sc = data
    
     '''
     sns.set_style('darkgrid')
     sns.set_context('paper')
        
     #cut out data starting with start time that will be plotted, for better scaling
     startind=np.where(sc.time > start)[0][0]  
     sc=sc[startind:]            
        
        
     
     fsize=10

     fig=plt.figure(figsize=(9,5), dpi=150)
    
   
    
     
     ax1 = plt.subplot(311) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)     
     plt.ylabel('B [nT] '+coord,fontsize=fsize)
     plt.legend(loc=3,ncol=4,fontsize=fsize-2)
     ax1.set_xlim(start,end)
     ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
    
    
     plt.ylim(-np.nanmax(sc.bt)-3,np.nanmax(sc.bt)+3)
     #plt.setp(ax1.get_xticklabels(), fontsize=6)
     plt.setp(ax1.get_xticklabels(), visible=False)

     plt.title(sc_label+' data       start: '+start.strftime("%Y-%b-%d %H:%M")+' UT       end: '+end.strftime("%Y-%b-%d %H:%M")+' UT',fontsize=fsize)


     ax2 = plt.subplot(312,sharex=ax1) 
     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=1)
     plt.ylabel('V [km s$^{-1}$]',fontsize=fsize)
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     ax2.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(100))
     plt.ylim(200,800)
     plt.setp(ax2.get_xticklabels(), visible=False) 
    
    
     ax3 = plt.subplot(313,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=1)
     plt.ylabel('N [ccm$^{-3}$]',fontsize=fsize)
     ax3.set_xlim(start,end)
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     ax3.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
     ax3.set_ylim((0, np.nanmax(sc.np)+10))
     #plt.setp(ax3.get_xticklabels(), visible=False)
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %Hh') )

     #ax4 = plt.subplot(414,sharex=ax1) 
     #ax4.plot_date(sc.time,sc.tp,'-k',label='Tp',linewidth=0.7)
     #plt.ylabel('T [MK]',fontsize=fsize)
     #ax4.set_xlim(start,end)
     #ax4.set_ylim((0, 1))
     #ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
        

     plt.figtext(0.04,0.01,'Austrian Space Weather Office   GeoSphere Austria', color='black', ha='left',fontsize=fsize-4, style='italic')

     plt.figtext(0.97,0.01,'helioforecast.space', color='black', ha='right',fontsize=fsize-4, style='italic')
     
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

     #if now2 exists as keyword, save as the file with just now2 in filename:     
     if 'now2' in kwargs:
        plotfile=path+sc_label+'_now2.png'
        plt.savefig(plotfile)
        print('saved as ',plotfile)                
       
        plotfile=path+sc_label+'_now2.pdf'
        plt.savefig(plotfile)
        print('saved as ',plotfile)
        

        
      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
      
 

def plot_insitu(sc, start, end, sc_label, path, **kwargs):
     '''
     sc = data
    
     '''
     
     start=parse_time(start).datetime
     end=parse_time(end).datetime
     #print(start)
     #print(end)
     
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
     ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )
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
     plt.ylim((0, 100))
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
   
     
       
  
            
            
            
            
            
################################## SIRCAT  ############################################################
  
    
    
    
    

def plot_sircat_events(sc,sci,ic,name,icplotsdir):

  
    fileind='sircat/indices_sircat/SIRCAT_indices_'+name+'.p'
    
    #get indices of events for this spacecraft
    [sir_start_ind,hss_start_ind, sir_end_ind,hss_end_ind]=pickle.load(open(fileind, 'rb'))     
    
    
    if name=='STEREO-A' or name=='STEREO-B':
    
        for i in np.arange(np.size(sci)):            

            plot_insitu_sircat_mag_plasma(sc[sir_start_ind[i]-60*2*24:sir_end_ind[i]+3*60*24],\
                                 ic.sir_start_time[sci[i]]-datetime.timedelta(days=2), \
                                 ic.sir_end_time[sci[i]]+datetime.timedelta(days=3),name, icplotsdir,ic,sci[i])
            plt.close('all')       

    if name=='MAVEN':
    
        for i in np.arange(np.size(sci)):            

            plot_insitu_sircat_MAVEN(sc[sir_start_ind[i]-6*7:sir_end_ind[i]+6*7],\
                                 ic.sir_start_time[sci[i]]-datetime.timedelta(days=3), \
                                 ic.sir_end_time[sci[i]]+datetime.timedelta(days=5),name, icplotsdir,ic,sci[i])
            plt.close('all')       
            
            
            
    if name=='Wind':
    
        for i in np.arange(np.size(sci)):            

            plot_insitu_sircat_mag_plasma(sc[hss_start_ind[i]-60*2*24:hss_end_ind[i]+3*60*24],\
                                 ic.hss_start_time[sci[i]]-datetime.timedelta(days=2), \
                                 ic.hss_end_time[sci[i]]+datetime.timedelta(days=3),name, icplotsdir,ic,sci[i])
            plt.close('all')       

            
    if name=='PSP':
    
        for i in np.arange(np.size(sci)):            

            plot_insitu_sircat_mag_plasma(sc[sir_start_ind[i]-60*2*24:sir_end_ind[i]+3*60*24],\
                                 ic.sir_start_time[sci[i]]-datetime.timedelta(days=2), \
                                 ic.sir_end_time[sci[i]]+datetime.timedelta(days=3),name, icplotsdir,ic,sci[i])
            plt.close('all') 
            
            
            
            
            


def plot_insitu_sircat_MAVEN(sc, start, end, sc_label, path, ic,i, **kwargs):
     '''
     sc ... data
    
     '''
     
     start=parse_time(start).datetime
     end=parse_time(end).datetime
     #print(start)
     #print(end)
     
     sns.set_style('darkgrid')
     sns.set_context('paper')
                
     color_sir='black'        
     color_vtmax='tomato'
        

     fig=plt.figure(figsize=(9,6), dpi=150)
     
     #sharex means that zooming in works with all subplots
     ax1 = plt.subplot(411) 
     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
     #plot vertical lines    
     ax1.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[-500,500],linestyle='-',linewidth=1,color=color_sir)            
     ax1.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[-500,500],linestyle='--',linewidth=1,color=color_sir)            
     ax1.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[-500,500],color=color_sir,linestyle='-',linewidth=1)            
     #ax1.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[-500,500],color=color_vtmax,linestyle='-',linewidth=1,marker='')            
     #ax1.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[-500,500],'-k',linewidth=1)    
     plt.ylabel('B [nT]')
     plt.legend(loc=1,ncol=4,fontsize=6)
     ax1.set_xlim(start,end)
     if np.isnan(np.nanmin(sc.bt))==False:
         ax1.set_ylim(-np.nanmax(sc.bt)-5,np.nanmax(sc.bt)+5)   
     ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )
     #plt.ylim((-20, 20))
     #ax1.set_xticklabels([]) does not work with sharex
     #plt.setp(ax1.get_xticklabels(), fontsize=6)
     plt.setp(ax1.get_xticklabels(), visible=False)
     plt.title(sc_label+' data, start: '+start.strftime("%Y-%b-%d %H:%M")+'  end: '+end.strftime("%Y-%b-%d %H:%M"))

    
     ax2 = plt.subplot(412,sharex=ax1) 
     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=0.7)
     #plot vertical lines
     ax2.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[0,3000],linestyle='-',linewidth=1,color=color_sir)            
     ax2.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[0,3000],linestyle='--',linewidth=1,color=color_sir)            
     ax2.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[0,3000],color=color_sir,linestyle='-',linewidth=1)            
     #ax2.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[0,3000],color=color_vtmax,linestyle='-',linewidth=1,marker='')            
     #ax2.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[0,3000],'-k',linewidth=1)    
     plt.ylabel('V [km/s]')
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     if np.isnan(np.nanmin(sc.vt))==False:
        ax2.set_ylim(np.nanmin(sc.vt)-20,np.nanmax(sc.vt)+100)  
     #plt.ylim((250, 800))
     #ax2.set_xticklabels([])
     plt.setp(ax2.get_xticklabels(), visible=False)

    
     ax3 = plt.subplot(413,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)
     #plot vertical lines
     ax3.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[-10,1000],color=color_sir,linestyle='-',linewidth=1)       
     ax3.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[-10,1000],color=color_sir,linestyle='--',linewidth=1)            
     ax3.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[-10,1000],color=color_sir,linestyle='-',linewidth=1)            
     #ax3.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[-10,1000],color=color_vtmax,linestyle='-',linewidth=1.5,marker='')            
     #ax3.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[0,1000],'-k',linewidth=1)    
     plt.ylabel('N [ccm-3]')
     ax3.set_xlim(start,end)
     if np.isnan(np.nanmin(sc.np))==False:
         ax3.set_ylim(0,np.nanmax(sc.np)+10)   
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     #plt.ylim((0, 50))
     #ax3.set_xticklabels([])
     plt.setp(ax3.get_xticklabels(), visible=False)


     ax4 = plt.subplot(414,sharex=ax1) 
     ax4.plot_date(sc.time,sc.tp/1e6,'-k',linewidth=0.7)    
     #plot vertical lines
     ax4.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[0,10],'-k',linewidth=1,label='sir_start_time / sir_end_time') 
     ax4.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[0,10],'--k',linewidth=1,label='hss_start_time = stream interface time') 
     ax4.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[0,10],color=color_sir,linewidth=1,linestyle='-',label='sir_end_time')         
     #ax4.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[0,10],color=color_vtmax,linewidth=1,linestyle='-',label='hss_vtmax_time',marker='')
     #ax4.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[0,10],'-k',linewidth=1)    
     plt.ylabel('T [MK]')
     ax4.set_xlim(start,end)
     if np.isnan(np.nanmin(sc.tp))==False:
         ax4.set_ylim(0,np.nanmax(sc.tp/1e6)+0.1)   
     ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.legend(loc=1,ncol=1,fontsize=6)
     #plt.ylim((0, 0.5))
     
     plt.tight_layout()
     #plt.show()

     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     plotfile=path+ic.sircat_id[i]+'.png'  

     plt.savefig(plotfile)
     print('saved as ',plotfile)
   
                 
            
            
            

def plot_insitu_sircat_mag_plasma(sc, start, end, sc_label, path, ic,i):
     '''
     sc ... data
    
     '''
     
     start=parse_time(start).datetime
     end=parse_time(end).datetime
     #print(start)
     #print(end)
     
     sns.set_style('darkgrid')
     sns.set_context('paper')
        
        
     color_sir='black'         
     color_vtmax='mediumseagreen'
        

     fig=plt.figure(figsize=(9,6), dpi=150)
     
     #sharex means that zooming in works with all subplots
     ax1 = plt.subplot(411) 
     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',linewidth=0.5)
  
     ax1.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[-500,500],'-k',linewidth=1)            
     ax1.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[-500,500],'--k',linewidth=1)            
     ax1.plot_date([ic.sir_vtmax_time[i],ic.sir_vtmax_time[i]],[-500,500],color=color_vtmax,linestyle='-',linewidth=1.5)            
     ax1.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[-500,500],color=color_sir,linestyle='-',linewidth=1)            
     ax1.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[-500,500],color=color_vtmax,linestyle='--',linewidth=1)            
     ax1.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[-500,500],color=color_sir,linestyle='--',linewidth=1)    
     plt.ylabel('B [nT]')
     plt.legend(loc=3,ncol=4,fontsize=8)
     ax1.set_xlim(start,end)
     #if np.isnan(np.nanmin(sc.bt))==False:
     #ax1.set_ylim(-np.nanmax(sc.bt)-5,np.nanmax(sc.bt)+5)   
     ax1.set_ylim(-np.nanmax(sc.bt)-10,np.nanmax(sc.bt)+10)   
       
     ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )
     #plt.ylim((-20, 20))
     #ax1.set_xticklabels([]) does not work with sharex
     #plt.setp(ax1.get_xticklabels(), fontsize=6)
     plt.setp(ax1.get_xticklabels(), visible=False)
     plt.title(sc_label+' data, start: '+start.strftime("%Y-%b-%d %H:%M")+'  end: '+end.strftime("%Y-%b-%d %H:%M"))

    
     ax2 = plt.subplot(412,sharex=ax1) 
     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=0.7)
     #plot vertical lines
     ax2.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[0,3000],'-k',linewidth=1)            
     ax2.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[0,3000],'--k',linewidth=1)                    
     ax2.plot_date([ic.sir_vtmax_time[i],ic.sir_vtmax_time[i]],[0,3000],color=color_vtmax,linewidth=1.5, linestyle='-')                 
     ax2.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[0,3000],color=color_sir,linestyle='-',linewidth=1)            
     ax2.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[0,3000],color=color_vtmax,linestyle='--',linewidth=1)            
     ax2.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[0,3000],'--k',linewidth=1)    
     ax2.set_ylabel('V [km s$^{-1}]$')
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     ax2.set_ylim(250,800)  
     try: 
        ax2.set_ylim(np.nanmin(sc.vt)-20,np.nanmax(sc.vt)+100)  
     except:
        ax2.set_ylim(200,1000)  
        
     plt.setp(ax2.get_xticklabels(), visible=False)

    
     ax3 = plt.subplot(413,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)
     #plot vertical lines
     ax3.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[-100,10000],'-k',linewidth=1)            
     ax3.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[-100,10000],'--k',linewidth=1)            
     ax3.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[-100,10000],color=color_sir,linestyle='-',linewidth=1)            
     ax3.plot_date([ic.sir_vtmax_time[i],ic.sir_vtmax_time[i]],[-100,10000],color=color_vtmax,linewidth=1.5, linestyle='-')                 
     ax3.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[-100,10000],color=color_vtmax,linestyle='-',linewidth=1)            
     ax3.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[-100,10000],'--k',linewidth=1)    
     plt.ylabel('N [ccm-3]')
     ax3.set_xlim(start,end)
     ax3.set_ylim(0,np.nanmax(sc.np)+10)   
     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     #plt.ylim((0, 50))
     #ax3.set_xticklabels([])
     plt.setp(ax3.get_xticklabels(), visible=False)


     ax4 = plt.subplot(414,sharex=ax1) 
     ax4.plot_date(sc.time,sc.tp/1e6,'-k',linewidth=0.7)    
     #plot vertical lines
     ax4.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[-10,10],'-k',linewidth=1,label='sir_start_time / sir_end_time') 
     ax4.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[-10,10],'--k',linewidth=1,label='hss_start_time / hss_end_time') 
     ax4.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[-10,10],color=color_sir,linewidth=1, linestyle='-')              
     ax4.plot_date([ic.sir_vtmax_time[i],ic.sir_vtmax_time[i]],[-100,100],color=color_vtmax,linewidth=1.5,linestyle='-', label='sir_vtmax_time')    
     ax4.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[-100,100],color=color_vtmax,linewidth=1.5,linestyle='--',label='hss_vtmax_time')
     ax4.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[-10,10],'--k',linewidth=1)    
     plt.ylabel('T [MK]')
     ax4.set_xlim(start,end)
     ax4.set_ylim(0,np.nanmax(sc.tp/1e6)+0.2)   
     ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     legend=ax4.legend(loc=1,ncol=1,fontsize=6)
    
     for lh in legend.legendHandles:
        lh.set_markerfacecolor('none')
        lh.set_markeredgecolor('none')
        #plt.ylim((0, 0.5))
     
     plt.tight_layout()
     #plt.show()

     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     plotfile=path+ic.sircat_id[i]+'.png'  

     plt.savefig(plotfile)
     print('saved as ',plotfile)
   
     
     
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
            
            
################################## ICMECAT  ############################################################
  
def plot_stereo_hi_fov(pos, time_num, timeind,ax,sc):    
    
    #plots the STA FOV HI1 HI2
    
    #STB never flipped the camera:    
    if sc=='B': 
        ang1d=-4
        ang2d=-24
        ang3d=-18
        ang4d=-88
        lcolor='blue'
    
    if sc=='A': 
        ang1d=4
        ang2d=24
        ang3d=18
        ang4d=88
        lcolor='red'

        #STA flipped during conjunction and then again reversing when crossing the Sun Earth line in 2023
        
        if np.logical_and(time_num>mdates.date2num(datetime.datetime(2015,11,1)),time_num<mdates.date2num(datetime.datetime(2023,8,12))   ):  
            ang1d=-4
            ang2d=-24
            ang3d=-18
            ang4d=-88
                        

    #calculate endpoints
    
    #sta position
    x0=pos.x[timeind]
    y0=pos.y[timeind]
    z0=0
    
    #sta position 180° rotated    
    x1=-pos.x[timeind]
    y1=-pos.y[timeind]
    z1=0
    
    #rotate by 4 deg for HI1 FOV
    ang=np.deg2rad(ang1d)
    rot=np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])    
    [x2,y2]=np.dot(rot,[x1,y1])
    z2=0    
    #add to sta position
    x2f=x0+x2
    y2f=y0+y2    
    z2f=0
    
    #rotate 2
    ang2=np.deg2rad(ang2d)
    rot2=np.array([[np.cos(ang2), -np.sin(ang2)], [np.sin(ang2), np.cos(ang2)]])    
    [x3,y3]=np.dot(rot2,[x1,y1])
    z3=0    
    #add to sta position
    x3f=x0+x3
    y3f=y0+y3    
    z3f=0
    
    #rotate 3
    ang3=np.deg2rad(ang3d)
    rot3=np.array([[np.cos(ang3), -np.sin(ang3)], [np.sin(ang3), np.cos(ang3)]])    
    [x4,y4]=np.dot(rot3,[x1,y1])
    z4=0    
    #add to sta position
    x4f=x0+x4
    y4f=y0+y4    
    z4f=0    

    #rotate 4
    ang4=np.deg2rad(ang4d)
    rot4=np.array([[np.cos(ang4), -np.sin(ang4)], [np.sin(ang4), np.cos(ang4)]])    
    [x5,y5]=np.dot(rot4,[x1,y1])
    z5=0    
    #add to sta position
    x5f=x0+x5
    y5f=y0+y5    
    z5f=0    

    
    #convert to polar coordinates and plot
    [r0,t0,lon0]=hd.cart2sphere(x0,y0,z0)    
    #[r1,t1,lon1]=hd.cart2sphere(x1,y1,z1)    
    [r2,t2,lon2]=hd.cart2sphere(x2f,y2f,z2f)    
    [r3,t3,lon3]=hd.cart2sphere(x3f,y3f,z3f)    
    [r4,t4,lon4]=hd.cart2sphere(x4f,y4f,z4f)    
    [r5,t5,lon5]=hd.cart2sphere(x5f,y5f,z5f)    
    
    #ax.plot([lon0,lon1],[r0,r1],'--r',alpha=0.5)
    ax.plot([lon0,lon2],[r0,r2],linestyle='-',color=lcolor,alpha=0.3)
    ax.plot([lon0,lon3],[r0,r3],linestyle='-',color=lcolor,alpha=0.3)
    ax.plot([lon0,lon4],[r0,r4],linestyle='--',color=lcolor,alpha=0.3)
    ax.plot([lon0,lon5],[r0,r5],linestyle='--',color=lcolor,alpha=0.3)

    
    
    
    
    
    
    
    
def plot_stereo_hi_fov_old(pos, time_num, timeind,ax,sc):    
    
    
    #plots the STA FOV HI1 HI2
    
    #STB never flipped the camera:    
    if sc=='B': 
        ang1d=-4
        ang2d=-24
        ang3d=-18
        ang4d=-88
        lcolor='blue'
    
    if sc=='A': 
        ang1d=4
        ang2d=24
        ang3d=18
        ang4d=88
        lcolor='red'

        #STA flipped during conjunction
        if time_num>mdates.date2num(datetime.datetime(2015,11,1)):  
            ang1d=-4
            ang2d=-24
            ang3d=-18
            ang4d=-88



    #calculate endpoints
    
    #sta position
    x0=pos.x[timeind]
    y0=pos.y[timeind]
    z0=pos.z[timeind]
    
    #sta position 180° rotated    
    x1=-pos.x[timeind]
    y1=-pos.y[timeind]
    z1=pos.z[timeind]
    
    #rotate by 4 deg for HI1 FOV
    ang=np.deg2rad(ang1d)
    rot=np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])    
    [x2,y2]=np.dot(rot,[x1,y1])
    z2=0    
    #add to sta position
    x2f=x0+x2
    y2f=y0+y2    
    z2f=z0+z2
    
    #rotate 2
    ang2=np.deg2rad(ang2d)
    rot2=np.array([[np.cos(ang2), -np.sin(ang2)], [np.sin(ang2), np.cos(ang2)]])    
    [x3,y3]=np.dot(rot2,[x1,y1])
    z3=0    
    #add to sta position
    x3f=x0+x3
    y3f=y0+y3    
    z3f=z0+z3    
    
    #rotate 3
    ang3=np.deg2rad(ang3d)
    rot3=np.array([[np.cos(ang3), -np.sin(ang3)], [np.sin(ang3), np.cos(ang3)]])    
    [x4,y4]=np.dot(rot3,[x1,y1])
    z4=0    
    #add to sta position
    x4f=x0+x4
    y4f=y0+y4    
    z4f=z0+z4    

    #rotate 4
    ang4=np.deg2rad(ang4d)
    rot4=np.array([[np.cos(ang4), -np.sin(ang4)], [np.sin(ang4), np.cos(ang4)]])    
    [x5,y5]=np.dot(rot4,[x1,y1])
    z5=0    
    #add to sta position
    x5f=x0+x5
    y5f=y0+y5    
    z5f=z0+z5    

    
    #convert to polar coordinates and plot
    [r0,t0,lon0]=hd.cart2sphere(x0,y0,z0)    
    #[r1,t1,lon1]=hd.cart2sphere(x1,y1,z1)    
    [r2,t2,lon2]=hd.cart2sphere(x2f,y2f,z2f)    
    [r3,t3,lon3]=hd.cart2sphere(x3f,y3f,z3f)    
    [r4,t4,lon4]=hd.cart2sphere(x4f,y4f,z4f)    
    [r5,t5,lon5]=hd.cart2sphere(x5f,y5f,z5f)    
    
    #ax.plot([lon0,lon1],[r0,r1],'--r',alpha=0.5)
    ax.plot([lon0,lon2],[r0,r2],linestyle='-',color=lcolor,alpha=0.3)
    ax.plot([lon0,lon3],[r0,r3],linestyle='-',color=lcolor,alpha=0.3)
    ax.plot([lon0,lon4],[r0,r4],linestyle='--',color=lcolor,alpha=0.3)
    ax.plot([lon0,lon5],[r0,r5],linestyle='--',color=lcolor,alpha=0.3)

    
    
    
    
    
    
    
    
    
    
    
    
    
  

def plot_icmecat_positions_mag(time_date1,frame,ax,pos):
    '''
    sc = data
    '''

    plot_orbit=False
    plot_parker=True

    
    sns.set_style('darkgrid')
    sns.set_context('paper')    
    
    time1=mdates.date2num(time_date1)
        
    frame='HEEQ'
    sun_rot=26.24
   
    AUkm=149597870.7   

    #for parker spiral   
    theta=np.arange(0,np.deg2rad(180),0.01)
    res_in_days=1/24
    k=0
    
    fadeind=int(100/res_in_days)
    fsize=17 
    symsize_planet=140
    symsize_spacecraft=100
    
    #order in pos array
    #[p_psp, p_solo, p_sta, p_bepi, p_l1, p_stb_new, p_uly_new, p_mes_new, p_earth, p_mercury, p_venus, p_mars, p_jupiter, p_saturn, p_uranus, p_neptune])
    #open( 'results/positions/positions_HEEQ_1hr.p', "rb" ) )

    psp=pos[0]
    solo=pos[1]
    sta=pos[2]
    stb=pos[3]
    bepi=pos[4]
    l1=pos[5]
    uly=pos[6]
    mes=pos[7]
    earth=pos[8]
    mercury=pos[9]
    venus=pos[10]
    mars=pos[11]
    
    #find index for psp
    dct=time1-psp.time
    psp_timeind=np.argmin(abs(dct))

    dct=time1-bepi.time
    bepi_timeind=np.argmin(abs(dct))

    dct=time1-solo.time
    solo_timeind=np.argmin(abs(dct))

    dct=time1-earth.time
    earth_timeind=np.argmin(abs(dct))
    
    dct=time1-mars.time
    mars_timeind=np.argmin(abs(dct))

    dct=time1-venus.time
    venus_timeind=np.argmin(abs(dct))
    
    dct=time1-sta.time
    sta_timeind=np.argmin(abs(dct))
    
    dct=time1-stb.time
    stb_timeind=np.argmin(abs(dct))

    dct=time1-mes.time
    mes_timeind=np.argmin(abs(dct))
    
    dct=time1-mercury.time
    mercury_timeind=np.argmin(abs(dct))




    
    backcolor='black'
    psp_color='black'
    bepi_color='blue'
    solo_color='green'
        

    ax.scatter(venus.lon[venus_timeind], venus.r[venus_timeind]*np.cos(venus.lat[earth_timeind]), s=symsize_planet, c='orange', alpha=1,lw=0,zorder=3)
    ax.scatter(mercury.lon[earth_timeind], mercury.r[earth_timeind]*np.cos(mercury.lat[earth_timeind]), s=symsize_planet, c='dimgrey', alpha=1,lw=0,zorder=3)
    ax.scatter(earth.lon[earth_timeind], earth.r[earth_timeind]*np.cos(earth.lat[earth_timeind]), s=symsize_planet, c='mediumseagreen', alpha=1,lw=0,zorder=3)
    ax.scatter(sta.lon[sta_timeind], sta.r[sta_timeind]*np.cos(sta.lat[sta_timeind]), s=symsize_spacecraft, c='red', marker='s', alpha=1,lw=0,zorder=3)
    ax.scatter(mars.lon[mars_timeind], mars.r[mars_timeind]*np.cos(mars.lat[mars_timeind]), s=symsize_planet, c='tomato', alpha=1,lw=0,zorder=3)

    
    #text thats always there on the plot
    ax.text(sta.lon[sta_timeind]-0.15,sta.r[sta_timeind],'STEREO-A', color='red', ha='center',fontsize=fsize-4,verticalalignment='top')
    ax.text(0,0,'Sun', color='black', ha='center',fontsize=fsize-5,verticalalignment='top')
    ax.text(0,earth.r[earth_timeind]+0.12,'Earth', color='mediumseagreen', ha='center',fontsize=fsize-5,verticalalignment='center')
    
    
    
    
    #plot stereo hi fov
    plot_stereo_hi_fov(sta,time1, sta_timeind, ax,'A')
    
    if time1<mdates.date2num(datetime.datetime(2014,9,26)):  
        plot_stereo_hi_fov(stb,time1, stb_timeind, ax,'B')
    
    ########### Sun
    ax.scatter(0,0,s=100,c='yellow',alpha=1, edgecolors='black', linewidth=0.3)


    #parker spiral
    if plot_parker:
        for q in np.arange(0,12):
            #parker spiral
            #sidereal rotation
            omega=2*np.pi/(sun_rot*60*60*24) #solar rotation in seconds
            v=400/AUkm #km/s
            r0=695000/AUkm
            r=v/omega*theta+r0*7
            ax.plot(-theta+np.deg2rad(0+(360/24.47)*res_in_days*k+360/12*q), r, alpha=0.4, lw=0.5,color='grey',zorder=2)
    
    

    yset=0.43
    ydiff=0.04
    plt.figtext(0.01,yset,'              R     lon     lat', fontsize=fsize+2, ha='left',color=backcolor)

    
    if psp_timeind > 0:
        #plot trajectory
        ax.scatter(psp.lon[psp_timeind], psp.r[psp_timeind]*np.cos(psp.lat[psp_timeind]), s=symsize_spacecraft, c=psp_color, marker='s', alpha=1,lw=0,zorder=3)
        #plot position as text
        psp_text='PSP:   '+str(f'{psp.r[psp_timeind]:6.2f}')+str(f'{np.rad2deg(psp.lon[psp_timeind]):8.1f}')+str(f'{np.rad2deg(psp.lat[psp_timeind]):8.1f}')
        f5=plt.figtext(0.01,yset-ydiff*4,psp_text, fontsize=fsize, ha='left',color=psp_color)
        plt.text(psp.lon[psp_timeind]-0.15,psp.r[psp_timeind],'Parker Solar Probe', color='black', ha='center',fontsize=fsize-4,verticalalignment='top')

        if plot_orbit: 
            ax.plot(psp.lon[psp_timeind:psp_timeind+fadeind], psp.r[psp_timeind:psp_timeind+fadeind]*np.cos(psp.lat[psp_timeind:psp_timeind+fadeind]), c=psp_color, alpha=0.6,lw=1,zorder=3)

           
            
            
    if bepi_timeind > 0:
        ax.scatter(bepi.lon[bepi_timeind], bepi.r[bepi_timeind]*np.cos(bepi.lat[bepi_timeind]), s=symsize_spacecraft, c=bepi_color, marker='s', alpha=1,lw=0,zorder=3)
        bepi_text='Bepi:   '+str(f'{bepi.r[bepi_timeind]:6.2f}')+str(f'{np.rad2deg(bepi.lon[bepi_timeind]):8.1f}')+str(f'{np.rad2deg(bepi.lat[bepi_timeind]):8.1f}')
        f6=plt.figtext(0.01,yset-ydiff*5,bepi_text, fontsize=fsize, ha='left',color=bepi_color)
        plt.text(bepi.lon[bepi_timeind]-0.15,bepi.r[bepi_timeind],'Bepi Colombo', color='blue', ha='center',fontsize=fsize-4,verticalalignment='top')
  
        if plot_orbit: 
            ax.plot(bepi.lon[bepi_timeind:bepi_timeind+fadeind], bepi.r[bepi_timeind:bepi_timeind+fadeind]*np.cos(bepi.lat[bepi_timeind:bepi_timeind+fadeind]), c=bepi_color, alpha=0.6,lw=1,zorder=3)



    if solo_timeind > 0:
        ax.scatter(solo.lon[solo_timeind], solo.r[solo_timeind]*np.cos(solo.lat[solo_timeind]), s=symsize_spacecraft, c=solo_color, marker='s', alpha=1,lw=0,zorder=3)
        solo_text='SolO:  '+str(f'{solo.r[solo_timeind]:6.2f}')+str(f'{np.rad2deg(solo.lon[solo_timeind]):8.1f}')+str(f'{np.rad2deg(solo.lat[solo_timeind]):8.1f}')
        f7=plt.figtext(0.01,yset-ydiff*6,solo_text, fontsize=fsize, ha='left',color=solo_color)
        plt.text(solo.lon[solo_timeind]-0.15,solo.r[solo_timeind],'Solar Orbiter', color='green', ha='center',fontsize=fsize-4,verticalalignment='top')

        if plot_orbit: 
            ax.plot(solo.lon[solo_timeind:solo_timeind+fadeind], solo.r[solo_timeind:solo_timeind+fadeind]*np.cos(solo.lat[solo_timeind:solo_timeind+fadeind]), c=solo_color, alpha=0.6,lw=1,zorder=3)

    if plot_orbit: 
        ax.plot(sta.lon[sta_timeind:sta_timeind+fadeind], sta.r[sta_timeind:earth_timeind+fadeind]*np.cos(sta.lat[sta_timeind:sta_timeind+fadeind]), c='red', alpha=0.6,lw=1,zorder=3)

    
        
    #STEREO-B    
    if time1<mdates.date2num(datetime.datetime(2014,9,26)):        

        #marker and text on plot - stb has similar indices to earth_timeind
        ax.scatter(stb.lon[stb_timeind], stb.r[stb_timeind]*np.cos(stb.lat[stb_timeind]), s=symsize_spacecraft, c='blue', marker='s',alpha=1,lw=0,zorder=3)        
        #trajectory
        ax.plot(stb.lon[stb_timeind:stb_timeind+fadeind], stb.r[stb_timeind:stb_timeind+fadeind]*np.cos(stb.lat[stb_timeind:stb_timeind+fadeind]), c='blue', alpha=0.6,lw=1,zorder=3)        
        plt.text(stb.lon[stb_timeind],stb.r[stb_timeind]+0.12,'STEREO-B', color='blue', ha='center',fontsize=fsize-5,verticalalignment='center')

        stb_text='STB:   '+str(f'{stb.r[stb_timeind]:6.2f}')+str(f'{np.rad2deg(stb.lon[stb_timeind]):8.1f}')+str(f'{np.rad2deg(stb.lat[stb_timeind]):8.1f}')
        f13=plt.figtext(0.01,yset-ydiff*4,stb_text, fontsize=fsize, ha='left',color='blue')
        
        
        
    #VEX    
    if time1<mdates.date2num(datetime.datetime(2014,11,26)):        
        vex_text='VEX:   '+str(f'{venus.r[venus_timeind]:6.2f}')+str(f'{np.rad2deg(venus.lon[venus_timeind]):8.1f}')+str(f'{np.rad2deg(venus.lat[venus_timeind]):8.1f}')
        f11=plt.figtext(0.01,yset-ydiff*5,vex_text, fontsize=fsize, ha='left',color='orange')
        plt.text(venus.lon[venus_timeind],venus.r[venus_timeind]+0.12,'VEX', color='orange', ha='center',fontsize=fsize-5,verticalalignment='center')

    #MESSENGER    
    if time1<mdates.date2num(datetime.datetime(2011,3,18)):        
     
        #marker and text on plot
        ax.scatter(mes.lon[mes_timeind], mes.r[mes_timeind]*np.cos(mes.lat[mes_timeind]), s=symsize_planet, c='darkgrey', marker='s',alpha=1,lw=0,zorder=3)
        plt.text(mes.lon[mes_timeind],mes.r[mes_timeind]+0.12,'MESSENGER', color='darkgrey', ha='center',fontsize=fsize-5,verticalalignment='center')

        #side legend
        mes_text='MES:   '+str(f'{mes.r[mes_timeind]:6.2f}')+str(f'{np.rad2deg(mes.lon[mes_timeind]):8.1f}')+str(f'{np.rad2deg(mes.lat[mes_timeind]):8.1f}')       
        f12=plt.figtext(0.01,yset-ydiff*6,mes_text, fontsize=fsize, ha='left',color='darkgrey')
        
        
    if np.logical_and(time1>mdates.date2num(datetime.datetime(2011,3,18)), time1<mdates.date2num(datetime.datetime(2015,4,30))):        
        mes_text='MES:   '+str(f'{mercury.r[mercury_timeind]:6.2f}')+str(f'{np.rad2deg(mercury.lon[mercury_timeind]):8.1f}')+str(f'{np.rad2deg(mercury.lat[earth_timeind]):8.1f}')
        f12=plt.figtext(0.01,yset-ydiff*6,mes_text, fontsize=fsize, ha='left',color='darkgrey')
        plt.text(mercury.lon[mercury_timeind],mercury.r[mercury_timeind]+0.12,'MESSENGER', color='darkgrey', ha='center',fontsize=fsize-5,verticalalignment='center')



        
        
        
        
        
        
        
        
        
        
        



    if frame=='HEEQ': earth_text='Earth: '+str(f'{earth.r[earth_timeind]:6.2f}')+str(f'{0.0:8.1f}')+str(f'{np.rad2deg(earth.lat[earth_timeind]):8.1f}')
    else: earth_text='Earth: '+str(f'{earth.r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(earth.lon[earth_timeind]):8.1f}')+str(f'{np.rad2deg(earth.lat[earth_timeind]):8.1f}')

    mars_text='Mars:  '+str(f'{mars.r[mars_timeind]:6.2f}')+str(f'{np.rad2deg(mars.lon[mars_timeind]):8.1f}')+str(f'{np.rad2deg(mars.lat[mars_timeind]):8.1f}')
    sta_text='STA:   '+str(f'{sta.r[sta_timeind]:6.2f}')+str(f'{np.rad2deg(sta.lon[sta_timeind]):8.1f}')+str(f'{np.rad2deg(sta.lat[sta_timeind]):8.1f}')

    
 
    f10=plt.figtext(0.01,yset-ydiff*1,earth_text, fontsize=fsize, ha='left',color='mediumseagreen')
    f9=plt.figtext(0.01,yset-ydiff*2,mars_text, fontsize=fsize, ha='left',c='orangered')
    f8=plt.figtext(0.01,yset-ydiff*3,sta_text, fontsize=fsize, ha='left',c='red')

    
    
    
    ########### time
    plt.figtext(0.68,0.45,time_date1.strftime("%Y %B %d  %H:%M"),fontsize=fsize+6, ha='left',c='black')
    
    ########## legend    
    plt.figtext(0.99,0.01,'helioforecast.space     Austrian Space Weather Office     GeoSphere Austria', color='black', ha='right',fontsize=fsize-8)
    #plt.figtext(0.85,0.05,'――― 100 days future trajectory', color='black', ha='center',fontsize=fsize-3)

     

    #set axes
    ax.set_theta_zero_location('E')
 
    #plt.rgrids((0.10,0.39,0.72,1.00,1.52),('0.10','0.39','0.72','1.0','1.52 AU'),angle=125, fontsize=fsize,alpha=0.9, color=backcolor)
    plt.rgrids((0.1,0.3,0.5,0.7,1.0,1.3,1.6),('0.1','0.3','0.5','0.7','1.0','1.3','1.6 AU'),angle=125, fontsize=fsize-2,alpha=0.4, color=backcolor)

    #ax.set_ylim(0, 1.75) with Mars
    ax.set_ylim(0, 1.25) 

    plt.thetagrids(range(0,360,45),(u'0\u00b0 '+frame+' longitude',u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'- 135\u00b0',u'- 90\u00b0',u'- 45\u00b0'), fmt='%d',ha='left',fontsize=fsize,color=backcolor, zorder=5, alpha=0.9)

    plt.tight_layout()

     
     
     
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


def plot_icmecat_positions_mag_plasma(time_date1,frame,ax,pos,name):
    
    sns.set_style('darkgrid')
    sns.set_context('paper')    
    
    #convert old matplotlib so consistent with the positions file (need to change this)
    #time1=mdates.date2num(time_date1)-mdates.date2num(np.datetime64('0000-12-31'))
    time1=mdates.date2num(time_date1)
    
    frame='HEEQ'
    sun_rot=26.24
   
    #sidereal solar rotation rate
    #if frame=='HCI': sun_rot=24.47
    #synodic
    #if frame=='HEEQ': 
    AUkm=149597870.7   

    #for parker spiral   
    theta=np.arange(0,np.deg2rad(180),0.01)
    res_in_days=1/24
    k=0
    
    plot_orbit=False
    plot_parker=True
    fadeind=int(100/res_in_days)
    fsize=17 
    symsize_planet=140
    symsize_spacecraft=100

    
    #order in pos array
   


    psp=pos[0]
    solo=pos[1]
    sta=pos[2]
    stb=pos[3]
    bepi=pos[4]
    l1=pos[5]
    juno=pos[6]
    juice=pos[7]
    uly=pos[8]
    mes=pos[9]
    earth=pos[10]
    mercury=pos[11]
    venus=pos[12]
    mars=pos[13]
    
    #find index for psp
    dct=time1-psp.time
    psp_timeind=np.argmin(abs(dct))

    dct=time1-bepi.time
    bepi_timeind=np.argmin(abs(dct))

    dct=time1-solo.time
    solo_timeind=np.argmin(abs(dct))

    dct=time1-juno.time
    juno_timeind=np.argmin(abs(dct))
    
    dct=time1-uly.time
    uly_timeind=np.argmin(abs(dct))

    dct=time1-earth.time
    earth_timeind=np.argmin(abs(dct))
    
    dct=time1-venus.time
    venus_timeind=np.argmin(abs(dct))
    
    dct=time1-mars.time
    mars_timeind=np.argmin(abs(dct))

    dct=time1-mercury.time
    mercury_timeind=np.argmin(abs(dct))

    dct=time1-sta.time
    sta_timeind=np.argmin(abs(dct))
    
    dct=time1-stb.time
    stb_timeind=np.argmin(abs(dct))

    dct=time1-mes.time
    mes_timeind=np.argmin(abs(dct))


    
    backcolor='black'
    psp_color='black'
    bepi_color='blue'
    solo_color='green'
        

    ax.scatter(venus.lon[venus_timeind], venus.r[venus_timeind]*np.cos(venus.lat[venus_timeind]), s=symsize_planet, c='gold', alpha=1,lw=0,zorder=3)
    ax.scatter(mercury.lon[mercury_timeind], mercury.r[mercury_timeind]*np.cos(mercury.lat[mercury_timeind]), s=symsize_planet, c='dimgrey', alpha=1,lw=0,zorder=3)
    ax.scatter(earth.lon[earth_timeind], earth.r[earth_timeind]*np.cos(earth.lat[earth_timeind]), s=symsize_planet, c='mediumseagreen', alpha=1,lw=0,zorder=3)    
    ax.scatter(mars.lon[mars_timeind], mars.r[mars_timeind]*np.cos(mars.lat[mars_timeind]), s=symsize_planet, c='orangered', alpha=1,lw=0,zorder=3)

        
    #text thats always there on the plot
    ax.text(0,0,'Sun', color='black', ha='center',fontsize=fsize-5,verticalalignment='top')
    ax.text(0,earth.r[earth_timeind]+0.12,'Earth', color='mediumseagreen', ha='center',fontsize=fsize-5,verticalalignment='center')
    
    # Sun
    ax.scatter(0,0,s=100,c='yellow',alpha=1, edgecolors='black', linewidth=0.3)

    
    #parker spiral
    if plot_parker:
        for q in np.arange(0,12):
            #parker spiral
            #sidereal rotation
            omega=2*np.pi/(sun_rot*60*60*24) #solar rotation in seconds
            v=400/AUkm #km/s
            r0=695000/AUkm
            r=v/omega*theta+r0*7
            ax.plot(-theta+np.deg2rad(0+(360/24.47)*res_in_days*k+360/12*q), r*2, alpha=0.4, lw=0.5,color='grey',zorder=2)
    
    
    xset1=0.49
    xset2=0.80 #0.77
    yset=0.97
    ydiff=0.04
    plt.figtext(xset1,yset,'              R     lon     lat', fontsize=fsize+2, ha='left',color=backcolor)


    earth_text='Earth: '+str(f'{earth.r[earth_timeind]:6.2f}')+str(f'{0.0:8.1f}')+str(f'{np.rad2deg(earth.lat[earth_timeind]):8.1f}')
    f10=plt.figtext(xset1,yset-ydiff*1,earth_text, fontsize=fsize, ha='left',color='mediumseagreen')
       
    
    #make a ring around if its the active spacecraft
    if name=='Wind':
               ax.scatter(earth.lon[earth_timeind], earth.r[earth_timeind]*np.cos(earth.lat[earth_timeind]), marker='o',s=350, zorder=4,edgecolors='black',facecolors='none')

                 
            
            
    
    mars_text='Mars:  '+str(f'{mars.r[mars_timeind]:6.2f}')+str(f'{np.rad2deg(mars.lon[mars_timeind]):8.1f}')+str(f'{np.rad2deg(mars.lat[mars_timeind]):8.1f}')
    f9=plt.figtext(xset1,yset-ydiff*2,mars_text, fontsize=fsize, ha='left',c='tomato')

    #stereo-a    
    if time1>mdates.date2num(datetime.datetime(2007,1,1)):        

        
        sta_text='STA:   '+str(f'{sta.r[sta_timeind]:6.2f}')+str(f'{np.rad2deg(sta.lon[sta_timeind]):8.1f}')+str(f'{np.rad2deg(sta.lat[sta_timeind]):8.1f}')


        ax.text(sta.lon[sta_timeind]-0.15,sta.r[sta_timeind],'STEREO-A', color='red', ha='center',fontsize=fsize-4,verticalalignment='top')
        ax.scatter(sta.lon[sta_timeind], sta.r[sta_timeind]*np.cos(sta.lat[sta_timeind]), s=symsize_spacecraft, c='red', marker='s', alpha=1,lw=0,zorder=3)
        f8=plt.figtext(xset1,yset-ydiff*3,sta_text, fontsize=fsize, ha='left',c='red')

        #plot stereo hi fov
        plot_stereo_hi_fov(sta,time1, sta_timeind, ax,'A')
        
        #make a ring around if its the active spacecraft
        if name=='STEREO-A':
               ax.scatter(sta.lon[sta_timeind], sta.r[sta_timeind]*np.cos(sta.lat[sta_timeind]), marker='o',s=350, zorder=4,edgecolors='black',facecolors='none')

        

        if plot_orbit: 
            ax.plot(sta.lon[earth_timeind:earth_timeind+fadeind], sta.r[earth_timeind:earth_timeind+fadeind]*np.cos(sta.lat[earth_timeind:earth_timeind+fadeind]), c='red', alpha=0.6,lw=1,zorder=3)

    
    
    if psp_timeind > 0:
        #plot trajectory
        ax.scatter(psp.lon[psp_timeind], psp.r[psp_timeind]*np.cos(psp.lat[psp_timeind]), s=symsize_spacecraft, c=psp_color, marker='s', alpha=1,lw=0,zorder=3)
        #plot position as text
        psp_text='PSP:   '+str(f'{psp.r[psp_timeind]:6.2f}')+str(f'{np.rad2deg(psp.lon[psp_timeind]):8.1f}')+str(f'{np.rad2deg(psp.lat[psp_timeind]):8.1f}')
        f5=plt.figtext(xset2,yset-ydiff*1,psp_text, fontsize=fsize, ha='left',color=psp_color)
        plt.text(psp.lon[psp_timeind]-0.15,psp.r[psp_timeind],'Parker Solar Probe', color='black', ha='center',fontsize=fsize-4,verticalalignment='top')

        #make a ring around if its the active spacecraft
        if name=='PSP':
               ax.scatter(psp.lon[psp_timeind], psp.r[psp_timeind]*np.cos(psp.lat[psp_timeind]), marker='o',s=350, zorder=4,edgecolors='black',facecolors='none')

        
        
        if plot_orbit: 
            ax.plot(psp.lon[psp_timeind:psp_timeind+fadeind], psp.r[psp_timeind:psp_timeind+fadeind]*np.cos(psp.lat[psp_timeind:psp_timeind+fadeind]), c=psp_color, alpha=0.6,lw=1,zorder=3)

            
            
    if bepi_timeind > 0:
        ax.scatter(bepi.lon[bepi_timeind], bepi.r[bepi_timeind]*np.cos(bepi.lat[bepi_timeind]), s=symsize_spacecraft, c=bepi_color, marker='s', alpha=1,lw=0,zorder=3)
        bepi_text='Bepi:   '+str(f'{bepi.r[bepi_timeind]:6.2f}')+str(f'{np.rad2deg(bepi.lon[bepi_timeind]):8.1f}')+str(f'{np.rad2deg(bepi.lat[bepi_timeind]):8.1f}')
        f6=plt.figtext(xset2,yset-ydiff*2,bepi_text, fontsize=fsize, ha='left',color=bepi_color)
        plt.text(bepi.lon[bepi_timeind]-0.15,bepi.r[bepi_timeind],'BepiColombo', color='blue', ha='center',fontsize=fsize-4,verticalalignment='top')
        
        
        #make a ring around if its the active spacecraft
        if name=='BepiColombo':
               ax.scatter(bepi.lon[bepi_timeind], bepi.r[bepi_timeind]*np.cos(bepi.lat[bepi_timeind]), marker='o',s=350, zorder=4,edgecolors='black',facecolors='none')

        
  
        if plot_orbit: 
            ax.plot(bepi.lon[bepi_timeind:bepi_timeind+fadeind], bepi.r[bepi_timeind:bepi_timeind+fadeind]*np.cos(bepi.lat[bepi_timeind:bepi_timeind+fadeind]), c=bepi_color, alpha=0.6,lw=1,zorder=3)



    if solo_timeind > 0:
        ax.scatter(solo.lon[solo_timeind], solo.r[solo_timeind]*np.cos(solo.lat[solo_timeind]), s=symsize_spacecraft, c=solo_color, marker='s', alpha=1,lw=0,zorder=3)
        solo_text='SolO:  '+str(f'{solo.r[solo_timeind]:6.2f}')+str(f'{np.rad2deg(solo.lon[solo_timeind]):8.1f}')+str(f'{np.rad2deg(solo.lat[solo_timeind]):8.1f}')
        f7=plt.figtext(xset2,yset-ydiff*3,solo_text, fontsize=fsize, ha='left',color=solo_color)
        plt.text(solo.lon[solo_timeind]-0.15,solo.r[solo_timeind],'Solar Orbiter', color='green', ha='center',fontsize=fsize-4,verticalalignment='top')
        
        #make a ring around if its the active spacecraft
        if name=='SolarOrbiter':
               ax.scatter(solo.lon[solo_timeind], solo.r[solo_timeind]*np.cos(solo.lat[solo_timeind]), marker='o',s=350, zorder=4,edgecolors='black',facecolors='none')


        if plot_orbit: 
            ax.plot(solo.lon[solo_timeind:solo_timeind+fadeind], solo.r[solo_timeind:solo_timeind+fadeind]*np.cos(solo.lat[solo_timeind:solo_timeind+fadeind]), c=solo_color, alpha=0.6,lw=1,zorder=3)

 
    
        
    #STEREO-B    
    if time1>mdates.date2num(datetime.datetime(2007,1,1)) and time1<mdates.date2num(datetime.datetime(2014,9,26)):        

        ax.scatter(stb.lon[stb_timeind], stb.r[stb_timeind]*np.cos(stb.lat[stb_timeind]), s=symsize_spacecraft, c='blue', marker='s',alpha=1,lw=0,zorder=3)        
        #ax.plot(stb.lon[stb_timeind:stb_timeind+fadeind], stb.r[stb_timeind:stb_timeind+fadeind]*np.cos(stb.lat[stb_timeind:stb_timeind+fadeind]), c='blue', alpha=0.6,lw=1,zorder=3)        
        plt.text(stb.lon[stb_timeind],stb.r[stb_timeind]+0.12,'STEREO-B', color='blue', ha='center',fontsize=fsize-5,verticalalignment='center')

        stb_text='STB:   '+str(f'{stb.r[stb_timeind]:6.2f}')+str(f'{np.rad2deg(stb.lon[stb_timeind]):8.1f}')+str(f'{np.rad2deg(stb.lat[stb_timeind]):8.1f}')
        f13=plt.figtext(xset2,yset-ydiff*1,stb_text, fontsize=fsize, ha='left',color='blue')
        
        plot_stereo_hi_fov(stb,time1, stb_timeind, ax,'B')


        
        
    #VEX    
    if time1>mdates.date2num(datetime.datetime(2007,1,1)) and time1<mdates.date2num(datetime.datetime(2014,11,26)):        
        vex_text='VEX:   '+str(f'{venus.r[venus_timeind]:6.2f}')+str(f'{np.rad2deg(venus.lon[venus_timeind]):8.1f}')+str(f'{np.rad2deg(venus.lat[venus_timeind]):8.1f}')
        f11=plt.figtext(xset2,yset-ydiff*2,vex_text, fontsize=fsize, ha='left',color='orange')
        plt.text(venus.lon[earth_timeind],venus.r[earth_timeind]+0.12,'VEX', color='orange', ha='center',fontsize=fsize-5,verticalalignment='center')

    #MESSENGER    
    if time1>mdates.date2num(datetime.datetime(2007,1,1)) and time1<mdates.date2num(datetime.datetime(2011,3,18)):        
     
        #marker and text on plot
        ax.scatter(mes.lon[mes_timeind], mes.r[mes_timeind]*np.cos(mes.lat[mes_timeind]), s=symsize_planet, c='darkgrey', marker='s',alpha=1,lw=0,zorder=3)
        plt.text(mes.lon[mes_timeind],mes.r[mes_timeind]+0.12,'MESSENGER', color='darkgrey', ha='center',fontsize=fsize-5,verticalalignment='center')

        #side legend
        mes_text='MES:   '+str(f'{mes.r[mes_timeind]:6.2f}')+str(f'{np.rad2deg(mes.lon[mes_timeind]):8.1f}')+str(f'{np.rad2deg(mes.lat[mes_timeind]):8.1f}')       
        f12=plt.figtext(xset2,yset-ydiff*3,mes_text, fontsize=fsize, ha='left',color='darkgrey')
        
        
    if np.logical_and(time1>mdates.date2num(datetime.datetime(2011,3,18)), time1<mdates.date2num(datetime.datetime(2015,4,30))):        
        mes_text='MES:   '+str(f'{mercury.r[mercury_timeind]:6.2f}')+str(f'{np.rad2deg(mercury.lon[mercury_timeind]):8.1f}')+str(f'{np.rad2deg(mercury.lat[earth_timeind]):8.1f}')
        f12=plt.figtext(xset2,yset-ydiff*3,mes_text, fontsize=fsize, ha='left',color='darkgrey')
        plt.text(mercury.lon[mercury_timeind],mercury.r[mercury_timeind]+0.12,'MESSENGER', color='darkgrey', ha='center',fontsize=fsize-5,verticalalignment='center')


    
    #Juno when active    
    if time1>mdates.date2num(datetime.datetime(2011, 8, 25, 15, 19)) and time1<mdates.date2num(datetime.datetime(2016, 6, 28, 23, 59)):        

        ax.scatter(juno.lon[juno_timeind], juno.r[juno_timeind]*np.cos(juno.lat[juno_timeind]), s=symsize_spacecraft, c='yellowgreen', marker='s', alpha=1,lw=0,zorder=3)
        
        
        #add annotation only when Juno is the active spacecraft; when not it must be < 1.2 AU, so plots for all other spacecraft are fine
        if name=='Juno':
            plt.text(juno.lon[juno_timeind],juno.r[juno_timeind]+0.12,'Juno', color='yellowgreen', ha='center',fontsize=fsize-5,verticalalignment='center')
        else:
            if juno.r[juno_timeind]< 1.2: 
                 plt.text(juno.lon[juno_timeind],juno.r[juno_timeind]+0.12,'Juno', color='yellowgreen', ha='center',fontsize=fsize-5,verticalalignment='center')

                

        #position text always
        juno_text='Juno:   '+str(f'{juno.r[juno_timeind]:6.2f}')+str(f'{np.rad2deg(juno.lon[juno_timeind]):8.1f}')+str(f'{np.rad2deg(juno.lat[juno_timeind]):8.1f}')
        f11=plt.figtext(xset2,yset-ydiff*4,juno_text, fontsize=fsize, ha='left',color='yellowgreen')

    #ring when Juno active
    if name=='Juno':    
        ax.scatter(juno.lon[juno_timeind], juno.r[juno_timeind]*np.cos(juno.lat[juno_timeind]), marker='o',s=350, zorder=4, edgecolors='black', facecolors='none')


        
    ############Ulysses when active    
    if time1>mdates.date2num(datetime.datetime(1990, 10, 7, 0, 0)) and time1<mdates.date2num(datetime.datetime(2009, 12, 31, 22, 0)):        

        ax.scatter(uly.lon[uly_timeind], uly.r[uly_timeind]*np.cos(uly.lat[uly_timeind]), s=symsize_spacecraft, c='chocolate', marker='s', alpha=1,lw=0,zorder=3)
       
        
        #add annotation only when Juno is the active spacecraft; when not it must be < 1.2 AU, so plots for all other spacecraft are fine
        if name=='ULYSSES':
            plt.text(uly.lon[uly_timeind],uly.r[uly_timeind]*np.cos(uly.lat[uly_timeind])+0.12,'Ulysses', color='chocolate', ha='center',fontsize=fsize-5,verticalalignment='center')
        else:     
            if uly.r[uly_timeind]< 1.2: 
                plt.text(uly.lon[uly_timeind],uly.r[uly_timeind]*np.cos(uly.lat[uly_timeind])+0.12,'Ulysses', color='chocolate', ha='center',fontsize=fsize-5,verticalalignment='center')
        
        
        uly_text='Ulysses: '+str(f'{uly.r[uly_timeind]:6.2f}')+str(f'{np.rad2deg(uly.lon[uly_timeind]):8.1f}')+str(f'{np.rad2deg(uly.lat[uly_timeind]):8.1f}')
        f11=plt.figtext(xset2,yset-ydiff*4,uly_text, fontsize=fsize, ha='left',color='chocolate')        
        
    #ring when Ulysses active
    if name=='ULYSSES':    
        ax.scatter(uly.lon[uly_timeind], uly.r[uly_timeind]*np.cos(uly.lat[uly_timeind]), marker='o',s=350, zorder=4, edgecolors='black', facecolors='none')
        
    
    
    ########### time
    #this is the icme_start_time
    #plt.figtext(0.50,0.08,time_date1.strftime("%Y %B %d  %H:%M")+' UT',fontsize=fsize, ha='left',c='black')

    #this is the time thats used from the positions file for Earth
    plt.figtext(0.50,0.08,mdates.num2date(earth.time[earth_timeind]).strftime("%Y %B %d  %H:%M")+' UT',fontsize=fsize, ha='left',c='black')

    
    ########## legend    
    plt.figtext(0.99,0.01,'helioforecast.space/icmecat   Austrian Space Weather Office   GeoSphere Austria', color='black', ha='right',fontsize=fsize-9)
    #plt.figtext(0.87,0.05,'――― 100 days future trajectory', color='black', ha='center',fontsize=fsize-3)

     
    #set axes
    ax.set_theta_zero_location('E')
 
    #plt.rgrids((0.10,0.39,0.72,1.00,1.52),('0.10','0.39','0.72','1.0','1.52 AU'),angle=125, fontsize=fsize,alpha=0.9, color=backcolor)
    plt.rgrids((0.1,0.3,0.5,0.7,1.0,1.3,1.6,2.0,2.5),('0.1','0.3','0.5','0.7','1.0','1.3','1.6 AU','2.0','2.5'),angle=125, fontsize=fsize-2,alpha=0.4, color=backcolor)

    #ax.set_ylim(0, 1.75) with Mars
    ax.set_ylim(0, 1.25) 
    
    ############# Juno and Ulysses different limits
    
    if name=='Juno':
        ax.set_ylim(0,juno.r[juno_timeind]+0.2)       
        
        
        #for later times make grid different
        if time1>mdates.date2num(datetime.datetime(2014, 4, 1)):
             plt.rgrids((0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5),('0.5','1.0','1.5','2.0','2.5','3.0','3.5 AU','4.0','4.5','5.0','5.5'),angle=125, fontsize=fsize-2,alpha=0.4, color=backcolor)
             ax.set_ylim(0,juno.r[juno_timeind]+0.2)
   
                
    if name=='ULYSSES':
        ax.set_ylim(0,uly.r[uly_timeind]+0.2)
        if uly.r[uly_timeind] > 2.0:
            plt.rgrids((0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5),('0.5','1.0','1.5','2.0','2.5','3.0','3.5 AU','4.0','4.5','5.0','5.5'),angle=125, fontsize=fsize-2,alpha=0.4, color=backcolor)
            ax.set_ylim(0,uly.r[uly_timeind]+0.2)
            
                            
            

    plt.thetagrids(range(0,360,45),(u'0\u00b0 '+frame,u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'- 135\u00b0',u'- 90\u00b0',u'- 45\u00b0'), fmt='%d',ha='left',fontsize=fsize,color=backcolor, zorder=5, alpha=0.9)

    plt.tight_layout()

     
     
    
    
    
    
    
    
    
    
    
   
  
     

def plot_insitu_icmecat_mag_plasma(sc, start, end, sc_label, path, ic,i, pos):
     '''
     sc ... data    
     '''
     
     start=parse_time(start).datetime
     end=parse_time(end).datetime
     #print(start)
     #print(end)
     
     sns.set_style('darkgrid')
     sns.set_context('paper')
        
     fs=13

     fig=plt.figure(1,figsize=(16,10), dpi=150)    
       
     ax1 = plt.subplot2grid((4,2), (0, 0))
     
     plt.title(sc_label+'    ICME start: '+\
               ic.icme_start_time[i].strftime("%Y-%b-%d %H:%M")+ \
               '   MO start: '+ic.mo_start_time[i].strftime("%Y-%b-%d %H:%M")+ \
               '   MO end: '+ic.mo_end_time[i].strftime("%Y-%b-%d %H:%M"),fontsize=10)
     

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
        
     #plot vertical lines
     ax1.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[-2000,2000],'-k',linewidth=1)            
     ax1.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[-2000,2000],'-k',linewidth=1)            
     ax1.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[-2000,2000],'-k',linewidth=1)            
       
        
     coord_string='RTN'   
    
     if sc_label=='STEREO-B':     coord_string='SCEQ'   
     if sc_label=='MESSENGER':    coord_string='SCEQ'   
     if sc_label=='VEX':          coord_string='SCEQ'   


     ax1.set_ylabel('B [nT] '+coord_string, fontsize=fs)
     leg=ax1.legend(loc='lower right',ncol=4,fontsize=fs-3)
     for line in leg.get_lines():
         line.set_linewidth(2.0)        
       
        
     ax1.set_xlim(start,end)

     bscale=np.nanmax(sc.bt)
    
     if np.isfinite(np.nanmin(sc.bt)):
            ax1.set_ylim(-bscale-5,bscale+5)   


     #ax1.text(sc.time[30],  bscale, 'west', fontsize=fs-4, color='green')
     #ax1.text(sc.time[30], -bscale, 'east', fontsize=fs-4, color='green')
            
     #ax1.text(sc.time[300],  bscale, 'north', fontsize=fs-4, color='blue')
     #ax1.text(sc.time[300], -bscale, 'south', fontsize=fs-4, color='blue')

     #add components with position on full figure
     ax1.text(0.08,0.94,'west', fontsize=fs-4, color='green',ha='center', va='center', transform=fig.transFigure)
     ax1.text(0.08,0.77,'east', fontsize=fs-4, color='green',ha='center', va='center', transform=fig.transFigure)
              
     ax1.text(0.11,0.94,'north', fontsize=fs-4, color='blue',ha='center', va='center', transform=fig.transFigure)
     ax1.text(0.11,0.77,'south', fontsize=fs-4, color='blue',ha='center', va='center', transform=fig.transFigure)
    
     ax1.tick_params(axis='x', labelsize=fs)  
     ax1.tick_params(axis='y', labelsize=fs)  
    
     plt.setp(ax1.get_xticklabels(), visible=False)
        

     ax2 = plt.subplot2grid((4,2), (1, 0),sharex=ax1)
     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=0.7)    
     
     #PSP has a lot of values that are lots of nans with a few data points in between
     if sc_label=='PSP': ax2.plot_date(sc.time,sc.vt,'ko',markersize=1)    
     

     #plot vertical lines
     ax2.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[0,3000],'-k',linewidth=1)            
     ax2.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[0,3000],'-k',linewidth=1)            
     ax2.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[0,3000],'-k',linewidth=1)            
     ax2.set_ylabel('V [km s$^{-1}]$',fontsize=fs)
     ax2.set_xlim(start,end)
     ax2.set_ylim(250,800)   
    
     ax2.tick_params(axis='x', labelsize=fs)
     ax2.tick_params(axis='y', labelsize=fs)  
    
     #check if plasma data exists
     if np.isfinite(np.nanmin(sc.vt)):
         ax2.set_ylim(np.nanmin(sc.vt)-20,np.nanmax(sc.vt)+100)   

     plt.setp(ax2.get_xticklabels(), visible=False)

    
     ax3 = plt.subplot2grid((4,2), (2, 0),sharex=ax1)
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7) 
    
     if sc_label=='PSP': ax3.plot_date(sc.time,sc.np,'ko',markersize=1)   
        
     #plot vertical lines
     ax3.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[0,10000],'-k',linewidth=1)
     ax3.plot_date([ic.mo_start_time[i],  ic.mo_start_time[i]],  [0,10000],'-k',linewidth=1)
     ax3.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],        [0,10000],'-k',linewidth=1)            
     ax3.set_ylabel('N [cm$^{-3}$]',fontsize=fs)
     ax3.set_xlim(start,end)
     ax3.set_ylim(0,100)   
     if np.isfinite(np.nanmin(sc.np)):
         ax3.set_ylim(0,np.nanmax(sc.np)+10)   
     #plt.ylim((0, 50))
     #ax3.set_xticklabels([])
    
     ax3.tick_params(axis='x', labelsize=fs)  
     ax3.tick_params(axis='y', labelsize=fs)  
    
     plt.setp(ax3.get_xticklabels(), visible=False)


     ax4 = plt.subplot2grid((4,2), (3, 0),sharex=ax1)
     ax4.plot_date(sc.time,sc.tp/1e6,'-k',label='Tp',linewidth=0.7)
     if sc_label=='PSP': ax4.plot_date(sc.time,sc.tp/1e6,'ko',markersize=1)   
        

     #plot vertical lines
     ax4.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[0,100],'-k',linewidth=1)            
     ax4.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[0,100],'-k',linewidth=1)            
     ax4.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[0,100],'-k',linewidth=1)            
     ax4.set_ylabel('T [MK]', fontsize=fs)
     ax4.set_xlim(start,end)
     ax4.set_ylim(0,1)   
     if np.isfinite(np.nanmin(sc.tp)):
         ax4.set_ylim(0,np.nanmax(sc.tp/1e6)+0.2)   
     ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b %d') )
     ax4.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        
     ax4.tick_params(axis='x', labelsize=fs-2)  
     ax4.tick_params(axis='y', labelsize=fs)  
        
    
     ################positions plot
     ax5 = plt.subplot2grid((4,2), (0, 1), rowspan=4, projection='polar')

     plot_icmecat_positions_mag_plasma(ic.icme_start_time[i],'HEEQ',ax5,pos,sc_label)
       
         
        
     if sc_label=='PSP': 
        plt.figtext(0.5,0.03,'Data source: Parker Solar Probe (FIELDS, UCBerkeley; SWEAP, UMich)',fontsize=11, ha='left',color='black')
       
     if sc_label=='STEREO-A': 
        plt.figtext(0.5,0.03,'Data source: STEREO-Ahead (IMPACT, UCBerkeley; PLASTIC, UNH)',fontsize=11, ha='left',color='black')

     if sc_label=='STEREO-B': 
        plt.figtext(0.5,0.03,'Data source: STEREO-Behind (IMPACT, UCBerkeley; PLASTIC, UNH)',fontsize=11, ha='left',color='black')
        
     if sc_label=='Wind': 
        plt.figtext(0.5,0.03,'Data source: Wind (MFI, SWE, NASA/Goddard)',fontsize=11, ha='left',color='black')

     if sc_label=='SolarOrbiter': 
        plt.figtext(0.5,0.03,'Data source: Solar Orbiter (MAG, Imperial College; SWA, University College London)',fontsize=11, ha='left',color='black')
        
     if sc_label=='BepiColombo': 
        plt.figtext(0.5,0.03,'Data source: BepiColombo MPO-MAG ib (IGEP/IWF/ISAS/IC)',fontsize=11, ha='left',color='black')
  
     if sc_label=='Juno': 
        plt.figtext(0.5,0.03,'Data source: Juno (MAG, NASA Goddard)',fontsize=11, ha='left',color='black')

     if sc_label=='ULYSSES': 
        plt.figtext(0.5,0.03,'Data source: Ulysses merged mag/plasma from CDAweb (Imperial College, UK)',fontsize=11, ha='left',color='black')
        
        
    
     logo = plt.imread('logo/GSA_Basislogo_Positiv_RGB_XXS.png')
     newax = fig.add_axes([0.93,0.93,0.06,0.06], anchor='NE', zorder=1)
     newax.imshow(logo)
     newax.axis('off')


        
     plt.tight_layout()
     #plt.show()

     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'

     plotfile=path+ic.icmecat_id[i]+'.png'

     plt.savefig(plotfile)
     #print('saved as ',plotfile)
     plt.close(1)
   
     
     
     
     
     
     
     

def plot_insitu_icmecat_mag(sc, start, end, sc_label, path, ic, i, pos):
     '''
     sc = data    
     '''
     
     start=parse_time(start).datetime
     end=parse_time(end).datetime
     #print(start)
     #print(end)

     sns.set_style('darkgrid')
     sns.set_context('paper')

     fig=plt.figure(figsize=(15,15), dpi=100)
     
     ax1 = plt.subplot(211) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
    
     #plot vertical lines
     ax1.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[-500,500],'-k',linewidth=1)            
     ax1.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[-500,500],'-k',linewidth=1)            
     ax1.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[-500,500],'-k',linewidth=1)            

     
     coord_string='RTN'   
            
     coord_string='RTN'   
    
     if sc_label=='STEREO-B':     coord_string='SCEQ'   
     if sc_label=='MESSENGER':    coord_string='SCEQ'   
     if sc_label=='VEX':          coord_string='SCEQ'   



     ax1.set_ylabel('B [nT] '+coord_string,fontsize=15)

     leg=ax1.legend(loc='lower right',ncol=4,fontsize=15)
     for line in leg.get_lines():
         line.set_linewidth(2.0)   

        
     ax1.set_xlim(start,end)
     bscale=np.nanmax(sc.bt)        
     ax1.set_ylim(-bscale-5,bscale+5)   

     ax1.text(sc.time[30],  bscale, 'west    north', fontsize=10, color='black')
     ax1.text(sc.time[30], -bscale, 'east    south', fontsize=10, color='black')
            
        
        
     ax1.tick_params(axis='x', labelsize=15)  
     ax1.tick_params(axis='y', labelsize=15)    

     ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d %H') )

     plt.title(sc_label+'    ICME start: '+\
               ic.icme_start_time[i].strftime("%Y-%b-%d %H:%M")+ \
               '   MO start: '+ic.mo_start_time[i].strftime("%Y-%b-%d %H:%M")+ \
               '   MO end: '+ic.mo_end_time[i].strftime("%Y-%b-%d %H:%M"),fontsize=15)
        
     if sc_label=='BepiColombo': 
        plt.figtext(0.01,0.01,'Data source: BepiColombo MPO-MAG (IGEP/IWF/ISAS/IC)',fontsize=13, ha='left',color='black')
  
        
     if sc_label=='VEX': 
        plt.figtext(0.01,0.01,'Data source: Venus Express (MAG, IWF/OEAW)',fontsize=13, ha='left',color='black')  

     if sc_label=='MESSENGER': 
        plt.figtext(0.01,0.01,'Data source: MESSENGER (MAG, JHU/APL)',fontsize=13, ha='left',color='black')



     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'

     plt.tight_layout()
    
     ################positions plot
     ax2 = plt.subplot(212,projection='polar') 
     plot_icmecat_positions_mag(ic.icme_start_time[i],'HEEQ',ax2,pos)
    
    
     plotfile=path+ic.icmecat_id[i]+'.png'
       
     plt.savefig(plotfile)
     #print('saved as ',plotfile)
     plt.close(1)
     
     

        
        
def plot_icmecat_events_multi(sc,sci,ic,name,icplotsdir,data_path,pos):

  
    fileind='icmecat/indices_icmecat/ICMECAT_indices_'+name+'.p'
   
    #get indices of events for this spacecrat
    [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))  
    
    #plasma available?
    if name=='Wind': plasma=True
    if name=='STEREO-A': plasma=True
    if name=='STEREO-B': plasma=True
    if name=='ULYSSES': plasma=True
    if name=='MAVEN': 
        plasma=True
        #load MSL rad data
        rad=hd.load_msl_rad(data_path)   
        #load HI ARRCAT
        file='arrcat/HELCATS_ARRCAT_v20_pandas.p'
        [arrcat,arrcat_header]=pickle.load( open(file, 'rb'))           
        #load WSA HUX
        w1=hd.load_mars_wsa_hux()        
        #load Huang SIRCAT
        msir=hd.load_maven_sir_huang()

    if name=='PSP': plasma=True
    if name=='VEX': plasma=False
    if name=='MESSENGER': plasma=False
    if name=='BepiColombo': plasma=True
    if name=='SolarOrbiter': plasma=True              
        
    
    for i in np.arange(np.size(sci)):    
        
        beginind=90*24
    
        if plasma == True:
            
            if name!='MAVEN':
                plot_insitu_icmecat_mag_plasma(sc[icme_start_ind[i]-beginind:mo_end_ind[i]+90*24],\
                             ic.icme_start_time[sci[i]]-datetime.timedelta(days=1.5), \
                             ic.mo_end_time[sci[i]]+datetime.timedelta(days=1.5),name, icplotsdir,ic,sci[i],pos)
            if name == 'MAVEN':                 
                plot_insitu_icmecat_maven(sc[icme_start_ind[i]-7*6:mo_end_ind[i]+7*6],\
                             ic.icme_start_time[sci[i]]-datetime.timedelta(days=7), \
                             ic.mo_end_time[sci[i]]+datetime.timedelta(days=7),\
                             name,icplotsdir,ic,sci[i],rad,arrcat,msir,w1)
                
            plt.close('all')
        else:
            plot_insitu_icmecat_mag(sc[icme_start_ind[i]-beginind:mo_end_ind[i]+90*24], \
                                    ic.icme_start_time[sci[i]]-datetime.timedelta(days=1.5), \
                                    ic.mo_end_time[sci[i]]+datetime.timedelta(days=1.5),name, icplotsdir,ic,sci[i],pos)
            plt.close('all')      
        
        
        
        
        
        

        
        
        
        
        
     
    

def plot_icmecat_events(sc,sci,ic,name,icplotsdir,data_path,pos):

  
    fileind='icmecat/indices_icmecat/ICMECAT_indices_'+name+'.p'
   
    #get indices of events for this spacecrat
    [icme_start_ind, mo_start_ind,mo_end_ind]=pickle.load(open(fileind, 'rb'))  
    
    
        
    #plasma available?
    if name=='Wind': plasma=True
    if name=='STEREO-A': plasma=True
    if name=='STEREO-B': plasma=True
  
    
    if name=='MAVEN': 
        plasma=True
        #load MSL rad data
        rad=hd.load_msl_rad(data_path)   
        #load HI ARRCAT
        file='arrcat/HELCATS_ARRCAT_v20_pandas.p'
        [arrcat,arrcat_header]=pickle.load( open(file, 'rb'))           
        #load WSA HUX
        w1=hd.load_mars_wsa_hux()        
        #load Huang SIRCAT
        msir=hd.load_maven_sir_huang()

    if name=='PSP': plasma=True
    if name=='VEX': plasma=False
    if name=='MESSENGER': plasma=False
    if name=='BepiColombo': plasma=False
    if name=='SolarOrbiter': plasma=True   
  

    if name=='ULYSSES': plasma=True
    if name=='Juno': plasma=True      #so it has the similar style
    
        
    
    for i in np.arange(np.size(sci)):    
        
    
        if plasma == True:                
                
            if name == 'MAVEN':                 
                plot_insitu_icmecat_maven(sc[icme_start_ind[i]-7*6:mo_end_ind[i]+7*6],\
                             ic.icme_start_time[sci[i]]-datetime.timedelta(days=7), \
                             ic.mo_end_time[sci[i]]+datetime.timedelta(days=7),\
                             name,icplotsdir,ic,sci[i],rad,arrcat,msir,w1)
                
            elif name == 'ULYSSES': #different time resolution for ULYSSES            
                plot_insitu_icmecat_mag_plasma(sc[icme_start_ind[i]-300:mo_end_ind[i]+300],\
                             ic.icme_start_time[sci[i]]-datetime.timedelta(days=2), \
                             ic.mo_end_time[sci[i]]+datetime.timedelta(days=2),name, icplotsdir,ic,sci[i],pos)
                
            else:             
                beginind=90*24
                plot_insitu_icmecat_mag_plasma(sc[icme_start_ind[i]-beginind:mo_end_ind[i]+90*24],\
                             ic.icme_start_time[sci[i]]-datetime.timedelta(days=1.5), \
                             ic.mo_end_time[sci[i]]+datetime.timedelta(days=1.5),name, icplotsdir,ic,sci[i],pos)

                
                
                
            plt.close('all')
        else:
            plot_insitu_icmecat_mag(sc[icme_start_ind[i]-beginind:mo_end_ind[i]+90*24], \
                                    ic.icme_start_time[sci[i]]-datetime.timedelta(days=1.5), \
                                    ic.mo_end_time[sci[i]]+datetime.timedelta(days=1.5),name, icplotsdir,ic,sci[i],pos)
            plt.close('all')
    
   
   
     
    
      
     
     
     

def plot_insitu_icmecat_maven(sc, start, end, sc_label, path, ic,i, rad, arrcat, msir, wsa, **kwargs):
    '''
    sc ... data
    '''

    #print(rad[0])
    #print(arrcat_mars[0])
     

    start=parse_time(start).datetime
    end=parse_time(end).datetime
    #print(start)
    #print(end)
        
    #slice msl rad data
    msl_start_ind=np.where(rad.time > start)[0][0]
    msl_end_ind=np.where(rad.time > end)[0][0]
    rad_slice=rad[msl_start_ind:msl_end_ind]
    
    #slice wsa data
    #wsa_start_ind=np.where(wsa.time > start)[0][0]
    #wsa_end_ind=np.where(wsa.time > end)[0][0]
    #wsa_slice=wsa[wsa_start_ind:wsa_end_ind]
            
    #locate all mars events in HI arrcat 
    ac_mars=arrcat.loc[arrcat['target_name'] == 'Mars']
    #get arrival times
    arr=parse_time(ac_mars['target_arrival_time']).plot_date
    #get error in HI arrival time derived from HI speed in hours, convert to days
    err=ac_mars['target_arrival_time_err'].values/24

    

    sns.set_style('darkgrid')
    sns.set_context('paper')

    fig=plt.figure(figsize=(9,8), dpi=150)
    
 

    #sharex means that zooming in works with all subplots
    ax1 = plt.subplot(411) 
    
    plt.title(sc_label+'    ICME start: '+\
               ic.icme_start_time[i].strftime("%Y-%b-%d %H:%M")+ \
               '   MO start: '+ic.mo_start_time[i].strftime("%Y-%b-%d %H:%M")+ \
               '   MO end: '+ic.mo_end_time[i].strftime("%Y-%b-%d %H:%M"))
    
    
    ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
    ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
    ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
    ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.7)
    #plot vertical lines
    ax1.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[-500,500],'-k',linewidth=1)            
    ax1.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[-500,500],'-k',linewidth=1)            
    ax1.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[-500,500],'-k',linewidth=1)  
    
    ax1.plot_date([arr,arr],[-1000,1000],'-b',linewidth=1) 
    #go through all events to plot error bar as shade
    for k  in np.arange(ac_mars.shape[0]):
        ax1.axvspan(arr[k]-err[k],arr[k]+err[k], alpha=0.1, color='blue')
    
    plt.ylabel('B [nT, MSO]')
    leg=ax1.legend(loc='lower right',ncol=4,fontsize=8)
    for line in leg.get_lines():
         line.set_linewidth(2.0)   
    ax1.set_xlim(start,end)
    #if np.isnan(np.nanmin(sc.bt))==False:
    ax1.set_ylim(-np.nanmax(sc.bt)-3,np.nanmax(sc.bt)+3)   
    ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )
    #plt.ylim((-20, 20))
    #ax1.set_xticklabels([]) does not work with sharex
    #plt.setp(ax1.get_xticklabels(), fontsize=6)
    plt.setp(ax1.get_xticklabels(), visible=False)


    
    
    
    ax2 = plt.subplot(412,sharex=ax1) 
    ax2.plot_date(sc.time,sc.vt,'-k',label='MAVEN',linewidth=0.8)
    #plot vertical lines
    ax2.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[0,3000],'-k',linewidth=1)            
    ax2.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[0,3000],'-k',linewidth=1)            
    ax2.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[0,3000],'-k',linewidth=1)                
    
    #over plot wsa model from Martin Reiss
    ax2.plot_date(wsa.time,wsa.vt,'-g',label='WSA/HUX',linewidth=0.8)

    #plot ARRCAT
    ax2.plot_date([arr,arr],[-1000,1000],'-b',linewidth=1) 
    for k  in np.arange(ac_mars.shape[0]):
        ax2.axvspan(arr[k]-err[k],arr[k]+err[k], alpha=0.1, color='blue')

#     #plot Huang SIR MAVEN catalog start end times
#     for i  in np.arange(msir.shape[0]):
#         ax2.axvspan(msir.start[i],msir.end[i], alpha=0.2, color='coral')        
#     #plot Huang SIR MAVEN catalog stream interface
#     for i  in np.arange(msir.shape[0]):
#         ax2.plot([msir.si[i],msir.si[i]],[0,1000], alpha=0.2, color='black')


    plt.legend(loc=3,ncol=2,fontsize=7)        
    ax2.set_ylabel('V [km s$^{-1}]$')
    ax2.set_xlim(start,end)
    #check plasma data exists
    #if np.isnan(np.nanmin(sc.vt))==False:
    ax2.set_ylim(np.nanmin(sc.vt)-50,np.nanmax(sc.vt)+100)   
    ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
    #ax2.set_xticklabels([])
    plt.setp(ax2.get_xticklabels(), visible=False)


    ax3 = plt.subplot(413,sharex=ax1) 
    ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.8)
    ax3.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[0,1000],'-k',linewidth=1)            
    ax3.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[0,1000],'-k',linewidth=1)            
    ax3.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[0,1000],'-k',linewidth=1)            
    ax3.plot_date([arr,arr],[-1000,1000],'-b',linewidth=1) 
    for k  in np.arange(ac_mars.shape[0]):
        ax3.axvspan(arr[k]-err[k],arr[k]+err[k], alpha=0.1, color='blue')

    ax3.set_ylabel('N [cm$^{-3}$]')
    ax3.set_xlim(start,end)
    #if np.isnan(np.nanmin(sc.np))==False:
    ax3.set_ylim(0,np.nanmax(sc.np)+5)   
    ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
    #plt.ylim((0, 50))
    #ax3.set_xticklabels([])
    plt.setp(ax3.get_xticklabels(), visible=False)


    #     ax4 = plt.subplot(514,sharex=ax1) 
    #     ax4.plot_date(sc.time,sc.tp/1e6,'-k',label='Tp',linewidth=0.8)
    #     #plot vertical lines
    #     ax4.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[0,10],'-k',linewidth=1)            
    #     ax4.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[0,10],'-k',linewidth=1)            
    #     ax4.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[0,10],'-k',linewidth=1)    
    #     ax4.plot_date([arr,arr],[-1000,1000],'-b',linewidth=1) 
    #     for k  in np.arange(ac_mars.shape[0]):
    #         ax4.axvspan(arr[k]-err[k],arr[k]+err[k], alpha=0.1, color='blue')
    #     plt.ylabel('T [MK]')
    #     ax4.set_xlim(start,end)
    #     #if np.isnan(np.nanmin(sc.tp))==False:   
    #     ax4.set_ylim(0,np.nanmax(sc.tp/1e6)+0.2)   
    #     plt.setp(ax4.get_xticklabels(), visible=False)

    #plot MSL RAD
    ax5 = plt.subplot(414,sharex=ax1) 
    ax5.plot_date(rad_slice.time,rad_slice.dose_sol,'-k',label='MSL/RAD dose_sol',linewidth=0.8)
    #ax5.plot_date(rad_slice.time,rad_slice.dose_hour,'-r',label='MSL/RAD dose_hour',linewidth=0.7)
    #plot vertical lines
    ax5.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[0,1000],'-k',linewidth=1)            
    ax5.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[0,1000],'-k',linewidth=1)            
    ax5.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[0,1000],'-k',linewidth=1)
    ax5.plot_date([arr,arr],[-1000,1000],'-b',linewidth=1) 
    for k  in np.arange(ac_mars.shape[0]):
        ax5.axvspan(arr[k]-err[k],arr[k]+err[k], alpha=0.1, color='blue')

    #use this fake data for one event for labeling CME arrivals as blue as otherwise there are too many blue vertical lines 
    #from plotting all arrcat events in the legend
    ax5.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[-100,-90],'-b',linewidth=1,label='HI CME arrivals')    
    
    ax5.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
    ax5.set_ylim(np.min(rad_slice.dose_sol)-15,np.max(rad_slice.dose_sol)+15)
    plt.legend(loc=3,ncol=2,fontsize=7)
    plt.ylabel('RAD dose [uGy/day]')

    
    

    plt.tight_layout()
    #plt.show()

    #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'

    plotfile=path+ic.icmecat_id[i]+'.png'


    plt.savefig(plotfile)
    print('saved as ',plotfile)










     
     
     
     
     
     
     
     
     
     
     
     
     

def plot_insitu_measure(sc, start, end, sc_label, path, **kwargs):
     '''
     sc = data
    
     '''
     
     start=parse_time(start).datetime
     end=parse_time(end).datetime
     #print(start)
     #print(end)
     
     #cut out data from full array for faster plotting    
     startind=np.where(sc.time > start)[0][0]   
     endind=np.where(sc.time > end)[0][0]    
        
    
     sc=sc[startind:endind]
        
    
     #t1=parse_time('2010-10-30T09:16').datetime
     #t2=parse_time('2010-10-30T16:00').datetime
     #t3=parse_time('2010-10-30T20:33').datetime
    
    
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
     ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )
     #plt.ylim((-20, 20))
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
     
  
     def onclick(event):

        #global stepout
        #print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
           #   (event.button, event.x, event.y, event.xdata, event.ydata))
        #print('time',str(mdates.num2date(event.xdata)))
        nearestm1=np.argmin(abs(parse_time(sc.time).plot_date-event.xdata))
        seltime=str(parse_time(sc.time[nearestm1]).iso)
        print(seltime[0:10]+'T'+seltime[11:16]+'Z',event.ydata )  
        #if event.button == 3: stepout=True


  
     cid = fig.canvas.mpl_connect('button_press_event', onclick)

    
     

     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     #plt.savefig(plotfile)
     #print('saved as ',plotfile)
   
 

def plot_insitu_measure_mag(sc, start1, end1, sc_label, path, **kwargs):
     '''
     sc = data
    
     '''
     
     start1=parse_time(start1).datetime
     end1=parse_time(end1).datetime
     #make timezone aware
        
     start1=start1.replace(tzinfo=datetime.timezone.utc)
     end1=end1.replace(tzinfo=datetime.timezone.utc)
     print(start1)
     print(end1)


     #cut out data from full array for faster plotting    
     startind=np.where(sc.time > start1)[0][0]   
        
     if sc.time[-1]> end1: 
        endind=np.where(sc.time > end1)[0][0] 
     else: 
        endind=len(sc)
    
    
     sc=sc[startind:endind]
        
    
     #t1=parse_time('2010-10-30T09:16').datetime
     #t2=parse_time('2010-10-30T16:00').datetime
     #t3=parse_time('2010-10-30T20:33').datetime
    
    
     sns.set_style('darkgrid')
     sns.set_context('paper')

     fig=plt.figure(figsize=(9,6), dpi=150)
     
     #sharex means that zooming in works with all subplots
     ax1 = plt.subplot(111) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
     
     plt.ylabel('B [nT]')
     plt.legend(loc=3,ncol=4,fontsize=8)
     ax1.set_xlim(start1,end1)
     ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )
     #plt.ylim((-20, 20))
     #ax1.set_xticklabels([]) does not work with sharex
     #plt.setp(ax1.get_xticklabels(), fontsize=6)

     plt.title(sc_label+' data, start: '+start1.strftime("%Y-%b-%d %H:%M")+'  end: '+end1.strftime("%Y-%b-%d %H:%M"))
     


     
     plt.tight_layout()
     #plt.show()
     
  
     def onclick(event):

        #global stepout
        #print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
           #   (event.button, event.x, event.y, event.xdata, event.ydata))
        #print('time',str(mdates.num2date(event.xdata)))
        nearestm1=np.argmin(abs(parse_time(sc.time).plot_date-event.xdata))
        seltime=str(parse_time(sc.time[nearestm1]).iso)
        print(seltime[0:10]+'T'+seltime[11:16]+'Z',event.ydata )  
        #if event.button == 3: stepout=True


  
     cid = fig.canvas.mpl_connect('button_press_event', onclick)

    
     

     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     #plt.savefig(plotfile)
     #print('saved as ',plotfile)
   
    
   
 
def plot_insitu_measure_mag_notz(sc, start1, end1, sc_label, path, **kwargs):
     '''
     sc = data
    
     '''
     
     start1=parse_time(start1).datetime
     end1=parse_time(end1).datetime
       


     #cut out data from full array for faster plotting    
     startind=np.where(sc.time > start1)[0][0]   
        
     if sc.time[-1]> end1: 
        endind=np.where(sc.time > end1)[0][0] 
     else: 
        endind=len(sc)
    
    
     sc=sc[startind:endind]
        
    
     #t1=parse_time('2010-10-30T09:16').datetime
     #t2=parse_time('2010-10-30T16:00').datetime
     #t3=parse_time('2010-10-30T20:33').datetime
    
    
     sns.set_style('darkgrid')
     sns.set_context('paper')

     fig=plt.figure(figsize=(9,6), dpi=150)
     
     #sharex means that zooming in works with all subplots
     ax1 = plt.subplot(111) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
     
     plt.ylabel('B [nT]')
     plt.legend(loc=3,ncol=4,fontsize=8)
     ax1.set_xlim(start1,end1)
     ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )
     #plt.ylim((-20, 20))
     #ax1.set_xticklabels([]) does not work with sharex
     #plt.setp(ax1.get_xticklabels(), fontsize=6)

     plt.title(sc_label+' data, start: '+start1.strftime("%Y-%b-%d %H:%M")+'  end: '+end1.strftime("%Y-%b-%d %H:%M"))
     


     
     plt.tight_layout()
     #plt.show()
     
  
     def onclick(event):

        #global stepout
        #print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
           #   (event.button, event.x, event.y, event.xdata, event.ydata))
        #print('time',str(mdates.num2date(event.xdata)))
        nearestm1=np.argmin(abs(parse_time(sc.time).plot_date-event.xdata))
        seltime=str(parse_time(sc.time[nearestm1]).iso)
        print(seltime[0:10]+'T'+seltime[11:16]+'Z',event.ydata )  
        #if event.button == 3: stepout=True


  
     cid = fig.canvas.mpl_connect('button_press_event', onclick)

    
     

     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     #plt.savefig(plotfile)
     #print('saved as ',plotfile)
   
              
   
         

def plot_insitu_measure_maven(sc, start, end, sc_label, path, rad, arrcat,wsa,msir,**kwargs):
    '''
    sc ... data
    '''

    #print(rad[0])
    #print(arrcat_mars[0])
     

    start=parse_time(start).datetime
    end=parse_time(end).datetime
    #print(start)
    #print(end)
        
    #slice msl rad data 
    msl_start_ind=np.where(rad.time > start)[0][0]
    msl_end_ind=np.where(rad.time > end)[0][0]
    rad_slice=rad[msl_start_ind:msl_end_ind]
    
    #locate all mars events in HI arrcat 
    ac_mars=arrcat.loc[arrcat['target_name'] == 'Mars']
    #get arrival times
    arr=parse_time(ac_mars['target_arrival_time']).plot_date
    #get error in HI arrival time derived from HI speed in hours, convert to days
    err=ac_mars['target_arrival_time_err'].values/24

    
    #################### make figure

    sns.set_style('darkgrid')
    sns.set_context('paper')

    fig=plt.figure(figsize=(8,8), dpi=150)
    
    plt.title(sc_label+' data, start: '+start.strftime("%Y-%b-%d %H:%M")+'  end: '+end.strftime("%Y-%b-%d %H:%M"))


    #sharex means that zooming in works with all subplots
    ax1 = plt.subplot(511) 
    ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
    ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
    ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
    ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
    
    ax1.plot_date([arr,arr],[-100,100],'-b',linewidth=1) 
    #go through all events to plot error bar as shade
    for i  in np.arange(ac_mars.shape[0]):
        ax1.axvspan(arr[i]-err[i],arr[i]+err[i], alpha=0.3, color='blue')
        
        
        
    #plot Huang SIR MAVEN catalog start end times
    for i  in np.arange(msir.shape[0]):
        ax1.axvspan(msir.start[i],msir.end[i], alpha=0.2, color='coral')        
    #plot Huang SIR MAVEN catalog stream interface
    for i  in np.arange(msir.shape[0]):
        ax1.plot([msir.si[i],msir.si[i]],[-100,100], alpha=0.2, color='black')

    plt.ylabel('B [nT]')
    plt.legend(loc=3,ncol=4,fontsize=8)
    ax1.set_xlim(start,end)
    #if np.isnan(np.nanmin(sc.bt))==False:
    ax1.set_ylim(-np.nanmax(sc.bt)-3,np.nanmax(sc.bt)+3)   
    ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )
    plt.ylim((-20, 20))
    #ax1.set_xticklabels([]) does not work with sharex
    #plt.setp(ax1.get_xticklabels(), fontsize=6)
    plt.setp(ax1.get_xticklabels(), visible=False)
    
    


    ax2 = plt.subplot(512,sharex=ax1) 
    ax2.plot_date(sc.time,sc.vt,'-k',label='V MAVEN',linewidth=0.7)
    #over plot wsa model from Martin Reiss
    ax2.plot_date(wsa.time,wsa.vt,'-g',label='V WSA',linewidth=0.7)

    #plot ARRCAT
    ax2.plot_date([arr,arr],[-1000,1000],'-b',linewidth=1) 
    for i  in np.arange(ac_mars.shape[0]):
        ax2.axvspan(arr[i]-err[i],arr[i]+err[i], alpha=0.3, color='blue')
        
    #plot Huang SIR MAVEN catalog start end times
    for i  in np.arange(msir.shape[0]):
        ax2.axvspan(msir.start[i],msir.end[i], alpha=0.2, color='coral')        
    #plot Huang SIR MAVEN catalog stream interface
    for i  in np.arange(msir.shape[0]):
        ax2.plot([msir.si[i],msir.si[i]],[0,1000], alpha=0.2, color='black')

    plt.ylabel('V [km/s]')
    plt.legend(loc=3,ncol=2,fontsize=8)
    ax2.set_xlim(start,end)
    #check plasma data exists
    #if np.isnan(np.nanmin(sc.vt))==False:
    ax2.set_ylim(250,750)   
    ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
    #plt.ylim((250, 800))
    #ax2.set_xticklabels([])
    plt.setp(ax2.get_xticklabels(), visible=False)


    

    ax3 = plt.subplot(513,sharex=ax1) 
    ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)
    ax3.plot_date([arr,arr],[-1000,1000],'-b',linewidth=1) 
    for i  in np.arange(ac_mars.shape[0]):
        ax3.axvspan(arr[i]-err[i],arr[i]+err[i], alpha=0.3, color='blue')
    plt.ylabel('N [ccm-3]')
    ax3.set_xlim(start,end)
    #if np.isnan(np.nanmin(sc.np))==False:
    ax3.set_ylim(0,np.nanmax(sc.np)+5)   
    ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
    #plt.ylim((0, 50))
    #ax3.set_xticklabels([])

    #plot Huang SIR MAVEN catalog start end times
    for i  in np.arange(msir.shape[0]):
        ax3.axvspan(msir.start[i],msir.end[i], alpha=0.2, color='coral')        
    #plot Huang SIR MAVEN catalog stream interface
    for i  in np.arange(msir.shape[0]):
        ax3.plot([msir.si[i],msir.si[i]],[0,1000], alpha=0.2, color='black')

    plt.setp(ax3.get_xticklabels(), visible=False)


    ax4 = plt.subplot(514,sharex=ax1) 
    ax4.plot_date(sc.time,sc.tp/1e6,'-k',label='Tp',linewidth=0.7)
    ax4.plot_date([arr,arr],[-1000,1000],'-b',linewidth=1) 
    for i  in np.arange(ac_mars.shape[0]):
        ax4.axvspan(arr[i]-err[i],arr[i]+err[i], alpha=0.3, color='blue')
    plt.ylabel('T [MK]')
    ax4.set_xlim(start,end)
    #if np.isnan(np.nanmin(sc.tp))==False:   
    ax4.set_ylim(0,1)   
    plt.setp(ax4.get_xticklabels(), visible=False)
    
    #plot Huang SIR MAVEN catalog start end times
    for i  in np.arange(msir.shape[0]):
        ax4.axvspan(msir.start[i],msir.end[i], alpha=0.2, color='coral')        
    #plot Huang SIR MAVEN catalog stream interface
    for i  in np.arange(msir.shape[0]):
        ax4.plot([msir.si[i],msir.si[i]],[0,1000], alpha=0.2, color='black')
   
    
    
    
    

    #plot MSL RAD and ARRCAT
    ax5 = plt.subplot(515,sharex=ax1) 
    ax5.plot_date(rad.time,rad.dose_sol,'-k',label='MSL/RAD dose_sol',linewidth=0.7)
    #ax5.plot_date(rad_slice.time,rad_slice.dose_sol,'-k',label='MSL/RAD dose_sol',linewidth=0.7)
    
    ax5.plot_date([arr,arr],[-1000,1000],'-b',linewidth=1) 
    for i  in np.arange(ac_mars.shape[0]):
        ax5.axvspan(arr[i]-err[i],arr[i]+err[i], alpha=0.3, color='blue')
        
        
    #plot Huang SIR MAVEN catalog start end times
    for i  in np.arange(msir.shape[0]):
        ax5.axvspan(msir.start[i],msir.end[i], alpha=0.2, color='coral')        
    #plot Huang SIR MAVEN catalog stream interface
    for i  in np.arange(msir.shape[0]):
        ax5.plot([msir.si[i],msir.si[i]],[0,1000], alpha=0.2, color='black')


    
    #ax5.plot_date(rad_slice.time,rad_slice.dose_hour,'-r',label='MSL/RAD dose_hour',linewidth=0.7)
    plt.ylabel('RAD dose_sol [uGy/day]')
    ax5.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
    ax5.set_ylim(np.min(rad_slice.dose_sol)-15,np.max(rad_slice.dose_sol)+15)
    #plt.legend(loc=1,ncol=3,fontsize=6)
       
    

    plt.tight_layout()
    #plt.show()

    #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
    #plotfile=path+ic.icmecat_id[i]+'.png'
    #plt.savefig(plotfile)
    #print('saved as ',plotfile)



  
    def onclick(event):

        #global stepout
        #print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
           #   (event.button, event.x, event.y, event.xdata, event.ydata))
        #print('time',str(mdates.num2date(event.xdata)))
        nearestm1=np.argmin(abs(parse_time(sc.time).plot_date-event.xdata))
        seltime=str(parse_time(sc.time[nearestm1]).iso)
        print(seltime[0:10]+'T'+seltime[11:16]+'Z',event.ydata )  
        #if event.button == 3: stepout=True


  
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    
     

     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     #plt.savefig(plotfile)
     #print('saved as ',plotfile)
   
        
     
     
        
        
        
        
        
     

def plot_positions(time_date1, path,frame, **kwargs):
    '''
    '''
    
    sns.set_style('darkgrid')
    #sns.set_style('whitegrid')
    sns.set_context('paper')    
    
    time1=mdates.date2num(time_date1)

    #old version
    #[psp, solo, sta, stb, bepi, l1, earth, mercury, venus, mars, jupiter, saturn, uranus, neptune]=pickle.load( open( 'results/positions/positions_psp_solo_sta_bepi_wind_planets_HEEQ_10min_rad.p', "rb" ) )

    #new version
    #made with positions.ipynb 10 min, in rad, matplotlib datenumber
    [psp, bepi, solo, sta, juice, earth, mercury, venus, mars, jupiter, saturn, uranus, neptune]=pickle.load( open( 'results/positions/positions_2020_all_HEEQ_10min_rad_cm.p', "rb" ) )    
    res_in_days=1/(24*6) #10 min resolution of positions file
    
 
    #sidereal solar rotation rate
    if frame=='HCI': sun_rot=24.47
    #synodic
    if frame=='HEEQ': sun_rot=26.24
   
    AUkm=149597870.7   

    #for parker spiral   
    theta=np.arange(0,np.deg2rad(180),0.01)

    fadeind=int(80/res_in_days)

    #fadeind=int(100/res_in_days)
    k=0
    
    plot_orbit=True
    plot_parker=True
    fsize=17 
    symsize_planet=140
    symsize_spacecraft=100
    
    #find current indices
    dct=time1-psp.time
    psp_timeind=np.argmin(abs(dct))

    dct=time1-bepi.time
    bepi_timeind=np.argmin(abs(dct))

    dct=time1-solo.time
    solo_timeind=np.argmin(abs(dct))
    
    dct=time1-juice.time
    juice_timeind=np.argmin(abs(dct))


    dct=time1-earth.time
    earth_timeind=np.argmin(abs(dct))
    
    #dct=time1-l1.time
    #l1_timeind=np.argmin(abs(dct))
    
    dct=time1-venus.time
    venus_timeind=np.argmin(abs(dct))
    
    dct=time1-mars.time
    mars_timeind=np.argmin(abs(dct))

    dct=time1-jupiter.time
    jupiter_timeind=np.argmin(abs(dct))

    
    dct=time1-mercury.time
    mercury_timeind=np.argmin(abs(dct))

    dct=time1-sta.time
    sta_timeind=np.argmin(abs(dct))
    
    
    
    
    #dct=time1-stb.time
    #stb_timeind=np.argmin(abs(dct))


    
    ################## figure    
    fig=plt.figure(1, figsize=(15,10), dpi=150)
    ax = plt.subplot(111, projection='polar')
    backcolor='black'
    psp_color='black'
    bepi_color='blue'
    solo_color='green'
    juice_color='gold'

    ax.scatter(venus.lon[venus_timeind], venus.r[venus_timeind]*np.cos(venus.lat[venus_timeind]), s=symsize_planet, c='orange', alpha=1,lw=0,zorder=3)
    ax.scatter(mercury.lon[mercury_timeind], mercury.r[mercury_timeind]*np.cos(mercury.lat[mercury_timeind]), s=symsize_planet, c='dimgrey', alpha=1,lw=0,zorder=3)
    ax.scatter(earth.lon[earth_timeind], earth.r[earth_timeind]*np.cos(earth.lat[earth_timeind]), s=symsize_planet, c='mediumseagreen', alpha=1,lw=0,zorder=3)
    ax.scatter(sta.lon[sta_timeind], sta.r[sta_timeind]*np.cos(sta.lat[sta_timeind]), s=symsize_spacecraft, c='red', marker='s', alpha=1,lw=0,zorder=3)
    ax.scatter(mars.lon[mars_timeind], mars.r[mars_timeind]*np.cos(mars.lat[mars_timeind]), s=symsize_planet, c='orangered', alpha=1,lw=0,zorder=3)

    plt.text(sta.lon[sta_timeind]+0.1,sta.r[sta_timeind],'STEREO-A', color='red', ha='center',fontsize=fsize-4,verticalalignment='top')
  
    plt.text(0,0,'Sun', color='black', ha='center',fontsize=fsize-5,verticalalignment='top')
    plt.text(0,earth.r[earth_timeind]+0.08,'Earth', color='mediumseagreen', ha='center',fontsize=fsize-5,verticalalignment='center')
    
    plt.figtext(0.01,0.01,'Austrian Space Weather Office   GeoSphere Austria', color='black', ha='left',fontsize=fsize-4, style='italic')
    plt.figtext(0.99,0.01,'helioforecast.space', color='black', ha='right',fontsize=fsize-4, style='italic')

    plt.figtext(0.80,0.13,'─   80 days future trajectory', color='black', ha='left',fontsize=fsize-3)
    plt.figtext(0.80,0.1 ,'- -  80 days past trajectory', color='black', ha='left',fontsize=fsize-3)
  
    
    '''
    plt.figtext(0.95,0.5,'Wind', color='mediumseagreen', ha='center',fontsize=fsize+3)
    plt.figtext(0.95,0.25,'STEREO-A', color='red', ha='center',fontsize=fsize+3)
    plt.figtext(0.9,0.9,'Mercury', color='dimgrey', ha='center',fontsize=fsize+5)
    plt.figtext(0.9	,0.8,'Venus', color='orange', ha='center',fontsize=fsize+5)
    plt.figtext(0.9,0.7,'Earth', color='mediumseagreen', ha='center',fontsize=fsize+5)
    #plt.figtext(0.9,0.7,'Mars', color='orangered', ha='center',fontsize=fsize+5)
    plt.figtext(0.9,0.6,'STEREO-A', color='red', ha='center',fontsize=fsize+5)
    plt.figtext(0.9,0.5,'Parker Solar Probe', color='black', ha='center',fontsize=fsize+5)
    plt.figtext(0.9,0.4,'Bepi Colombo', color='blue', ha='center',fontsize=fsize+5)
    plt.figtext(0.9,0.3,'Solar Orbiter', color='green', ha='center',fontsize=fsize+5)

    '''


    f10=plt.figtext(0.01,0.93,'              R     lon     lat', fontsize=fsize+2, ha='left',color=backcolor)

    
    if psp_timeind > 0:
        #plot trajectory
        ax.scatter(psp.lon[psp_timeind], psp.r[psp_timeind]*np.cos(psp.lat[psp_timeind]), s=symsize_spacecraft, c=psp_color, marker='s', alpha=1,lw=0,zorder=3)
        #plot position as text
        psp_text='PSP:   '+str(f'{psp.r[psp_timeind]:6.2f}')+str(f'{np.rad2deg(psp.lon[psp_timeind]):8.1f}')+str(f'{np.rad2deg(psp.lat[psp_timeind]):8.1f}')
        f5=plt.figtext(0.01,0.78,psp_text, fontsize=fsize, ha='left',color=psp_color)
        plt.text(psp.lon[psp_timeind]-0.15,psp.r[psp_timeind],'Parker Solar Probe', color='black', ha='center',fontsize=fsize-4,verticalalignment='top')

        if plot_orbit: 
            ax.plot(psp.lon[psp_timeind:psp_timeind+fadeind], psp.r[psp_timeind:psp_timeind+fadeind]*np.cos(psp.lat[psp_timeind:psp_timeind+fadeind]), c=psp_color, alpha=0.6,lw=1,zorder=3)
            ax.plot(psp.lon[psp_timeind-fadeind:psp_timeind], psp.r[psp_timeind-fadeind:psp_timeind]*np.cos(psp.lat[psp_timeind-fadeind:psp_timeind]), c=psp_color, linestyle='--', alpha=0.4,lw=1,zorder=3)


    if bepi_timeind > 0:
        ax.scatter(bepi.lon[bepi_timeind], bepi.r[bepi_timeind]*np.cos(bepi.lat[bepi_timeind]), s=symsize_spacecraft, c=bepi_color, marker='s', alpha=1,lw=0,zorder=3)
        bepi_text='Bepi:   '+str(f'{bepi.r[bepi_timeind]:6.2f}')+str(f'{np.rad2deg(bepi.lon[bepi_timeind]):8.1f}')+str(f'{np.rad2deg(bepi.lat[bepi_timeind]):8.1f}')
        f6=plt.figtext(0.01,0.74,bepi_text, fontsize=fsize, ha='left',color=bepi_color)
        plt.text(bepi.lon[bepi_timeind]-0.15,bepi.r[bepi_timeind]+0.08,'Bepi Colombo', color='blue', ha='center',fontsize=fsize-4,verticalalignment='top')
  
        if plot_orbit: 
            ax.plot(bepi.lon[bepi_timeind:bepi_timeind+fadeind], bepi.r[bepi_timeind:bepi_timeind+fadeind]*np.cos(bepi.lat[bepi_timeind:bepi_timeind+fadeind]), c=bepi_color, alpha=0.6,lw=1,zorder=3)
            ax.plot(bepi.lon[bepi_timeind-fadeind:bepi_timeind], bepi.r[bepi_timeind-fadeind:bepi_timeind]*np.cos(bepi.lat[bepi_timeind-fadeind:bepi_timeind]), c=bepi_color, linestyle='--', alpha=0.4,lw=1,zorder=3)



    if solo_timeind > 0:
        ax.scatter(solo.lon[solo_timeind], solo.r[solo_timeind]*np.cos(solo.lat[solo_timeind]), s=symsize_spacecraft, c=solo_color, marker='s', alpha=1,lw=0,zorder=3)
        solo_text='SolO:  '+str(f'{solo.r[solo_timeind]:6.2f}')+str(f'{np.rad2deg(solo.lon[solo_timeind]):8.1f}')+str(f'{np.rad2deg(solo.lat[solo_timeind]):8.1f}')
        f7=plt.figtext(0.01,0.7,solo_text, fontsize=fsize, ha='left',color=solo_color)
        plt.text(solo.lon[solo_timeind]-0.15,solo.r[solo_timeind],'Solar Orbiter', color='green', ha='center',fontsize=fsize-4,verticalalignment='top')

        if plot_orbit: 
            ax.plot(solo.lon[solo_timeind:solo_timeind+fadeind], solo.r[solo_timeind:solo_timeind+fadeind]*np.cos(solo.lat[solo_timeind:solo_timeind+fadeind]), c=solo_color, alpha=0.6,lw=1,zorder=3)
            ax.plot(solo.lon[solo_timeind-fadeind:solo_timeind], solo.r[solo_timeind-fadeind:solo_timeind]*np.cos(solo.lat[solo_timeind-fadeind:solo_timeind]), c=solo_color, linestyle='--',alpha=0.5,lw=1,zorder=3)

            
           
            
    if juice_timeind > 0:
        ax.scatter(juice.lon[juice_timeind], juice.r[juice_timeind]*np.cos(juice.lat[juice_timeind]), s=symsize_spacecraft, c=juice_color, marker='s', alpha=1,lw=0,zorder=3)
        juice_text='JUICE:  '+str(f'{juice.r[juice_timeind]:6.2f}')+str(f'{np.rad2deg(juice.lon[juice_timeind]):8.1f}')+str(f'{np.rad2deg(juice.lat[juice_timeind]):8.1f}')
        f7=plt.figtext(0.01,0.66,juice_text, fontsize=fsize, ha='left',color=juice_color)
        if plot_orbit: 
            fadestart=juice_timeind-fadeind
            if  fadestart < 0: fadestart=0            
            ax.plot(juice.lon[fadestart:juice_timeind+fadeind], juice.r[fadestart:juice_timeind+fadeind]*np.cos(juice.lat[fadestart:juice_timeind+fadeind]), c=juice_color, alpha=0.6,lw=1,zorder=3)

             
            
            
            
            
            
            
            

    if plot_orbit: 
        ax.plot(sta.lon[sta_timeind:sta_timeind+fadeind], sta.r[sta_timeind:sta_timeind+fadeind]*np.cos(sta.lat[sta_timeind:sta_timeind+fadeind]), c='red', alpha=0.6,lw=1,zorder=3)
        ax.plot(sta.lon[sta_timeind-fadeind:sta_timeind], sta.r[sta_timeind-fadeind:sta_timeind]*np.cos(sta.lat[sta_timeind-fadeind:sta_timeind]), c='red', linestyle='--', alpha=0.5,lw=1,zorder=3)

        
        
        
    if time1<mdates.date2num(datetime.datetime(2014,11,26)):        
        vex_text='VEX:   '+str(f'{venus.r[venus_timeind]:6.2f}')+str(f'{np.rad2deg(venus.lon[venus_timeind]):8.1f}')+str(f'{np.rad2deg(venus.lat[venus_timeind]):8.1f}')
        f11=plt.figtext(0.01,0.78,vex_text, fontsize=fsize, ha='left',color='orange')
        plt.text(venus.lon[venus_timeind],venus.r[venus_timeind]+0.12,'VEX', color='orange', ha='center',fontsize=fsize-5,verticalalignment='center')

    #if time1<mdates.date2num(datetime.datetime(2011,3,18))):        
    #    mes_text='MES:   '+str(f'{mercury.r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(mercury.lon[earth_timeind]):8.1f}')+str(f'{np.rad2deg(mercury.lat[earth_timeind]):8.1f}')
    #    f12=plt.figtext(0.01,0.78,mes_text, fontsize=fsize, ha='left',color='darkgrey')
    #    plt.text(mercury.lon[earth_timeind],mercury.r[earth_timeind]+0.12,'MESSENGER', color='darkgrey', ha='center',fontsize=fsize-5,verticalalignment='center')
        
        
    if np.logical_and(time1>mdates.date2num(datetime.datetime(2011,3,18)), time1<mdates.date2num(datetime.datetime(2015,4,30))):        
        mes_text='MES:   '+str(f'{mercury.r[mercury_timeind]:6.2f}')+str(f'{np.rad2deg(mercury.lon[mercury_timeind]):8.1f}')+str(f'{np.rad2deg(mercury.lat[earth_timeind]):8.1f}')
        f12=plt.figtext(0.01,0.74,mes_text, fontsize=fsize, ha='left',color='darkgrey')
        plt.text(mercury.lon[mercury_timeind],mercury.r[mercury_timeind]+0.12,'MESSENGER', color='darkgrey', ha='center',fontsize=fsize-5,verticalalignment='center')






    if frame=='HEEQ': earth_text='Earth: '+str(f'{earth.r[earth_timeind]:6.2f}')+str(f'{0.0:8.1f}')+str(f'{np.rad2deg(earth.lat[earth_timeind]):8.1f}')
    else: earth_text='Earth: '+str(f'{earth.r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(earth.lon[earth_timeind]):8.1f}')+str(f'{np.rad2deg(earth.lat[earth_timeind]):8.1f}')

    mars_text='Mars:  '+str(f'{mars.r[mars_timeind]:6.2f}')+str(f'{np.rad2deg(mars.lon[mars_timeind]):8.1f}')+str(f'{np.rad2deg(mars.lat[mars_timeind]):8.1f}')
    sta_text='STA:   '+str(f'{sta.r[sta_timeind]:6.2f}')+str(f'{np.rad2deg(sta.lon[sta_timeind]):8.1f}')+str(f'{np.rad2deg(sta.lat[sta_timeind]):8.1f}')

    #Sun
    ax.scatter(0,0,s=100,c='yellow',alpha=1, edgecolors='black', linewidth=0.3)
    
    
    #add longitudes for Mars and Jupiter
    ax.scatter(mars.lon[mars_timeind], 1.1, s=symsize_spacecraft, c='red', marker=6, alpha=0.8,lw=0,zorder=3)
    plt.text(mars.lon[mars_timeind],1.0,'to Mars', color='red', ha='center',fontsize=fsize-4,verticalalignment='top',alpha=0.8)

    ax.scatter(jupiter.lon[jupiter_timeind], 1.1, s=symsize_spacecraft, c='darkgrey', marker=6, alpha=0.8,lw=0,zorder=3)
    plt.text(jupiter.lon[jupiter_timeind],1.0,'to Jupiter', color='darkgrey', ha='center',fontsize=fsize-4,verticalalignment='top')

    #ax.annotate('', xy=(jupiter.lon[jupiter_timeind], 1.1), xytext=(0, 0),  arrowprops=dict(arrowstyle='->', lw=0), color='red')

 
    f10=plt.figtext(0.01,0.9,earth_text, fontsize=fsize, ha='left',color='mediumseagreen')
    f9=plt.figtext(0.01,0.86,mars_text, fontsize=fsize, ha='left',c='orangered')
    f8=plt.figtext(0.01,0.82,sta_text, fontsize=fsize, ha='left',c='red')
    
    #time
    plt.figtext(0.95,0.9,time_date1.strftime("%Y %B %d  %H:%M UT"),fontsize=fsize+4, ha='right',c='black')

    #parker spiral
    if plot_parker:
        for q in np.arange(0,12):
            #parker spiral
            #sidereal rotation
            omega=2*np.pi/(sun_rot*60*60*24) #solar rotation in seconds
            v=400/AUkm #km/s
            r0=695000/AUkm
            r=v/omega*theta+r0*7
            ax.plot(-theta+np.deg2rad(0+(360/24.47)*res_in_days*k+360/12*q), r, alpha=0.4, lw=0.5,color='grey',zorder=2)
     

    #set axes

    ax.set_theta_zero_location('E')
 
    #plt.rgrids((0.10,0.39,0.72,1.00,1.52),('0.10','0.39','0.72','1.0','1.52 AU'),angle=125, fontsize=fsize,alpha=0.9, color=backcolor)
    plt.rgrids((0.1,0.3,0.5,0.7,1.0,1.3,1.6),('0.1','0.3','0.5','0.7','1.0','1.3','1.6 AU'),angle=125, fontsize=fsize-2,alpha=0.4, color=backcolor)

    #ax.set_ylim(0, 1.75) with Mars
    ax.set_ylim(0, 1.15) 

    plt.thetagrids(range(0,360,45),(u'0\u00b0 '+frame+' longitude',u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'- 135\u00b0',u'- 90\u00b0',u'- 45\u00b0'), fmt='%d',ha='left',fontsize=fsize,color=backcolor, zorder=5, alpha=0.9)

    plt.tight_layout()

    #plotfile=path+%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
    plotfile=path+'positions_'+time_date1.strftime("%Y_%b_%d")+'.png'

    plt.savefig(plotfile)
    print('saved as ',plotfile)
    
    
    #write txt file with positions 30 days prior to and later than current time
    #and pickle file
    
    #psp, solo, sta, stb, bepi, l1, earth, mercury, venus, mars, jupiter, saturn, uranus, neptune
    #10 min resolution
    psp_cut=psp[psp_timeind-6*24*30:psp_timeind+6*24*30]
    solo_cut=solo[solo_timeind-6*24*30:solo_timeind+6*24*30]
    sta_cut=sta[sta_timeind-6*24*30:sta_timeind+6*24*30]
    bepi_cut=bepi[bepi_timeind-6*24*30:bepi_timeind+6*24*30]
    #l1_cut=l1[l1_timeind-6*24*30:l1_timeind+6*24*30]
    earth_cut=psp[earth_timeind-6*24*30:earth_timeind+6*24*30]
    mercury_cut=mercury[mercury_timeind-6*24*30:mercury_timeind+6*24*30]
    venus_cut=venus[venus_timeind-6*24*30:venus_timeind+6*24*30]
    mars_cut=mars[mars_timeind-6*24*30:mars_timeind+6*24*30]
    jupiter_cut=jupiter[jupiter_timeind-6*24*30:jupiter_timeind+6*24*30]
    
    
    
    #with l1
    #pickle.dump([psp_cut,solo_cut,sta_cut,bepi_cut,l1_cut, earth_cut, mercury_cut, venus_cut, mars_cut, jupiter_cut], open( path+'positions_now.p', "wb" ) ) 
    pickle.dump([psp_cut,solo_cut,sta_cut,bepi_cut, earth_cut, mercury_cut, venus_cut, mars_cut, jupiter_cut], open( path+'positions_now.p', "wb" ) ) 

    print('saved as ',path+'positions_now.p')
    
    
    #print(np.concatenate((psp_cut, solo_cut, sta_cut, bepi_cut, l1_cut), axis=0))

    #make adjustments for txt file output
    output_format='%Y-%m-%dT%H:%MZ'
    psp_cut.time=[mdates.num2date(ts).strftime(output_format) for ts in psp_cut.time]
    psp_cut.lon=np.rad2deg(psp_cut.lon)
    psp_cut.lat=np.rad2deg(psp_cut.lat)
    
    solo_cut.time=[mdates.num2date(ts).strftime(output_format) for ts in solo_cut.time]
    solo_cut.lon=np.rad2deg(solo_cut.lon)
    solo_cut.lat=np.rad2deg(solo_cut.lat)
    
    bepi_cut.time=[mdates.num2date(ts).strftime(output_format) for ts in bepi_cut.time]
    bepi_cut.lon=np.rad2deg(bepi_cut.lon)
    bepi_cut.lat=np.rad2deg(bepi_cut.lat)

    sta_cut.time=[mdates.num2date(ts).strftime(output_format) for ts in sta_cut.time]
    sta_cut.lon=np.rad2deg(sta_cut.lon)
    sta_cut.lat=np.rad2deg(sta_cut.lat)

    #l1_cut.time=[mdates.num2date(ts).strftime(output_format) for ts in l1_cut.time]
    #l1_cut.lon=np.rad2deg(l1_cut.lon)
    #l1_cut.lat=np.rad2deg(l1_cut.lat)

    
    #with L1
    #np.savetxt(path+'positions_now.txt', np.concatenate((psp_cut,solo_cut,bepi_cut,sta_cut,l1_cut), axis=0), delimiter=' ', fmt='%s %s %f %f %f %f %f %f', header='spacecraft time     R [AU] lon [deg] lat [deg]   x [AU]   y [AU]  z [AU]  HEEQ coordinates / ASWO, GeoSphere Austria created '+str(datetime.datetime.utcnow())[0:16])

    
    np.savetxt(path+'positions_now.txt', np.concatenate((psp_cut,solo_cut,bepi_cut,sta_cut), axis=0), delimiter=' ', fmt='%s %s %f %f %f %f %f', header='spacecraft time     R [AU] lon [deg] lat [deg]   x [AU]   y [AU]  z [AU]  HEEQ coordinates / ASWO, GeoSphere Austria created '+str(datetime.datetime.utcnow())[0:16])

#    np.savetxt(path+'current_positions.txt', np.concatenate((psp_cut, solo_cut, sta_cut, bepi_cut, l1_cut), axis=0), delimiter=' ', fmt='%s %f %f %f %f %f %f %f ')
    print('saved as ',path+'positions_now.txt')

    

    #if now exists as keyword, save as the file with just now in filename:     
    if 'pdf' in kwargs:
        plotfile=path+'positions_pdf.pdf'
        plt.savefig(plotfile)
        print('saved as ',plotfile)

    

    #if now exists as keyword, save as the file with just now in filename:     
    if 'now' in kwargs:
        plotfile=path+'positions_now.png'
        plt.savefig(plotfile)
        print('saved as ',plotfile)

    if 'eps' in kwargs:
        plotfile=path+'positions_'+time_date1.strftime("%Y_%b_%d")+'.eps'
        plt.savefig(plotfile)
        print('saved as ',plotfile)

    plt.show()
    plt.close('all')

     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
    
def make_positions(time1,frame):
    'makes spacecraft positions for input time'
    'does not work yet'


    time1_num=mdates.date2num(time1)    

    spice.furnish(spicedata.get_kernel('psp_pred'))
    psp=spice.Trajectory('SPP')
    psp.generate_positions(time1_num,'Sun',frame)
    psp.change_units(astropy.units.AU)  
    [psp_r, psp_lat, psp_lon]=cart2sphere(psp.x,psp.y,psp.z)

    spice.furnish(spicedata.get_kernel('bepi_pred'))
    bepi=spice.Trajectory('BEPICOLOMBO MPO') # or BEPICOLOMBO MMO
    bepi.generate_positions(time1_num,'Sun',frame)
    bepi.change_units(astropy.units.AU)  
    [bepi_r, bepi_lat, bepi_lon]=cart2sphere(bepi.x,bepi.y,bepi.z)

    spice.furnish(spicedata.get_kernel('solo_2020'))
    solo=spice.Trajectory('Solar Orbiter')
    solo.generate_positions(time1_num, 'Sun',frame)
    solo.change_units(astropy.units.AU)
    [solo_r, solo_lat, solo_lon]=cart2sphere(solo.x,solo.y,solo.z)


    planet_kernel=spicedata.get_kernel('planet_trajectories')

    earth=spice.Trajectory('399')  #399 for Earth, not barycenter (because of moon)
    earth.generate_positions(time1_num,'Sun',frame)
    earth.change_units(astropy.units.AU)  
    [earth_r, earth_lat, earth_lon]=cart2sphere(earth.x,earth.y,earth.z)

    mercury=spice.Trajectory('1')  #barycenter
    mercury.generate_positions(time1_num,'Sun',frame)  
    mercury.change_units(astropy.units.AU)  
    [mercury_r, mercury_lat, mercury_lon]=cart2sphere(mercury.x,mercury.y,mercury.z)

    venus=spice.Trajectory('2')  
    venus.generate_positions(time1_num,'Sun',frame)  
    venus.change_units(astropy.units.AU)  
    [venus_r, venus_lat, venus_lon]=cart2sphere(venus.x,venus.y,venus.z)

    mars_time_num=earth_time_num
    mars=spice.Trajectory('4')  
    mars.generate_positions(time1_num,'Sun',frame)  
    mars.change_units(astropy.units.AU)  
    [mars_r, mars_lat, mars_lon]=cart2sphere(mars.x,mars.y,mars.z)

    spice.furnish(spicedata.get_kernel('stereo_a_pred'))
    sta=spice.Trajectory('-234')  
    sta.generate_positions(time1_num,'Sun',frame)  
    sta.change_units(astropy.units.AU)  
    [sta_r, sta_lat, sta_lon]=cart2sphere(sta.x,sta.y,sta.z)

    psp=np.rec.array([psp_time_num,psp_r,psp_lon,psp_lat, psp.x, psp.y,psp.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
    bepi=np.rec.array([bepi_time_num,bepi_r,bepi_lon,bepi_lat,bepi.x, bepi.y,bepi.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
    solo=np.rec.array([solo_time_num,solo_r,solo_lon,solo_lat,solo.x, solo.y,solo.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
    sta=np.rec.array([sta_time_num,sta_r,sta_lon,sta_lat,sta.x, sta.y,sta.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
    earth=np.rec.array([earth_time_num,earth_r,earth_lon,earth_lat, earth.x, earth.y,earth.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
    venus=np.rec.array([venus_time_num,venus_r,venus_lon,venus_lat, venus.x, venus.y,venus.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
    mars=np.rec.array([mars_time_num,mars_r,mars_lon,mars_lat, mars.x, mars.y,mars.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])
    mercury=np.rec.array([mercury_time_num,mercury_r,mercury_lon,mercury_lat,mercury.x, mercury.y,mercury.z],dtype=[('time','f8'),('r','f8'),('lon','f8'),('lat','f8'),('x','f8'),('y','f8'),('z','f8')])

    return [psp,bepi,solo,sta,earth,venus,mars,mercury]
     
     
     
     
     
     
     
     
     
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
     #ax1.plot_date([e2,e2],[-100,100],'--k',linewidth=1)
     #ax1.plot_date([e3,e3],[-100,100],'-b',linewidth=1)

     
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
     #ax2.plot_date([e2,e2],[0,5000],'--k',linewidth=1)
     #ax2.plot_date([e3,e3],[0,5000],'-b',linewidth=1)
     

     plt.ylabel('V [km/s]')
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((np.nanmin(sc.vt)-50,np.nanmax(sc.vt)+100))
     #ax2.set_xticklabels([])
     plt.setp(ax2.get_xticklabels(), visible=False)



     ax3 = plt.subplot(313,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)

     ax3.plot_date([e1,e1],[0,500],'-k',linewidth=1)     
     #ax3.plot_date([e2,e2],[0,500],'--k',linewidth=1)
     #ax3.plot_date([e3,e3],[0,500],'-b',linewidth=1)


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

     sc_time_num=mdates.date2num(sc.time)

     def onclick(event):

        global stepout
        #print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
           #   (event.button, event.x, event.y, event.xdata, event.ydata))
        #print('time',str(mdates.num2date(event.xdata)))
        nearestm1=np.argmin(abs(sc_time_num-event.xdata))
        seltime=str(mdates.num2date(sc_time_num[nearestm1]))
        print(seltime[0:10]+'T'+seltime[11:16]+'Z',event.ydata )  
        if event.button == 3: stepout=True

     #measure with mouse times
     #cid = fig.canvas.mpl_connect('button_press_event', onclick)

     #print(stepout) 
     #if stepout: break


     #plt.show(block=True)



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
   


     
     
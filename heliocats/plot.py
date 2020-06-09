#plots.py
#plotting routines for heliocats
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy
import copy
import matplotlib.dates as mdates
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
import heliosat
from numba import njit
from astropy.time import Time
import astropy
import importlib

import heliopy.data.spice as spicedata
import heliopy.spice as spice

from config import data_path

from heliocats import data as hd
importlib.reload(hd) #reload again while debugging


from heliocats import cats as hc
importlib.reload(hc) #reload again while debugging




    


####################################### 



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
     plt.ylim((-20, 20))
     #ax1.set_xticklabels([]) does not work with sharex
     #plt.setp(ax1.get_xticklabels(), fontsize=6)
     plt.setp(ax1.get_xticklabels(), visible=False)

     plt.title(sc_label+' data       start: '+start.strftime("%Y-%b-%d %H:%M")+'       end: '+end.strftime("%Y-%b-%d %H:%M"),fontsize=fsize)


     ax2 = plt.subplot(412,sharex=ax1) 
     ax2.plot_date(sc.time,sc.vt,'-k',label='V',linewidth=0.7)

     plt.ylabel('V [km/s]',fontsize=fsize)
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.ylim((250, 800))
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
     ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d') )
     plt.ylim((0, 0.5))
     
     plt.figtext(0.99,0.01,'Möstl, Bailey / Helio4Cast', color='black', ha='right',fontsize=fsize-3)

     
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
     ax1.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[-500,500],linestyle='-',linewidth=1,color=color_sir,marker='')            
     ax1.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[-500,500],linestyle='--',linewidth=1,color=color_sir,marker='')            
     ax1.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[-500,500],color=color_sir,linestyle='-',marker='',linewidth=1)            
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
     ax2.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[0,3000],'-k',linewidth=1,color=color_sir)            
     ax2.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[0,3000],'--k',linewidth=1,color=color_sir)            
     ax2.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[0,3000],color=color_sir,linestyle='-',linewidth=1,marker='')            
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
     ax3.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[-10,1000],color=color_sir,linestyle='-',linewidth=1,marker='')       
     ax3.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[-10,1000],color=color_sir,linestyle='--',linewidth=1,marker='')            
     ax3.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[-10,1000],color=color_sir,linestyle='-',linewidth=1,marker='')            
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
     ax4.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[0,10],color=color_sir,linewidth=1,linestyle='-',label='sir_end_time',marker='')         
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
   
                 
            
            
            

def plot_insitu_sircat_mag_plasma(sc, start, end, sc_label, path, ic,i, **kwargs):
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
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
     #plot vertical lines
     #ax1.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[-500,500],'-k',linewidth=1)            

     ax1.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[-500,500],'-k',linewidth=1)            
     ax1.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[-500,500],'--k',linewidth=1)            
     ax1.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[-500,500],color=color_sir,linestyle='-',linewidth=1,marker='')            
     ax1.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[-500,500],color=color_vtmax,linestyle='-',linewidth=1,marker='')            
     ax1.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[-500,500],color=color_sir,linestyle='--',linewidth=1,marker='')    
     plt.ylabel('B [nT]')
     plt.legend(loc=1,ncol=4,fontsize=6)
     ax1.set_xlim(start,end)
     #if np.isnan(np.nanmin(sc.bt))==False:
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
     ax2.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[0,3000],'-k',linewidth=1)            
     ax2.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[0,3000],'--k',linewidth=1)            
     ax2.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[0,3000],color=color_sir,linestyle='-',linewidth=1,marker='')            
     ax2.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[0,3000],color=color_vtmax,linestyle='-',linewidth=1,marker='')            
     ax2.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[0,3000],'--k',linewidth=1)    
     plt.ylabel('V [km/s]')
     ax2.set_xlim(start,end)
     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     ax2.set_ylim(np.nanmin(sc.vt)-20,np.nanmax(sc.vt)+100)  
     #plt.ylim((250, 800))
     #ax2.set_xticklabels([])
     plt.setp(ax2.get_xticklabels(), visible=False)

    
     ax3 = plt.subplot(413,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)
     #plot vertical lines
     ax3.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[0,1000],'-k',linewidth=1)            
     ax3.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[0,1000],'--k',linewidth=1)            
     ax3.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[0,1000],color=color_sir,linestyle='-',linewidth=1,marker='')            
     ax3.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[0,1000],color=color_vtmax,linestyle='-',linewidth=1,marker='')            
     ax3.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[0,1000],'--k',linewidth=1)    
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
     ax4.plot_date([ic.sir_start_time[i],ic.sir_start_time[i]],[0,10],'-k',linewidth=1,label='sir_start_time / sir_end_time') 
     ax4.plot_date([ic.hss_start_time[i],ic.hss_start_time[i]],[0,10],'--k',linewidth=1,label='hss_start_time / hss_end_time') 
     ax4.plot_date([ic.sir_end_time[i],ic.sir_end_time[i]],[0,10],color=color_sir,linewidth=1,linestyle='-',marker='')              
     ax4.plot_date([ic.hss_vtmax_time[i],ic.hss_vtmax_time[i]],[0,10],color=color_vtmax,linewidth=1.5,linestyle='-',label='hss_vtmax_time',marker='')
     ax4.plot_date([ic.hss_end_time[i],ic.hss_end_time[i]],[0,10],'--k',linewidth=1)    
     plt.ylabel('T [MK]')
     ax4.set_xlim(start,end)
     ax4.set_ylim(0,np.nanmax(sc.tp/1e6)+0.2)   
     ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     plt.legend(loc=1,ncol=1,fontsize=6)
     #plt.ylim((0, 0.5))
     
     plt.tight_layout()
     #plt.show()

     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     plotfile=path+ic.sircat_id[i]+'.png'  

     plt.savefig(plotfile)
     print('saved as ',plotfile)
   
     
     
     
  
            
            
################################## ICMECAT  ############################################################
  
      
   
  
     

def plot_insitu_icmecat_mag_plasma(sc, start, end, sc_label, path, ic,i, **kwargs):
     '''
     sc ... data
    
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
        
     #plot vertical lines
     ax1.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[-500,500],'-k',linewidth=1)            
     ax1.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[-500,500],'-k',linewidth=1)            
     ax1.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[-500,500],'-k',linewidth=1)            

    
     plt.ylabel('B [nT]')
     plt.legend(loc=3,ncol=4,fontsize=8)
     ax1.set_xlim(start,end)
     #if np.isnan(np.nanmin(sc.bt))==False:
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
     ax2.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[0,3000],'-k',linewidth=1)            
     ax2.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[0,3000],'-k',linewidth=1)            
     ax2.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[0,3000],'-k',linewidth=1)            



     plt.ylabel('V [km/s]')
     ax2.set_xlim(start,end)
     #check plasma data exists
     if np.isnan(np.nanmin(sc.vt))==False:
         ax2.set_ylim(np.nanmin(sc.vt)-20,np.nanmax(sc.vt)+100)   


     ax2.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     #plt.ylim((250, 800))
     #ax2.set_xticklabels([])
     plt.setp(ax2.get_xticklabels(), visible=False)


     ax3 = plt.subplot(413,sharex=ax1) 
     ax3.plot_date(sc.time,sc.np,'-k',label='Np',linewidth=0.7)
     
     #plot vertical lines
     ax3.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[0,1000],'-k',linewidth=1)            
     ax3.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[0,1000],'-k',linewidth=1)            
     ax3.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[0,1000],'-k',linewidth=1)            


     plt.ylabel('N [ccm-3]')
     ax3.set_xlim(start,end)
     if np.isnan(np.nanmin(sc.np))==False:
         ax3.set_ylim(0,np.nanmax(sc.np)+10)   


     ax3.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     #plt.ylim((0, 50))
     #ax3.set_xticklabels([])
     plt.setp(ax3.get_xticklabels(), visible=False)


     ax4 = plt.subplot(414,sharex=ax1) 
     ax4.plot_date(sc.time,sc.tp/1e6,'-k',label='Tp',linewidth=0.7)
    
     #plot vertical lines
     ax4.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[0,10],'-k',linewidth=1)            
     ax4.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[0,10],'-k',linewidth=1)            
     ax4.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[0,10],'-k',linewidth=1)            



     plt.ylabel('T [MK]')
     ax4.set_xlim(start,end)
     if np.isnan(np.nanmin(sc.tp))==False:
         ax4.set_ylim(0,np.nanmax(sc.tp/1e6)+0.2)   

     ax4.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%b-%d %H') )
     #plt.ylim((0, 0.5))
     
     plt.tight_layout()
     #plt.show()

     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'

     plotfile=path+ic.icmecat_id[i]+'.png'
  

     plt.savefig(plotfile)
     print('saved as ',plotfile)
   
     
     
     
     
     
     
     

def plot_insitu_icmecat_mag(sc, start, end, sc_label, path, ic, i, **kwargs):
     '''
     sc = data
    
     '''
     
     start=parse_time(start).datetime
     end=parse_time(end).datetime
     #print(start)
     #print(end)
     
     sns.set_style('darkgrid')
     sns.set_context('paper')

     fig=plt.figure(figsize=(9,3), dpi=150)
     
     ax1 = plt.subplot(111) 

     ax1.plot_date(sc.time,sc.bx,'-r',label='Bx',linewidth=0.5)
     ax1.plot_date(sc.time,sc.by,'-g',label='By',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bz,'-b',label='Bz',linewidth=0.5)
     ax1.plot_date(sc.time,sc.bt,'-k',label='Btotal',lw=0.5)
    
     #plot vertical lines
     ax1.plot_date([ic.icme_start_time[i],ic.icme_start_time[i]],[-500,500],'-k',linewidth=1)            
     ax1.plot_date([ic.mo_start_time[i],ic.mo_start_time[i]],[-500,500],'-k',linewidth=1)            
     ax1.plot_date([ic.mo_end_time[i],ic.mo_end_time[i]],[-500,500],'-k',linewidth=1)            

     
     plt.ylabel('B [nT]')
     plt.legend(loc=3,ncol=4,fontsize=8)
     ax1.set_xlim(start,end)
     ax1.set_ylim(-np.nanmax(sc.bt)-5,np.nanmax(sc.bt)+5)   

     ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b-%d') )

     plt.title(sc_label+' data, start: '+start.strftime("%Y-%b-%d %H:%M")+'  end: '+end.strftime("%Y-%b-%d %H:%M"))
     
    
     plt.tight_layout()

     #plotfile=path+sc_label+'_'+start.strftime("%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
     
     
     plotfile=path+ic.icmecat_id[i]+'.png'
       
     plt.savefig(plotfile)
     print('saved as ',plotfile)
   
     
     
     
    

def plot_icmecat_events(sc,sci,ic,name,icplotsdir):

  
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
        rad=hd.load_msl_rad()   
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
                
    
    for i in np.arange(np.size(sci)):    
    #for i in np.arange(1):     
        if plasma == True:
            print(i)
            
            if name!='MAVEN':
                plot_insitu_icmecat_mag_plasma(sc[icme_start_ind[i]-90*24:mo_end_ind[i]+90*24],\
                             ic.icme_start_time[sci[i]]-datetime.timedelta(days=1.5), \
                             ic.mo_end_time[sci[i]]+datetime.timedelta(days=1.5),name, icplotsdir,ic,sci[i])
            if name == 'MAVEN':                 
                plot_insitu_icmecat_maven(sc[icme_start_ind[i]-7*6:mo_end_ind[i]+7*6],\
                             ic.icme_start_time[sci[i]]-datetime.timedelta(days=7), \
                             ic.mo_end_time[sci[i]]+datetime.timedelta(days=7),\
                             name,icplotsdir,ic,sci[i],rad,arrcat,msir,w1)
                
            plt.close('all')
        else:
            plot_insitu_icmecat_mag(sc[icme_start_ind[i]-90*24:mo_end_ind[i]+90*24], \
                                    ic.icme_start_time[sci[i]]-datetime.timedelta(days=1.5), \
                                    ic.mo_end_time[sci[i]]+datetime.timedelta(days=1.5),name, icplotsdir,ic,sci[i])
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
    
    plt.title(sc_label+' data, start: '+start.strftime("%Y-%b-%d %H:%M")+'  end: '+end.strftime("%Y-%b-%d %H:%M"))


    #sharex means that zooming in works with all subplots
    ax1 = plt.subplot(411) 
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
    plt.legend(loc=3,ncol=4,fontsize=7)
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
    plt.ylabel('V [km/s]')
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

    plt.ylabel('N [ccm-3]')
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
    sc = data

    '''
    sns.set_style('darkgrid')
    sns.set_context('paper')    
    
    time1=mdates.date2num(time_date1)

    #made with sc_positions_for_vr
    [psp, bepi, solo, sta, earth, venus, mars, mercury,jupiter, saturn, uranus, neptune,frame]=pickle.load( open( 'results/positions_HEEQ_1hr.p', "rb" ) )
    #[psp,bepi,solo,sta,earth,venus,mars,mercury]=make_positions(time1,frame)
    

    #sidereal solar rotation rate
    if frame=='HCI': sun_rot=24.47
    #synodic
    if frame=='HEEQ': sun_rot=26.24
   
    AUkm=149597870.7   

    #for parker spiral   
    theta=np.arange(0,np.deg2rad(180),0.01)
    res_in_days=1/24
    k=0
    
    plot_orbit=True
    plot_parker=True
    fadeind=int(100/res_in_days)
    fsize=17 
    symsize_planet=140
    symsize_spacecraft=100
    
    #find index for psp
    dct=time1-psp.time
    psp_timeind=np.argmin(abs(dct))

    dct=time1-bepi.time
    bepi_timeind=np.argmin(abs(dct))

    dct=time1-solo.time
    solo_timeind=np.argmin(abs(dct))

    #all others same time as Earth
    dct=time1-earth.time
    earth_timeind=np.argmin(abs(dct))
    
    
    ################## figure    
    fig=plt.figure(1, figsize=(15,10), dpi=150)
    ax = plt.subplot(111, projection='polar')
    backcolor='black'
    psp_color='black'
    bepi_color='blue'
    solo_color='green'

    ax.scatter(venus.lon[earth_timeind], venus.r[earth_timeind]*np.cos(venus.lat[earth_timeind]), s=symsize_planet, c='orange', alpha=1,lw=0,zorder=3)
    ax.scatter(mercury.lon[earth_timeind], mercury.r[earth_timeind]*np.cos(mercury.lat[earth_timeind]), s=symsize_planet, c='dimgrey', alpha=1,lw=0,zorder=3)
    ax.scatter(earth.lon[earth_timeind], earth.r[earth_timeind]*np.cos(earth.lat[earth_timeind]), s=symsize_planet, c='mediumseagreen', alpha=1,lw=0,zorder=3)
    ax.scatter(sta.lon[earth_timeind], sta.r[earth_timeind]*np.cos(sta.lat[earth_timeind]), s=symsize_spacecraft, c='red', marker='s', alpha=1,lw=0,zorder=3)
    ax.scatter(mars.lon[earth_timeind], mars.r[earth_timeind]*np.cos(mars.lat[earth_timeind]), s=symsize_planet, c='orangered', alpha=1,lw=0,zorder=3)

    plt.text(sta.lon[earth_timeind]-0.15,sta.r[earth_timeind],'STEREO-A', color='red', ha='center',fontsize=fsize-4,verticalalignment='top')
  
    plt.text(0,0,'Sun', color='black', ha='center',fontsize=fsize-5,verticalalignment='top')
    plt.text(0,earth.r[earth_timeind]+0.12,'Earth', color='mediumseagreen', ha='center',fontsize=fsize-5,verticalalignment='center')
    
    plt.figtext(0.99,0.01,'C. Möstl / Helio4Cast', color='black', ha='right',fontsize=fsize-8)

    plt.figtext(0.85,0.1,'――― 100 days future trajectory', color='black', ha='center',fontsize=fsize-3)

    
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


    if bepi_timeind > 0:
        ax.scatter(bepi.lon[bepi_timeind], bepi.r[bepi_timeind]*np.cos(bepi.lat[bepi_timeind]), s=symsize_spacecraft, c=bepi_color, marker='s', alpha=1,lw=0,zorder=3)
        bepi_text='Bepi:   '+str(f'{bepi.r[bepi_timeind]:6.2f}')+str(f'{np.rad2deg(bepi.lon[bepi_timeind]):8.1f}')+str(f'{np.rad2deg(bepi.lat[bepi_timeind]):8.1f}')
        f6=plt.figtext(0.01,0.74,bepi_text, fontsize=fsize, ha='left',color=bepi_color)
        plt.text(bepi.lon[bepi_timeind]-0.15,bepi.r[bepi_timeind],'Bepi Colombo', color='blue', ha='center',fontsize=fsize-4,verticalalignment='top')
  
        if plot_orbit: 
            ax.plot(bepi.lon[bepi_timeind:bepi_timeind+fadeind], bepi.r[bepi_timeind:bepi_timeind+fadeind]*np.cos(bepi.lat[bepi_timeind:bepi_timeind+fadeind]), c=bepi_color, alpha=0.6,lw=1,zorder=3)



    if solo_timeind > 0:
        ax.scatter(solo.lon[solo_timeind], solo.r[solo_timeind]*np.cos(solo.lat[solo_timeind]), s=symsize_spacecraft, c=solo_color, marker='s', alpha=1,lw=0,zorder=3)
        solo_text='SolO:  '+str(f'{solo.r[solo_timeind]:6.2f}')+str(f'{np.rad2deg(solo.lon[solo_timeind]):8.1f}')+str(f'{np.rad2deg(solo.lat[solo_timeind]):8.1f}')
        f7=plt.figtext(0.01,0.7,solo_text, fontsize=fsize, ha='left',color=solo_color)
        plt.text(solo.lon[solo_timeind]-0.15,solo.r[solo_timeind],'Solar Orbiter', color='green', ha='center',fontsize=fsize-4,verticalalignment='top')

        if plot_orbit: 
            ax.plot(solo.lon[solo_timeind:solo_timeind+fadeind], solo.r[solo_timeind:solo_timeind+fadeind]*np.cos(solo.lat[solo_timeind:solo_timeind+fadeind]), c=solo_color, alpha=0.6,lw=1,zorder=3)

    if plot_orbit: 
        ax.plot(sta.lon[earth_timeind:earth_timeind+fadeind], sta.r[earth_timeind:earth_timeind+fadeind]*np.cos(sta.lat[earth_timeind:earth_timeind+fadeind]), c='red', alpha=0.6,lw=1,zorder=3)




    if frame=='HEEQ': earth_text='Earth: '+str(f'{earth.r[earth_timeind]:6.2f}')+str(f'{0.0:8.1f}')+str(f'{np.rad2deg(earth.lat[earth_timeind]):8.1f}')
    else: earth_text='Earth: '+str(f'{earth.r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(earth.lon[earth_timeind]):8.1f}')+str(f'{np.rad2deg(earth.lat[earth_timeind]):8.1f}')

    mars_text='Mars:  '+str(f'{mars.r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(mars.lon[earth_timeind]):8.1f}')+str(f'{np.rad2deg(mars.lat[earth_timeind]):8.1f}')
    sta_text='STA:   '+str(f'{sta.r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(sta.lon[earth_timeind]):8.1f}')+str(f'{np.rad2deg(sta.lat[earth_timeind]):8.1f}')

    #Sun
    ax.scatter(0,0,s=100,c='yellow',alpha=1, edgecolors='black', linewidth=0.3)


 
    f10=plt.figtext(0.01,0.9,earth_text, fontsize=fsize, ha='left',color='mediumseagreen')
    f9=plt.figtext(0.01,0.86,mars_text, fontsize=fsize, ha='left',c='orangered')
    f8=plt.figtext(0.01,0.82,sta_text, fontsize=fsize, ha='left',c='red')
    
    #time
    plt.figtext(0.65,0.9,time_date1.strftime("%Y %B %d  %H:%M"),fontsize=fsize+6, ha='left',c='black')

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
    ax.set_ylim(0, 1.25) 

    plt.thetagrids(range(0,360,45),(u'0\u00b0 '+frame+' longitude',u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'- 135\u00b0',u'- 90\u00b0',u'- 45\u00b0'), fmt='%d',ha='left',fontsize=fsize,color=backcolor, zorder=5, alpha=0.9)

    plt.tight_layout()

    #plotfile=path+%Y_%b_%d")+'_'+end.strftime("%Y_%b_%d")+'.png'
    plotfile=path+'positions_'+time_date1.strftime("%Y_%b_%d")+'.png'

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
   


     
     
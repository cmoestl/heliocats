# for installation of a conda environment to run this notebook, see instructions in README.md
# conda dependencies are listed under environment.yml, and pip in requirements.txt

from scipy import stats
import scipy.io
from matplotlib import cm
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import astropy.constants as const
import sunpy.time
import time
import pickle
import seaborn as sns
import os
import urllib
import json


#matplotlib.use('TkAgg')

#warnings.filterwarnings('ignore') # some numpy mean-of-empty-slice runtime warnings



plt.close('all')
print()
print('cme_stats.py main program.')
print('ICME parameters at all 4 terrestrial planets.')
print('Christian Moestl, IWF Graz, Austria')



resdir='results'
if os.path.isdir(resdir) == False: os.mkdir(resdir)

outputdirectory='results/plots_stats'
#create plot directory
if os.path.isdir(outputdirectory) == False: os.mkdir(outputdirectory)

#load ICMECAT
file='icmecat/HELCATS_ICMECAT_v20.p'
print('loaded ', file)
ic=pickle.load(open(file, "rb" ) )  
ic.keys() 


wini=np.where(ic.sc_insitu == 'Wind')[:][0] 
vexi=np.where(ic.sc_insitu == 'VEX')[:][0]  
mesi=np.where(ic.sc_insitu == 'MESSENGER')[:][0]   
stai=np.where(ic.sc_insitu == 'STEREO-A')[:][0]    
stbi=np.where(ic.sc_insitu == 'STEREO-B')[:][0]    
mavi=np.where(ic.sc_insitu == 'MAVEN')[:][0]    
ulyi=np.where(ic.sc_insitu == 'ULYSSES')[:][0]   

#print(ic.icme_start_time)

imercind=np.where(np.logical_and(ic.sc_insitu =='MESSENGER', ic.icme_start_time > sunpy.time.parse_time('2011-03-18').datetime))
print(imercind)


#pspi=np.where(ic.sc_insitu == 'ParkerSolarProbe')[:][0]    


#get parameters
ic.keys()  

########################### get all parameters from ICMECAT for easier handling
# id for each event
#iid=ic.icmecat_id.to_numpy()

# observing spacecraft
#isc=ic.sc_insitu.to_numpy()  
#parameters
#mo_bmax=ic.mo_bmax.to_numpy()
#sc_heliodistance=ic.mo_sc_heliodistance.to_numpy()



###############################################################
############### set limits of solar minimum, rising/declining phase and solar maximum

minstart=sunpy.time.parse_time('2007-01-01').datetime
minend=sunpy.time.parse_time('2009-12-31').datetime


minstart_num=mdates.date2num(sunpy.time.parse_time('2007-01-01').datetime)
minend_num=mdates.date2num(sunpy.time.parse_time('2009-12-31').datetime)

risestart=sunpy.time.parse_time('2010-01-01').datetime
riseend=sunpy.time.parse_time('2011-06-30').datetime


risestart_num=mdates.date2num(sunpy.time.parse_time('2010-01-01').datetime)
riseend_num=mdates.date2num(sunpy.time.parse_time('2011-06-30').datetime)


maxstart=sunpy.time.parse_time('2011-07-01').datetime
maxend=sunpy.time.parse_time('2014-12-31').datetime

maxstart_num=mdates.date2num(sunpy.time.parse_time('2011-07-01').datetime)
maxend_num=mdates.date2num(sunpy.time.parse_time('2014-12-31').datetime)


declstart=sunpy.time.parse_time('2015-01-01').datetime
declend=sunpy.time.parse_time('2017-12-31').datetime

declstart_num=mdates.date2num(sunpy.time.parse_time('2015-01-01').datetime)
declend_num=mdates.date2num(sunpy.time.parse_time('2017-12-31').datetime)


############### extract events by limits of solar min, 
############### rising, max, too few events for MAVEN and Ulysses

iallind_min=np.where(np.logical_and(ic.icme_start_time > minstart,ic.icme_start_time < minend))[0]
iallind_rise=np.where(np.logical_and(ic.icme_start_time > risestart,ic.icme_start_time < riseend))[0]
iallind_max=np.where(np.logical_and(ic.icme_start_time > maxstart,ic.icme_start_time < maxend))[0]

iwinind_min=iallind_min[np.where(ic.sc_insitu[iallind_min]=='Wind')]
iwinind_rise=iallind_rise[np.where(ic.sc_insitu[iallind_rise]=='Wind')]
iwinind_max=iallind_max[np.where(ic.sc_insitu[iallind_max]=='Wind')]

ivexind_min=iallind_min[np.where(ic.sc_insitu[iallind_min]=='VEX')]
ivexind_rise=iallind_rise[np.where(ic.sc_insitu[iallind_rise]=='VEX')]
ivexind_max=iallind_max[np.where(ic.sc_insitu[iallind_max]=='VEX')]

imesind_min=iallind_min[np.where(ic.sc_insitu[iallind_min]=='MESSENGER')]
imesind_rise=iallind_rise[np.where(ic.sc_insitu[iallind_rise]=='MESSENGER')]
imesind_max=iallind_max[np.where(ic.sc_insitu[iallind_max]=='MESSENGER')]

istaind_min=iallind_min[np.where(ic.sc_insitu[iallind_min]=='STEREO-A')]
istaind_rise=iallind_rise[np.where(ic.sc_insitu[iallind_rise]=='STEREO-A')]
istaind_max=iallind_max[np.where(ic.sc_insitu[iallind_max]=='STEREO-A')]

istbind_min=iallind_min[np.where(ic.sc_insitu[iallind_min]=='STEREO-B')]
istbind_rise=iallind_rise[np.where(ic.sc_insitu[iallind_rise]=='STEREO-B')]
istbind_max=iallind_max[np.where(ic.sc_insitu[iallind_max]=='STEREO-B')]


# select the events at Mercury extra after orbit insertion
# no events available for solar minimum!
imercind_min=iallind_min[np.where(np.logical_and(ic.sc_insitu[iallind_min] =='MESSENGER',ic.icme_start_time[iallind_min] > sunpy.time.parse_time('2011-03-18').datetime))]
imercind_rise=iallind_rise[np.where(np.logical_and(ic.sc_insitu[iallind_rise] =='MESSENGER',ic.icme_start_time[iallind_rise] > sunpy.time.parse_time('2011-03-18').datetime))]
imercind_max=iallind_max[np.where(np.logical_and(ic.sc_insitu[iallind_max] =='MESSENGER',ic.icme_start_time[iallind_max] > sunpy.time.parse_time('2011-03-18').datetime))]




##################### (1a) ICME DURATION VS DISTANCE and linear fit  #####################

sns.set_context("talk")     #sns.set_style("darkgrid")  
sns.set_style("ticks",{'grid.linestyle': '--'})

print('-------------------------------------------------')
print()
print('1a ICME DURATION VS DISTANCE')
print()

fig=plt.figure(1,figsize=(12,11	))
fsize=15
ax1 = plt.subplot2grid((2,1), (0, 0))
#x axis
xfit=np.linspace(0,2,1000)

#force through origin, fit with y=kx
scx=ic.mo_sc_heliodistance[:,np.newaxis]
durfit_f, _, _, _ =np.linalg.lstsq(scx,ic.icme_duration, rcond=None)

scxmin=ic.mo_sc_heliodistance[iallind_min][:,np.newaxis]
durfitmin_f, _, _, _ =np.linalg.lstsq(scxmin,ic.icme_duration[iallind_min],rcond=None)

scxrise=ic.mo_sc_heliodistance[iallind_rise][:,np.newaxis]
durfitrise_f, _, _, _ =np.linalg.lstsq(scxrise,ic.icme_duration[iallind_rise],rcond=None)

scxmax=ic.mo_sc_heliodistance[iallind_max][:,np.newaxis]
durfitmax_f, _, _, _ =np.linalg.lstsq(scxmax,ic.icme_duration[iallind_max],rcond=None)

#make the y axis for the fits forced through the origin
ydurfitall_f=durfit_f*xfit
ydurfitmin_f=durfitmin_f*xfit
ydurfitrise_f=durfitrise_f*xfit
ydurfitmax_f=durfitmax_f*xfit

plt.plot(ic.mo_sc_heliodistance,ic.icme_duration,'o',color='blue',markersize=4, alpha=0.4)
#for plotting min/rise/max differently
#plt.plot(sc_heliodistance[iallind_min],icme_durations[iallind_min],'o',color='dimgre',markersize=3, alpha=0.4,label='D min')
#plt.plot(sc_heliodistance[iallind_rise],icme_durations[iallind_rise],'o',color='grey',markersize=3, alpha=0.7,label='D rise')
#plt.plot(sc_heliodistance[iallind_max],icme_durations[iallind_max],'o',color='black',markersize=3, alpha=0.8,label='D max')

#plot fits
plt.plot(xfit,ydurfitall_f,'-',color='blue', lw=3, alpha=0.9,label='fit')
plt.plot(xfit,ydurfitmin_f,'--',color='black', lw=2, alpha=0.9,label='min fit')
plt.plot(xfit,ydurfitrise_f,'-.',color='black', lw=2, alpha=0.9,label='rise fit')
plt.plot(xfit,ydurfitmax_f,'-',color='black', lw=2, alpha=0.9,label='max fit')

print()

print('linear fit results, hours vs AU')	
print('overall:    D[h]={:.2f} R[AU] '.format(durfit_f[0]))
print('minimum:    D[h]={:.2f} R[AU] '.format(durfitmin_f[0]))
print('rise phase: D[h]={:.2f} R[AU]'.format(durfitrise_f[0]))
print('maximum:    D[h]={:.2f} R[AU]'.format(durfitmax_f[0]))

label_level=85
label_level=30

#plt.annotate('overall:',xy=(0.06,label_level),fontsize=11)
#plt.annotate('D[h]={:.2f} R[AU]'.format(durfit_f[0]),xy=(0.15,label_level),fontsize=11)
plt.annotate('min:',xy=(0.06,label_level-4),fontsize=11) 
plt.annotate('D[h]={:.2f} R[AU] '.format(durfitmin_f[0]),xy=(0.12,label_level-4),fontsize=11)
plt.annotate('rise:',xy=(0.06,label_level-8),fontsize=11) 
plt.annotate('D[h]={:.2f} R[AU]'.format(durfitrise_f[0]),xy=(0.12,label_level-8),fontsize=11)
plt.annotate('max:',xy=(0.06,label_level-12),fontsize=11)    
plt.annotate('D[h]={:.2f} R[AU]'.format(durfitmax_f[0]),xy=(0.12,label_level-12),fontsize=11)

'''
#planet limits
plt.axvspan(np.min(pos.mars[0]),np.max(pos.mars[0]), color='orangered', alpha=0.2)
plt.axvspan(np.min(pos.mercury[0]),np.max(pos.mercury[0]), color='darkgrey', alpha=0.2)
plt.axvspan(np.min(pos.venus[0]),np.max(pos.venus[0]), color='orange', alpha=0.2)
plt.axvspan(np.min(pos.earth[0]),np.max(pos.earth[0]), color='mediumseagreen', alpha=0.2)
#plt.axvspan(np.min(pos.sta[0]),np.max(pos.sta[0]), color='red', alpha=0.2)  #STEREO-A
#plt.axvspan(np.min(pos.stb[0]),np.max(pos.stb[0]), color='blue', alpha=0.2)  #STEREO-B
#Parker Probe minimum
plt.plot([0.046,0.046],[0,110], color='black', linestyle='--', linewidth=1)

#label_level=100
label_level=38

plt.annotate('Mars', xy=(1.5,label_level), ha='center',fontsize=fsize)
plt.annotate('Mercury', xy=(0.38,label_level), ha='center',fontsize=fsize)
plt.annotate('Venus', xy=(0.72,label_level), ha='center',fontsize=fsize)
plt.annotate('Earth', xy=(1,label_level), ha='center',fontsize=fsize)

plt.annotate('PSP', xy=(0.05,label_level), ha='left',fontsize=fsize)



#for PSP zoom
plt.xlim(0,1.2)
plt.ylim(0,40)
'''


ax1.set_xticks(np.arange(0,2,0.1))
#plt.xlim(0,max(sc_heliodistance)+0.3)
#plt.ylim(0,110)

plt.legend(loc=4,fontsize=fsize-1)
plt.xlabel('Heliocentric distance R [AU]',fontsize=fsize)
plt.ylabel('ICME duration D [hours]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 
#plt.grid()





##################### (1b) ICME DURATION VS TIME ########################################

print('')
print('1b ICME DURATION VS TIME')

#tests for gaussians and Hathaway's function for solar cycle, not used

#Wind
#tfit=mdates.date2num(sunpy.time.parse_time('2009-04-01'))+np.arange(0,365*10)
#t0=mdates.date2num(sunpy.time.parse_time('2009-01-01'))

#Gaussian
#sigma=1000
#bfitmax=30
#mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
#ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='mediumseagreen',linestyle='-',markersize=0, label='Earth fit')

#Hathaway 2015 equation 6 page 40
#average cycle sunspot number 
#A=100 #amplitude ##195 for sunspot
#b=100*12 #56*12 for months to days
#c=0.8
#4 free parameters A, b, c, t0

#Fwind=A*(((tfit-t0)/b)**3) * 1/(np.exp((((tfit-t0)/b)**2))-c)
#plt.plot_date(tfit, Fwind,'o',color='mediumseagreen',linestyle='-',markersize=0, label='Earth fit')

#xaxis: 10 years, daily data point
#xfit2=mdates.date2num(sunpy.time.parse_time('2007-01-01'))+np.arange(0,365*10)
#MESSENGER
#sigma=1000
#bfitmax=10
#mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
#ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit2-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='darkgrey',linestyle='-',markersize=0, label='Mercury fit')

#VEX
#inital guess
#sigma=1000
#bfitmax=20
#mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
#ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit2-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit2, ygauss*1/max(ygauss)*bfitmax,'o',color='orange',linestyle='-',markersize=0, label='Venus fit')

#for Mars: reconstruct likely parameters if sigma is quite similar for all fits, take mean of those sigmas and adjust bfitmax as function of distance with power law)
#plot reconstructed function for Mars
#bfitmax=40
#plt.plot_date(xfit2, Fwind,'o',color='steelblue',linestyle='--',markersize=0, label='Mars reconstr.')



################## plot 

ax2 = plt.subplot2grid((2,1), (1, 0))
markers=6
linew=0

#plot durations for all planets
ax2.plot_date(ic.icme_start_time[mesi],ic.icme_duration[mesi], \
    'o',color='darkgrey',markersize=markers,linestyle='-',linewidth=linew,label='MESSENGER')
ax2.plot_date(ic.icme_start_time[vexi],ic.icme_duration[vexi], \
    'o',color='orange',markersize=markers,linestyle='-',linewidth=linew, label='Venus')
ax2.plot_date(ic.icme_start_time[wini],ic.icme_duration[wini], \
    'o',color='mediumseagreen',markersize=markers, linestyle='-', linewidth=linew, label='Earth')
ax2.plot_date(ic.icme_start_time[mavi],ic.icme_duration[mavi], \
    'o',color='steelblue',markersize=markers,linestyle='-',linewidth=linew, label='Mars')


#limits solar min/rise/maxax2.set_ylim(0,80)
vlevel=130
spanalpha=0.05

plt.axvspan(minstart,minend, color='green', alpha=spanalpha)
plt.annotate('solar minimum',xy=(minstart_num+(minend_num-minstart_num)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(minstart_num+10,vlevel),ha='left')
plt.annotate('>',xy=(minend_num-10,vlevel),ha='right')

plt.axvspan(risestart,riseend, color='yellow', alpha=spanalpha)
plt.annotate('rising phase',xy=(risestart_num+(riseend_num-risestart_num)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(risestart_num+10,vlevel),ha='left')
plt.annotate('>',xy=(riseend_num-10,vlevel),ha='right')

plt.axvspan(maxstart,maxend, color='red', alpha=spanalpha)
plt.annotate('solar maximum',xy=(maxstart_num+(maxend_num-maxstart_num)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(maxstart_num+10,vlevel),ha='left')
plt.annotate('>',xy=(maxend_num,vlevel),ha='right')


#plot means as horizontal lines for each sub interval
plt.plot_date( [minstart,minend], [np.mean(ic.icme_duration[iwinind_min]),np.mean(ic.icme_duration[iwinind_min])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [minstart,minend], [np.mean(ic.icme_duration[ivexind_min]),np.mean(ic.icme_duration[ivexind_min])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [minstart,minend], [np.mean(ic.icme_duration[imesind_min]),np.mean(ic.icme_duration[imesind_min])], color='darkgrey', linestyle='-', markersize=0) 

plt.plot_date( [risestart,riseend], [np.mean(ic.icme_duration[iwinind_rise]),np.mean(ic.icme_duration[iwinind_rise])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [risestart,riseend], [np.mean(ic.icme_duration[ivexind_rise]),np.mean(ic.icme_duration[ivexind_rise])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [risestart,riseend], [np.mean(ic.icme_duration[imesind_rise]),np.mean(ic.icme_duration[imesind_rise])], color='darkgrey', linestyle='-', markersize=0) 

plt.plot_date( [maxstart,maxend], [np.mean(ic.icme_duration[iwinind_max]),np.mean(ic.icme_duration[iwinind_max])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [maxstart,maxend], [np.mean(ic.icme_duration[ivexind_max]),np.mean(ic.icme_duration[ivexind_max])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [maxstart,maxend], [np.mean(ic.icme_duration[imesind_max]),np.mean(ic.icme_duration[imesind_max])], color='darkgrey', linestyle='-', markersize=0) 

plt.xlim(sunpy.time.parse_time('2007-01-01').datetime, sunpy.time.parse_time('2016-12-31').datetime)
plt.ylabel('ICME duration D [hours]',fontsize=fsize)
plt.xlabel('year',fontsize=fsize)
plt.tight_layout()
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 
plt.legend(loc=1,fontsize=fsize-1)

#panel labels
plt.figtext(0.01,0.98,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.01,0.485,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')

plt.show()
plt.savefig('results/plots_stats/ic.icme_duration_distance_time_paper.pdf', dpi=300)
plt.savefig('results/plots_stats/ic.icme_duration_distance_time_paper.png', dpi=300)






















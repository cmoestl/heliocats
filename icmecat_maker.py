'''
 icmecat_maker.py

 makes the ICMECATv2.0

 Author: C. Moestl, IWF Graz, Austria
 twitter @chrisoutofspace, https://github.com/cmoestl/heliocats
 last update January 2020

 python > 3.7 
 
 install a conda environment to run this code: https://github.com/cmoestl/heliocats

 needs file /heliocats/data.py
 saves under /data and /results and /icmecat

 current status:
 work in progress
 
 

to do:

- despike sta stb wind all
- go through all ICMEs and extract data
- (new B and V for STA, Wind and PSP converted to SCEQ components, plasma correct for new PSP, wind, sta)



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



from scipy import stats
import scipy.io
from matplotlib import cm
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import  DateFormatter
import numpy as np
import astropy.constants as const
import time
import pickle
import seaborn as sns
import os
import urllib
import json
import importlib
import heliopy.data.spice as spicedata
import heliopy.spice as spice
import astropy
import pandas as pd
import openpyxl
import heliosat
from sunpy.time import parse_time
import datetime
import seaborn as sns
import copy
from numba import njit


from heliocats import data as hd
importlib.reload(hd) #reload again while debugging


#where the final data are located
data_path='/nas/helio/data/insitu_python/'



###################################### ICMECAT operations ################################


def load_helcats_icmecat_master_from_excel(file):

    print('load HELCATS ICMECAT from file:', file)
    ic=pd.read_excel(file)

    #convert times to datetime objects
    for i in np.arange(0,ic.shape[0]):    
    
        a=str(ic.icme_start_time[i]).strip() #remove leading and ending blank spaces if any
        ic.icme_start_time.loc[i]=parse_time(a).datetime

        a=str(ic.mo_start_time[i]).strip() #remove leading and ending blank spaces if any
        ic.mo_start_time.loc[i]=parse_time(a).datetime 

        a=str(ic.mo_end_time[i]).strip() #remove leading and ending blank spaces if any
        ic.mo_end_time.loc[i]=parse_time(a).datetime 

        
        a=str(ic.icme_end_time[i]).strip() #remove leading and ending blank spaces if any
        if a!= '9999-99-99T99:99Z':
            ic.icme_end_time.loc[i]=parse_time(a).datetime 
        else: ic.icme_end_time.loc[i]=np.nan

    return ic




    
##########################################################################################
######################################## MAIN PROGRAM ####################################
##########################################################################################

##################################### (1) load new data with HelioSat and heliocats.data


load_data=0


if load_data >0:

    print('load new Wind, STEREO-A, MAVEN, and ParkerProbe data')

    #filemav='maven_2014_2018.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )

    #filemav='maven_2014_2018_removed.p'
    #[mav,hmav]=pickle.load(open(filemav, 'rb' ) )
    
    filemav='maven_2014_2018_removed_smoothed.p'
    [mav,hmav]=pickle.load(open(data_path+filemav, 'rb' ) )


    #filewin2="data/wind_jan2018_nov2019_GSE_HEEQ.p"
    #win2=pickle.load(open(filewin2, "rb" ) )  
    
    filesta2='sta_2018_now_beacon.p'
    sta2=pickle.load(open(data_path+filesta2, "rb" ) )  


    filepsp='psp_2018_2019.p'
    [psp,hpsp]=pickle.load(open(data_path+filepsp, "rb" ) )  


    # ADD BepiColombo  
    # ADD Solar Orbiter



    ##################################### (2) load HELCATS DATACAT

    #download if you need this file and change the path, url for this file is: ###########********* TO DO
    #add mes vex orbit position
    [vex,win,mes,sta,stb,uly,hvex,hwin,hmes,hsta,hstb,huly]=hd.load_helcats_datacat(data_path+'helcats_all_data_removed.p') 




################################ (3) make ICMECAT 



ic=load_helcats_icmecat_master_from_excel('icmecat/HELCATS_ICMECAT_v20_master.xlsx')

#get indices for all spacecraft
wini=np.where(ic.sc_insitu == 'Wind')[:][0] 
vexi=np.where(ic.sc_insitu == 'VEX')[:][0]  
mesi=np.where(ic.sc_insitu == 'MESSENGER')[:][0]   
stai=np.where(ic.sc_insitu == 'STEREO-A')[:][0]    
stbi=np.where(ic.sc_insitu == 'STEREO-B')[:][0]    
mavi=np.where(ic.sc_insitu == 'MAVEN')[:][0]    
pspi=np.where(ic.sc_insitu == 'ParkerSolarProbe')[:][0]    



#get parameters

win_istart=mdates.date2num(ic.icme_start_time[wini])   
win_iend=mdates.date2num(ic.icme_end_time[wini])   
ic.icme_duration.loc[wini]=np.round((win_iend-win_istart)*24,2)



sta_mstart=mdates.date2num(ic.mo_start_time[stai])   
sta_mend=mdates.date2num(ic.mo_end_time[stai])   
ic.mo_duration.loc[stai]=np.round((sta_mend-sta_mstart)*24,2)




################################ (4) save ICMECAT #################################

ic_date=copy.deepcopy(ic)  

#pickle, excel, json, csv, txt (cdf? votable?)

#save as pickle with datetime
file='icmecat/HELCATS_ICMECAT_v20.p'
pickle.dump(ic, open(file, 'wb'))



#use date and time format from master table
ic2=pd.read_excel('icmecat/HELCATS_ICMECAT_v20_master.xlsx')
ic.icme_start_time=ic2.icme_start_time
ic.mo_start_time=ic2.mo_start_time
ic.mo_end_time=ic2.mo_end_time
ic.icme_end_time=ic2.icme_end_time
del(ic2)

#save as Excel
file='icmecat/HELCATS_ICMECAT_v20.xlsx'
ic.to_excel(file,sheet_name='ICMECATv2.0')

#save as json
file='icmecat/HELCATS_ICMECAT_v20.json'
ic.to_json(file)

#save as csv
file='icmecat/HELCATS_ICMECAT_v20.csv'
ic.to_csv(file)


#save as hdf needs pip install tables
#file='icmecat/HELCATS_ICMECAT_v20.hdf'
#ic.to_hdf(file,key='icmecat')


#save as .mat does not work yet
#ile='icmecat/HELCATS_ICMECAT_v20.mat'
#icdict=ic.to_dict()
#scipy.io.savemat(file,ic.values)


#save as txt
file='icmecat/HELCATS_ICMECAT_v20.txt'
np.savetxt(file, ic.values.astype(str), fmt='%s' )




#icl=pickle.load(open(file, 'rb' ) )

######################################################################################
################################### END MAIN #########################################
######################################################################################


sys.exit()
















'''
sns.set_style('darkgrid')
fig1=plt.figure(figsize=[25, 12],dpi=100)
 

plt.suptitle('PSP first 2 orbits')

ax1=plt.subplot(421)
plt.plot_date(t1,bx1,'-r',label='BR',linewidth=0.5)
plt.plot_date(t1,by1,'-g',label='BT',linewidth=0.5)
plt.plot_date(t1,bz1,'-b',label='BN',linewidth=0.5)
plt.plot_date(t1,bt1,'-k',label='Btotal')
ax1.set_ylabel('B [nT]')
plt.legend(loc=2)
ax1.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d %H') )
ax1.set_xlim(t1[0],t1[-1])

ax2=plt.subplot(422)
plt.plot_date(t2,bx2,'-r',label='BR', linewidth=0.5)
plt.plot_date(t2,by2,'-g',label='BT', linewidth=0.5)
plt.plot_date(t2,bz2,'-b',label='BN', linewidth=0.5)
plt.plot_date(t2,bt2,'-k',label='Btotal')
ax2.set_ylabel('B [nT]')
plt.legend(loc=2)
ax2.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d %H') )
ax2.set_xlim(t2[0],t2[-1])



ax3=plt.subplot(423, sharex=ax1)
plt.plot_date(t_swe1,v1,'-k',label='V', linewidth=0.8)
ax3.set_ylabel('V [km/s]')
ax3.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d %H') )
ax3.set_xlim(t1[0],t1[-1])



ax4=plt.subplot(424,sharex=ax2)
plt.plot_date(t_swe2,v2,'-k',label='V', linewidth=0.8)
ax4.set_xlim(t2[0],t2[-1])
ax4.set_ylabel('V [km/s]')
ax4.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d %H') )


ax5=plt.subplot(425,sharex=ax1)
plt.plot_date(t1,r1,'-k')
ax5.set_ylabel('R [AU]')
ax5.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d') )
ax5.set_xlim(t1[0],t1[-1])


ax6=plt.subplot(426,sharex=ax2)
plt.plot_date(t2,r2,'-k')
ax6.set_ylabel('R [AU]')
ax6.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d') )
ax6.set_xlim(t2[0],t2[-1])


ax7=plt.subplot(427,sharex=ax1)
plt.plot_date(t1,np.rad2deg(lon1),'-r')
ax7.set_xlim(t1[0],t1[-1])
ax7.set_ylabel('HEEQ longitude [deg]')
ax7.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d') )


ax8=plt.subplot(428,sharex=ax2)
plt.plot_date(t2,np.rad2deg(lon2),'-r')
ax8.set_xlim(t2[0],t2[-1])
ax8.set_ylabel('HEEQ longitude [deg]')
ax8.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d') )


plt.tight_layout()



plt.show()



plt.savefig('results/psp_orbit12.jpg')

sys.exit()

#hd.convert_MAVEN_mat_to_pickle()

mav=hd.load_MAVEN()

plt.plot_date(mav['timeD'],mav['BT'],'-k') 






#smooth with median filter?

print('done')
'''











'''
#ignore warnings
#import warnings
#warnings.filterwarnings('ignore')

#use one time to control plots
#plt.show(block=True)


#CHECK for MAVEN program  for MARS EVENTS at MAVEN -> HI / MSL LIST FOR VERIFICATION

#url_hicat='https://www.helcats-fp7.eu/catalogues/data/HCME_WP3_V05.json'
#download, see URLLIB https://docs.python.org/3/howto/urllib2.html
#with urllib.request.urlopen(url_plasma) as url:
#   pr = json.loads	(url.read().decode())
#hicat_file='cats/HCME_WP3_V05.json'
#h = json.loads(open(hicat_file).read())
##kill first row which stems from the description part
#hmf_speed=np.zeros(len(h['data']))
#fpf_speed=np.zeros(len(h['data']))
#for k in np.arange(0,len(h['data']),1):
#  hmf_speed[k]=h['data'][k][25]
#  fpf_speed[k]=h['data'][k][9]


#for line in open('cats/MSL_event_table_forstner_2019.txt'):
#   line = line.split() # to deal with blank 
#   print(line)
   
#http://docs.astropy.org/en/stable/io/ascii/
#alternative
#from astropy.io import ascii
#data = ascii.read('cats/MSL_event_table_forstner_2019.txt')  


### ORDER: 1 durations, 2 B field, 3 time inside, 4 ICME frequency


plt.close('all')
print()
print('Start cme_stats.py main program.')
print('ICME parameters at all 4 terrestrial planets.')
print('Christian Moestl, IWF Graz, Austria, last update: November 2018.')

######################## get cats


#solar radius
Rs_in_AU=float(const.R_sun/const.au)


filename_icmecat='cats/HELCATS_ICMECAT_v20_SCEQ.sav'
i=getcat(filename_icmecat)

#now this is a scipy recarray  
#access each element of the array see http://docs.scipy.org/doc/numpy/user/basics.rec.html
#access variables
#i.icmecat['id']
#look at contained variables
#print(i.icmecat.dtype)

#get spacecraft and planet positions

pos=getcat('cats/positions_2007_2023_HEEQ_6hours.sav')
pos_time_num=time_to_num_cat(pos.time)[0]




print()
print()

print('Spacecraft positions available:')
print(pos.keys())
print()





########## LOAD OMNI2 data for sunspot data
load_url_current_directory('omni2_all_years.dat', \
        'ftp://nssdcftp.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat')
        
#if omni2 hourly data is not yet converted and saved as pickle, do it:
if not os.path.exists('omni2_all_years_pickle.p'):
 #load OMNI2 dataset from .dat file with a function from dst_module.py
 o=get_omni2_data()
 #contains: o. time,day,hour,btot,bx,by,bz,bygsm,bzgsm,speed,speedx,den,pdyn,dst,kp
 #save for faster loading later
 pickle.dump(o, open('omni2_all_years_pickle.p', 'wb') )
else:  o=pickle.load(open('omni2_all_years_pickle.p', 'rb') )
print('loaded OMNI2 data')





######## load a file with the time and B total of the original in situ data, 
#to be used in part (3)  - has 355 MB to download
#STA, STB, Wind, MES, VEX, MAVEN
load_url_current_directory('insitu_data_time_btot_moestl_2019_paper.p', \
                     'https://oeawcloud.oeaw.ac.at/index.php/s/z7nZVbH8qBjaTf5/download')

#info for myself: file insitu_data_time_btot_moestl_2019_paper.p
#was done with program extract_data in process_data package (not on github)

[win_time,win_btot,sta_time,sta_btot,stb_time,stb_btot, \
mav_time, mav_btot, vex_time, vex_btot, mes_time, mes_btot]= \
pickle.load(open( "insitu_data_time_btot_moestl_2019_paper.p", "rb" ) )

print('loaded time/Btotal data for STEREO-A/B, Wind, VEX, MESSENGER, MAVEN')




########################### get all parameters from ICMECAT for easier handling
# id for each event
iid=i.icmecat['id']
# need to decode all strings 
iid=decode_array(iid)

# observing spacecraft
isc=i.icmecat['sc_insitu'] #string
isc=decode_array(isc)

# all times need to be converted from the IDL format to matplotlib format, 
# also make strings for each date
icme_start_time=i.icmecat['ICME_START_TIME']
[icme_start_time_num,icme_start_time_str]=time_to_num_cat(icme_start_time)

mo_start_time=i.icmecat['MO_START_TIME']
[mo_start_time_num,mo_start_time_str]=time_to_num_cat(mo_start_time)

mo_end_time=i.icmecat['MO_END_TIME']
[mo_end_time_num,mo_end_time_str]=time_to_num_cat(mo_end_time)

# this time exists only for Wind
icme_end_time=i.icmecat['ICME_END_TIME']
[icme_end_time_num,icme_end_time_str]=time_to_num_cat(icme_end_time)

# get variables from ICMECAT
sc_heliodistance=i.icmecat['SC_HELIODISTANCE']
sc_long_heeq=i.icmecat['SC_LONG_HEEQ']
sc_lat_heeq=i.icmecat['SC_LAT_HEEQ']
mo_bmax=i.icmecat['MO_BMAX']
mo_bmean=i.icmecat['MO_BMEAN']
mo_bstd=i.icmecat['MO_BSTD']
mo_bzmean=i.icmecat['MO_BZMEAN']
mo_bzmin=i.icmecat['MO_BZMIN']
mo_duration=i.icmecat['MO_DURATION']
mo_mva_axis_long=i.icmecat['MO_MVA_AXIS_LONG']
mo_mva_axis_lat=i.icmecat['MO_MVA_AXIS_LAT']
mo_mva_ratio=i.icmecat['MO_MVA_RATIO']
sheath_speed=i.icmecat['SHEATH_SPEED']
sheath_speed_std=i.icmecat['SHEATH_SPEED_STD']
mo_speed=i.icmecat['MO_SPEED']
mo_speed_st=i.icmecat['MO_SPEED_STD']
sheath_density=i.icmecat['SHEATH_DENSITY']
sheath_density_std=i.icmecat['SHEATH_DENSITY_STD']
mo_density=i.icmecat['MO_DENSITY']
mo_density_std=i.icmecat['MO_DENSITY_STD']
sheath_temperature=i.icmecat['SHEATH_TEMPERATURE']
sheath_temperature_std=i.icmecat['SHEATH_TEMPERATURE_STD']
mo_temperature=i.icmecat['MO_TEMPERATURE']
mo_temperature_std=i.icmecat['MO_TEMPERATURE_STD']
sheath_pdyn=i.icmecat['SHEATH_PDYN']
sheath_pdyn_std=i.icmecat['SHEATH_PDYN_STD']
mo_pdyn=i.icmecat['MO_PDYN']
mo_pdyn_std=i.icmecat['MO_PDYN_STD']


#calculate icme_durations in hours
icme_durations=(mo_end_time_num-icme_start_time_num)*24 #hours
#same for magnetic obstacle (MO)
mo_durations=(mo_end_time_num-mo_start_time_num)*24 #hours


#get indices of events by different spacecraft
ivexind=np.where(isc == 'VEX')
istaind=np.where(isc == 'STEREO-A')
istbind=np.where(isc == 'STEREO-B')
iwinind=np.where(isc == 'Wind')
imesind=np.where(isc == 'MESSENGER')
iulyind=np.where(isc == 'ULYSSES')
imavind=np.where(isc == 'MAVEN')


#take MESSENGER only at Mercury, only events after orbit insertion
imercind=np.where(np.logical_and(isc =='MESSENGER', \
          icme_start_time_num > mdates.date2num(sunpy.time.parse_time('2011-03-18'))))



############### set limits of solar minimum, rising/declining phase and solar maximum

minstart=mdates.date2num(sunpy.time.parse_time('2007-01-01'))
minend=mdates.date2num(sunpy.time.parse_time('2009-12-31'))

risestart=mdates.date2num(sunpy.time.parse_time('2010-01-01'))
riseend=mdates.date2num(sunpy.time.parse_time('2011-06-30'))

maxstart=mdates.date2num(sunpy.time.parse_time('2011-07-01'))
maxend=mdates.date2num(sunpy.time.parse_time('2014-12-31'))

declstart=mdates.date2num(sunpy.time.parse_time('2015-01-01'))
declend=mdates.date2num(sunpy.time.parse_time('2017-12-31'))


############### extract events by limits of solar min, 
############### rising, max, too few events for MAVEN and Ulysses

iallind_min=np.where(np.logical_and(icme_start_time_num > minstart,icme_start_time_num < minend))[0]
iallind_rise=np.where(np.logical_and(icme_start_time_num > risestart,icme_start_time_num < riseend))[0]
iallind_max=np.where(np.logical_and(icme_start_time_num > maxstart,icme_start_time_num < maxend))[0]

iwinind_min=iallind_min[np.where(isc[iallind_min]=='Wind')]
iwinind_rise=iallind_rise[np.where(isc[iallind_rise]=='Wind')]
iwinind_max=iallind_max[np.where(isc[iallind_max]=='Wind')]

ivexind_min=iallind_min[np.where(isc[iallind_min]=='VEX')]
ivexind_rise=iallind_rise[np.where(isc[iallind_rise]=='VEX')]
ivexind_max=iallind_max[np.where(isc[iallind_max]=='VEX')]

imesind_min=iallind_min[np.where(isc[iallind_min]=='MESSENGER')]
imesind_rise=iallind_rise[np.where(isc[iallind_rise]=='MESSENGER')]
imesind_max=iallind_max[np.where(isc[iallind_max]=='MESSENGER')]

istaind_min=iallind_min[np.where(isc[iallind_min]=='STEREO-A')]
istaind_rise=iallind_rise[np.where(isc[iallind_rise]=='STEREO-A')]
istaind_max=iallind_max[np.where(isc[iallind_max]=='STEREO-A')]

istbind_min=iallind_min[np.where(isc[iallind_min]=='STEREO-B')]
istbind_rise=iallind_rise[np.where(isc[iallind_rise]=='STEREO-B')]
istbind_max=iallind_max[np.where(isc[iallind_max]=='STEREO-B')]


# select the events at Mercury extra after orbit insertion
# no events available for solar minimum!
imercind_min=iallind_min[np.where(np.logical_and(isc[iallind_min] =='MESSENGER',icme_start_time_num[iallind_min] > mdates.date2num(sunpy.time.parse_time('2011-03-18'))))]
imercind_rise=iallind_rise[np.where(np.logical_and(isc[iallind_rise] =='MESSENGER',icme_start_time_num[iallind_rise] > mdates.date2num(sunpy.time.parse_time('2011-03-18'))))]
imercind_max=iallind_max[np.where(np.logical_and(isc[iallind_max] =='MESSENGER',icme_start_time_num[iallind_max] > mdates.date2num(sunpy.time.parse_time('2011-03-18'))))]












print()
print()






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
scx=sc_heliodistance[:,np.newaxis]
durfit_f, _, _, _ =np.linalg.lstsq(scx,icme_durations, rcond=None)

scxmin=sc_heliodistance[iallind_min][:,np.newaxis]
durfitmin_f, _, _, _ =np.linalg.lstsq(scxmin,icme_durations[iallind_min],rcond=None)

scxrise=sc_heliodistance[iallind_rise][:,np.newaxis]
durfitrise_f, _, _, _ =np.linalg.lstsq(scxrise,icme_durations[iallind_rise],rcond=None)

scxmax=sc_heliodistance[iallind_max][:,np.newaxis]
durfitmax_f, _, _, _ =np.linalg.lstsq(scxmax,icme_durations[iallind_max],rcond=None)

#make the y axis for the fits forced through the origin
ydurfitall_f=durfit_f*xfit
ydurfitmin_f=durfitmin_f*xfit
ydurfitrise_f=durfitrise_f*xfit
ydurfitmax_f=durfitmax_f*xfit

plt.plot(sc_heliodistance,icme_durations,'o',color='blue',markersize=4, alpha=0.4)
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

ax1.set_xticks(np.arange(0,2,0.1))
#plt.xlim(0,max(sc_heliodistance)+0.3)
#plt.ylim(0,110)


#for PSP zoom
plt.xlim(0,1.2)
plt.ylim(0,40)


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
ax2.plot_date(icme_start_time_num[imesind],icme_durations[imesind], \
    'o',color='darkgrey',markersize=markers,linestyle='-',linewidth=linew,label='MESSENGER')
ax2.plot_date(icme_start_time_num[ivexind],icme_durations[ivexind], \
    'o',color='orange',markersize=markers,linestyle='-',linewidth=linew, label='Venus')
ax2.plot_date(icme_start_time_num[iwinind],icme_durations[iwinind], \
    'o',color='mediumseagreen',markersize=markers, linestyle='-', linewidth=linew, label='Earth')
ax2.plot_date(icme_start_time_num[imavind],icme_durations[imavind], \
    'o',color='steelblue',markersize=markers,linestyle='-',linewidth=linew, label='Mars')


#limits solar min/rise/maxax2.set_ylim(0,80)
vlevel=130
spanalpha=0.05

plt.axvspan(minstart,minend, color='green', alpha=spanalpha)
plt.annotate('solar minimum',xy=(minstart+(minend-minstart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(minstart+10,vlevel),ha='left')
plt.annotate('>',xy=(minend-10,vlevel),ha='right')

plt.axvspan(risestart,riseend, color='yellow', alpha=spanalpha)
plt.annotate('rising phase',xy=(risestart+(riseend-risestart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(risestart+10,vlevel),ha='left')
plt.annotate('>',xy=(riseend-10,vlevel),ha='right')

plt.axvspan(maxstart,maxend, color='red', alpha=spanalpha)
plt.annotate('solar maximum',xy=(maxstart+(maxend-maxstart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(maxstart+10,vlevel),ha='left')
plt.annotate('>',xy=(maxend,vlevel),ha='right')


#plot means as horizontal lines for each sub interval
plt.plot_date( [minstart,minend], [np.mean(icme_durations[iwinind_min]),np.mean(icme_durations[iwinind_min])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [minstart,minend], [np.mean(icme_durations[ivexind_min]),np.mean(icme_durations[ivexind_min])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [minstart,minend], [np.mean(icme_durations[imesind_min]),np.mean(icme_durations[imesind_min])], color='darkgrey', linestyle='-', markersize=0) 

plt.plot_date( [risestart,riseend], [np.mean(icme_durations[iwinind_rise]),np.mean(icme_durations[iwinind_rise])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [risestart,riseend], [np.mean(icme_durations[ivexind_rise]),np.mean(icme_durations[ivexind_rise])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [risestart,riseend], [np.mean(icme_durations[imesind_rise]),np.mean(icme_durations[imesind_rise])], color='darkgrey', linestyle='-', markersize=0) 

plt.plot_date( [maxstart,maxend], [np.mean(icme_durations[iwinind_max]),np.mean(icme_durations[iwinind_max])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [maxstart,maxend], [np.mean(icme_durations[ivexind_max]),np.mean(icme_durations[ivexind_max])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [maxstart,maxend], [np.mean(icme_durations[imesind_max]),np.mean(icme_durations[imesind_max])], color='darkgrey', linestyle='-', markersize=0) 

plt.xlim(mdates.date2num(sunpy.time.parse_time('2007-01-01')), mdates.date2num(sunpy.time.parse_time('2016-12-31')))
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
plt.savefig('plots_stats/icme_durations_distance_time_paper.pdf', dpi=300)
plt.savefig('plots_stats/icme_durations_distance_time_paper.png', dpi=300)



#Sunspot number
#ax3 = ax2.twinx()
#ax3.plot_date(o.time, o.spot, '-', color='black', alpha=0.5,linewidth=0.2)
#ax3.set_ylabel('Sunspot number')
#ax3.set_ylim([0,220])


#results on durations


 
print()
print()

print('DURATION results, mean +/- std [hours]')

print()
print('Mercury ', round(np.mean(icme_durations[imercind]),1),' +/- ', round(np.std(icme_durations[imercind]),1))
#print('min     ', round(np.mean(icme_durations[imercind_min]),1), ' +/- ', round(np.std(icme_durations[imercind_min]),1))
print('min     no events')
print('rise    ', round(np.mean(icme_durations[imercind_rise]),1), ' +/- ', round(np.std(icme_durations[imercind_rise]),1))
print('max     ', round(np.mean(icme_durations[imercind_max]),1), ' +/- ', round(np.std(icme_durations[imercind_max]),1))

print()
print('Venus   ', round(np.mean(icme_durations[ivexind]),1),' +/- ', round(np.std(icme_durations[ivexind]),1))
print('min     ', round(np.mean(icme_durations[ivexind_min]),1), ' +/- ', round(np.std(icme_durations[ivexind_min]),1))
print('rise    ', round(np.mean(icme_durations[ivexind_rise]),1), ' +/- ', round(np.std(icme_durations[ivexind_rise]),1))
print('max     ', round(np.mean(icme_durations[ivexind_max]),1), ' +/- ', round(np.std(icme_durations[ivexind_max]),1))


print()
print('Earth   ', round(np.mean(icme_durations[iwinind]),1),' +/- ', round(np.std(icme_durations[iwinind]),1))
print('min     ', round(np.mean(icme_durations[iwinind_min]),1), ' +/- ', round(np.std(icme_durations[iwinind_min]),1))
print('rise    ', round(np.mean(icme_durations[iwinind_rise]),1), ' +/- ', round(np.std(icme_durations[iwinind_rise]),1))
print('max     ', round(np.mean(icme_durations[iwinind_max]),1), ' +/- ', round(np.std(icme_durations[iwinind_max]),1))

print()
print('MAVEN only declining phase')
print('decl    ',round(np.mean(icme_durations[imavind]),1),' +/- ', round(np.std(icme_durations[imavind]),1))







print()
print()















#################################### (2) Bfield plot ICMECAT  ############################

print('-------------------------------------------------')
print('2 MO FIELD VS DISTANCE and time, 3 figure panels')
print()


sns.set_context("talk")     
sns.set_style("ticks",{'grid.linestyle': '--'})
fig=plt.figure(2,figsize=(12,12	))
fsize=15
ax1 = plt.subplot2grid((2,2), (0, 0))

# xfit starts here at 1 Rs  because there should not be a 0 in xfit ** bmaxfit[0] later	
xfit=np.linspace(Rs_in_AU,2,1000)


print('Fit results for B in Form: y=B0*x^k')

####### power law fits for all events
bmaxfit=np.polyfit(np.log10(sc_heliodistance),np.log10(mo_bmax),1)
b=10**bmaxfit[1]
bmaxfitfun=b*(xfit**bmaxfit[0])
print('bmax:       ',round(10**bmaxfit[1],2),' x ^', round(bmaxfit[0],2))

bmeanfit=np.polyfit(np.log10(sc_heliodistance),np.log10(mo_bmean),1)
b=10**bmeanfit[1]
bmeanfitfun=b*(xfit**bmeanfit[0])
print('bmean:      ', round(10**bmeanfit[1],2),' x ^',round(bmeanfit[0],2))


##fit with only minimum events
bmeanfit_min=np.polyfit(np.log10(sc_heliodistance[iallind_min]),np.log10(mo_bmean[iallind_min]),1)
bmeanfitfun_min=(10**bmeanfit_min[1])*(xfit**bmeanfit_min[0])
print('bmean_min:  ', round(10**bmeanfit_min[1],2),' x ^', round(bmeanfit_min[0],2))

##fit with only rising events
bmeanfit_rise=np.polyfit(np.log10(sc_heliodistance[iallind_rise]),np.log10(mo_bmean[iallind_rise]),1)
bmeanfitfun_rise=(10**bmeanfit_rise[1])*(xfit**bmeanfit_rise[0])
print('bmean_rise: ', round(10**bmeanfit_rise[1],2),' x ^', round(bmeanfit_rise[0],2))

##fit with only maximum events
bmeanfit_max=np.polyfit(np.log10(sc_heliodistance[iallind_max]),np.log10(mo_bmean[iallind_max]),1)
bmeanfitfun_max=(10**bmeanfit_max[1])*(xfit**bmeanfit_max[0])
print('bmean_max:  ', round(10**bmeanfit_max[1],2),' x ^',round(bmeanfit_max[0],2))




################# plots 2a
plt.plot(sc_heliodistance,mo_bmean,'o',color='black',markersize=5, alpha=0.7,label='$\mathregular{<B>}$')
plt.plot(xfit,bmeanfitfun,'-',color='black', lw=2, alpha=0.7,label='$\mathregular{<B> \\ fit}$')

plt.plot(sc_heliodistance,mo_bmax,'o',color='dodgerblue',markersize=5, alpha=0.7,label='$\mathregular{B_{max}}$')
plt.plot(xfit,bmaxfitfun,'-',color='dodgerblue', lw=2, alpha=0.7,label='$\mathregular{B_{max} \\ fit}$')

#plt.text(1.1,120,'$\mathregular{<B> [nT]= {:.2f} R[AU]^{{:.2f}}}$'.format(10**bmeanfit[1],bmeanfit[0]), fontsize=10)
plt.text(1.1,60,r'<B> [nT]= {:.2f} R[AU]^{:.2f}'.format(10**bmeanfit[1],bmeanfit[0]), fontsize=10)
plt.text(1.1,100,r'Bmax [nT]= {:.2f} R[AU]^{:.2f}'.format(10**bmaxfit[1],bmaxfit[0]), fontsize=10,color='dodgerblue')


#mars limits
plt.axvspan(np.min(pos.mars[0]),np.max(pos.mars[0]), color='orangered', alpha=0.2)
#plt.figtext(0.8,0.8,'Mars',color='orangered')

plt.axvspan(np.min(pos.mercury[0]),np.max(pos.mercury[0]), color='darkgrey', alpha=0.2)
#plt.figtext(0.25,0.8,'Mercury',color='darkgrey')

plt.axvspan(np.min(pos.venus[0]),np.max(pos.venus[0]), color='orange', alpha=0.2)
#plt.figtext(0.42,0.8,'Venus',color='orange')

plt.axvspan(np.min(pos.earth[0]),np.max(pos.earth[0]), color='mediumseagreen', alpha=0.2)
#plt.figtext(0.6,0.8,'Earth',color='mediumseagreen')


#plt.figtext(0.65,0.2,' D[h]={:.2f} R[AU] + {:.2f}'.format(durfit[0],durfit[1]))
plt.xlim(0,1.8)
plt.ylim(0,max(mo_bmax)+20)


plt.legend(loc=1,fontsize=fsize)

plt.xlabel('Heliocentric distance R [AU]',fontsize=fsize)
plt.ylabel('Magnetic field in MO B [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 
#plt.grid()




######################## 2b logarithmic plot with Sun

#for the bmean fit, append one value for the coronal field 
#patsourakos georgoulis 2016: 0.03 G for 10 Rs #10^5 nT is 1 Gauss
mo_bmean_sun=np.append(mo_bmean,10**5*0.03) 
mo_bmax_sun=np.append(mo_bmax,10**5*0.03) 
sc_heliodistance_sun=np.append(sc_heliodistance,10*Rs_in_AU)

print()
bmeanfit_sun=np.polyfit(np.log10(sc_heliodistance_sun),np.log10(mo_bmean_sun),1)
b=10**bmeanfit_sun[1]
bmeanfitfun_sun=b*(xfit**bmeanfit_sun[0])
print('bmean_sun:  ', round(10**bmeanfit_sun[1],2),' x ^',round(bmeanfit_sun[0],2))

bmaxfit_sun=np.polyfit(np.log10(sc_heliodistance_sun),np.log10(mo_bmax_sun),1)
b=10**bmaxfit_sun[1]
bmaxfitfun_sun=b*(xfit**bmaxfit_sun[0])
print('bmax_sun:   ', round(10**bmaxfit_sun[1],2),' x ^',round(bmaxfit_sun[0],2))


ax3 = plt.subplot2grid((2,2), (0, 1))

plt.plot(sc_heliodistance_sun,np.log10(mo_bmean_sun),'o',color='black',markersize=2, alpha=0.7,label='$\mathregular{<B>}$')
plt.plot(sc_heliodistance,np.log10(mo_bmax),'o',color='dodgerblue',markersize=2, alpha=0.7,label='$\mathregular{B_{max}}$')
plt.plot(xfit,np.log10(bmeanfitfun_sun),'-',color='black', lw=2, alpha=0.7,label='$\mathregular{<B> fit}$')
plt.plot(xfit,np.log10(bmaxfitfun_sun),'-',color='dodgerblue', lw=2, alpha=0.7,label='$\mathregular{B_{max} fit}$')
plt.ylim(0.1,4)


#Planet labels and shades
ax3.annotate('Mars', xy=(1.5,2.8), ha='center',fontsize=fsize-2)
ax3.annotate('Mercury', xy=(0.38,2.8), ha='center',fontsize=fsize-2)
ax3.annotate('Venus', xy=(0.72,2.8), ha='center',fontsize=fsize-2)
ax3.annotate('Earth', xy=(1,2.8), ha='center',fontsize=fsize-2)
plt.plot([0.046,0.046],[0,10], color='black', linestyle='--', linewidth=1)
plt.axvspan(np.min(pos.mars[0]),np.max(pos.mars[0]), color='orangered', alpha=0.2)
plt.axvspan(np.min(pos.mercury[0]),np.max(pos.mercury[0]), color='darkgrey', alpha=0.2)
plt.axvspan(np.min(pos.venus[0]),np.max(pos.venus[0]), color='orange', alpha=0.2)
plt.axvspan(np.min(pos.earth[0]),np.max(pos.earth[0]), color='mediumseagreen', alpha=0.2)
#plt.xlim(0,1.8)
#PSP
plt.xlim(0,1.2)

plt.text(0.1,0.2,r'<B> [nT]= {:.2f} R[AU]^{:.2f}'.format(10**bmeanfit[1],bmeanfit[0]), fontsize=10)
plt.text(0.1,0.6,r'Bmax [nT]= {:.2f} R[AU]^{:.2f}'.format(10**bmaxfit[1],bmaxfit[0]), fontsize=10,color='dodgerblue')



plt.legend(loc=1,fontsize=fsize-5)

plt.xlabel('Heliocentric distance R [AU]',fontsize=fsize)
plt.ylabel('MFR Magnetic field log(B) [nT]',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


#panel labels
plt.figtext(0.03,0.96,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.515,0.96,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.03,0.49,'c',color='black', fontsize=fsize, ha='left',fontweight='bold')


plt.text(1.1,60,r'<B> [nT]= {:.2f} R[AU]^{:.2f}'.format(10**bmeanfit[1],bmeanfit[0]), fontsize=10)
plt.text(1.1,100,r'Bmax [nT]= {:.2f} R[AU]^{:.2f}'.format(10**bmaxfit[1],bmaxfit[0]), fontsize=10,color='dodgerblue')



################################# 2c MO B vs. time
ax2 = plt.subplot2grid((2,2), (1, 0), colspan=2)
markers=6
linew=0
plt.plot_date(icme_start_time_num[imesind],mo_bmean[imesind],'o',color='darkgrey',markersize=markers,linestyle='-',linewidth=linew,label='MESSENGER')
plt.plot_date(icme_start_time_num[ivexind],mo_bmean[ivexind],'o',color='orange',markersize=markers,linestyle='-',linewidth=linew, label='Venus')
plt.plot_date(icme_start_time_num[iwinind],mo_bmean[iwinind],'o',color='mediumseagreen',markersize=markers, linestyle='-', linewidth=linew, label='Earth')
plt.plot_date(icme_start_time_num[imavind],mo_bmean[imavind],'o',color='steelblue',markersize=markers,linestyle='-',linewidth=linew, label='Mars')


#NOT USED:
#add gaussian fits for MESSENGER, VEX, Wind (MAVEN too few data points)
#instead of gaussian, fit solar cycle functions in Hathaway 2015 solar cycle living reviews equation 6
#xaxis: 10 years, daily data point
#xfit=mdates.date2num(sunpy.time.parse_time('2007-01-01'))+np.arange(0,365*10)
#MESSENGER
#sigma=1000
#bfitmax=80
#mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
#ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='darkgrey',linestyle='-',markersize=0, label='Mercury fit')
#VEX
#inital guess
#sigma=1000
#bfitmax=40
#mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
#ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='orange',linestyle='-',markersize=0, label='Venus fit')
#Wind
#sigma=1000
#bfitmax=10
#mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
#ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit-mu)**2)/(2*sigma**2) )
#normalize with 1/max(ygauss)
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='mediumseagreen',linestyle='-',markersize=0, label='Earth fit')
#for Mars: reconstruct likely parameters if sigma is quite similar for all fits, take mean of those sigmas and adjust bfitmax as function of distance with power law)
#plot reconstructed function for Mars
#bfitmax=6
#plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='steelblue',linestyle='--',markersize=0, label='Mars reconstr.')


plt.legend(loc=1,fontsize=fsize-2)

#limits solar min/rise/max

vlevel=150

plt.axvspan(minstart,minend, color='green', alpha=0.1)
plt.annotate('solar minimum',xy=(minstart+(minend-minstart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(minstart+10,vlevel),ha='left')
plt.annotate('>',xy=(minend-10,vlevel),ha='right')

plt.axvspan(risestart,riseend, color='yellow', alpha=0.1)
plt.annotate('rising phase',xy=(risestart+(riseend-risestart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(risestart+10,vlevel),ha='left')
plt.annotate('>',xy=(riseend-10,vlevel),ha='right')

plt.axvspan(maxstart,maxend, color='red', alpha=0.1)
plt.annotate('solar maximum',xy=(maxstart+(maxend-maxstart)/2,vlevel),color='black', ha='center')
plt.annotate('<',xy=(maxstart+10,vlevel),ha='left')
plt.annotate('>',xy=(maxend,vlevel),ha='right')


plt.plot_date( [minstart,minend], [np.mean(mo_bmean[iwinind_min]), \
               np.mean(mo_bmean[iwinind_min])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [minstart,minend], [np.mean(mo_bmean[ivexind_min]), \
               np.mean(mo_bmean[ivexind_min])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [minstart,minend], [np.mean(mo_bmean[imesind_min]), \
               np.mean(mo_bmean[imesind_min])], color='darkgrey', linestyle='-', markersize=0) 

plt.plot_date( [risestart,riseend], [np.mean(mo_bmean[iwinind_rise]), \
               np.mean(mo_bmean[iwinind_rise])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [risestart,riseend], [np.mean(mo_bmean[ivexind_rise]), \
               np.mean(mo_bmean[ivexind_rise])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [risestart,riseend], [np.mean(mo_bmean[imesind_rise]), \
               np.mean(mo_bmean[imesind_rise])], color='darkgrey', linestyle='-', markersize=0) 

plt.plot_date( [maxstart,maxend], [np.mean(mo_bmean[iwinind_max]), \
               np.mean(mo_bmean[iwinind_max])], color='mediumseagreen', linestyle='-',markersize=0 ) 
plt.plot_date( [maxstart,maxend], [np.mean(mo_bmean[ivexind_max]), \
               np.mean(mo_bmean[ivexind_max])], color='orange', linestyle='-', markersize=0) 
plt.plot_date( [maxstart,maxend], [np.mean(mo_bmean[imesind_max]), \
               np.mean(mo_bmean[imesind_max])], color='darkgrey', linestyle='-', markersize=0) 


plt.ylabel('Magnetic field in MO [nT]', fontsize=fsize)
plt.xlabel('Year', fontsize=fsize)

fsize=14

plt.xlim(mdates.date2num(sunpy.time.parse_time('2007-01-01')), mdates.date2num(sunpy.time.parse_time('2016-12-31')))

plt.ylabel('Magnetic field in MO [nT]', fontsize=fsize)
plt.xlabel('Year', fontsize=fsize)
plt.tight_layout()
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 
plt.legend(loc=1,fontsize=fsize-1)



plt.tight_layout()
plt.savefig('plots_stats/icme_total_field_distance_time_paper.pdf', dpi=300)
plt.savefig('plots_stats/icme_total_field_distance_time_paper.png', dpi=300)



########### RESULTS

print()
print()

print('Magnetic field B MO_BMEAN results, mean +/- std [nT]')

print()
print('Mercury ', round(np.mean(mo_bmean[imercind]),1),' +/- ', round(np.std(mo_bmean[imercind]),1))
#print('min     ', round(np.mean(mo_bmean[imercind_min]),1), ' +/- ', round(np.std(mo_bmean[imercind_min]),1))
print('min      no events')
print('rise    ', round(np.mean(mo_bmean[imercind_rise]),1), ' +/- ', round(np.std(mo_bmean[imercind_rise]),1))
print('max     ', round(np.mean(mo_bmean[imercind_max]),1), ' +/- ', round(np.std(mo_bmean[imercind_max]),1))

print()
print('Venus   ', round(np.mean(mo_bmean[ivexind]),1),' +/- ', round(np.std(mo_bmean[ivexind]),1))
print('min     ', round(np.mean(mo_bmean[ivexind_min]),1), ' +/- ', round(np.std(mo_bmean[ivexind_min]),1))
print('rise    ', round(np.mean(mo_bmean[ivexind_rise]),1), ' +/- ', round(np.std(mo_bmean[ivexind_rise]),1))
print('max     ', round(np.mean(mo_bmean[ivexind_max]),1), ' +/- ', round(np.std(mo_bmean[ivexind_max]),1))


print()
print('Earth   ', round(np.mean(mo_bmean[iwinind]),1),' +/- ', round(np.std(mo_bmean[iwinind]),1))
print('min     ', round(np.mean(mo_bmean[iwinind_min]),1), ' +/- ', round(np.std(mo_bmean[iwinind_min]),1))
print('rise    ', round(np.mean(mo_bmean[iwinind_rise]),1), ' +/- ', round(np.std(mo_bmean[iwinind_rise]),1))
print('max     ', round(np.mean(mo_bmean[iwinind_max]),1), ' +/- ', round(np.std(mo_bmean[iwinind_max]),1))

print()
print('MAVEN only declining phase')
print('decl    ',round(np.mean(mo_bmean[imavind]),1),' +/- ', round(np.std(mo_bmean[imavind]),1))





















































############### (3) TIME SPENT INSIDE CMEs in %, as function of solar cycle ##############
print()
print('-------------------------------------------------')
print('3 Time spent inside CMEs for each planet, 3 figure panels')
print()


#make bin for each year for yearly histograms

yearly_start_times=[mdates.date2num(sunpy.time.parse_time('2007-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2008-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2009-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2010-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2011-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2012-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2013-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2014-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2015-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2016-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2017-01-01'))]
                  
yearly_end_times=[mdates.date2num(sunpy.time.parse_time('2007-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2008-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2009-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2010-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2011-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2012-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2013-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2014-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2015-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2016-12-31')),
                  mdates.date2num(sunpy.time.parse_time('2017-12-31'))]

yearly_mid_times=[mdates.date2num(sunpy.time.parse_time('2007-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2008-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2009-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2010-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2011-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2012-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2013-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2014-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2015-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2016-07-01')),
                  mdates.date2num(sunpy.time.parse_time('2017-07-01'))]


########### 1a Determine SOLAR WIND ONLY MISSIONS STA, STB, Wind
########### all are available in 1 minute time resolution 

#define arrays and fill with nan
total_data_days_yearly_win=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_win.fill(np.nan)

total_data_days_yearly_sta=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_sta.fill(np.nan)

total_data_days_yearly_stb=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_stb.fill(np.nan)

total_data_days_yearly_mes=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_mes.fill(np.nan)

total_data_days_yearly_merc=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_merc.fill(np.nan)

total_data_days_yearly_vex=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_vex.fill(np.nan)

total_data_days_yearly_mav=np.zeros(np.size(yearly_mid_times))
total_data_days_yearly_mav.fill(np.nan)


#go through each year and search for available data
#time is available for all dates, so no NaNs in time
#all NaNs in Btotal variable

for i in range(np.size(yearly_mid_times)):

  #Wind index for this year  
  thisyear=np.where(np.logical_and((win_time > yearly_start_times[i]),(win_time < yearly_end_times[i])))
  #get np.size of available data for each year
  datas=np.size(np.where(np.isnan(win_btot[thisyear])==False))
  #wind is  in 1 minute resolution
  min_in_days=1/(60*24)
  #calculate available days from number of datapoints (each 1 minute) 
  #divided by number of minutes in 1 days
  #this should only be the case if data is available this year, otherwise set to NaN
  if datas >0: total_data_days_yearly_win[i]=datas*min_in_days

  #same for STEREO-A
  thisyear=np.where(np.logical_and((sta_time > yearly_start_times[i]),(sta_time < yearly_end_times[i])))
  datas=np.size(np.where(np.isnan(sta_btot[thisyear])==False))
  if datas >0: total_data_days_yearly_sta[i]=datas*min_in_days

  #same for STEREO-B
  thisyear=np.where(np.logical_and((stb_time > yearly_start_times[i]),(stb_time < yearly_end_times[i])))
  datas=np.size(np.where(np.isnan(stb_btot[thisyear])==False))
  if datas >0: total_data_days_yearly_stb[i]=datas*min_in_days

  #same for MESSENGER
  thisyear=np.where(np.logical_and((mes_time > yearly_start_times[i]),(mes_time < yearly_end_times[i])))
  datas=np.size(np.where(np.isnan(mes_btot[thisyear])==False))
  if datas >0: total_data_days_yearly_mes[i]=datas*min_in_days

  #same for Mercury alone
  #start with 2011
  if i == 4: 
   thisyear=np.where(np.logical_and((mes_time > mdates.date2num(sunpy.time.parse_time('2011-03-18'))),(mes_time < yearly_end_times[i])))
   datas=np.size(np.where(np.isnan(mes_btot[thisyear])==False))
   if datas >0: total_data_days_yearly_merc[i]=datas*min_in_days
  #2012 onwards 
  if i > 4: 
   thisyear=np.where(np.logical_and((mes_time > yearly_start_times[i]),(mes_time < yearly_end_times[i])))
   datas=np.size(np.where(np.isnan(mes_btot[thisyear])==False))
   if datas >0: total_data_days_yearly_merc[i]=datas*min_in_days
   
   


  #same for VEX
  thisyear=np.where(np.logical_and((vex_time > yearly_start_times[i]),(vex_time < yearly_end_times[i])))
  datas=np.size(np.where(np.isnan(vex_btot[thisyear])==False))
  if datas >0: total_data_days_yearly_vex[i]=datas*min_in_days

  #for MAVEN different time resolution
  thisyear=np.where(np.logical_and((mav_time > yearly_start_times[i]),(mav_time < yearly_end_times[i])))
  datas=np.size(np.where(np.isnan(mav_btot[thisyear])==False))
  datas_ind=np.where(np.isnan(mav_btot[thisyear])==False)
  #sum all time intervals for existing data points, but avoid counting gaps where diff is > 1 orbit (0.25 days)
  alldiff=np.diff(mav_time[datas_ind])
  smalldiff_ind=np.where(alldiff <0.25)  
  if datas >0: total_data_days_yearly_mav[i]=np.sum(alldiff[smalldiff_ind])

print('Data days each year:')

print()
print('Wind')
print(np.round(total_data_days_yearly_win,1))
print('STA')
print(np.round(total_data_days_yearly_sta,1))
print('STB')
print(np.round(total_data_days_yearly_stb,1))
print('MERCURY')
print(np.round(total_data_days_yearly_merc,1))
print('MES')
print(np.round(total_data_days_yearly_mes,1))
print('VEX')
print(np.round(total_data_days_yearly_vex,1))
print('MAV')
print(np.round(total_data_days_yearly_mav,1))



















############################### 3a get time inside ICME percentage for full time range


###############
if not os.path.exists('data_icme_indices_moestl_2019_paper.p'):

 ###### check for each spacecraft and all data points if it is inside an ICME
 # save results as pickle because it takes a minute to calculate

 #get all win_btot correct data indices -> 
 win_data_ind=np.where(np.isnan(win_btot)==False)
 win_icme_ind=np.int32(0) #needs to be integer because these will be array indices
 #Wind: go through each icme
 for i in np.arange(np.size(icme_start_time_num[iwinind])): 
    this_icme_ind=np.where(np.logical_and( (win_time[win_data_ind] > icme_start_time_num[iwinind][i]),\
                                            win_time[win_data_ind] < icme_end_time_num[iwinind][i] ))[0]
    win_icme_ind=np.append(win_icme_ind,this_icme_ind)	


 sta_data_ind=np.where(np.isnan(sta_btot)==False)
 sta_icme_ind=np.int32(0)
 for i in np.arange(np.size(icme_start_time_num[istaind])): 
    this_icme_ind=np.where(np.logical_and( (sta_time[sta_data_ind] > icme_start_time_num[istaind][i]),\
                                            sta_time[sta_data_ind] < mo_end_time_num[istaind][i] ))[0]
    sta_icme_ind=np.append(sta_icme_ind,this_icme_ind)	

 stb_data_ind=np.where(np.isnan(stb_btot)==False)
 stb_icme_ind=np.int32(0)
 for i in np.arange(np.size(icme_start_time_num[istbind])): 
    this_icme_ind=np.where(np.logical_and( (stb_time[stb_data_ind] > icme_start_time_num[istbind][i]),\
                                            stb_time[stb_data_ind] < mo_end_time_num[istbind][i] ))[0]
    stb_icme_ind=np.append(stb_icme_ind,this_icme_ind)	

 vex_data_ind=np.where(np.isnan(vex_btot)==False)
 vex_icme_ind=np.int32(0)
 for i in np.arange(np.size(icme_start_time_num[ivexind])): 
    this_icme_ind=np.where(np.logical_and( (vex_time[vex_data_ind] > icme_start_time_num[ivexind][i]),\
                                            vex_time[vex_data_ind] < mo_end_time_num[ivexind][i] ))[0]
    vex_icme_ind=np.append(vex_icme_ind,this_icme_ind)	

 mes_data_ind=np.where(np.isnan(mes_btot)==False)
 mes_icme_ind=np.int32(0)
 for i in np.arange(np.size(icme_start_time_num[imesind])): 
    this_icme_ind=np.where(np.logical_and( (mes_time[mes_data_ind] > icme_start_time_num[imesind][i]),\
                                            mes_time[mes_data_ind] < mo_end_time_num[imesind][i] ))[0]
    mes_icme_ind=np.append(mes_icme_ind,this_icme_ind)	

 mav_data_ind=np.where(np.isnan(mav_btot)==False)
 mav_icme_ind=np.int32(0)
 for i in np.arange(np.size(icme_start_time_num[imavind])): 
    this_icme_ind=np.where(np.logical_and( (mav_time[mav_data_ind] > icme_start_time_num[imavind][i]),\
                                            mav_time[mav_data_ind] < mo_end_time_num[imavind][i] ))[0]
    mav_icme_ind=np.append(mav_icme_ind,this_icme_ind)	


 merc_data_ind=np.where(np.logical_and(np.isnan(mes_btot)==False,mes_time > \
                       mdates.date2num(sunpy.time.parse_time('2011-03-18'))))
 merc_icme_ind=np.int32(0)
 for i in np.arange(np.size(icme_start_time_num[imercind])): 
    this_icme_ind=np.where(np.logical_and( (mes_time[merc_data_ind] > icme_start_time_num[imercind][i]),\
                                            mes_time[merc_data_ind] < mo_end_time_num[imercind][i] ))[0]
    merc_icme_ind=np.append(merc_icme_ind,this_icme_ind)	
 

 pickle.dump([win_icme_ind,win_data_ind, sta_icme_ind,sta_data_ind, stb_icme_ind,stb_data_ind,  vex_icme_ind,vex_data_ind, mes_icme_ind,mes_data_ind,
 merc_icme_ind,merc_data_ind, mav_icme_ind, mav_data_ind],open( "data_icme_indices_moestl_2019_paper.p", "wb" ) )
################

[win_icme_ind,win_data_ind, sta_icme_ind,sta_data_ind, stb_icme_ind,stb_data_ind,  vex_icme_ind,vex_data_ind, mes_icme_ind,mes_data_ind,\
 merc_icme_ind,merc_data_ind,mav_icme_ind, mav_data_ind]= pickle.load(open( "data_icme_indices_moestl_2019_paper.p", "rb" ) )

print()
print()
print('Percentage of time inside ICMEs average over full time range')    
print()
print('Mercury MESSENGER:******* ',np.round((np.size(merc_icme_ind)/np.size(merc_data_ind)*100),1))
print('Venus VEX:         ',np.round((np.size(vex_icme_ind)/np.size(vex_data_ind)*100),1))
print('Earth Wind:        ',np.round((np.size(win_icme_ind)/np.size(win_data_ind)*100),1))
print('Mars MAVEN:        ',np.round((np.size(mes_icme_ind)/np.size(mes_data_ind)*100),1))
print()
print('MESSENGER:         ',np.round((np.size(mes_icme_ind)/np.size(mes_data_ind)*100),1))
print('STB:               ',np.round((np.size(stb_icme_ind)/np.size(stb_data_ind)*100),1))
print('STA:               ',np.round((np.size(sta_icme_ind)/np.size(sta_data_ind)*100),1))









################## 3b get time inside ICME percentage for yearly time range, using results for arrays from 3a



################################# make array for time inside percentages

inside_win_perc=np.zeros(np.size(yearly_mid_times))
inside_win_perc.fill(np.nan)

inside_sta_perc=np.zeros(np.size(yearly_mid_times))
inside_sta_perc.fill(np.nan)

inside_stb_perc=np.zeros(np.size(yearly_mid_times))
inside_stb_perc.fill(np.nan)

inside_mes_perc=np.zeros(np.size(yearly_mid_times))
inside_mes_perc.fill(np.nan)

inside_merc_perc=np.zeros(np.size(yearly_mid_times))
inside_merc_perc.fill(np.nan)

inside_vex_perc=np.zeros(np.size(yearly_mid_times))
inside_vex_perc.fill(np.nan)

inside_mav_perc=np.zeros(np.size(yearly_mid_times))
inside_mav_perc.fill(np.nan)



#count ratio of datapoint inside ICME to all available datapoints for each year, Wind, STA, STB, Mercury, MESSENGER, VEX, MAVEN

for i in range(np.size(yearly_mid_times)):

    thisyear_icme=np.where(np.logical_and((win_time[win_icme_ind] > yearly_start_times[i]),(win_time[win_icme_ind] < yearly_end_times[i])))
    thisyear_data=np.where(np.logical_and((win_time[win_data_ind] > yearly_start_times[i]),(win_time[win_data_ind] < yearly_end_times[i])))
    if np.size(thisyear_data) >0:inside_win_perc[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((sta_time[sta_icme_ind] > yearly_start_times[i]),(sta_time[sta_icme_ind] < yearly_end_times[i])))
    thisyear_data=np.where(np.logical_and((sta_time[sta_data_ind] > yearly_start_times[i]),(sta_time[sta_data_ind] < yearly_end_times[i])))
    if np.size(thisyear_data) >0:inside_sta_perc[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((stb_time[stb_icme_ind] > yearly_start_times[i]),(stb_time[stb_icme_ind] < yearly_end_times[i])))
    thisyear_data=np.where(np.logical_and((stb_time[stb_data_ind] > yearly_start_times[i]),(stb_time[stb_data_ind] < yearly_end_times[i])))
    if np.size(thisyear_data) >0:inside_stb_perc[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((mes_time[mes_icme_ind] > yearly_start_times[i]),(mes_time[mes_icme_ind] < yearly_end_times[i])))
    thisyear_data=np.where(np.logical_and((mes_time[mes_data_ind] > yearly_start_times[i]),(mes_time[mes_data_ind] < yearly_end_times[i])))
    if np.size(thisyear_data) >0:inside_mes_perc[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((vex_time[vex_icme_ind] > yearly_start_times[i]),(vex_time[vex_icme_ind] < yearly_end_times[i])))
    thisyear_data=np.where(np.logical_and((vex_time[vex_data_ind] > yearly_start_times[i]),(vex_time[vex_data_ind] < yearly_end_times[i])))
    if np.size(thisyear_data) >0:inside_vex_perc[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((mes_time[merc_icme_ind] > yearly_start_times[i]),(mes_time[merc_icme_ind] < yearly_end_times[i])))
    thisyear_data=np.where(np.logical_and((mes_time[merc_data_ind] > yearly_start_times[i]),(mes_time[merc_data_ind] < yearly_end_times[i])))
    if np.size(thisyear_data) >0: inside_merc_perc[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((mav_time[mav_icme_ind] > yearly_start_times[i]),(mav_time[mav_icme_ind] < yearly_end_times[i])))
    thisyear_data=np.where(np.logical_and((mav_time[mav_data_ind] > yearly_start_times[i]),(mav_time[mav_data_ind] < yearly_end_times[i])))
    if np.size(thisyear_data) >0: inside_mav_perc[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)


print()
print()
print('Percentage of time inside ICMEs for each year')    
print()
print('Mercury MESSENGER: ',inside_merc_perc)
print('Venus VEX:         ',inside_vex_perc)
print('Earth Wind:        ',inside_win_perc)
print('Mars MAVEN:        ',inside_mav_perc)
print()
print('MESSENGER:         ',inside_mes_perc)
print('STB:               ',inside_stb_perc)
print('STA:               ',inside_sta_perc)








####################make the same thing not yearly, but for the 3 solar cycle phases


cycle_start_times=[minstart, risestart, maxstart]
cycle_end_times=[minend, riseend, maxend]

############################################


inside_win_cycle=np.zeros(np.size(cycle_start_times))
inside_win_cycle.fill(np.nan)

inside_sta_cycle=np.zeros(np.size(cycle_start_times))
inside_sta_cycle.fill(np.nan)

inside_stb_cycle=np.zeros(np.size(cycle_start_times))
inside_stb_cycle.fill(np.nan)

inside_mes_cycle=np.zeros(np.size(cycle_start_times))
inside_mes_cycle.fill(np.nan)

inside_merc_cycle=np.zeros(np.size(cycle_start_times))
inside_merc_cycle.fill(np.nan)

inside_vex_cycle=np.zeros(np.size(cycle_start_times))
inside_vex_cycle.fill(np.nan)

inside_mav_cycle=np.zeros(np.size(cycle_start_times))
inside_mav_cycle.fill(np.nan)


#count ratio of datapoint inside ICME to all available datapoints for each year, Wind, STA, STB, Mercury, MESSENGER, VEX, MAVEN

for i in range(np.size(cycle_start_times)):

    thisyear_icme=np.where(np.logical_and((win_time[win_icme_ind] > cycle_start_times[i]),(win_time[win_icme_ind] < cycle_end_times[i])))
    thisyear_data=np.where(np.logical_and((win_time[win_data_ind] > cycle_start_times[i]),(win_time[win_data_ind] < cycle_end_times[i])))
    if np.size(thisyear_data) >0:inside_win_cycle[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((sta_time[sta_icme_ind] > cycle_start_times[i]),(sta_time[sta_icme_ind] < cycle_end_times[i])))
    thisyear_data=np.where(np.logical_and((sta_time[sta_data_ind] > cycle_start_times[i]),(sta_time[sta_data_ind] < cycle_end_times[i])))
    if np.size(thisyear_data) >0:inside_sta_cycle[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((stb_time[stb_icme_ind] > cycle_start_times[i]),(stb_time[stb_icme_ind] < cycle_end_times[i])))
    thisyear_data=np.where(np.logical_and((stb_time[stb_data_ind] > cycle_start_times[i]),(stb_time[stb_data_ind] < cycle_end_times[i])))
    if np.size(thisyear_data) >0:inside_stb_cycle[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((mes_time[mes_icme_ind] > cycle_start_times[i]),(mes_time[mes_icme_ind] < cycle_end_times[i])))
    thisyear_data=np.where(np.logical_and((mes_time[mes_data_ind] > cycle_start_times[i]),(mes_time[mes_data_ind] < cycle_end_times[i])))
    if np.size(thisyear_data) >0:inside_mes_cycle[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((vex_time[vex_icme_ind] > cycle_start_times[i]),(vex_time[vex_icme_ind] < cycle_end_times[i])))
    thisyear_data=np.where(np.logical_and((vex_time[vex_data_ind] > cycle_start_times[i]),(vex_time[vex_data_ind] < cycle_end_times[i])))
    if np.size(thisyear_data) >0:inside_vex_cycle[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((mes_time[merc_icme_ind] > cycle_start_times[i]),(mes_time[merc_icme_ind] < cycle_end_times[i])))
    thisyear_data=np.where(np.logical_and((mes_time[merc_data_ind] > cycle_start_times[i]),(mes_time[merc_data_ind] < cycle_end_times[i])))
    if np.size(thisyear_data) >0: inside_merc_cycle[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)

    thisyear_icme=np.where(np.logical_and((mav_time[mav_icme_ind] > cycle_start_times[i]),(mav_time[mav_icme_ind] < cycle_end_times[i])))
    thisyear_data=np.where(np.logical_and((mav_time[mav_data_ind] > cycle_start_times[i]),(mav_time[mav_data_ind] < cycle_end_times[i])))
    if np.size(thisyear_data) >0: inside_mav_cycle[i]=round(np.size(thisyear_icme)/np.size(thisyear_data)*100,1)


print()
print()
print('Percentage of time inside ICMEs for min, rise, max:')    
print()
print('Mercury MESSENGER: ',inside_merc_cycle)
print('Venus VEX:         ',inside_vex_cycle)
print('Earth Wind:        ',inside_win_cycle)
print('Mars MAVEN:        ',inside_mav_cycle)
print()
print('MESSENGER:         ',inside_mes_cycle)
print('STB:               ',inside_stb_cycle)
print('STA:               ',inside_sta_cycle)






sys.exit()
  
  

### fix that VEX MESSENGER impact frequency is less than 1 AU by multiplying with a factor of 1.5
#check exact values with frequency plot

#inside_vex_perc=inside_vex_perc*1.5
#inside_mes_perc=inside_mes_perc*1.5


sns.set_context("talk")     
sns.set_style("ticks",{'grid.linestyle': '--'})

fig=plt.figure(5,figsize=(10,10	))

ax1 = plt.subplot(211) 

plt.plot_date(yearly_mid_times,inside_win_perc,'o',color='mediumseagreen',markersize=8, linestyle='-')
plt.plot_date(yearly_mid_times,inside_merc_perc,'o',color='darkgrey',markersize=8,linestyle='-')
plt.plot_date(yearly_mid_times,inside_vex_perc,'o',color='orange',markersize=8,linestyle='-')
plt.plot_date(yearly_mid_times,inside_stb_perc,'o',color='royalblue',markersize=8,linestyle='-')
plt.plot_date(yearly_mid_times,inside_sta_perc,'o',color='red',markersize=8,linestyle='-')
plt.plot_date(yearly_mid_times,inside_mav_perc,'o',color='steelblue',markersize=8,linestyle='-')

plt.ylabel('Time inside ICME [%]')

#plt.xlim(yearly_bin_edges[0],yearly_bin_edges[10])
ax1.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)

#sets planet / spacecraft labels
xoff=0.15
yoff=0.85
fsize=14

plt.figtext(xoff,yoff-0.03*0,'Mercury',color='darkgrey', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*1,'Venus',color='orange', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*2,'Earth',color='mediumseagreen', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*3,'Mars',color='steelblue', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*4,'STEREO-A',color='red', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*5,'STEREO-B',color='royalblue', fontsize=fsize, ha='left')
#panel labels
plt.figtext(0.02,0.98,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.02,0.48,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')



#limits solar min/rise/max

vlevel=22
fsize=11

plt.axvspan(minstart,minend, color='green', alpha=0.1)
plt.annotate('solar minimum',xy=(minstart+(minend-minstart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(minstart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(minend-10,vlevel),ha='right', fontsize=fsize)


plt.axvspan(risestart,riseend, color='yellow', alpha=0.1)
plt.annotate('rising phase',xy=(risestart+(riseend-risestart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(risestart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(riseend-10,vlevel),ha='right', fontsize=fsize)

plt.axvspan(maxstart,maxend, color='red', alpha=0.1)
plt.annotate('solar maximum',xy=(maxstart+(maxend-maxstart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(maxstart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(maxend,vlevel),ha='right', fontsize=fsize)


plt.ylim((0,25))
fsize=15
plt.ylabel('Time inside ICME [%]')
plt.xlabel('Year',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


plt.tight_layout()










#plt.ylim(0,45)
plt.xlim(yearly_start_times[0],yearly_end_times[9])

#sns.despine()




#### plot time inside vs. heliocentric distance

pos_wind_perc=np.zeros(np.size(yearly_mid_times))
pos_wind_perc.fill(np.nan)
pos_wind_perc_std=np.zeros(np.size(yearly_mid_times))
pos_wind_perc_std.fill(np.nan)

pos_sta_perc=np.zeros(np.size(yearly_mid_times))
pos_sta_perc.fill(np.nan)
pos_sta_perc_std=np.zeros(np.size(yearly_mid_times))
pos_sta_perc_std.fill(np.nan)

pos_stb_perc=np.zeros(np.size(yearly_mid_times))
pos_stb_perc.fill(np.nan)
pos_stb_perc_std=np.zeros(np.size(yearly_mid_times))
pos_stb_perc_std.fill(np.nan)

#pos_mes_perc=np.zeros(np.size(yearly_mid_times))
#pos_mes_perc.fill(np.nan)
#pos_mes_perc_std=np.zeros(np.size(yearly_mid_times))
#pos_mes_perc_std.fill(np.nan)


pos_merc_perc=np.zeros(np.size(yearly_mid_times))
pos_merc_perc.fill(np.nan)
pos_merc_perc_std=np.zeros(np.size(yearly_mid_times))
pos_merc_perc_std.fill(np.nan)

pos_vex_perc=np.zeros(np.size(yearly_mid_times))
pos_vex_perc.fill(np.nan)
pos_vex_perc_std=np.zeros(np.size(yearly_mid_times))
pos_vex_perc_std.fill(np.nan)

pos_mav_perc=np.zeros(np.size(yearly_mid_times))
pos_mav_perc.fill(np.nan)
pos_mav_perc_std=np.zeros(np.size(yearly_mid_times))
pos_mav_perc_std.fill(np.nan)


allpositions=np.zeros([np.size(yearly_mid_times), 6])
allinside=np.zeros([np.size(yearly_mid_times), 6])

#calculate average distance +/- std for each year
#go through each year 
for i in range(np.size(yearly_mid_times)):
  
  #select those positions that are inside the current year
  thisyear=np.where(np.logical_and((pos_time_num > yearly_start_times[i]),(pos_time_num < yearly_end_times[i])))
  
  #pos_mes_perc[i]=np.mean(pos.messenger[0][thisyear])
  #pos_mes_perc_std[i]=np.std(pos.messenger[0][thisyear])
  pos_merc_perc[i]=np.mean(pos.mercury[0][thisyear])
  pos_merc_perc_std[i]=np.std(pos.mercury[0][thisyear])
  

  pos_mav_perc[i]=np.mean(pos.mars[0][thisyear])
  pos_mav_perc_std[i]=np.std(pos.mars[0][thisyear])

  pos_vex_perc[i]=np.mean(pos.venus[0][thisyear])
  pos_vex_perc_std[i]=np.std(pos.venus[0][thisyear])

  pos_wind_perc[i]=np.mean(pos.earth_l1[0][thisyear])
  pos_wind_perc_std[i]=np.std(pos.earth_l1[0][thisyear])

  pos_sta_perc[i]=np.mean(pos.sta[0][thisyear])
  pos_sta_perc_std[i]=np.std(pos.sta[0][thisyear])

  pos_stb_perc[i]=np.mean(pos.stb[0][thisyear])
  pos_stb_perc_std[i]=np.std(pos.stb[0][thisyear])
  
  allpositions[i][:]=(pos_merc_perc[i], pos_mav_perc[i], pos_vex_perc[i],pos_wind_perc[i],pos_sta_perc[i],pos_stb_perc[i])
  allinside[i][:]=(inside_merc_perc[i], inside_mav_perc[i], inside_vex_perc[i],inside_win_perc[i],inside_sta_perc[i],inside_stb_perc[i])
  
 
  



#***make alpha variable for each year?

ax3 = plt.subplot(212) 


#for every year linear fit **check if power law works better

#for fit plotting
xfit=np.linspace(0,2,1000)

#allpositions[i] and allinside[i] are the data for each year
#no fit for 2016 as only MAVEN data is available


for i in range(np.size(yearly_mid_times)-2):
 #make linear fits ignoring NaN
 notnan=np.where(np.isfinite(allinside[i]) > 0)
 durfit=np.polyfit(allpositions[i][notnan],allinside[i][notnan],1)
 #this is similar to D=durfit[0]*xfit+durfit[1]
 durfitfun=np.poly1d(durfit)
 print('year',i+2007)
 print('time inside linear fit: D[hours]={:.2f}r[AU]+{:.2f}'.format(durfit[0],durfit[1]))
 plt.plot(xfit,durfitfun(xfit),'-',color='black', lw=2, alpha=i/10+0.2)#,label='fit')
 
 plt.errorbar(pos_merc_perc[i], inside_merc_perc[i], xerr=pos_merc_perc_std[i],yerr=0,fmt='o',color='darkgrey',markersize=8,linestyle=' ',alpha=i/10+0.2)
 plt.errorbar(pos_mav_perc[i], inside_mav_perc[i],xerr=pos_mav_perc_std[i],fmt='o',color='steelblue',markersize=8,linestyle=' ',alpha=i/10+0.2)
 plt.errorbar(pos_sta_perc[i], inside_sta_perc[i],xerr=pos_sta_perc_std[i],fmt='o',color='red',markersize=8,linestyle=' ',alpha=i/10+0.2)
 plt.errorbar(pos_stb_perc[i], inside_stb_perc[i],xerr=pos_stb_perc_std[i],fmt='o',color='royalblue',markersize=8,linestyle=' ',alpha=i/10+0.2)
 plt.errorbar(pos_wind_perc[i], inside_win_perc[i],xerr=pos_wind_perc_std[i],fmt='o',color='mediumseagreen',markersize=8, linestyle=' ',alpha=i/10+0.2)
 plt.errorbar(pos_vex_perc[i], inside_vex_perc[i],xerr=pos_vex_perc_std[i],fmt='o',color='orange',markersize=8,linestyle=' ',alpha=i/10+0.2)
 
 plt.annotate(str(i+2007), xy=(0.1,5+2.5*i), alpha=i/10+0.2)
 
 
 #reconstruct Mars time inside from linear fits but not for years 2015 /2016
 if i < 8: inside_mav_perc[i]=durfitfun(pos_mav_perc[i])


#mars limits
plt.axvspan(np.min(pos.mars[0]),np.max(pos.mars[0]), color='orangered', alpha=0.2)
#plt.figtext(0.8,0.8,'Mars',color='orangered')
plt.axvspan(np.min(pos.mercury[0]),np.max(pos.mercury[0]), color='darkgrey', alpha=0.2)
#plt.figtext(0.25,0.8,'Mercury',color='darkgrey')
plt.axvspan(np.min(pos.venus[0]),np.max(pos.venus[0]), color='orange', alpha=0.2)
#plt.figtext(0.42,0.8,'Venus',color='orange')
plt.axvspan(np.min(pos.earth[0]),np.max(pos.earth[0]), color='mediumseagreen', alpha=0.2)
#plt.figtext(0.6,0.8,'Earth',color='mediumseagreen')
plt.xlim(0,1.8)

#solar probe plus 10 to 36 Rs close approaches

#plt.axvspan(Rs_in_AU*10,Rs_in_AU*36,color='magenta', alpha=0.2)

plt.ylabel('Time inside ICME [%]')
plt.xlabel('Heliocentric distance [AU]')
ax3.set_xticks(np.arange(0,2,0.2))

#add reconstructed Mars time inside on previous plot
ax1.plot_date(yearly_mid_times,inside_mav_perc,'o',color='steelblue',markersize=8,linestyle='--')


plt.ylim((0,25))

plt.tight_layout()

plt.show()
plt.savefig('plots_stats/time_inside_CMEs_paper.pdf', dpi=300)
plt.savefig('plots_stats/time_inside_CMEs_paper.png', dpi=300)





#TI

print()
print()
print('--------------------------------------------------')
print()
print()

print('Time Inside')

print()
print('Mercury +/-')
print(round(np.nanmean(inside_merc_cycle),1))
print(round(np.nanstd(inside_merc_cycle),1))
#print('min')
#print(round(inside_merc_cycle[0],1))
print('rise')
print(round(inside_merc_cycle[1],1))
print('max')
print(round(inside_merc_cycle[2],1))


print()
print('Venus +/-')
print(round(np.nanmean(inside_vex_cycle),1))
print(round(np.nanstd(inside_vex_cycle),1))
print('min')
print(round(inside_vex_cycle[0],1))
print('rise')
print(round(inside_vex_cycle[1],1))
print('max')
print(round(inside_vex_cycle[2],1))


print()
print('Earth +/-')
print(round(np.nanmean(inside_win_cycle),1))
print(round(np.nanstd(inside_win_cycle),1))
print('min')
print(round(inside_win_cycle[0],1))
print('rise')
print(round(inside_win_cycle[1],1))
print('max')
print(round(inside_win_cycle[2],1))


#only declining phase
print('MAVEN')
print(round(inside_mav_cycle[1],1))




###################### MAVEN

#from processing program
#all in days
totaldays=385 
total_icme_duration_maven=np.sum(icme_durations[imavind])/24

print()
print()
print()


print('MAVEN results from 385 days of data, Dec 2014-Feb 2016, with gaps where no solar wind is available')
print('MAVEN total days of observations with solar wind data:')
print(totaldays)
print('MAVEN total ICME durations:')
print(total_icme_duration_maven)
print('Mars is in percent of time inside ICMEs, for intervals in 2014-2016 (declining phase):')
print(total_icme_duration_maven/totaldays*100)
print('on average, Mars is hit by a CME every ... days')
print(totaldays/np.size(imavind))
print('The ICME average duration is, in hours')
print(np.mean(icme_durations[imavind]))






































sys.exit()





##################### (4) arrival frequencies in ICMECAT  ################################



#################################### (3) time spent inside ICMEs, in % ############################
print()
print('-------------------------------------------------')
print('4 ARRIVAL FREQUENCY ICMECAT, 2 panels')
print()


yearly_bin_edges=[mdates.date2num(sunpy.time.parse_time('2007-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2008-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2009-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2010-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2011-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2012-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2013-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2014-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2015-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2016-01-01')),
                  mdates.date2num(sunpy.time.parse_time('2017-01-01'))]

#bin width in days         
binweite=360/8


sns.set_context("talk")     
sns.set_style("ticks",{'grid.linestyle': '--'})

fig=plt.figure(4,figsize=(12,10	))

fsize=15
ax1 = plt.subplot(211) 

plt.plot_date(icme_start_time_num[iwinind],sc_heliodistance[iwinind],fmt='o',color='mediumseagreen',markersize=5)
plt.plot_date(icme_start_time_num[imesind],sc_heliodistance[imesind],fmt='o',color='darkgrey',markersize=5)
plt.plot_date(icme_start_time_num[ivexind],sc_heliodistance[ivexind],fmt='o',color='orange',markersize=5)
plt.plot_date(icme_start_time_num[istbind],sc_heliodistance[istbind],fmt='o',color='royalblue',markersize=5)
plt.plot_date(icme_start_time_num[istaind],sc_heliodistance[istaind],fmt='o',color='red',markersize=5)
plt.plot_date(icme_start_time_num[iulyind],sc_heliodistance[iulyind],fmt='o',color='brown',markersize=5)
plt.plot_date(icme_start_time_num[imavind],sc_heliodistance[imavind],fmt='o',color='steelblue',markersize=5)




fsize=15
plt.ylabel('Heliocentric distance R [AU]',fontsize=fsize)
plt.xlabel('Year',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


plt.xlim(yearly_bin_edges[0],yearly_bin_edges[10])
ax1.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myformat)



##############

ax2 = plt.subplot(212) 

(histwin, bin_edgeswin) = np.histogram(icme_start_time_num[iwinind], yearly_bin_edges)
(histvex, bin_edgesvex) = np.histogram(icme_start_time_num[ivexind], yearly_bin_edges)
(histmes, bin_edgesmes) = np.histogram(icme_start_time_num[imesind], yearly_bin_edges)
(histstb, bin_edgesstb) = np.histogram(icme_start_time_num[istbind], yearly_bin_edges)
(histsta, bin_edgessta) = np.histogram(icme_start_time_num[istaind], yearly_bin_edges)
(histmav, bin_edgesmav) = np.histogram(icme_start_time_num[imavind], yearly_bin_edges)

#********
#recalculate number of ICMEs as events per month or day, including data gaps


cycle_bin_edges=[minstart, minend, riseend, maxend]

(histwincyc, bin_edgescyc) = np.histogram(icme_start_time_num[iwinind], cycle_bin_edges)
(histvexcyc, bin_edgescyc) = np.histogram(icme_start_time_num[ivexind], cycle_bin_edges)
(histmescyc, bin_edgescyc) = np.histogram(icme_start_time_num[imesind], cycle_bin_edges)
(histstbcyc, bin_edgescyc) = np.histogram(icme_start_time_num[istbind], cycle_bin_edges)
(histstacyc, bin_edgescyc) = np.histogram(icme_start_time_num[istaind], cycle_bin_edges)
(histmavcyc, bin_edgescyc) = np.histogram(icme_start_time_num[imavind], cycle_bin_edges)

#use total_data_days_vex etc. from previous plot 
histwincyc=histwincyc/total_data_days_win_cycle*365
histvexcyc=histvexcyc/total_data_days_vex_cycle*365
histmescyc=histmescyc/total_data_days_mes_cycle*365
histstbcyc=histstbcyc/total_data_days_stb_cycle*365
histstacyc=histstacyc/total_data_days_sta_cycle*365
histmavcyc=histmavcyc/total_data_days_mav_cycle*365


#normalize each dataset for data gaps

histwin=histwin/total_data_days_wind*365
histvex=histvex/total_data_days_vex*365
histmes=histmes/total_data_days_mes*365
histsta=histsta/total_data_days_sta*365
histstb=histstb/total_data_days_stb*365
histmav=histmav/total_data_days_mav*365

binedges=bin_edgeswin
pickle.dump([binedges,histwin,histvex,histmes,histsta,histstb,histmav], open( "plots_stats/icme_frequency.p", "wb" ), protocol=2 )
#[binedges,histwin,histvex,histmes,histsta,histstb,histmav]=pickle.load( open( "plots_stats/stats/icme_frequency.p", "rb" ) )

#binweite=45
ax2.bar(bin_edgeswin[:-1]+30,histwin, width=binweite,color='mediumseagreen', alpha=0.5)
ax2.bar(bin_edgesvex[:-1]+30+binweite,histvex, width=binweite,color='orange', alpha=0.5)
ax2.bar(bin_edgesmes[:-1]+30+ binweite*2,histmes, width=binweite,color='darkgrey', alpha=0.5)
ax2.bar(bin_edgesstb[:-1]+30+binweite*3,histstb, width=binweite,color='royalblue', alpha=0.5)
ax2.bar(bin_edgessta[:-1]+30+binweite*4,histsta, width=binweite,color='red', alpha=0.5)
#ax2.bar(bin_edgessta[:-1]+30+binweite*5,histuly, width=binweite,color='brown', alpha=0.5)
ax2.bar(bin_edgesmav[:-1]+30+binweite*6,histmav, width=binweite,color='steelblue', alpha=0.5)

plt.xlim(yearly_bin_edges[0],yearly_bin_edges[10])
ax2.xaxis_date()
myformat = mdates.DateFormatter('%Y')
ax2.xaxis.set_major_formatter(myformat)
#sets planet / spacecraft labels
xoff=0.85
yoff=0.45
fsize=12
plt.figtext(xoff,yoff,'Earth L1',color='mediumseagreen', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*1,'VEX',color='orange', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*2,'MESSENGER',color='dimgrey', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*3,'STEREO-A',color='red', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*4,'STEREO-B',color='royalblue', fontsize=fsize, ha='left')
#plt.figtext(xoff,yoff-0.03*5,'Ulysses',color='brown', fontsize=fsize, ha='left')
plt.figtext(xoff,yoff-0.03*5,'MAVEN',color='steelblue', fontsize=fsize, ha='left')
#panel labels
plt.figtext(0.02,0.98,'a',color='black', fontsize=fsize, ha='left',fontweight='bold')
plt.figtext(0.02,0.48,'b',color='black', fontsize=fsize, ha='left',fontweight='bold')

plt.ylim(0,48)


#limits solar min/rise/max

vlevel=44
fsize=13

plt.axvspan(minstart,minend, color='green', alpha=0.1)
plt.annotate('solar minimum',xy=(minstart+(minend-minstart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(minstart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(minend-10,vlevel),ha='right', fontsize=fsize)


plt.axvspan(risestart,riseend, color='yellow', alpha=0.1)
plt.annotate('rising phase',xy=(risestart+(riseend-risestart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(risestart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(riseend-10,vlevel),ha='right', fontsize=fsize)

plt.axvspan(maxstart,maxend, color='red', alpha=0.1)
plt.annotate('solar maximum',xy=(maxstart+(maxend-maxstart)/2,vlevel),color='black', ha='center', fontsize=fsize)
plt.annotate('<',xy=(maxstart+10,vlevel),ha='left', fontsize=fsize)
plt.annotate('>',xy=(maxend,vlevel),ha='right', fontsize=fsize)


fsize=15
plt.ylabel('ICMEs per year',fontsize=fsize)
plt.xlabel('Year',fontsize=fsize)
plt.yticks(fontsize=fsize) 
plt.xticks(fontsize=fsize) 


plt.tight_layout()

#sns.despine()
plt.show()
plt.savefig('plots_stats/ICME_frequency_paper.pdf', dpi=300)
plt.savefig('plots_stats/ICME_frequency_paper.png', dpi=300)







print('for solar min 2007-2009 average ICME per year rate:')
mean07=np.mean([histwin[0],histvex[0],histsta[0],histstb[0],histmes[0]])
mean08=np.mean([histwin[1],histvex[1],histsta[1],histstb[1],histmes[1]])
mean09=np.mean([histwin[2],histvex[2],histsta[2],histstb[2],histmes[2]])
print(np.nanmean([mean07,mean08,mean09]))

print('for 2010 2011')
mean10=np.mean([histwin[3],histvex[3],histsta[3],histstb[3],histmes[3]])
mean11=np.mean([histwin[4],histvex[4],histsta[4],histstb[4],histmes[4]])
print(np.mean([mean10,mean11]))


print('for 2012 2013 2014')
mean12=np.mean([histwin[5],histvex[5],histsta[5],histstb[5],histmes[5]])
mean13=np.mean([histwin[6],histvex[6],histsta[6],histstb[6],histmes[6]])
mean14=np.mean([histwin[7],histvex[7],histsta[7],histstb[7],histmes[7]])

print(np.mean([mean12,mean13,mean14]))


'''









































































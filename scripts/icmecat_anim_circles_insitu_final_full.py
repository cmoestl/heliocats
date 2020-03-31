#animation ICMECAT rose scatter plot with HI CME circles and in situ data	

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


#for reading catalogues  
def getcat(filename):
  print('reading CAT '+filename)
  cat=scipy.io.readsav(filename, verbose='true')  
  print('done reading CAT')
  return cat  
  
  
def getpositions(filename):  
  print( 'reading positions in '+filename)
  pos=scipy.io.readsav(filename, verbose='true')  
  print( 'done reading positions')
  return pos

  
def time_to_num_cat(time_in):  

  #for time conversion from catalogue .sav to numerical time
  #this for 1-minute data or lower time resolution

  #for all catalogues
  #time_in is the time in format: 2007-11-17T07:20:00 or 2007-11-17T07:20Z
  #for times help see: 
  #http://docs.sunpy.org/en/latest/guide/time.html
  #http://matplotlib.org/examples/pylab_examples/date_demo2.html
  
  j=0
  #time_str=np.empty(np.size(time_in),dtype='S19')
  time_str= ['' for x in range(len(time_in))]
  #=np.chararray(np.size(time_in),itemsize=19)
  time_num=np.zeros(np.size(time_in))
  
  for i in time_in:

   #convert from bytes (output of scipy.readsav) to string
   time_str[j]=time_in[j][0:16].decode()+':00'
   year=int(time_str[j][0:4])
   time_str[j]
   #convert time to sunpy friendly time and to matplotlibdatetime
   #only for valid times so 9999 in year is not converted
   #pdb.set_trace()
   if year < 2100:
    	  time_num[j]=mdates.date2num(sunpy.time.parse_time(time_str[j]))
   j=j+1  
   #the date format in matplotlib is e.g. 735202.67569444
   #this is time in days since 0001-01-01 UTC, plus 1.
   
   #return time_num which is already an array and convert the list of strings to an array
  return time_num, np.array(time_str)


def decode_array(bytearrin):
 #for decoding the strings from the IDL .sav file to a list of python strings, not bytes 
 #make list of python lists with arbitrary length
 bytearrout= ['' for x in range(len(bytearrin))]
 for i in range(0,len(bytearrin)-1):
  bytearrout[i]=bytearrin[i].decode()
 #has to be np array so to be used with numpy "where"
 bytearrout=np.array(bytearrout)
 return bytearrout  

  
  
def IDL_time_to_num(time_in):  
 #convert IDL time to matplotlib datetime
 time_num=np.zeros(np.size(time_in))
 for ii in np.arange(0,np.size(time_in)):
   time_num[ii]=mdates.date2num(sunpy.time.parse_time(time_in[ii]))
   
 return time_num 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    

######################################################
#main program


plt.close('all')



sns.set_context("talk")     
sns.set_style("darkgrid")  
################## CONTROLS

#how much time is between frames
dayjump=0.25

#either keep or fade detections
fade=1
keep=0

#if keep is selected, the alpha for plotting each dot
keepalpha=0.7

#how long an ARRIVAL stays visible in fade mode
fadedays=30

#how big the circles are on the plot
bscale=4

#half width of the circles
lamda=30



################################


print( 'start icmecat animation program.')

#get ICMECAT
filename_icmecat='ALLCATS/HELCATS_ICMECAT_v10_SCEQ.sav'
i=getcat(filename_icmecat)

#get parameters
bmean=i.icmecat['MO_BMEAN']*bscale #bscale makes circles larger in movie
long=i.icmecat['SC_LONG_HEEQ']*np.pi/180 #heeq longitude converted to radians
rdist=i.icmecat['sc_heliodistance'] #AU
sc=i.icmecat['sc_insitu'] #string
sc=decode_array(sc)


#get indices of events in different spacecraft
vexind=np.where(sc == 'VEX')
staind=np.where(sc == 'STEREO-A')
stbind=np.where(sc == 'STEREO-B')
winind=np.where(sc == 'Wind')
mesind=np.where(sc == 'MESSENGER')
ulyind=np.where(sc == 'ULYSSES')


##################################### read in situ

print( 'read MESSENGER')
#get insitu data
mes= pickle.load( open( "DATACAT/MES_2007to2015_SCEQ_removed.p", "rb" ) )
#time conversion
#mes_time=IDL_time_to_num(mes.time)
print( 'read MESSENGER done.')



print ('read VEX')
#get insitu data
vex= pickle.load( open( "DATACAT/VEX_2007to2014_SCEQ_removed.p", "rb" ) )
#time conversion
#vex_time=IDL_time_to_num(vex.time)
print( 'read VEX done.')



print( 'read Wind')
#get insitu data
wind= pickle.load( open( "DATACAT/WIND_2007to2016_HEEQ.p", "rb" ) )
#time conversion
#wind_time=IDL_time_to_num(wind.time)
print( 'read Wind done.')




print( 'read STEREO-A')
#get insitu data
sta= pickle.load( open( "DATACAT/STA_2007to2015_SCEQ.p", "rb" ) )
#time conversion
#sta_time=IDL_time_to_num(sta.time)
print( 'read STA done.')




print( 'read STEREO-B')
#get insitu data
stb= pickle.load( open( "DATACAT/STB_2007to2014_SCEQ.p", "rb" ) )

#time conversion
#stb_time=IDL_time_to_num(stb.time)
print( 'read STB done.')

#save times
#pickle.dump([vex_time,wind_time,sta_time,stb_time,mes_time], open( "DATACAT/Insitu_times_mdates_2.p", "wb" ) )

#quicker when just reloading times
[vex_time,wind_time,sta_time,stb_time,mes_time]=pickle.load( open( "DATACAT/Insitu_times_mdates_2.p", "rb" ) )
#print 'loaded in situ times'
######################################




#get positions
pos=getpositions('DATACAT/positions_2007_2018_HEEQ_6hours.sav')
[pos_time_num,pos_time_str]=time_to_num_cat(pos.time)
#available as pos.mercury etc.

#get cme apex positions
h=getcat('ALLCATS/hicat_v3_cat_behind_visual.sav')
[h_time_num,h_time_str]=time_to_num_cat(h.all_apex_t_str)

all_apex_s=decode_array(h.all_apex_s)


#make time conversion for all icme_start_time variables
#save it as string
icme_start_time_str=i.icmecat['icme_start_time']
#save it as matplotlib date number
[icme_start_time_num,icme_start_time_str]=time_to_num_cat(icme_start_time_str)


#for each spacecraft, make a zeros array 
active_icme_vex=np.zeros(np.size(icme_start_time_num))
active_icme_stb=np.zeros(np.size(icme_start_time_num))
active_icme_sta=np.zeros(np.size(icme_start_time_num))
active_icme_win=np.zeros(np.size(icme_start_time_num))
active_icme_mes=np.zeros(np.size(icme_start_time_num))
active_icme_uly=np.zeros(np.size(icme_start_time_num))


#initiate plot
plt.figure(1, figsize=(12, 6), dpi=100, facecolor='w', edgecolor='w')


#full  movie April 2014 Jan 1 until end of November 2014
frame_time_num=mdates.date2num(sunpy.time.parse_time('2007-Apr-1'))



################################### plot over all frames
for k in np.arange(12680/4,(12680+120)/4,dayjump):  
#3169 is time in days
  
  
   
 start=time.time()
 #to current frame time, the days need to be added, so +k is done 
 #save frame time as string to write on plot
 frame_time_str=str(mdates.num2date(frame_time_num+k))
 print( 'current frame_time_num+k', frame_time_str)
 
 
 #for each frame time, check active ICMEs by looking into the full catalogue:
 
 for m in range(0,len(icme_start_time_num)):
 
  #calculate difference in icme_start_time to current frame
  icme_diff_to_frame=(frame_time_num+k)-icme_start_time_num[m]
 
  #for all icme_start_times that are later than the current frame, 
  #make them active for 30 days (fading) or infinite (keeping).
  
  
  #**********************for fading
  if fade > 0:
   if  icme_diff_to_frame > 0 and icme_diff_to_frame < fadedays:
     #check if this active icme belongs to a spacecraft
     #in1d compares to arrays; true or 1 if m is contained in vexind
      if np.in1d(m,vexind) == 1:
         active_icme_vex[m]=icme_diff_to_frame
     #same for the other spacecraft    
      if np.in1d(m,stbind) == 1:
         active_icme_stb[m]=icme_diff_to_frame
      if np.in1d(m,staind) == 1:
         active_icme_sta[m]=icme_diff_to_frame
      if np.in1d(m,winind) == 1:
         active_icme_win[m]=icme_diff_to_frame
      if np.in1d(m,mesind) == 1:
         active_icme_mes[m]=icme_diff_to_frame
      if np.in1d(m,ulyind) == 1:
         active_icme_uly[m]=icme_diff_to_frame

   else:
     #if no detection, set the index to 0
     active_icme_vex[m]=0
     active_icme_stb[m]=0
     active_icme_sta[m]=0
     active_icme_win[m]=0
     active_icme_mes[m]=0
     active_icme_uly[m]=0
     
#************************** for keeping
  if keep > 0:
   if  icme_diff_to_frame > 0:
          #check if this active icme belongs to a spacecraft
       #in1d compares to arrays; true or 1 if m is contained in vexind
       if np.in1d(m,vexind) == 1:
           active_icme_vex[m]=icme_diff_to_frame
       #same for the other spacecraft    
       if np.in1d(m,stbind) == 1:
           active_icme_stb[m]=icme_diff_to_frame
       if np.in1d(m,staind) == 1:
           active_icme_sta[m]=icme_diff_to_frame
       if np.in1d(m,winind) == 1:
           active_icme_win[m]=icme_diff_to_frame
       if np.in1d(m,mesind) == 1:
           active_icme_mes[m]=icme_diff_to_frame
   else:
      #if no detection, set the index to 0
      active_icme_vex[m]=0
      active_icme_stb[m]=0
      active_icme_sta[m]=0
      active_icme_win[m]=0
      active_icme_mes[m]=0
 
   
 
 #look which ICMEs are active
 active_index_vex=np.where(active_icme_vex > 0)
 active_index_stb=np.where(active_icme_stb > 0)
 active_index_sta=np.where(active_icme_sta > 0)
 active_index_win=np.where(active_icme_win > 0)
 active_index_mes=np.where(active_icme_mes > 0)
 active_index_uly=np.where(active_icme_uly > 0)


 #print 'active icme indices are:', active_index_vex
 
 print (' ')
 
  
 #check for active CME indices from HICAT (with the lists produced in IDL for the apex positions)
 #check where time is identical to frame time
 cmeind=np.where(h_time_num == frame_time_num+k)
 
 ############make plot 
 
 #                       rows - columns, starts with 0
 ax = plt.subplot2grid((5,2), (0, 0), rowspan=5, projection='polar')
 #ax = plt.subplot(121,projection='polar')
 


 ######################## 1 plot all active CME circles
 #ax.scatter(h.all_apex_long[cmeind]*np.pi/180,h.all_apex_r[cmeind], s=10, c='black', alpha=1, marker='s')
 
 #plot all active CME circles
 #if np.size(cmeind) >0:
 for p in range(0,np.size(cmeind)):
   #print p, h.all_apex_long[cmeind[0][p]], h.all_apex_r[cmeind[0][p]]
   #central d
   dir=np.array([np.cos(h.all_apex_long[cmeind[0][p]]*np.pi/180),np.sin(h.all_apex_long[cmeind[0][p]]*np.pi/180)])*h.all_apex_r[cmeind[0][p]]
   
   #points on circle, correct for longitude
   circ_ang = ((np.arange(111)*2-20)*np.pi/180)-(h.all_apex_long[cmeind[0][p]]*np.pi/180)
   
 
   #these equations are from moestl and davies 2013
   xc = 0+dir[0]/(1+np.sin(lamda*np.pi/180)) + (h.all_apex_r[cmeind[0][p]]*np.sin(lamda*np.pi/180)/(1+np.sin(lamda*np.pi/180)))*np.sin(circ_ang)
   yc = 0+dir[1]/(1+np.sin(lamda*np.pi/180)) + (h.all_apex_r[cmeind[0][p]]*np.sin(lamda*np.pi/180)/(1+np.sin(lamda*np.pi/180)))*np.cos(circ_ang)
   #now convert to polar coordinates
   rcirc=np.sqrt(xc**2+yc**2)
   longcirc=np.arctan2(yc,xc)
   #plot in correct color
   if all_apex_s[cmeind[0][p]] == 'A':    
    #make alpha dependent on distance to solar equatorial plane - maximum latitude is -40/+40 - 
    #so to make also the -/+40 latitude CME visible, divide by 50 so alpha > 0 for these events
    ax.plot(longcirc,rcirc, c='red', alpha=1-abs(h.all_apex_lat[cmeind[0][p]]/50), lw=1.5) 
   if all_apex_s[cmeind[0][p]] == 'B':
    ax.plot(longcirc,rcirc, c='royalblue', alpha=1-abs(h.all_apex_lat[cmeind[0][p]]/50), lw=1.5) 
  



 ####################### 3 plot ICME detections 
 #fader style plot alpha dependent on time difference - for this loop over each element:
 if fade >0:
 
  for y in range(0,np.size(active_index_vex)):
   z=active_index_vex[0][y] #access elements in tuple that is produced by where
   fadealpha=1-active_icme_vex[z]/(fadedays)  #fadedays is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='orange', alpha=fadealpha)

  for y in range(0,np.size(active_index_sta)):
   z=active_index_sta[0][y]
   fadealpha=1-active_icme_sta[z]/(fadedays)  #30 days is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='red', alpha=fadealpha)

  for y in range(0,np.size(active_index_stb)):
   z=active_index_stb[0][y]
   fadealpha=1-active_icme_stb[z]/(fadedays)  #30 days is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='royalblue', alpha=fadealpha)

  for y in range(0,np.size(active_index_win)):
   z=active_index_win[0][y]
   fadealpha=1-active_icme_win[z]/(fadedays)  #30 days is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='mediumseagreen', alpha=fadealpha)

  for y in range(0,np.size(active_index_mes)): 
   z=active_index_mes[0][y]
   fadealpha=1-active_icme_mes[z]/(fadedays)  #30 days is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='dimgrey', alpha=fadealpha)

  for y in range(0,np.size(active_index_uly)): 
   z=active_index_uly[0][y]
   fadealpha=1-active_icme_uly[z]/(fadedays)  #30 days is maximum difference in time, and alpha from 0 to 1
   ax.scatter(long[z], rdist[z], s=bmean[z], c='darkolivegreen', alpha=fadealpha)


 if keep >0: 
  ax.scatter(long[active_index_vex], rdist[active_index_vex], s=bmean[active_index_vex], c='orange', alpha=keepalpha)
  ax.scatter(long[active_index_sta], rdist[active_index_sta], s=bmean[active_index_sta], c='red', alpha=keepalpha)
  ax.scatter(long[active_index_stb], rdist[active_index_stb], s=bmean[active_index_stb], c='royalblue', alpha=keepalpha)
  ax.scatter(long[active_index_win], rdist[active_index_win], s=bmean[active_index_win], c='mediumseagreen', alpha=keepalpha)
  ax.scatter(long[active_index_mes], rdist[active_index_mes], s=bmean[active_index_mes], c='dimgrey', alpha=keepalpha)
 

 


 
 plt.suptitle('STEREO/HI modeled CMEs (SSEF30) + in situ ICME detections and data    HELCATS - HIGEOCAT ICMECAT DATACAT', fontsize=12)	
  
 #Sun
 ax.scatter(0,0,s=100,c='yellow',alpha=0.8, edgecolors='yellow')
 plt.figtext(0.30,0.5,'Sun', fontsize=10, ha='center')
 
 #Earth
 plt.figtext(0.30,0.25,'Earth', fontsize=10, ha='center')
 
 #units
 #plt.figtext(0.525,0.0735,'HEEQ longitude', fontsize=10, ha='left')
 #plt.figtext(0.655,0.164,'AU', fontsize=10, ha='center')

 

 #----------------- legend
 plt.figtext(0.05,0.02,'Mercury', color='dimgrey', ha='center', fontsize=12)
 plt.figtext(0.15,0.02,'MESSENGER', color='dimgrey', ha='center', fontsize=10)
 plt.figtext(0.25	,0.02,'Venus', color='orange', ha='center',fontsize=12)
 plt.figtext(0.35,0.02,'STEREO-A', color='red', ha='center',fontsize=12)
 plt.figtext(0.48,0.02,'STEREO-B', color='royalblue', ha='center',fontsize=12)
 plt.figtext(0.58,0.02,'Earth', color='mediumseagreen', ha='center',fontsize=12)
 plt.figtext(0.65,0.02,'Mars', color='orangered', ha='center',fontsize=10)
 plt.figtext(0.71,0.02,'MSL', color='magenta', ha='center', fontsize=10)
 plt.figtext(0.76,0.02,'Maven', color='steelblue', ha='center', fontsize=10)
 plt.figtext(0.83,0.02,'Ulysses', color='darkolivegreen', ha='center', fontsize=10)
 plt.figtext(0.90,0.02,'Rosetta', color='black', ha='center', fontsize=10)
 
 
 #add legend for bmean
 bleg=np.array([10,50,100])*bscale
 blegstr=['10 nT','50','100']

 blegr=np.zeros(len(bleg))+1.6
 blegt=np.radians(range(170,195,10))
 ax.scatter(blegt, blegr,s=bleg,c='violet', edgecolor='violet')

 for p in range(0,len(bleg)):
   ax.annotate(blegstr[p],xy=(blegt[p],blegr[p]-0.2), ha='center', va='center', fontsize=8)
  
  
 
 
 ############################## plot positions
 #check which index is closest in positions to current time
 #frame_time_num+k vs. pos_time_num
 
 timeind=np.where(frame_time_num+k-pos_time_num == min(abs((frame_time_num+k)-pos_time_num)))
 
 #index 1 is longitude, 0 is rdist
 ax.scatter(pos.venus[1,timeind], pos.venus[0,timeind], s=50, c='orange', alpha=1, lw=0)
 ax.scatter(pos.mercury[1,timeind], pos.mercury[0,timeind], s=50, c='dimgrey', alpha=1,lw=0)
 ax.scatter(pos.messenger[1,timeind], pos.messenger[0,timeind], s=25, c='dimgrey', alpha=1,lw=0,marker='s')
 ax.scatter(pos.sta[1,timeind], pos.sta[0,timeind], s=25, c='red', alpha=1,lw=0, marker='s')
 ax.scatter(pos.stb[1,timeind], pos.stb[0,timeind], s=25, c='royalblue', alpha=1,lw=0, marker='s')
 ax.scatter(pos.earth[1,timeind], pos.earth[0,timeind], s=50, c='mediumseagreen', alpha=1,lw=0)
 ax.scatter(pos.mars[1,timeind], pos.mars[0,timeind], s=50, c='orangered', alpha=1,lw=0)
 ax.scatter(pos.ulysses[1,timeind], pos.ulysses[0,timeind], s=25, c='darkolivegreen', alpha=1,lw=0,marker='s')
 ax.scatter(pos.msl[1,timeind], pos.msl[0,timeind], s=25, c='magenta', alpha=1,lw=0,marker='s')
 ax.scatter(pos.maven[1,timeind], pos.maven[0,timeind], s=25, c='steelblue', alpha=1,lw=0, marker='s')
 ax.scatter(pos.rosetta[1,timeind], pos.rosetta[0,timeind], s=25, c='black', alpha=1,lw=0, marker='s')
 
  

 
 #set axes
 plt.thetagrids(range(0,360,45),(u'0\u00b0 HEEQ longitude',u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'-135\u00b0',u'-90\u00b0',u'-45\u00b0'), fmt='%d', frac = 1.05,fontsize=10)
 ax.set_theta_zero_location('S')
 ax.set_ylim(0, 1.8)
 plt.rgrids((0.4,0.7,1.0,1.3,1.6),('0.4','0.7','1.0','1.3','1.6 AU'),fontsize=10)
 
 
 #plot text for date extra so it does not move 
 #year
 plt.figtext(0.47-0.22,0.9,frame_time_str[0:4], fontsize=13, ha='center')
 #month
 plt.figtext(0.51-0.22,0.9,frame_time_str[5:7], fontsize=13, ha='center')
 #day
 plt.figtext(0.54-0.22,0.9,frame_time_str[8:10], fontsize=13, ha='center')
 #hours
 plt.figtext(0.57-0.22,0.9,frame_time_str[11:13], fontsize=13, ha='center')

 #mysignature
 plt.figtext(0.96,0.01,r'$C. M\ddot{o}stl$', fontsize=7, ha='center')
 
 


 ############# 5 in situ data plots 

 plotstartdate=mdates.num2date(frame_time_num+k-3)
 plotenddate=mdates.num2date(frame_time_num+k+3)
 
 
 
 #slicing
 
 #take only those indices where the difference to frame_time_num+k is less than 3
 mes_ind_plot=np.where(abs(mes_time-(frame_time_num+k)) < 3)
 vex_ind_plot=np.where(abs(vex_time-(frame_time_num+k)) < 3)
 stb_ind_plot=np.where(abs(stb_time-(frame_time_num+k)) < 3)
 sta_ind_plot=np.where(abs(sta_time-(frame_time_num+k)) < 3)
 wind_ind_plot=np.where(abs(wind_time-(frame_time_num+k)) < 3)
  
 
 #rows - columns
 #MESSENGER
 ax2 = plt.subplot2grid((5,2), (0, 1))
 ax2.plot_date(mes_time[mes_ind_plot],mes.btot[mes_ind_plot],'-k', lw=0.3)
 ax2.plot_date(mes_time[mes_ind_plot],mes.bx[mes_ind_plot], '-r',lw=0.3)
 ax2.plot_date(mes_time[mes_ind_plot],mes.by[mes_ind_plot],'-g',lw=0.3)
 ax2.plot_date(mes_time[mes_ind_plot],mes.bz[mes_ind_plot],'-b',lw=0.3)
 plt.tick_params( axis='x', labelbottom='off')
 #current time
 plt.yticks(fontsize=9)
 plt.ylabel('B SCEQ [nT]', fontsize=9)
 ax2.plot_date([mdates.num2date(frame_time_num+k),mdates.num2date(frame_time_num+k)], [-120,120],'-k', lw=0.5, alpha=0.8)
 plt.xlim((plotstartdate, plotenddate))
 plt.ylim((-120, 120))


 
 #VEX
 ax3 = plt.subplot2grid((5,2), (1, 1))
 ax3.plot_date(vex_time[vex_ind_plot],vex.btot[vex_ind_plot],'-k', lw=0.3)
 ax3.plot_date(vex_time[vex_ind_plot],vex.bx[vex_ind_plot], '-r',lw=0.3)
 ax3.plot_date(vex_time[vex_ind_plot],vex.by[vex_ind_plot],'-g',lw=0.3)
 ax3.plot_date(vex_time[vex_ind_plot],vex.bz[vex_ind_plot],'-b',lw=0.3)
 plt.tick_params( axis='x', labelbottom='off')
 plt.yticks(fontsize=9)
 plt.ylabel('B SCEQ [nT]', fontsize=9)
 ax3.plot_date([mdates.num2date(frame_time_num+k),mdates.num2date(frame_time_num+k)], [-50,50],'-k', lw=0.5, alpha=0.8)
 plt.xlim((plotstartdate, plotenddate))
 plt.ylim((-50, 50))
 

 
 #Earth
 ax4 = plt.subplot2grid((5,2), (2, 1))
 ax4.plot_date(wind_time[wind_ind_plot],wind.btot[wind_ind_plot],'-k', lw=0.3)
 ax4.plot_date(wind_time[wind_ind_plot],wind.bx[wind_ind_plot], '-r',lw=0.3)
 ax4.plot_date(wind_time[wind_ind_plot],wind.by[wind_ind_plot],'-g',lw=0.3)
 ax4.plot_date(wind_time[wind_ind_plot],wind.bz[wind_ind_plot],'-b',lw=0.3)
 plt.tick_params( axis='x', labelbottom='off')
 plt.yticks(fontsize=9)
 plt.ylabel('B SCEQ [nT]', fontsize=9)
 ax4.plot_date([mdates.num2date(frame_time_num+k),mdates.num2date(frame_time_num+k)], [-50,50],'-k', lw=0.5, alpha=0.8)
 plt.xlim((plotstartdate, plotenddate))
 plt.ylim((-35, 35))
 
 #STA
 ax5 = plt.subplot2grid((5,2), (3, 1))
 ax5.plot_date(sta_time[sta_ind_plot],sta.btot[sta_ind_plot],'-k', lw=0.3)
 ax5.plot_date(sta_time[sta_ind_plot],sta.bx[sta_ind_plot], '-r',lw=0.3)
 ax5.plot_date(sta_time[sta_ind_plot],sta.by[sta_ind_plot],'-g',lw=0.3)
 ax5.plot_date(sta_time[sta_ind_plot],sta.bz[sta_ind_plot],'-b',lw=0.3)
 plt.tick_params( axis='x', labelbottom='off')
 plt.yticks(fontsize=9)
 plt.ylabel('B SCEQ [nT]', fontsize=9)
 ax5.plot_date([mdates.num2date(frame_time_num+k),mdates.num2date(frame_time_num+k)], [-50,50],'-k', lw=0.5, alpha=0.8)
 plt.xlim((plotstartdate, plotenddate))
 plt.ylim((-35, 35))


 #STB
 ax6 = plt.subplot2grid((5,2), (4, 1))
 ax6.plot_date(stb_time[stb_ind_plot],stb.btot[stb_ind_plot],'-k', lw=0.3)
 ax6.plot_date(stb_time[stb_ind_plot],stb.bx[stb_ind_plot], '-r',lw=0.3)
 ax6.plot_date(stb_time[stb_ind_plot],stb.by[stb_ind_plot],'-g',lw=0.3)
 ax6.plot_date(stb_time[stb_ind_plot],stb.bz[stb_ind_plot],'-b',lw=0.3)  
 plt.xlim((plotstartdate, plotenddate))
   
 myformat = mdates.DateFormatter('%m-%d')
 ax6.xaxis.set_major_formatter(myformat)  
 plt.yticks(fontsize=9)
 plt.ylabel('B SCEQ [nT]', fontsize=9)
 ax6.plot_date([mdates.num2date(frame_time_num+k),mdates.num2date(frame_time_num+k)], [-50,50],'-k', lw=0.5, alpha=0.8)
 plt.ylim((-35, 35))
 plt.xticks(fontsize=10)
 
 
 #labeling of spacecraft and longitude in HEEQ
 plt.figtext(0.92,0.82,'MESSENGER',color='dimgrey', fontsize=10, ha='left')
 plt.figtext(0.94,0.77,"%d" % (pos.messenger[1,timeind]*180/np.pi),color='black', fontsize=10, ha='right')
 
 plt.figtext(0.92,0.82-0.165,'VEX',color='orange', fontsize=10, ha='left')
 plt.figtext(0.94,0.77-0.165,"%d" % (pos.venus[1,timeind]*180/np.pi),color='black', fontsize=10, ha='right')

 plt.figtext(0.92,0.82-0.165*2,'Wind',color='mediumseagreen', fontsize=10, ha='left')
 plt.figtext(0.94,0.77-0.165*2,"%d" % (pos.earth[1,timeind]*180/np.pi),color='black', fontsize=10, ha='right')

 plt.figtext(0.92,0.82-0.165*3,'STEREO-A',color='red', fontsize=10, ha='left')
 plt.figtext(0.94,0.77-0.165*3,"%d" % (pos.sta[1,timeind]*180/np.pi),color='black', fontsize=10, ha='right')

 plt.figtext(0.92,0.82-0.165*4,'STEREO-B',color='royalblue', fontsize=10, ha='left')
 plt.figtext(0.94,0.77-0.165*4,"%d" % (pos.stb[1,timeind]*180/np.pi),color='black', fontsize=10, ha='right')

 #labeling in situ components
 plt.figtext(0.75,0.92,'Bx',color='red', fontsize=10, ha='left')
 plt.figtext(0.8,0.92,'By',color='green', fontsize=10, ha='left')
 plt.figtext(0.85,0.92,'Bz',color='blue', fontsize=10, ha='left')


 
 
 
 
 #save figure for frame - this starts with zero at the start time
 framestr = '%04i' % (k*4)  
 #framenr=framenr+1
 print( 'frame nr.', framestr) 
 #plt.show()
 if fade >0:
  plt.savefig('animations/animation_icmecat_6hour_fade_circ_insitu_final_full/icmecat_'+framestr+'.png',  dpi=300)
  #plt.savefig('animations/animation_icmecat_6hour_fade_circ_insitu_final_full/icmecat_'+framestr+'.jpg',  dpi=300)

# if keep >0: 
#  plt.savefig('animations/animation_icmecat_6hour_keep_circ_insitu_final_full/icmecat_'+framestr+'.jpg', format='jpg', dpi=300)


 end=time.time()
 print( 'took time in seconds:', (end-start) ,'for this frame')



 #clears plot window
 plt.clf()


############end of cycle



#make animation convert with automator into jpg before

#os.system('/Users/chris/movie/ffmpeg -r 15 -i /Users/chris/python/catpy/animations/animation_icmecat_6hour_fade_circ_insitu_final_full_jpg/icmecat_%04d.jpg -b 5000k -r 15 animations/icmecat_anim_6hour_fade_circ_insitu_final_full.mp4 -y')

print( 'made movie')




print( 'end icmecat animation program.')


#/Users/chris/movie/ffmpeg -r 15 -i /Users/chris/python/catpy/animations/animation_icmecat_6hour_fade_circ_insitu_all/icmecat_%04d.jpg -b 5000k -r 15 animations/icmecat_anim_6hour_fade_circ_insitu_all.mp4 -y









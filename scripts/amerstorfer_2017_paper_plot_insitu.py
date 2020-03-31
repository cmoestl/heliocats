

from scipy import stats
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import matplotlib as mpl
import numpy as np
import sunpy.time
import time
import pickle
import seaborn as sns

#initialize
#get current directory
#os.system('pwd')
#closes all plots
plt.close('all')



def IDL_time_to_num(time_in):  
 #convert IDL time to matplotlib datetime
 time_num=np.zeros(np.size(time_in))
 for ii in np.arange(0,np.size(time_in)):
   time_num[ii]=mdates.date2num(sunpy.time.parse_time(time_in[ii]))
   
 return time_num 
  
  

sns.set_context("talk")     
#sns.set_style("whitegrid")  
sns.set_style("ticks")  




#plot +/- 5 days
dur=4
frame_time_num=mdates.date2num(sunpy.time.parse_time('2010-Nov-8'))
 
read_data=0

if read_data >0:

 print('read STEREO-B')
 #get insitu data
 stb= pickle.load( open( "DATACAT/STB_2007to2014_SCEQ.p", "rb" ) )
 #time conversion
 stb_time=IDL_time_to_num(stb.time)
 print('read STB done.')


 print('read MESSENGER')	
 #get insitu data
 mes= pickle.load( open( "DATACAT/MES_2007to2014_SCEQ.p", "rb" ) )
 #time conversion
 mes_time=IDL_time_to_num(mes.time)
 print('read MESSENGER done.')

 #take only those indices of the plot time
 mes_ind_plot=np.where(abs(mes_time-(frame_time_num)) < dur)
 stb_ind_plot=np.where(abs(stb_time-(frame_time_num)) < dur)



 #
 pickle.dump((mes[mes_ind_plot],stb[stb_ind_plot], stb_time[stb_ind_plot], mes_time[mes_ind_plot]), open( "DATACAT/mes_stb_tanja_plot_2017paper.p", "wb" ) )
 


[mes,stb,stb_time,mes_time]=pickle.load( open( "DATACAT/mes_stb_tanja_plot_2017paper.p", "rb" ) )





############# 5 in situ data plots 

fig=plt.figure(2,figsize=(10,12))



plotstartdate=mdates.num2date(frame_time_num-dur)
plotenddate=mdates.num2date(frame_time_num+dur)
  
 
  

#linewidth  
weite=1  
 
#fontsize
fsize=15
 
 
#rows - columns
#MESSENGER
ax1 = plt.subplot2grid((5,1), (0, 0))
ax1.plot_date([mdates.num2date(frame_time_num-dur),mdates.num2date(frame_time_num+dur)], [0,0],'--k', lw=0.5, alpha=0.8)
ax1.plot_date(mes_time,mes.btot,'-k',lw=weite, label='$\mathregular{|B|}$')
ax1.plot_date(mes_time,mes.bx, '-',color='magenta',lw=weite,label='$\mathregular{B_x}$')
ax1.plot_date(mes_time,mes.by,'-g',lw=weite,label='$\mathregular{B_y}$')
ax1.plot_date(mes_time,mes.bz,'-b',lw=weite,label='$\mathregular{B_z}$' )
plt.tick_params( axis='x', labelbottom='off')
plt.ylabel('B SCEQ [nT]', fontsize=fsize)
plt.xlim((plotstartdate, plotenddate))

plt.legend(loc=0,ncol=4,fontsize=fsize)

arr_mes=mdates.date2num(sunpy.time.parse_time('2010-Nov-5 11:46'))
#arrival time MESSENGER
ax1.plot_date([mdates.num2date(arr_mes),mdates.num2date(arr_mes)], [-60,60],'-k', lw=1, alpha=0.8)




#STB
ax2 = plt.subplot2grid((5,1), (1, 0))
ax2.plot_date([mdates.num2date(frame_time_num-dur),mdates.num2date(frame_time_num+dur)], [0,0],'--k', lw=0.5, alpha=0.8)
ax2.plot_date(stb_time,stb.btot,'-k', lw=weite, label='$\mathregular{|B|}$')
ax2.plot_date(stb_time,stb.bx, '-',color='magenta',lw=weite,label='$\mathregular{B_x}$')
ax2.plot_date(stb_time,stb.by,'-g',lw=weite,label='$\mathregular{B_y}$')
ax2.plot_date(stb_time,stb.bz,'-b',lw=weite,label='$\mathregular{B_z}$')  
ax2.plot_date([mdates.num2date(frame_time_num-dur),mdates.num2date(frame_time_num+dur)], [0,0],'-k', lw=0.5, alpha=0.8)
plt.xlim((plotstartdate, plotenddate))
plt.ylabel('B SCEQ [nT]', fontsize=fsize)
plt.tick_params( axis='x', labelbottom='off')
plt.ylim((-20,20))


arr_stb=mdates.date2num(sunpy.time.parse_time('2010-Nov-7 19:05'))
#arrival time STB
ax2.plot_date([mdates.num2date(arr_stb),mdates.num2date(arr_stb)], [-25,25],'-k', lw=1.5, alpha=0.8)
   

   
#V
ax3 = plt.subplot2grid((5,1), (2, 0))
ax3.plot_date(stb_time,stb.vtot,'-k', lw=weite)
plt.xlim((plotstartdate, plotenddate))
plt.ylabel('V $\mathregular{[km \\ s^{-1}]}$', fontsize=fsize)
plt.tick_params( axis='x', labelbottom='off')
ax3.plot_date([mdates.num2date(arr_stb),mdates.num2date(arr_stb)], [250,700],'-k', lw=1.5, alpha=0.8)
   


#T
ax4 = plt.subplot2grid((5,1), (3, 0))
ax4.plot_date(stb_time,stb.temperature*1e-6,'-k', lw=weite)
plt.xlim((plotstartdate, plotenddate))
plt.ylabel('T [MK]', fontsize=fsize)
plt.tick_params( axis='x', labelbottom='off')
ax4.plot_date([mdates.num2date(arr_stb),mdates.num2date(arr_stb)], [0,1.5],'-k', lw=1.5, alpha=0.8)
plt.ylim((0,1.5))   

#N
ax5 = plt.subplot2grid((5,1), (4, 0))
ax5.plot_date(stb_time,stb.density,'-k', lw=weite)
plt.ylabel('N $\mathregular{[ccm^{-3}]}$', fontsize=fsize)
plt.xlim((plotstartdate, plotenddate))
ax5.plot_date([mdates.num2date(arr_stb),mdates.num2date(arr_stb)], [0,35],'-k', lw=1.5, alpha=0.8)
   
      
ax5.xaxis.set_major_formatter(mdates.DateFormatter('%b %d\n%Y')) 


plt.tight_layout()

plt.figtext(0.12,0.95,'MESSENGER',color='black', fontsize=15, ha='left')
plt.figtext(0.12,0.76,'STEREO-B',color='blue', fontsize=15, ha='left')


plt.savefig('plots/amerstorfer_2017_paper/insitu_MES_STB_Nov_2010_new.eps', format='eps',dpi=300)
plt.savefig('plots/amerstorfer_2017_paper/insitu_MES_STB_Nov_2010_new.png', format='png', dpi=300)



 

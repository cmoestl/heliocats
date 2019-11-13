#data.py
#load data for heliocats
#https://github.com/cmoestl/heliocats

import numpy as np
import scipy
import copy
import sunpy
import matplotlib.dates as mdates
import datetime
import urllib
import json
import os
import pdb
import scipy.io
import pickle
import sys
import cdflib


######################################################### MAVEN ####################################
'''
Hi Christian,

I just did a run on the MAVEN datasets that you had with my new technique (and using the Gruesbeck model). The mat file is here:

https://oeawcloud.oeaw.ac.at/index.php/s/Hjc8BVmtQT3k0ub

Password is:      maven2019

The save command in matlab I have used is the following:
save('Data-MAVEN-SolarWind.mat', 'timeD','np','Vx','Vy','Vz','VT','Tp','Bx','By','Bz','BT','Xsc','Ysc','Zsc');

Note that Xsc, Ysc, and Zsc are given in units of Mars radius (Rp = 3389.5 km). All other ones are the original units of the cdf file you gave me. I have also gone through the data and set to NaN all negative values of density np, temperature Tp and total magnetic field BT.

It would be interesting to compare your old results with these ones, also as a double check that the data was filtered correctly, although the 3D model of Gruesbeck+ 2018 does not assume an aberration angle of 4 degrees like the polar model of Edberg+ 2008 -- which I am not sure you took into account originally when processing the data.

Cheers,

Cyril

P.S. I'm writing a compendium for all the details of the technique I used. Hopefully this should be finished soon.

Hallo Christian,

Super, happy that it looks fine! :) Yes, agreed with the median filtering, this should take care of the spikes. 
There were also many strange spikes in Xsc, Ysc and Zsc (very large values >2.9e27, probably due to an issue with the SPICE kernel, 
about 1685 points, i.e., 0.087% of the data), so I set all of these anomalous data points to NaN too (all variables including Xsc, Ysc, Zsc).

Cheers,

Cyril


Mehr anzeigen von Christian Möstl



'''
def convert_MAVEN_mat_to_pickle():

    print('load MAVEN from MAT')
    file='data/MAVEN_2014to2018_removed_cyril.mat'
    mav = scipy.io.loadmat(file)
    print('save MAVEN as pickle')
    pickle.dump(mav, open("data/MAVEN_2014to2018_removed_cyril.p", "wb"))

def load_MAVEN():
    
    print('load MAVEN from pickle')
    file='data/MAVEN_2014to2018_removed_cyril.p'
    mav=pickle.load( open( file, 'rb' ) )
    return mav
    









'''
# use
# import importlib
# importlib.reload(cme_stats_module)
# to update module while working in command line 

#set breakpoints with pdb.set_trace



# LIST OF FUNCTIONS

# load_url_current_directory
# getpositions
# getcat
# gaussian
# dynamic_pressure
# decode_array
# time_to_num_cat
# get_omni2_data


def load_url_current_directory(filename,url):
#loads a file from any url to the current directory
#I use owncloud for the direct url links, 
#also works for dropbox when changing the last 0 to 1 in the url-> gives a direct link to files

 if not os.path.exists(filename):
  print('download file ', filename, ' from')
  print(url)
  try: 
    urllib.request.urlretrieve(url, filename)
    print('done')
  except urllib.error.URLError as e:
    print(' ', data_url,' ',e.reason)




def getpositions(filename):  
    pos=scipy.io.readsav(filename)  
    print
    print('positions file:', filename) 
    return pos


def getcat(filename):
   cat=scipy.io.readsav(filename)#, verbose='true')  
   return cat  
  
  
def gaussian(x, amp, mu, sig):
   return amp * exp(-(x-cen)**2 /wid)



def dynamic_pressure(density, speed):
   # make dynamic pressure from density and speed
   #assume pdyn is only due to protons
   #pdyn=np.zeros(len([density])) #in nano Pascals
   protonmass=1.6726219*1e-27  #kg
   pdyn=np.multiply(np.square(speed*1e3),density)*1e6*protonmass*1e9  #in nanoPascal
   return pdyn
  
def decode_array(bytearrin):
  #for decoding the strings from the IDL .sav file to a list of python strings, not bytes 
  #make list of python lists with arbitrary length
  bytearrout= ['' for x in range(len(bytearrin))]
  for i in range(0,len(bytearrin)-1):
    bytearrout[i]=bytearrin[i].decode()
  #has to be np array so to be used with numpy "where"
  bytearrout=np.array(bytearrout)
  return bytearrout  
  
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



def get_omni2_data():

  #FORMAT(2I4,I3,I5,2I3,2I4,14F6.1,F9.0,F6.1,F6.0,2F6.1,F6.3,F6.2, F9.0,F6.1,F6.0,2F6.1,F6.3,2F7.2,F6.1,I3,I4,I6,I5,F10.2,5F9.2,I3,I4,2F6.1,2I6,F5.1)
 #1963   1  0 1771 99 99 999 999 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 9999999. 999.9 9999. 999.9 999.9 9.999 99.99 9999999. 999.9 9999. 999.9 999.9 9.999 999.99 999.99 999.9  7  23    -6  119 999999.99 99999.99 99999.99 99999.99 99999.99 99999.99  0   3 999.9 999.9 99999 99999 99.9

#define variables from OMNI2 dataset
 #see http://omniweb.gsfc.nasa.gov/html/ow_data.html

 #omni2_url='ftp://nssdcftp.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
 
 #check how many rows exist in this file
 f=open('omni2_all_years.dat')
 dataset= len(f.readlines())
 #print(dataset)
 #global Variables
 spot=np.zeros(dataset) 
 btot=np.zeros(dataset) #floating points
 bx=np.zeros(dataset) #floating points
 by=np.zeros(dataset) #floating points
 bz=np.zeros(dataset) #floating points
 bzgsm=np.zeros(dataset) #floating points
 bygsm=np.zeros(dataset) #floating points

 speed=np.zeros(dataset) #floating points
 speedx=np.zeros(dataset) #floating points
 speed_phi=np.zeros(dataset) #floating points
 speed_theta=np.zeros(dataset) #floating points

 dst=np.zeros(dataset) #float
 kp=np.zeros(dataset) #float

 den=np.zeros(dataset) #float
 pdyn=np.zeros(dataset) #float
 year=np.zeros(dataset)
 day=np.zeros(dataset)
 hour=np.zeros(dataset)
 t=np.zeros(dataset) #index time
 
 
 j=0
 print('Read OMNI2 data ...')
 with open('omni2_all_years.dat') as f:
  for line in f:
   line = line.split() # to deal with blank 
   #print line #41 is Dst index, in nT
   dst[j]=line[40]
   kp[j]=line[38]
   
   if dst[j] == 99999: dst[j]=np.NaN
   #40 is sunspot number
   spot[j]=line[39]
   if spot[j] == 999: spot[j]=NaN

   #25 is bulkspeed F6.0, in km/s
   speed[j]=line[24]
   if speed[j] == 9999: speed[j]=np.NaN
 
   #get speed angles F6.1
   speed_phi[j]=line[25]
   if speed_phi[j] == 999.9: speed_phi[j]=np.NaN

   speed_theta[j]=line[26]
   if speed_theta[j] == 999.9: speed_theta[j]=np.NaN
   #convert speed to GSE x see OMNI website footnote
   speedx[j] = - speed[j] * np.cos(np.radians(speed_theta[j])) * np.cos(np.radians(speed_phi[j]))



   #9 is total B  F6.1 also fill ist 999.9, in nT
   btot[j]=line[9]
   if btot[j] == 999.9: btot[j]=np.NaN

   #GSE components from 13 to 15, so 12 to 14 index, in nT
   bx[j]=line[12]
   if bx[j] == 999.9: bx[j]=np.NaN
   by[j]=line[13]
   if by[j] == 999.9: by[j]=np.NaN
   bz[j]=line[14]
   if bz[j] == 999.9: bz[j]=np.NaN
 
   #GSM
   bygsm[j]=line[15]
   if bygsm[j] == 999.9: bygsm[j]=np.NaN
 
   bzgsm[j]=line[16]
   if bzgsm[j] == 999.9: bzgsm[j]=np.NaN 	
 
 
   #24 in file, index 23 proton density /ccm
   den[j]=line[23]
   if den[j] == 999.9: den[j]=np.NaN
 
   #29 in file, index 28 Pdyn, F6.2, fill values sind 99.99, in nPa
   pdyn[j]=line[28]
   if pdyn[j] == 99.99: pdyn[j]=np.NaN 		
 
   year[j]=line[0]
   day[j]=line[1]
   hour[j]=line[2]
   j=j+1     
   

 #convert time to matplotlib format
 #http://docs.sunpy.org/en/latest/guide/time.html
 #http://matplotlib.org/examples/pylab_examples/date_demo2.html

 times1=np.zeros(len(year)) #datetime time
 print('convert time start')
 for index in range(0,len(year)):
      #first to datetimeobject 
      timedum=datetime.datetime(int(year[index]), 1, 1) + datetime.timedelta(day[index] - 1) +datetime.timedelta(hours=hour[index])
      #then to matlibplot dateformat:
      times1[index] = mdates.date2num(timedum)
 print('convert time done')   #for time conversion

 print('all done.')
 print(j, ' datapoints')   #for reading data from OMNI file
 
 #make structured array of data
 omni_data=np.rec.array([times1,btot,bx,by,bz,bygsm,bzgsm,speed,speedx,den,pdyn,dst,kp,spot], \
 dtype=[('time','f8'),('btot','f8'),('bx','f8'),('by','f8'),('bz','f8'),\
 ('bygsm','f8'),('bzgsm','f8'),('speed','f8'),('speedx','f8'),('den','f8'),('pdyn','f8'),('dst','f8'),('kp','f8'), ('spot','f8')])
 
 return omni_data
 
 



import pickle
import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import sunpy.time




def IDL_time_to_num(time_in):  
 #convert IDL time to matplotlib datetime
 time_num=np.zeros(np.size(time_in))
 for ii in np.arange(0,np.size(time_in)):
   time_num[ii]=mdates.date2num(sunpy.time.parse_time(time_in[ii]))
   
 return time_num 
  

print ('read VEX')
#get insitu data
vex= pickle.load( open( "../catpy/DATACAT/VEX_2007to2014_SCEQ_removed.p", "rb" ) )

#time conversion
vex_time=IDL_time_to_num(vex.time)
print( 'read VEX done.')



plt.figure(1)
plt.plot_date(vex_time,vex.btot,'-k')
plt.show()

'''
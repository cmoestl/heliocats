#stats.py
#statistics stuff for heliocats
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy.io
import urllib
import os
from config import data_path
from sunpy.time import parse_time


'''
import copy
import matplotlib.dates as mdates
import matplotlib
import seaborn as sns
import datetime

import json
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

import heliopy.data.spice as spicedata
import heliopy.spice as spice
'''


####################################### 



def is_float(value):
  try:
    float(value)
    return True
  except:
    return False



def running_mean(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]


def expon(x, a, k, b):
    return a*np.exp(k*x) + b




def powerlaw(x, a, b):
    return a*x**b


def gaussian_nox0(x, a,sigma):
    return a*np.exp(-(x)**2/(2*sigma**2))


def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))



def hathaway_fit(x,x0, a, b, c):
    #for numerical fitting of hathaway function
    #x1 is input in days, convert to months
    x1=(x-x0)/30.42
    #print(x1[0])
    hatfunc=a*(((x1)/b)**3) * 1/(np.exp((((x1)/b)**2))-c)
    return hatfunc


def hathaway(x,x0, a, b, c):
    #Hathaway 2015 equation 6 page 40
    #average cycle sunspot number 
    #A=100 #amplitude ##195 for sunspot
    #b=100*12 #56*12 for months to days
    #c=0.8
    #4 free parameters A, b, c, t0
    
    #b=b*12 #months to days    
    x1=(parse_time(x).plot_date-parse_time(x0).plot_date)/30.42 #days
    #print(len(x1),len(x))
    
    hatfunc=a*(((x1)/b)**3) * 1/(np.exp((((x1)/b)**2))-c)
        
    return hatfunc
    



def dynamic_pressure(density, speed):
    '''
    make dynamic pressure from density and speed
    assume pdyn is only due to protons
    '''
    protonmass=1.6726219*1e-27  #kg
    pdyn=np.multiply(np.square(speed*1e3),density)*1e6*protonmass*1e9  #in nanoPascal
    
    return pdyn
  


def load_url_current_directory(filename,url):
    '''
    loads a file from any url to the current directory
    I use owncloud for the direct url links, 
    also works for dropbox when changing the last 0 to 1 in the url-> gives a direct link to files
    '''
    
    if not os.path.exists(filename):
        print('download file ', filename, ' from')
        print(url)
        try: 
            urllib.request.urlretrieve(url, filename)
            print('done')
        except urllib.error.URLError as e:
            print(' ', data_url,' ',e.reason)


def getcat(filename):
    cat = scipy.io.readsav(filename, verbose=False)
    return cat


def decode_array(bytearrin):
    '''
    for decoding the strings from the IDL .sav file to a list of python strings, not bytes
    make list of python lists with arbitrary length
    '''
    bytearrout = ['' for x in range(len(bytearrin))]
    for i in range(0, len(bytearrin)):
        bytearrout[i] = bytearrin[i].decode()
    # has to be np array so to be used with numpy "where"
    bytearrout = np.array(bytearrout)
    
    return bytearrout

    




'''
def hathaway(x, a, b, c, x0):
    return a*(((x-x0)/b)**3) * 1/(np.exp((((x-x0)/b)**2))-c)

paramh = scipy.optimize.curve_fit(hathaway, wt,wd)



#gaussian fits for daily new cases
newcases=np.gradient(cases)
paramg = scipy.optimize.curve_fit(gaussian, dates1, newcases )
ygfit=gaussian(dates_fut1,paramg[0][0],paramg[0][1],paramg[0][2])
#plt.plot_date(dates,newcases)
#plt.plot_date(dates_fut,ygfit,'-k')
print('Gaussian fit parameters:',paramg[0])

tests for gaussians and Hathaway's function for solar cycle, not used
Wind
tfit=mdates.date2num(sunpy.time.parse_time('2009-04-01'))+np.arange(0,365*10)
t0=mdates.date2num(sunpy.time.parse_time('2009-01-01'))

Gaussian
sigma=1000
bfitmax=30
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit-mu)**2)/(2*sigma**2) )
normalize with 1/max(ygauss)
plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='mediumseagreen',linestyle='-',markersize=0, label='Earth fit')

Hathaway 2015 equation 6 page 40
average cycle sunspot number 
A=100 amplitude 195 for sunspot
b=100*12 56*12 for months to days
c=0.8
4 free parameters A, b, c, t0

Fwind=A*(((tfit-t0)/b)**3) * 1/(np.exp((((tfit-t0)/b)**2))-c)
plt.plot_date(tfit, Fwind,'o',color='mediumseagreen',linestyle='-',markersize=0, label='Earth fit')

xaxis: 10 years, daily data point
xfit2=mdates.date2num(sunpy.time.parse_time('2007-01-01'))+np.arange(0,365*10)
MESSENGER
sigma=1000
bfitmax=10
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit2-mu)**2)/(2*sigma**2) )
normalize with 1/max(ygauss)
plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='darkgrey',linestyle='-',markersize=0, label='Mercury fit')

VEX
inital guess
sigma=1000
bfitmax=20
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit2-mu)**2)/(2*sigma**2) )
normalize with 1/max(ygauss)
plt.plot_date(xfit2, ygauss*1/max(ygauss)*bfitmax,'o',color='orange',linestyle='-',markersize=0, label='Venus fit')

for Mars: reconstruct likely parameters if sigma is quite similar for all fits, take mean of those sigmas and adjust bfitmax as function of distance with power law)
plot reconstructed function for Mars
bfitmax=40
plt.plot_date(xfit2, Fwind,'o',color='steelblue',linestyle='--',markersize=0, label='Mars reconstr.')
'''

'''
def hathaway(x, a, b, c, x0):
    return a*(((x-x0)/b)**3) * 1/(np.exp((((x-x0)/b)**2))-c)

paramh = scipy.optimize.curve_fit(hathaway, wt,wd)



#gaussian fits for daily new cases
newcases=np.gradient(cases)
paramg = scipy.optimize.curve_fit(gaussian, dates1, newcases )
ygfit=gaussian(dates_fut1,paramg[0][0],paramg[0][1],paramg[0][2])
#plt.plot_date(dates,newcases)
#plt.plot_date(dates_fut,ygfit,'-k')
print('Gaussian fit parameters:',paramg[0])

tests for gaussians and Hathaway's function for solar cycle, not used
Wind
tfit=mdates.date2num(sunpy.time.parse_time('2009-04-01'))+np.arange(0,365*10)
t0=mdates.date2num(sunpy.time.parse_time('2009-01-01'))

Gaussian
sigma=1000
bfitmax=30
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit-mu)**2)/(2*sigma**2) )
normalize with 1/max(ygauss)
plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='mediumseagreen',linestyle='-',markersize=0, label='Earth fit')

Hathaway 2015 equation 6 page 40
average cycle sunspot number 
A=100 amplitude 195 for sunspot
b=100*12 56*12 for months to days
c=0.8
4 free parameters A, b, c, t0

Fwind=A*(((tfit-t0)/b)**3) * 1/(np.exp((((tfit-t0)/b)**2))-c)
plt.plot_date(tfit, Fwind,'o',color='mediumseagreen',linestyle='-',markersize=0, label='Earth fit')

xaxis: 10 years, daily data point
xfit2=mdates.date2num(sunpy.time.parse_time('2007-01-01'))+np.arange(0,365*10)
MESSENGER
sigma=1000
bfitmax=10
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit2-mu)**2)/(2*sigma**2) )
normalize with 1/max(ygauss)
plt.plot_date(xfit, ygauss*1/max(ygauss)*bfitmax,'o',color='darkgrey',linestyle='-',markersize=0, label='Mercury fit')

VEX
inital guess
sigma=1000
bfitmax=20
mu=mdates.date2num(sunpy.time.parse_time('2013-01-01'))
ygauss=1/(sigma*np.sqrt(2*np.pi))*np.exp(-((xfit2-mu)**2)/(2*sigma**2) )
normalize with 1/max(ygauss)
plt.plot_date(xfit2, ygauss*1/max(ygauss)*bfitmax,'o',color='orange',linestyle='-',markersize=0, label='Venus fit')

for Mars: reconstruct likely parameters if sigma is quite similar for all fits, take mean of those sigmas and adjust bfitmax as function of distance with power law)
plot reconstructed function for Mars
bfitmax=40
plt.plot_date(xfit2, Fwind,'o',color='steelblue',linestyle='--',markersize=0, label='Mars reconstr.')




#old version of B fits
# print('Fit results for B in Form: y=B0*x^k')
# ####### power law fits for all events
# bmaxfit=np.polyfit(np.log10(ic.mo_sc_heliodistance),np.log10(ic.mo_bmax),1)
# b=10**bmaxfit[1]
# bmaxfitfun=b*(xfit**bmaxfit[0])
# print('bmax:       ',round(10**bmaxfit[1],2),' x ^', round(bmaxfit[0],2))

# bmeanfit=np.polyfit(np.log10(ic.mo_sc_heliodistance),np.log10(ic.mo_bmean),1)
# b=10**bmeanfit[1]
# bmeanfitfun=b*(xfit**bmeanfit[0])
# print('bmean:      ', round(10**bmeanfit[1],2),' x ^',round(bmeanfit[0],2))

# ##fit with only minimum events
# bmeanfit_min=np.polyfit(np.log10(ic.mo_sc_heliodistance[iall_min]),np.log10(ic.mo_bmean[iall_min]),1)
# bmeanfitfun_min=(10**bmeanfit_min[1])*(xfit**bmeanfit_min[0])
# print('bmean_min:  ', round(10**bmeanfit_min[1],2),' x ^', round(bmeanfit_min[0],2))

# ##fit with only rising events
# bmeanfit_rise=np.polyfit(np.log10(ic.mo_sc_heliodistance[iall_rise]),np.log10(ic.mo_bmean[iall_rise]),1)
# bmeanfitfun_rise=(10**bmeanfit_rise[1])*(xfit**bmeanfit_rise[0])
# print('bmean_rise: ', round(10**bmeanfit_rise[1],2),' x ^', round(bmeanfit_rise[0],2))

# ##fit with only maximum events
# bmeanfit_max=np.polyfit(np.log10(ic.mo_sc_heliodistance[iall_max]),np.log10(ic.mo_bmean[iall_max]),1)
# bmeanfitfun_max=(10**bmeanfit_max[1])*(xfit**bmeanfit_max[0])
# print('bmean_max:  ', round(10**bmeanfit_max[1],2),' x ^',round(bmeanfit_max[0],2))



'''
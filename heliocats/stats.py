#stats.py
#statistics stuff for heliocats
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy.io
import urllib
import os
from input import *


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


def expon(x, a, k, b):
    return a*np.exp(k*x) + b


def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def hathaway(x, amp, mu, sig):
   return amp * exp(-(x-cen)**2 /wid)



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

    





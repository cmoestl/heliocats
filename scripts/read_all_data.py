#for updating data every day for Wind and STEREO-A
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy
import copy
import matplotlib.dates as mdates
import datetime
import urllib
import json
import os
import pdb
from sunpy.time import parse_time
import scipy.io
import pickle
import sys
import cdflib
import matplotlib.pyplot as plt
import heliosat
from numba import njit
from astropy.time import Time
import heliopy.data.cassini as cassinidata
import heliopy.data.helios as heliosdata
import heliopy.spice as spice
import astropy


data_path='/nas/helio/data/insitu_python/'


filesta="sta_2018_now.p" 
filewin="wind_2018_now.p" 
[win,hwin]=pickle.load(open(data_path+filewin, "rb" ) )  
[sta,hsta]=pickle.load(open(data_path+filesta, "rb" ) ) 






#stats.py
#statistics stuff for heliocats
#https://github.com/cmoestl/heliocats

import numpy as np
import pandas as pd
import scipy
import copy
import matplotlib.dates as mdates
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

import heliopy.data.spice as spicedata
import heliopy.spice as spice

data_path='/nas/helio/data/insitu_python/'


####################################### 

#def (sc, start, end, sc_label, path, **kwargs):

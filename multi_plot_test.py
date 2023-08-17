#!/usr/bin/env python
# coding: utf-8

# In[1]:


#demonstrator for plotting with matplotlib with multiprocessing

#https://docs.python.org/3/library/multiprocessing.html

import multiprocessing as mp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import time
import sys
import os

matplotlib.use('Agg')

#### check for system type
#server
if sys.platform == 'linux': 
    print('system is linux')
    used=50
    
#mac
if sys.platform =='darwin':  
    print('system is mac')
    print(os.system('sysctl -n hw.physicalcpu'))
    print(os.system('sysctl -n hw.logicalcpu'))
    used=8
    
print()

plotsdir='plots_multi' 

if os.path.isdir(plotsdir) == False: os.mkdir(plotsdir)


# In[9]:


def make_plot(i):
    

    plt.figure(1,figsize=(20,10),dpi=100)
    plt.plot(data,'ok',markersize=0.1,alpha=0.1)
    plt.title(str(i)+title)    
    plt.savefig(plotsdir+'/'+str(i)+'.png',dpi=100)
    plt.close(1)

    
counter = np.arange(0,10,1)  # pool counter

#global variables
title ='Hello'    
data=  np.random.rand(1000000) #1 mio random numbers

print(data)
print(counter)
#print(len(np.arange(0,100,1e-4)))
#p=zip([parameters,data])


# In[10]:


t0 = time.time()


#define pool using fork and number of processes
#using fork works on both mac and linux
pool=mp.get_context('fork').Pool(processes=used)
print('Using multiprocessing, nr of cores',mp.cpu_count(), ', nr of processes used: ',used)

# Map the worker function onto the parameters    
pool.map(make_plot, counter)
pool.close()
pool.join()     

t1 = time.time()
multi_time=np.round(t1-t0,2)

print('done')
print('plotting takes', np.round(multi_time,2), 'seconds')    


# In[11]:


t0 = time.time()

for i in counter:
    make_plot(i)
    
t1 = time.time()
single_time=np.round(t1-t0,2)
print('plotting takes', np.round(single_time,2), 'seconds')      

print('multiprocessing is a factor,',np.round(single_time/multi_time,2), ' faster')


# In[ ]:


##example for writing in the same array with multiple processes



####


# In[ ]:





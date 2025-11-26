#!/usr/bin/env python
# coding: utf-8

# ## positions plot 3D
# 
# 
# Christian Möstl, Emma Davies, Eva Weiler
# 
# last update: November 2025
# 
# - Issues:
# 
# 
# 

# In[1]:


#switches
#debug_mode=0
#always turn off debug mode when deploying!

import pickle
import importlib
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
import sys
import numpy as np
import datetime
import os   
import time
import sunpy

import plotly.graph_objects as go
from plotly.offline import iplot, init_notebook_mode
from plotly.subplots import make_subplots
import plotly.io as pio
import plotly.express as px
pio.renderers.default = 'browser'

import astropy.constants as const
import astropy.units as u

from heliocats import data as hd

##### check for system type
#server
if sys.platform == 'linux': 
    print('system is linux')
    matplotlib.use('Agg') 
#mac
if sys.platform =='darwin':  
    print('system is mac')
    #for testing
    get_ipython().run_line_magic('matplotlib', 'inline')
    #matplotlib.use('Agg') 

au=const.au.value*1e-3 #au in km

################################################ CHECK  ##############################################

#make sure to convert the current notebook to a script
os.system('jupyter nbconvert --to script positions_3d.ipynb')   

####################################################################################################################

#test execution times
clock_start = time.time()


# ### load position files

# In[2]:


[psp, bepi, solo, sta, juice, earth, mercury, venus, mars, jupiter, saturn, uranus, neptune,l4,l5]=pickle.load( open( 'results/positions/positions_2020_all_HEEQ_1h_rad_cm.p', "rb" ) )   


# ### prepare data

# In[3]:


time1=mdates.date2num(datetime.datetime.utcnow())

#override current date
#time1=mdates.date2num(datetime.datetime(2028,1,1))

fadeind=24*80 #80 days for 1 hour resolution

#find current indices
dct=time1-psp.time
psp_timeind=np.argmin(abs(dct))

dct=time1-bepi.time
bepi_timeind=np.argmin(abs(dct))

dct=time1-solo.time
solo_timeind=np.argmin(abs(dct))

dct=time1-juice.time
juice_timeind=np.argmin(abs(dct))

dct=time1-earth.time
earth_timeind=np.argmin(abs(dct))

dct=time1-venus.time
venus_timeind=np.argmin(abs(dct))

dct=time1-mars.time
mars_timeind=np.argmin(abs(dct))

dct=time1-jupiter.time
jupiter_timeind=np.argmin(abs(dct))

dct=time1-mercury.time
mercury_timeind=np.argmin(abs(dct))

dct=time1-sta.time
sta_timeind=np.argmin(abs(dct))

dct=time1-l4.time
l4_timeind=np.argmin(abs(dct))

dct=time1-l5.time
l5_timeind=np.argmin(abs(dct))

#dct=time1-l1.time
#l1_timeind=np.argmin(abs(dct))


#current frame time
frame_time=mdates.num2date(earth['time'][earth_timeind])

#convert for plotting the time along the past and future trajectories
solo_time=mdates.num2date(solo.time)
psp_time=mdates.num2date(psp.time)
bepi_time=mdates.num2date(bepi.time)
sta_time=mdates.num2date(sta.time)


####################### SOLAR ORBITER
#for future times
solo_hover=np.stack((np.round(solo['r'][solo_timeind:solo_timeind+fadeind],3), 
                     np.round(np.rad2deg(solo['lon'][solo_timeind:solo_timeind+fadeind]),1),
                     np.round(np.rad2deg(solo['lat'][solo_timeind:solo_timeind+fadeind]),1),solo_time[solo_timeind:solo_timeind+fadeind]),axis=-1)

#for past
solo_hover_past=np.stack((np.round(solo['r'][solo_timeind-fadeind:solo_timeind],3), 
                     np.round(np.rad2deg(solo['lon'][solo_timeind-fadeind:solo_timeind]),1),
                     np.round(np.rad2deg(solo['lat'][solo_timeind-fadeind:solo_timeind]),1),solo_time[solo_timeind-fadeind:solo_timeind]), axis=-1)
#for now
solo_hover_now=np.stack((np.round(solo['r'][solo_timeind],3), 
                     np.round(np.rad2deg(solo['lon'][solo_timeind]),1),
                     np.round(np.rad2deg(solo['lat'][solo_timeind]),1),solo_time[solo_timeind]), axis=-1)


###################### PSP
#for future times
psp_hover=np.stack((np.round(psp['r'][psp_timeind:psp_timeind+fadeind],3), 
                     np.round(np.rad2deg(psp['lon'][psp_timeind:psp_timeind+fadeind]),1),
                     np.round(np.rad2deg(psp['lat'][psp_timeind:psp_timeind+fadeind]),1),psp_time[psp_timeind:psp_timeind+fadeind]), axis=-1)

#for past
psp_hover_past=np.stack((np.round(psp['r'][psp_timeind-fadeind:psp_timeind],3), 
                     np.round(np.rad2deg(psp['lon'][psp_timeind-fadeind:psp_timeind]),1),
                     np.round(np.rad2deg(psp['lat'][psp_timeind-fadeind:psp_timeind]),1),psp_time[psp_timeind-fadeind:psp_timeind]), axis=-1)

#for now
psp_hover_now=np.stack((np.round(psp['r'][psp_timeind],3), 
                     np.round(np.rad2deg(psp['lon'][psp_timeind]),1),
                     np.round(np.rad2deg(psp['lat'][psp_timeind]),1),psp_time[psp_timeind]), axis=-1)


##################### STEREO-A
#for future times
sta_hover=np.stack((np.round(sta['r'][sta_timeind:sta_timeind+fadeind],3), 
                     np.round(np.rad2deg(sta['lon'][sta_timeind:sta_timeind+fadeind]),1),
                     np.round(np.rad2deg(sta['lat'][sta_timeind:sta_timeind+fadeind]),1),sta_time[sta_timeind:sta_timeind+fadeind]), axis=-1)

#for past
sta_hover_past=np.stack((np.round(sta['r'][sta_timeind-fadeind:sta_timeind],3), 
                     np.round(np.rad2deg(sta['lon'][sta_timeind-fadeind:sta_timeind]),1),
                     np.round(np.rad2deg(sta['lat'][sta_timeind-fadeind:sta_timeind]),1),sta_time[sta_timeind-fadeind:sta_timeind]), axis=-1)

#for now
sta_hover_now=np.stack((np.round(sta['r'][sta_timeind],3), 
                     np.round(np.rad2deg(sta['lon'][sta_timeind]),1),
                     np.round(np.rad2deg(sta['lat'][sta_timeind]),1),sta_time[sta_timeind]), axis=-1)


#################### BEPI COLOMBO
#for future times
bepi_hover=np.stack((np.round(bepi['r'][bepi_timeind:bepi_timeind+fadeind],3), 
                     np.round(np.rad2deg(bepi['lon'][bepi_timeind:bepi_timeind+fadeind]),1),
                     np.round(np.rad2deg(bepi['lat'][bepi_timeind:bepi_timeind+fadeind]),1),bepi_time[bepi_timeind:bepi_timeind+fadeind]), axis=-1)

#for past
bepi_hover_past=np.stack((np.round(bepi['r'][bepi_timeind-fadeind:bepi_timeind],3), 
                     np.round(np.rad2deg(bepi['lon'][bepi_timeind-fadeind:bepi_timeind]),1),
                     np.round(np.rad2deg(bepi['lat'][bepi_timeind-fadeind:bepi_timeind]),1),bepi_time[bepi_timeind-fadeind:bepi_timeind]), axis=-1)

#for now
bepi_hover_now=np.stack((np.round(bepi['r'][bepi_timeind],3), 
                     np.round(np.rad2deg(bepi['lon'][bepi_timeind]),1),
                     np.round(np.rad2deg(bepi['lat'][bepi_timeind]),1),bepi_time[bepi_timeind]), axis=-1)


#solo_hover2=np.stack((solo['r'], np.rad2deg(solo['lon']),np.rad2deg(solo['lat'])), axis=-1)
#print(solo_hover_now)
#solo_hover_now[2]


# ### Make figure

# In[4]:


fig = go.Figure()

mfac=15
sunsize=5 #solar radii
planetsize=3 #solar radii, Earth size

######################### SUN AND PLANETS

# Create data for a sphere
theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)
theta, phi = np.meshgrid(theta, phi)

r = (700*1e3)/(149.5*1e6)*sunsize

# Convert spherical coordinates to Cartesian coordinates
[x,y,z]=hd.sphere2cart(r, theta,phi)
# Create 3D surface plot

#Sun
fig.add_trace(go.Surface(x=x, y=y,z=z, colorscale='hot', showscale=False, name='5 R_Sun')) #5 solar radii

#Earth
r = (700*1e3)/(149.5*1e6)*planetsize # 3 solar radii
[x,y,z]=hd.sphere2cart(r, theta,phi)
fig.add_trace(go.Surface(x=x+earth['x'][earth_timeind]/au, y=y+earth['y'][earth_timeind]/au,z=z+earth['z'][earth_timeind]/au, colorscale='Blues', showscale=False, name='Earth'))

#mercury
r = (700*1e3)/(149.5*1e6)*planetsize*2440/6371  # 3 solar radii times factor relative to Earth
[x,y,z]=hd.sphere2cart(r, theta,phi)
fig.add_trace(go.Surface(x=x+mercury['x'][mercury_timeind]/au, y=y+mercury['y'][mercury_timeind]/au,z=z+mercury['z'][mercury_timeind]/au, colorscale='gray', showscale=False, name='Mercury'))

#venus
r = (700*1e3)/(149.5*1e6)*planetsize*6052/6371 # 3 solar radii times factor relative to Earth
[x,y,z]=hd.sphere2cart(r, theta,phi)
fig.add_trace(go.Surface(x=x+venus['x'][venus_timeind]/au, y=y+venus['y'][venus_timeind]/au,z=z+venus['z'][venus_timeind]/au, colorscale='solar', showscale=False, name='Venus'))

#mars
r = (700*1e3)/(149.5*1e6)*planetsize*3390/6371 # 3 solar radii times factor relative to Earth
[x,y,z]=hd.sphere2cart(r, theta,phi)
fig.add_trace(go.Surface(x=x+mars['x'][mars_timeind]/au, y=y+mars['y'][mars_timeind]/au,z=z+mars['z'][mars_timeind]/au, colorscale='magma', showscale=False, name='Mars'))

################### add solar system geometries

gridc='white'
num_points = 100
# Create theta values (angles) for the circle
theta_values = np.linspace(0, 2*np.pi, num_points)
r = 1  # radius of the circle
x_values = r * np.cos(theta_values)
y_values = r * np.sin(theta_values)

fig.add_trace(go.Scatter3d(
    x=x_values,
    y=y_values,
    z=np.zeros(num_points),  # Set z-values to zero for 2D appearance
    mode='lines', name='1 AU',
    line=dict(color=gridc, width=1)
))

fig.add_trace(go.Scatter3d(
    x=x_values*0.86,
    y=y_values*0.86,
    z=np.zeros(num_points),  # Set z-values to zero for 2D appearance
    mode='lines', name='0.86 AU',
    line=dict(color=gridc, width=1)
))

fig.add_trace(go.Scatter3d(
    x=x_values*0.72,
    y=y_values*0.72,
    z=np.zeros(num_points),  # Set z-values to zero for 2D appearance
    mode='lines', name='0.7 AU',
    line=dict(color=gridc, width=1)
))

fig.add_trace(go.Scatter3d(
    x=x_values*0.3,
    y=y_values*0.3,
    z=np.zeros(num_points),  # Set z-values to zero for 2D appearance
    mode='lines', name='0.3 AU',
    line=dict(color=gridc, width=1)
))

fig.add_trace(go.Scatter3d(
    x=x_values*0.1,
    y=y_values*0.1,
    z=np.zeros(num_points),  # Set z-values to zero for 2D appearance
    mode='lines', name='0.1 AU',
    line=dict(color=gridc, width=1)
))

#add Sun- 1 AU line for HEEQ latitude 0
fig.add_trace(go.Scatter3d(
    x=np.linspace(0,1,num_points),
    y=np.zeros(num_points),
    z=np.zeros(num_points), 
    mode='lines', name='Sun-1 AU line',
    line=dict(color=gridc, width=1)
))


#add L5
fig.add_trace(go.Scatter3d(
    x=[l5['x'][l5_timeind]/au],
    y=[l5['y'][l5_timeind]/au],
    z=[l5['z'][l5_timeind]/au], 
    mode='markers', name='L5',
    marker=dict(
        size=7,
        symbol='diamond', 
        color=['gold', 'yellow'],
        )
    ))


#add L4
fig.add_trace(go.Scatter3d(
    x=[l4['x'][l4_timeind]/au],
    y=[l4['y'][l4_timeind]/au],
    z=[l4['z'][l4_timeind]/au], 
    mode='markers', name='L4',
    marker=dict(
        size=7,
        symbol='diamond', 
        color=['purple', 'blue'],
        )
    ))   

######################### parker spiral
res_in_days=1/24
k=0
sun_rot=26.24
v=400/au #km/s
r0=695000/au

omega=2*np.pi/(sun_rot*60*60*24) #solar rotation in seconds

def polar2cart(r,lon):
    x = r * np.cos(lon)
    y = r * np.sin(lon)
    return (x, y)

plon=np.arange(0,np.deg2rad(90),0.01)
pr=np.zeros(len(plon))
plat=np.zeros(len(plon))

pz=np.zeros(len(plon)) #for 3D plot

for i in np.arange(0,12):
        pr=v/omega*(plon)+r0       
        [px,py]=polar2cart(pr,-plon+np.deg2rad(360/12*i))
        fig.add_trace(go.Scatter3d(x=px, y=py, z=pz,mode='lines',line=dict(color='white', width=0.5), showlegend=False))
   
   


msize=5  #markersize

################### Orbits future and past

# SOLO 
fig.add_trace(go.Scatter3d(x=solo['x'][solo_timeind:solo_timeind+fadeind]/au, y=solo['y'][solo_timeind:solo_timeind+fadeind]/au, 
                           z=solo['z'][solo_timeind:solo_timeind+fadeind]/au, name='Solar Orbiter',mode='lines',line=dict(color='lightgreen', width=3),
                           customdata=solo_hover,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]}',showlegend=False  ))

fig.add_trace(go.Scatter3d(x=solo['x'][solo_timeind-fadeind:solo_timeind]/au, y=solo['y'][solo_timeind-fadeind:solo_timeind]/au, 
                           z=solo['z'][solo_timeind-fadeind:solo_timeind]/au, name='Solar Orbiter',mode='lines',line=dict(color='lightgreen', dash='dash' , width=3), 
                           customdata=solo_hover_past,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]}', showlegend=False))

#PSP
fig.add_trace(go.Scatter3d(x=psp['x'][psp_timeind:psp_timeind+fadeind]/au, y=psp['y'][psp_timeind:solo_timeind+fadeind]/au, 
                           z=psp['z'][psp_timeind:psp_timeind+fadeind]/au, name='PSP',mode='lines',line=dict(color='white', width=3),
                           customdata=psp_hover,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]}',showlegend=False  ))

fig.add_trace(go.Scatter3d(x=psp['x'][psp_timeind-fadeind:psp_timeind]/au, y=psp['y'][psp_timeind-fadeind:psp_timeind]/au, 
                           z=psp['z'][psp_timeind-fadeind:psp_timeind]/au, name='PSP',mode='lines',line=dict(color='white', dash='dash' , width=3), 
                           customdata=psp_hover_past,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]}', showlegend=False))

#Bepi
fig.add_trace(go.Scatter3d(x=bepi['x'][bepi_timeind:bepi_timeind+fadeind]/au, y=bepi['y'][bepi_timeind:bepi_timeind+fadeind]/au, 
                           z=bepi['z'][bepi_timeind:bepi_timeind+fadeind]/au, name='BepiColombo',mode='lines',line=dict(color='lightskyblue', width=3),
                           customdata=bepi_hover,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]}',showlegend=False  ))

fig.add_trace(go.Scatter3d(x=bepi['x'][bepi_timeind-fadeind:bepi_timeind]/au, y=bepi['y'][bepi_timeind-fadeind:bepi_timeind]/au, 
                           z=bepi['z'][bepi_timeind-fadeind:bepi_timeind]/au, name='BepiColomb',mode='lines',line=dict(color='lightskyblue',dash='dash', width=3), 
                           customdata=bepi_hover_past,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]}',showlegend=False))

#STEREO-A
fig.add_trace(go.Scatter3d(x=sta['x'][sta_timeind:sta_timeind+fadeind]/au, y=sta['y'][sta_timeind:sta_timeind+fadeind]/au, 
                           z=sta['z'][sta_timeind:sta_timeind+fadeind]/au, name='STEREO-A',mode='lines',line=dict(color='tomato', width=3),
                           customdata=sta_hover,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]}',showlegend=False))

fig.add_trace(go.Scatter3d(x=sta['x'][sta_timeind-fadeind:sta_timeind]/au, y=sta['y'][sta_timeind-fadeind:sta_timeind]/au, 
                           z=sta['z'][sta_timeind-fadeind:sta_timeind]/au, name='STEREO-A',mode='lines',line=dict(color='tomato',dash='dash', width=3), 
                           customdata=sta_hover_past,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]}',showlegend=False))



#### spacecraft positions

fig.add_trace(go.Scatter3d(x=[solo['x'][solo_timeind]/au], y=[solo['y'][solo_timeind]/au], z=[solo['z'][solo_timeind]/au], name='Solar Orbiter',
        mode='markers',marker=dict(color='lightgreen', size=msize, symbol=['square']),
        customdata=[solo_hover_now],hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]} '+'<extra></extra>'))

fig.add_trace(go.Scatter3d(x=[psp['x'][psp_timeind]/au], y=[psp['y'][psp_timeind]/au], z=[psp['z'][psp_timeind]/au], name='PSP',
        mode='markers',marker=dict(color='white', size=msize, symbol=['square']),
        customdata=[psp_hover_now],hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]} '+'<extra></extra>'))

fig.add_trace(go.Scatter3d(x=[sta['x'][sta_timeind]/au], y=[sta['y'][sta_timeind]/au], z=[sta['z'][sta_timeind]/au], name='STEREO-A',
        mode='markers',marker=dict(color='tomato', size=msize, symbol=['square']),
        customdata=[sta_hover_now],hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]} '+'<extra></extra>'))

fig.add_trace(go.Scatter3d(x=[bepi['x'][bepi_timeind]/au], y=[bepi['y'][bepi_timeind]/au], z=[bepi['z'][bepi_timeind]/au], name='BepiColombo',
        mode='markers',marker=dict(color='lightskyblue', size=msize, symbol=['square']),
        customdata=[bepi_hover_now],hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°<br>%{customdata[3]} '+'<extra></extra>'))


######## set camera position
zoom=1.0
fig.update_layout(
    scene=dict( aspectmode='data',
        camera=dict(
            eye=dict(x=0, y=-zoom, z=zoom),  # Set the position of the camera
            center=dict(x=0, y=0, z=0),      # Set the point the camera is looking at
            up=dict(x=0, y=0, z=1),          # Set the up vector of the camera
    ))
)

########### Update layout for black background and styling
fig.update_layout(
    title='Spacecraft positions '+frame_time.strftime('%Y %b %d %H:00 UTC'),
    scene=dict(
        xaxis=dict(title='X [HEEQ]', gridcolor='#444', zerolinecolor='#666', backgroundcolor='black'),
        yaxis=dict(title='Y [HEEQ]', gridcolor='#444', zerolinecolor='#666', backgroundcolor='black'),
        zaxis=dict(title='Z [HEEQ]', gridcolor='#444', zerolinecolor='#666', backgroundcolor='black'),
        bgcolor='black'
    ),
    paper_bgcolor='black',
    plot_bgcolor='black',
    font=dict(color='white'),
    margin=dict(l=0, r=0, t=50, b=0)
)

fig.show()        
#save as image
#pio.write_image(fig, 'results/positions/position_3D.png',scale=1, width=2000, height=1200)

fig.write_html(f'results/positions/position_3D.html')

##pio.write_image(fig, 'results/positions/position_3D.png',scale=2, width=1500, height=850)


# In[5]:


clock_end = time.time()
print('done, all took ',np.round(clock_end - clock_start,2),' seconds.')

sys.exit()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ### tests

# In[ ]:


#######time slider example


import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Create figure
fig = go.Figure()

# Create frames for animation
frames = []
num_frames = 20

for frame_idx in range(num_frames):
    # Time parameter that changes with each frame
    t = frame_idx / (num_frames - 1)
    
    # Generate data for a 3D parametric curve
    n = 100
    theta = np.linspace(0, 4 * np.pi, n)
    
    # The shape evolves with time t
    z = np.sin(theta + t * 2 * np.pi) * 2
    r = 2 + np.cos(theta + t * np.pi)
    x = r * np.cos(3 * theta)
    y = r * np.sin(3 * theta)
    
    # Create colormap values
    colorscale = np.linspace(0, 1, n)
    
    # Create a frame for this time step
    frames.append(
        go.Frame(
            data=[go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode='markers',
                marker=dict(
                    size=5,
                    color=colorscale,
                    colorscale='Viridis',
                    opacity=0.8
                )
            )],
            name=f"frame{frame_idx}"
        )
    )

# Initial data (first frame)
t_initial = 0
theta = np.linspace(0, 4 * np.pi, n)
z_initial = np.sin(theta + t_initial * 2 * np.pi) * 2
r_initial = 2 + np.cos(theta + t_initial * np.pi)
x_initial = r_initial * np.cos(3 * theta)
y_initial = r_initial * np.sin(3 * theta)

# Add initial data
fig.add_trace(
    go.Scatter3d(
        x=x_initial,
        y=y_initial,
        z=z_initial,
        mode='markers',
        marker=dict(
            size=5,
            color=np.linspace(0, 1, n),
            colorscale='Viridis',
            opacity=0.8
        )
    )
)

# Add frames to the figure
fig.frames = frames

# Create slider
sliders = [{
    'active': 0,
    'yanchor': 'top',
    'xanchor': 'left',
    'currentvalue': {
        'font': {'size': 16, 'color': 'white'},
        'prefix': 'Time: ',
        'visible': True,
        'xanchor': 'right'
    },
    'transition': {'duration': 100, 'easing': 'cubic-in-out'},
    'pad': {'b': 10, 't': 50},
    'len': 0.9,
    'x': 0.1,
    'y': 0,
    'steps': [
        {
            'args': [
                [f'frame{k}'],
                {
                    'frame': {'duration': 100, 'redraw': True},
                    'mode': 'immediate',
                    'transition': {'duration': 100}
                }
            ],
            'label': f'{k/(num_frames-1):.2f}',
            'method': 'animate'
        }
        for k in range(num_frames)
    ]
}]

# Update layout for black background and styling
fig.update_layout(
    title='3D Parametric Curve with Time Slider',
    scene=dict(
        xaxis=dict(title='X Axis', gridcolor='#444', zerolinecolor='#666', backgroundcolor='black'),
        yaxis=dict(title='Y Axis', gridcolor='#444', zerolinecolor='#666', backgroundcolor='black'),
        zaxis=dict(title='Z Axis', gridcolor='#444', zerolinecolor='#666', backgroundcolor='black'),
        bgcolor='black',
        camera=dict(
            eye=dict(x=1.5, y=1.5, z=1.2)
        )
    ),
    paper_bgcolor='black',
    plot_bgcolor='black',
    font=dict(color='white'),
    margin=dict(l=0, r=0, t=50, b=100),
    sliders=sliders,
    updatemenus=[{
        'buttons': [
            {
                'args': [None, {'frame': {'duration': 500, 'redraw': True}, 
                                'fromcurrent': True, 
                                'transition': {'duration': 300, 'easing': 'quadratic-in-out'}}],
                'label': 'Play',
                'method': 'animate'
            },
            {
                'args': [[None], {'frame': {'duration': 0, 'redraw': True}, 
                                'mode': 'immediate',
                                'transition': {'duration': 0}}],
                'label': 'Pause',
                'method': 'animate'
            }
        ],
        'direction': 'left',
        'pad': {'r': 10, 't': 87},
        'showactive': False,
        'type': 'buttons',
        'x': 0.1,
        'xanchor': 'right',
        'y': 0,
        'yanchor': 'top'
    }]
)

# Make the plot responsive
fig.update_layout(
    autosize=True,
    width=1500,
    height=1000,
)

# Show the plot
fig.show()
fig.write_html(f'results/positions/slider_3d_test.html')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# ## positions plot 3D
# 
# 
# Christian Möstl, Emma Davies, Eva Weiler
# 
# April 2025
# 
# - Issues: add times
# 
# 
# 

# In[1]:


#switches
debug_mode=0
#always turn off debug mode when deploying!


import pickle
import importlib
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.dates as mdates
import sys
import numpy as np
import datetime
import scipy.signal
import urllib
import json
import os   
import time
import h5py
import pytz
import copy
import cdflib
import sunpy
import pickle


import plotly.graph_objects as go
from plotly.offline import iplot, init_notebook_mode
from plotly.subplots import make_subplots
import plotly.io as pio
import plotly.express as px
pio.renderers.default = 'browser'

import astropy.constants as const

import astropy.units as u
from heliocats import data as hd
from heliocats import plot as hp

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

au=const.au.value*1e-3

################################################ CHECK  ##############################################

#make sure to convert the current notebook to a script
os.system('jupyter nbconvert --to script positions_3d.ipynb')   


####################################################################################################################

#test execution times
clock_start = time.time()


# In[2]:


[psp, bepi, solo, sta, juice, earth, mercury, venus, mars, jupiter, saturn, uranus, neptune,l4,l5]=pickle.load( open( 'results/positions/positions_2020_all_HEEQ_1h_rad_cm.p', "rb" ) )   


# In[3]:


time1=mdates.date2num(datetime.datetime.utcnow())


#override current date
#time1=mdates.date2num(datetime.datetime(2028,1,1))


#psp['time']
#psp_timeind:psp_timeind+fadeind

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

#dct=time1-l1.time
#l1_timeind=np.argmin(abs(dct))

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
#solo['x'][solo_timeind]/au
#solo['y'][solo_timeind]/au
#solo['z'][solo_timeind]/au


#current frame time
frame_time=mdates.num2date(earth['time'][earth_timeind])



##need to make custom data for plotting the position for the position cutout converted to spherical coordinates

#for future times
solo_hover=np.stack((np.round(solo['r'][solo_timeind:solo_timeind+fadeind],3), 
                     np.round(np.rad2deg(solo['lon'][solo_timeind:solo_timeind+fadeind]),1),
                     np.round(np.rad2deg(solo['lat'][solo_timeind:solo_timeind+fadeind]),1)), axis=-1)

#for past
solo_hover_past=np.stack((np.round(solo['r'][solo_timeind-fadeind:solo_timeind],3), 
                     np.round(np.rad2deg(solo['lon'][solo_timeind-fadeind:solo_timeind]),1),
                     np.round(np.rad2deg(solo['lat'][solo_timeind-fadeind:solo_timeind]),1)), axis=-1)

#for now
solo_hover_now=np.stack((np.round(solo['r'][solo_timeind],3), 
                     np.round(np.rad2deg(solo['lon'][solo_timeind]),1),
                     np.round(np.rad2deg(solo['lat'][solo_timeind]),1)), axis=-1)


#for future times
psp_hover=np.stack((np.round(psp['r'][psp_timeind:psp_timeind+fadeind],3), 
                     np.round(np.rad2deg(psp['lon'][psp_timeind:psp_timeind+fadeind]),1),
                     np.round(np.rad2deg(psp['lat'][psp_timeind:psp_timeind+fadeind]),1)), axis=-1)

#for past
psp_hover_past=np.stack((np.round(psp['r'][psp_timeind-fadeind:psp_timeind],3), 
                     np.round(np.rad2deg(psp['lon'][psp_timeind-fadeind:psp_timeind]),1),
                     np.round(np.rad2deg(psp['lat'][psp_timeind-fadeind:psp_timeind]),1)), axis=-1)

#for now
psp_hover_now=np.stack((np.round(psp['r'][psp_timeind],3), 
                     np.round(np.rad2deg(psp['lon'][psp_timeind]),1),
                     np.round(np.rad2deg(psp['lat'][psp_timeind]),1)), axis=-1)



#for future times
sta_hover=np.stack((np.round(sta['r'][sta_timeind:sta_timeind+fadeind],3), 
                     np.round(np.rad2deg(sta['lon'][sta_timeind:sta_timeind+fadeind]),1),
                     np.round(np.rad2deg(sta['lat'][sta_timeind:sta_timeind+fadeind]),1)), axis=-1)

#for past
sta_hover_past=np.stack((np.round(sta['r'][sta_timeind-fadeind:sta_timeind],3), 
                     np.round(np.rad2deg(sta['lon'][sta_timeind-fadeind:sta_timeind]),1),
                     np.round(np.rad2deg(sta['lat'][sta_timeind-fadeind:sta_timeind]),1)), axis=-1)

#for now
sta_hover_now=np.stack((np.round(sta['r'][sta_timeind],3), 
                     np.round(np.rad2deg(sta['lon'][sta_timeind]),1),
                     np.round(np.rad2deg(sta['lat'][sta_timeind]),1)), axis=-1)




#for future times
bepi_hover=np.stack((np.round(bepi['r'][bepi_timeind:bepi_timeind+fadeind],3), 
                     np.round(np.rad2deg(bepi['lon'][bepi_timeind:bepi_timeind+fadeind]),1),
                     np.round(np.rad2deg(bepi['lat'][bepi_timeind:bepi_timeind+fadeind]),1)), axis=-1)

#for past
bepi_hover_past=np.stack((np.round(bepi['r'][bepi_timeind-fadeind:bepi_timeind],3), 
                     np.round(np.rad2deg(bepi['lon'][bepi_timeind-fadeind:bepi_timeind]),1),
                     np.round(np.rad2deg(bepi['lat'][bepi_timeind-fadeind:bepi_timeind]),1)), axis=-1)

#for now
bepi_hover_now=np.stack((np.round(bepi['r'][bepi_timeind],3), 
                     np.round(np.rad2deg(bepi['lon'][bepi_timeind]),1),
                     np.round(np.rad2deg(bepi['lat'][bepi_timeind]),1)), axis=-1)







#solo_hover2=np.stack((solo['r'], np.rad2deg(solo['lon']),np.rad2deg(solo['lat'])), axis=-1)


#print(solo_hover_now)
#solo_hover_now[2]


# In[4]:


fig = go.Figure()

mfac=15


sunsize=8
planetsize=3

############# add Sun
# Create data for a sphere
theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)
theta, phi = np.meshgrid(theta, phi)



r = (700*1e3)/(149.5*1e6)*sunsize

# Convert spherical coordinates to Cartesian coordinates
[x,y,z]=hd.sphere2cart(r, theta,phi)
# Create 3D surface plot

#Sun
fig.add_trace(go.Surface(x=x, y=y,z=z, colorscale='hot', showscale=False, name='5 R_Sun'))



#Earth
r = (700*1e3)/(149.5*1e6)*planetsize # 1 solar radii
[x,y,z]=hd.sphere2cart(r, theta,phi)
fig.add_trace(go.Surface(x=x+earth['x'][earth_timeind]/au, y=y+earth['y'][earth_timeind]/au,z=z+earth['z'][earth_timeind]/au, colorscale='blugrn', showscale=False, name='Earth'))

#mercury
r = (700*1e3)/(149.5*1e6)*planetsize  # 1 solar radii
[x,y,z]=hd.sphere2cart(r, theta,phi)
fig.add_trace(go.Surface(x=x+mercury['x'][mercury_timeind]/au, y=y+mercury['y'][mercury_timeind]/au,z=z+mercury['z'][mercury_timeind]/au, colorscale='algae', showscale=False, name='Mercury'))


#venus
r = (700*1e3)/(149.5*1e6)*planetsize  # 1 solar radii
[x,y,z]=hd.sphere2cart(r, theta,phi)
fig.add_trace(go.Surface(x=x+venus['x'][venus_timeind]/au, y=y+venus['y'][venus_timeind]/au,z=z+venus['z'][venus_timeind]/au, colorscale='solar', showscale=False, name='Venus'))


#mars
r = (700*1e3)/(149.5*1e6)*planetsize  # 1 solar radii
[x,y,z]=hd.sphere2cart(r, theta,phi)
fig.add_trace(go.Surface(x=x+mars['x'][mars_timeind]/au, y=y+mars['y'][mars_timeind]/au,z=z+mars['z'][mars_timeind]/au, colorscale='magma', showscale=False, name='Mars'))


#   ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
#             'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
#             'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
#             'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
#             'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
#             'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
#             'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
#             'orrd', 'oryel', 'oxy', 'peach', 'phase', 'picnic', 'pinkyl',
#             'piyg', 'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn',
#             'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu',
#             'rdgy', 'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar',
#             'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn',
#             'tealrose', 'tempo', 'temps', 'thermal', 'tropic', 'turbid',
#             'turbo', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr',
#             'ylorrd'].


################### add circle at 1 AU


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
    x=x_values*0.7,
    y=y_values*0.7,
    z=np.zeros(num_points),  # Set z-values to zero for 2D appearance
    mode='lines', name='0.7 AU',
    line=dict(color=gridc, width=1)
))


#fig.add_trace(go.Scatter3d(
#    x=x_values*0.5,
#    y=y_values*0.5,
#    z=np.zeros(num_points),  # Set z-values to zero for 2D appearance
#    mode='lines', name='0.5 AU',
#    line=dict(color=gridc, width=1)
#))

fig.add_trace(go.Scatter3d(
    x=x_values*0.3,
    y=y_values*0.3,
    z=np.zeros(num_points),  # Set z-values to zero for 2D appearance
    mode='lines', name='0.3 AU',
    line=dict(color=gridc, width=1)
))

#fig.add_trace(go.Scatter3d(
#    x=x_values*0.1,
#    y=y_values*0.1,
#    z=np.zeros(num_points),  # Set z-values to zero for 2D appearance
#    mode='lines', name='0.1 AU',
#    line=dict(color=gridc, width=1)
#))


#add Sun-Earth line for HEEQ latitude 0
fig.add_trace(go.Scatter3d(
    x=np.linspace(0,1,num_points),
    y=np.zeros(num_points),
    z=np.zeros(num_points), 
    mode='lines', name='Sun-Earth line',
    line=dict(color=gridc, width=1)
))



zoom=1.0

fig.update_layout(
    scene=dict( aspectmode='data',
        camera=dict(
            eye=dict(x=0, y=-zoom, z=zoom),  # Set the position of the camera
            center=dict(x=0, y=0, z=0),      # Set the point the camera is looking at
            up=dict(x=0, y=0, z=1),          # Set the up vector of the camera
    ))
)


# Update layout for black background and styling
fig.update_layout(
    title='Spacecraft positions '+frame_time.strftime('%Y %b %d %H:00'),
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



msize=5


##### Orbits

# SOLO 
fig.add_trace(go.Scatter3d(x=solo['x'][solo_timeind:solo_timeind+fadeind]/au, y=solo['y'][solo_timeind:solo_timeind+fadeind]/au, 
                           z=solo['z'][solo_timeind:solo_timeind+fadeind]/au, name='Solar Orbiter',mode='lines',line=dict(color='lightgreen', width=3),
                           customdata=solo_hover,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°'  ))

fig.add_trace(go.Scatter3d(x=solo['x'][solo_timeind-fadeind:solo_timeind]/au, y=solo['y'][solo_timeind-fadeind:solo_timeind]/au, 
                           z=solo['z'][solo_timeind-fadeind:solo_timeind]/au, name='Solar Orbiter',mode='lines',line=dict(color='lightgreen', dash='dash' , width=3), 
                           customdata=solo_hover_past,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°', showlegend=False))

#PSP
fig.add_trace(go.Scatter3d(x=psp['x'][psp_timeind:psp_timeind+fadeind]/au, y=psp['y'][psp_timeind:psp_timeind+fadeind]/au, 
                           z=psp['z'][psp_timeind:psp_timeind+fadeind]/au, name='Parker Solar Probe',mode='lines',line=dict(color='white', width=3),
                           customdata=psp_hover,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°'  ))

fig.add_trace(go.Scatter3d(x=psp['x'][psp_timeind-fadeind:psp_timeind]/au, y=psp['y'][psp_timeind-fadeind:psp_timeind]/au, 
                           z=psp['z'][psp_timeind-fadeind:psp_timeind]/au, name='Parker Solar Probe',mode='lines',line=dict(color='white',dash='dash', width=3),
                           customdata=psp_hover_past,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°', showlegend=False))


#Bepi
fig.add_trace(go.Scatter3d(x=bepi['x'][bepi_timeind:bepi_timeind+fadeind]/au, y=bepi['y'][bepi_timeind:bepi_timeind+fadeind]/au, 
                           z=bepi['z'][bepi_timeind:bepi_timeind+fadeind]/au, name='BepiColombo',mode='lines',line=dict(color='lightskyblue', width=3),
                           customdata=bepi_hover,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°'  ))


fig.add_trace(go.Scatter3d(x=bepi['x'][bepi_timeind-fadeind:bepi_timeind]/au, y=bepi['y'][bepi_timeind-fadeind:bepi_timeind]/au, 
                           z=bepi['z'][bepi_timeind-fadeind:bepi_timeind]/au, name='BepiColomb',mode='lines',line=dict(color='lightskyblue',dash='dash', width=3), 
                           customdata=bepi_hover_past,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°',showlegend=False))


#STEREO-A
fig.add_trace(go.Scatter3d(x=sta['x'][sta_timeind:sta_timeind+fadeind]/au, y=sta['y'][sta_timeind:sta_timeind+fadeind]/au, 
                           z=sta['z'][sta_timeind:sta_timeind+fadeind]/au, name='STEREO-A',mode='lines',line=dict(color='tomato', width=3),
                           customdata=sta_hover,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°',showlegend=False))


fig.add_trace(go.Scatter3d(x=sta['x'][sta_timeind-fadeind:sta_timeind]/au, y=sta['y'][sta_timeind-fadeind:sta_timeind]/au, 
                           z=sta['z'][sta_timeind-fadeind:sta_timeind]/au, name='STEREO-A',mode='lines',line=dict(color='tomato',dash='dash', width=3), 
                           customdata=sta_hover_past,hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}°',showlegend=False))




#parker spiral
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
   
   

#### spacecraft

fig.add_trace(go.Scatter3d(x=[solo['x'][solo_timeind]/au], y=[solo['y'][solo_timeind]/au], z=[solo['z'][solo_timeind]/au], name='Solar Orbiter',
        mode='markers',marker=dict(color='lightgreen', size=msize, symbol=['square']),
        customdata=[solo_hover_now],hovertemplate='R %{customdata[0]} au<br>lon %{customdata[1]}°<br>lat %{customdata[2]}° '+'<extra></extra>'))
                           
#                           customdata
#                   hovertemplate='<b>%{customdata[1]}</b><br>' +
#                  'Date: %{x}<br>' +
#                  'Temperature: %{y}°C<br>' +
#                  'Humidity: %{customdata[0]}%<br>' +
#                  'Weather: %{text}<br>' +
#                  '<extra></extra>'
                         
                           
                           
#                           hovertemplate='Solar Orbiter %{text}',
#              text=[np.round(solo['r'][solo_timeind],2)+' '+solo['lat'][solo_timeind]+'  '+solo['lon'][solo_timeind] ] ))


fig.add_trace(go.Scatter3d(x=[psp['x'][psp_timeind]/au], y=[psp['y'][psp_timeind]/au], z=[psp['z'][psp_timeind]/au], name='Parker Solar Probe',
    mode='markers',marker=dict(color='white', size=msize, symbol=['square'])))
              #hovertemplate='%{text}', text=str(psp['r'][psp_timeind])   ))


fig.add_trace(go.Scatter3d(x=[sta['x'][sta_timeind]/au], y=[sta['y'][sta_timeind]/au], z=[sta['z'][sta_timeind]/au], name='STEREO-A',
    mode='markers',marker=dict(color='tomato', size=msize, symbol=['square'])))


fig.add_trace(go.Scatter3d(x=[bepi['x'][bepi_timeind]/au], y=[bepi['y'][bepi_timeind]/au], z=[bepi['z'][bepi_timeind]/au], name='BepiColombo',
    mode='markers',marker=dict(color='lightskyblue', size=msize,symbol=['square'])))


#fig.add_trace(go.Scatter3d(x=x[iwin], y=y[iwin], z=z[iwin], name='Wind',mode='markers',marker=dict(color='mediumseagreen', size=msize) ))
#        hovertemplate='Wind<br>ID: %{text}', text=ic.icmecat_id[iwin] ))

#fig.add_trace(go.Scatter3d(x=x[ista], y=y[ista], z=z[ista], name='STEREO-A',mode='markers',marker=dict(color='red', size=msize) ))
#        hovertemplate='STEREO-A<br>ID: %{text}', text=ic.icmecat_id[ista] ))

    
    
fig.show()        

fig.write_html(f'results/positions/position_3D.html')

#save as image
pio.write_image(fig, 'results/positions/position_3D.png',scale=1, width=2000, height=1200)

##pio.write_image(fig, 'results/positions/position_3D.png',scale=2, width=1500, height=850)


# In[5]:


clock_end = time.time()
print('done, all took ',np.round(clock_end - clock_start,2),' seconds.')


# In[ ]:





# In[ ]:





# In[ ]:





# In[6]:


sys.exit()


# In[ ]:





# In[ ]:





# In[ ]:





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





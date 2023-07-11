#distribute_web.py
#for putting the real time results to helioforecast.space every hour

#pip install pysftp

import sys
import os
import pysftp
from web_data import *
#sys.path.insert(0, r'/home/cmoestl/.ssh')

#(1) ############### copy to helioforecast.space server


print('(1) copy files to https://helioforecast.space/static/sync/')
print()


sftp=pysftp.Connection(server, username=user, private_key=key,port=port)

sftp.chdir('www-helioforecast/helioforecast/public/sync/')  #change dir
print(sftp.pwd)
print()


#print('copy predstorm_real.png/txt to ',sftp.pwd) #show current dir
sftp.put(path+'Wind_now.png')  # upload file to public/ on remote
sftp.put(path+'STEREO-A_beacon_14_days_now.png')
sftp.put(path+'OMNI2_now.png')


sftp.put(path+'NOAA_RTSW_PREDSTORM_55days_now.png')
sftp.put(path+'NOAA_RTSW_PREDSTORM_14days_now.png')
sftp.put(path+'NOAA_RTSW_PREDSTORM_3days_now.png')
sftp.put(path+'OMNI2_and_NOAA_RTSW_now.png')


print()

sftp.chdir('icme_solar_cycle/')  #change dir
print(sftp.pwd)



#sftp.get('remote_file')         # get a remote file
sftp.close()

print()
print()



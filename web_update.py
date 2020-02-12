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

#print('copy predstorm_real.png/txt to ',sftp.pwd) #show current dir
sftp.put(path+'Wind_now.png')  # upload file to public/ on remote
sftp.put(path+'STEREO-A_beacon_now.png')
sftp.put(path+'OMNI2_now.png')
sftp.put(path+'NOAA_RTSW_now.png')
sftp.put(path+'NOAA_RTSW_now2.png')

sftp.put('/home/cmoestl/pycode/heliocats/icmecat/HELCATS_ICMECAT_v20.txt')
sftp.put('/home/cmoestl/pycode/heliocats/icmecat/HELCATS_ICMECAT_v20_header.txt')
sftp.put('/home/cmoestl/pycode/heliocats/icmecat/HELCATS_ICMECAT_v20.xlsx')




#sftp.get('remote_file')         # get a remote file
sftp.close()

print()
print()



#(2) ############### copy to  chris dropbox


drop=0
#import dropbox
#from dropbox.files import WriteMode

if drop > 0:

    print('(2) copy files to chris dropbox')
    print()


    print('upload predstorm_real.png/txt and 3 small aurora gifs to my dropbox ')  
   
    db = dropbox.Dropbox(dropkey)
    #print('Account infos: ', db.users_get_current_account())


    fname = path+'predstorm/predstorm_real.png'     # Name einer lokalen Datei
    dname = '/predstorm_real.png'  # Name der Datei in Dropbox
    f = open(fname, 'rb')
    response = db.files_upload(f.read(), dname,mode=WriteMode('overwrite'))

    fname = path+'predstorm/predstorm_real.txt'    # Name einer lokalen Datei
    dname = '/predstorm_real.txt'  # Name der Datei in Dropbox
    f = open(fname, 'rb')
    response = db.files_upload(f.read(), dname,mode=WriteMode('overwrite'))
    #print('uploaded:', response)
    
    
    fname = path+'auroramaps/results/aurora_realtime/prob_global_small.gif'     # Name einer lokalen Datei
    dname = '/prob_global_small.gif'  # Name der Datei in Dropbox
    f = open(fname, 'rb')
    response = db.files_upload(f.read(), dname,mode=WriteMode('overwrite'))
    

    fname = path+'auroramaps/results/aurora_realtime/prob_europe_small.gif'     # Name einer lokalen Datei
    dname = '/prob_europe_small.gif'  # Name der Datei in Dropbox
    f = open(fname, 'rb')
    response = db.files_upload(f.read(), dname,mode=WriteMode('overwrite'))
    
    
    fname = path+'auroramaps/results/aurora_realtime/prob_canada_small.gif'     # Name einer lokalen Datei
    dname = '/prob_canada_small.gif'  # Name der Datei in Dropbox
    f = open(fname, 'rb')
    response = db.files_upload(f.read(), dname,mode=WriteMode('overwrite'))
    
   
    
    
    f.close()

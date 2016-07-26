#!/usr/bin/env python

## This is call back script which will be executed after created the grib2 files.
## 
## Hycom Model Input requires analysis of 06, 09, 12, 15, 18, 21-hours from 
## yesterday and 00 & 03-hours from today date. All 3-hourly forecasts from
## today date.
## 
## While creating tar ball, all files must be in present directory, so that
## when user extract it, will produce only files instead of entire paths!
##
## And finally putting into their ftp server.
##
## Arulalan.T
## 04-Mar-2016.

import os, subprocess, datetime, getopt, sys, glob

pbzip2 = '/gpfs1/home/Libs/GNU/ZIPUTIL/pbzip2'
pigz = '/gpfs1/home/Libs/GNU/ZIPUTIL/pigz'

def createTarBalls(path, today, utc, stephr=3):    

    cdir = os.getcwd()
    os.chdir(path)
    anal_ftemp = 'anal*%s*%s*.grb2'
    tDay = datetime.datetime.strptime(today, "%Y%m%d")
    lag1 = datetime.timedelta(days=1)
    yDay = (tDay - lag1).strftime('%Y%m%d')
    lag4 = datetime.timedelta(days=4)
    y4Day = (tDay - lag4).strftime('%Y%m%d')
    # get past 11th day timestamp
    y11Day = (tDay - datetime.timedelta(days=11)).strftime('%Y%m%d')
    
    if not os.path.exists('../TarFiles'): os.makedirs('../TarFiles')
    print "currnet path : ", os.getcwd()
    # normal "$ tar cvjf fcst_20160223.tar.bz2 *fcst*grb2" cmd takes 6 minutes 43 seconds.
    #
    # where as in parallel bz2, "$ tar -c *fcst*grb2 | pbzip2 -v -c -f -p32 -m500 > fcst_20160223_parallel.tar.bz2" cmd takes only just 23 seconds alone, with 32 processors and 500MB RAM memory.
    #
    
    # Lets create tar ball for every day individually.
    for i, hr in enumerate(range(0, 241, 24)):
        # create forecast files tar file in parallel # -m500 need to be include for pbzip2
        cmd = "tar -c ./fcst_{%s..%d}h%s*.grb2 | %s  -v  -c -f -p32 > %s/fcst_day%s_solar_ind_0.25_%s.tar.gz" % (str(hr+1).zfill(2), hr+24, today, pigz, '../TarFiles', str(i+1).zfill(2), today)
        print cmd
        subprocess.call(cmd, shell=True)
    # end of for i, hr in enumerate(range(0, 241, 24)):
        
    # delete today's forecasts files, after tar ball has been created!    
    cmd = "rm -rf fcst*%s*.grb2" % today
    print cmd
    subprocess.call(cmd, shell=True)
    
    # remove before 3 days directory, neither today nor yesterday directory!!!
    y4DayPath = os.path.join(path, '../%s' % y4Day)
    if os.path.exists(y4DayPath):    
        cmd = "rm -rf %s" % y4DayPath
        print cmd
        subprocess.call(cmd, shell=True)
    # end of if os.path.exists(y4DayPath):            
        
    tarpath = os.path.abspath('../TarFiles')
    # do scp the tar files to ftp_server and nkn_server    
    cmd = 'ssh ncmlogin3 "scp -p %s/fcst_*%s.tar.gz  %s:/data/niwe/NCUM_SOLAR_ENERGY/0.25/"' % (tarpath, today, ftp_server)
    print cmd
    subprocess.call(cmd, shell=True)
    
    # remove past 11th day tar ball from ftp_server 
    cmd = 'ssh ncmlogin3 "ssh %s rm -rf /data/niwe/NCUM_SOLAR_ENERGY/0.25/*%s*tar.bz2"' % (ftp_server, y11Day)
    print cmd
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        print "past 11th day tar ball has been removed from ftp_server, already", e    
        
    os.chdir(cdir)  
# end of def createTarBalls(path, today, ...):

if __name__ == '__main__':

    ftp_server="prod@ftp"
    date = None
    outpath = None
    oftype = None
    utc = None
    helpmsg = 'solar_energy_create_tarball_g2files_put_into_ftp.py --date=20160302 --outpath=path --oftype=analysis --utc=00'
    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:o:t:z:", ["date=","outpath=", "oftype=", "utc="])
    except getopt.GetoptError:
        print helpmsg
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print helpmsg
            sys.exit()
        elif opt in ("-d", "--date"):
            date = arg
        elif opt in ("-o", "--outpath"):
            outpath = arg 
        elif opt in ("-t", "--oftype"):
            oftype = arg
        elif opt in ("-z", "--utc"):
            utc = arg
    # end of for opt, arg in opts:
    
    # create tar balls only if forecast & utc is 00, otherwise skip it!    
    if oftype == 'forecast' and utc == '00': 
        # pass the arg to function  
        createTarBalls(outpath, date, utc, stephr=3)    
    # end of if oftype == 'forecast' and utc == '00': 

#!/usr/bin/env python

## This is call back script which will be executed after created the grib2 files.
## 
## OSF Model Input requires analysis of 06, 12, 18 -hours from 
## yesterday and 00 hour from today date. All 6-hourly forecasts from
## today date.
## 
## While creating tar ball, all files must be in present directory, so that
## when user extract it, will produce only files instead of entire paths!
##
## And finally putting into their ftp server.
##
## Arulalan.T
## 04-Mar-2016.

import os, subprocess, datetime, getopt, sys

pbzip2 = '/gpfs1/home/Libs/GNU/ZIPUTIL/pbzip2'

def createTarBalls(path, today, utc, stephr=6):
    # create tar balls only if present utc is 00, otherwise skip it!    
    if utc != '00': return
    cdir = os.getcwd()
    os.chdir(path)
    anal_ftemp = 'anal*%s*%s.grb2'
    tDay = datetime.datetime.strptime(today, "%Y%m%d")
    lag = datetime.timedelta(days=1)
    yDay = (tDay - lag).strftime('%Y%m%d')
    # get yesterday's analysis files from 06hr onwards
    yanal_files = [anal_ftemp % (str(hr).zfill(2), yDay) for hr in range(6, 24, stephr)]
    # get today's analysis 00 and 03 hr
    tanal_files = [anal_ftemp % (str(hr).zfill(2), today) for hr in range(0, 6, stephr)]
    for yf in yanal_files:
        # move yesterday's analysis files to today's directory
        cmd = 'mv ../%s/%s .' % (yDay, yf)
        subprocess.call(cmd, shell=True)
    # end of for yf in yanal_files:
    if not os.path.exists('../TarFiles'): os.makedirs('../TarFiles')
    
    # normal "$ tar cvjf fcst_20160223.tar.bz2 *fcst*grb2" cmd takes 6 minutes 43 seconds.
    #
    # where as in parallel bz2, "$ tar -c *fcst*grb2 | pbzip2 -v -c -f -p32 -m500 > fcst_20160223_parallel.tar.bz2" cmd takes only just 23 seconds alone, with 32 processors and 500MB RAM memory.
    #
    # create analysis files tar file in parallel
    anal_files = '  '.join(yanal_files + tanal_files)  
    cmd = "tar -c  %s | %s -v  -c -f -p32 -m500 > %s/ncum_anal_%s.tar.bz2" % (anal_files, pbzip2, '../TarFiles', today)
    subprocess.call(cmd, shell=True)
    # create forecast files tar file in parallel
    cmd = "tar -c fcst*.grb2 | %s  -v  -c -f -p32 -m500 > %s/ncum_fcst_%s.tar.bz2" % (pbzip2, '../TarFiles', today)
    subprocess.call(cmd, shell=True)
    
    # delete analysis & forecasts files, after tar ball has been created!
    cmd = "rm %s" % anal_files
    subprocess.call(cmd, shell=True)
    cmd = "rm ncum_fcst_%s*.grb2" % today
    subprocess.call(cmd, shell=True)
     
    # remove yesterday's empty directory, not today directory!!!
    yDayPath = os.path.join(path, '../%s' % yDay)
    os.rmdir(yDayPath) 
    
    # do scp the tar files to ftp_server and nkn_server
    cmd = 'rsh ncmr0102 "scp -p ../TarFiles/ncum_anal_%s.tar.bz2  %s:/data/ftp/pub/outgoing/NCUM_INCOIS/OSF/"' % (today, ftp_server)
    subprocess.call(cmd, shell=True)
    cmd = 'rsh ncmr0102 "scp -p ../TarFiles/ncum_anal_%s.tar.bz2  %s:NCUM/osf/"' % (today, nkn_server)
    subprocess.call(cmd, shell=True)
    
    cmd = 'rsh ncmr0102 "scp -p ../TarFiles/ncum_fcst_%s.tar.bz2  %s:/data/ftp/pub/outgoing/NCUM_INCOIS/OSF/"' % (today, ftp_server)
    subprocess.call(cmd, shell=True)
    cmd = 'rsh ncmr0102 "scp -p ../TarFiles/ncum_fcst_%s.tar.bz2  %s:NCUM/osf/"' % (today, nkn_server)
    subprocess.call(cmd, shell=True)
    
    os.chdir(cdir)   
# end of def createTarBalls(path, today, ...):

if __name__ == '__main__':

    nkn_server="incois@nkn"
    ftp_server="prod@ftp"
    date = None
    outpath = None
    oftype = None
    utc = None
    helpmsg = 'osf_create_tarball_g2files_put_into_ftp.py --date=20160302 --outpath=path --oftype=analysis --utc=00'
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
    
    # pass the arg to function  
    createTarBalls(outpath, date, utc, stephr=6)

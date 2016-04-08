#!/usr/bin/env python

## This is call back script which will be executed after created the grib2 files.
## 
## OSF Model Input requires 6-hourly fcst files
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

def createTarBalls(path, oftype, today, utc, stephr=6):
    
    cdir = os.getcwd()
    os.chdir(path)
    tarpath = os.path.abspath('../TarFiles')
    if not os.path.exists(tarpath): os.makedirs(tarpath)
        
    if oftype == 'forecast':                
        # create fcst files tar file in parallel
        cmd = "tar -c ./fcst*%s*.grb2 | %s  -v  -c -f -p32 -m500 > %s/ncum_fcst_glb_0.5_%s.tar.bz2" % (today, pbzip2, '../TarFiles', today)
        print cmd
        subprocess.call(cmd, shell=True)
            
        # delete today's forecasts files, after tar ball has been created!    
        cmd = "rm -rf fcst*%s*.grb2" % today
        print cmd
        subprocess.call(cmd, shell=True)    
        
        # do scp the fcst tar files to ftp_server and nkn_server
        cmd = 'ssh ncmlogin3 "scp -p %s/ncum_fcst_glb_0.5_%s.tar.bz2  %s:/data/ftp/pub/outgoing/NCUM_INCOIS/OSF/0.5/"' % (tarpath, today, ftp_server)
        print cmd
        subprocess.call(cmd, shell=True)
        cmd = 'ssh ncmlogin3 "scp -p %s/ncum_fcst_glb_0.5_%s.tar.bz2  %s:NCUM/osf/0.5/"' % (tarpath, today, nkn_server)
        print cmd
        subprocess.call(cmd, shell=True)
    # end of if oftype == 'forecast':      
    
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
    
    # create tar balls only if utc is 00, otherwise skip it!    
    if oftype == 'forecast' and utc == '00': 
        # pass the arg to function  
        createTarBalls(outpath, oftype, date, utc, stephr=6)    
    # end of if oftype == 'forecast' and utc == '00': 

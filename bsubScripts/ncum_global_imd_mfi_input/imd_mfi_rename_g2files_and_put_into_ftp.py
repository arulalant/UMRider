#!/usr/bin/env python

## This is call back script which will be executed after created the grib2 files.
## 
## IMD MFI System require grib2 filenames must follows their own nomenclature
## In UMRider we can not create like that. So here we are just renaming the 
## grib2 files as per their nomenclature!
## 
## IBM people (GP Singh) will take care to put into ftp and further renaming it.

## Arulalan.T
## 15-Mar-2016.

import os, sys, getopt, subprocess


def renameFiles(outpath):
    cdir = os.getcwd()
    os.chdir(outpath)
    gfiles = [f for f in os.listdir(outpath) if f.endswith('.grb2')]
    
    mfi_name_structure = {'anal': 'YMXA00', 'fc186': 'YMXK00', 
    'fc180': 'YMXK00', 'fc240': 'YMXK00', 'fc102': 'YMXK00', 'fc108': 'YMXK00',
    'fc78': 'YMXK00', 'fc54': 'YMXX01', 'fc222': 'YMXK00', 'fc144': 'YMXK00',
    'fc204': 'YMXK00', 'fc120': 'YMXK00', 'fc90': 'YMXK00', 'fc96': 'YMXK00', 
    'fc168': 'YMXK00', 'fc228': 'YMXK00', 'fc126': 'YMXK00', 'fc12': 'YMXC00',
    'fc36': 'YMXG00', 'fc30': 'YMXF00', 'fc18': 'YMXD00', 'fc198': 'YMXK00', 
    'fc192': 'YMXK00', 'fc138': 'YMXK00', 'fc84': 'YMXK00', 'fc210': 'YMXK00',
    'fc114': 'YMXK00', 'fc72': 'YMXK00', 'fc60': 'YMXJ00', 'fc234': 'YMXK00',
    'fc174': 'YMXK00', 'fc42': 'YMXH00', 'fc216': 'YMXK00', 'fc132': 'YMXK00',
    'fc48': 'YMXI00', 'fc156': 'YMXK00', 'fc162': 'YMXK00', 'fc150': 'YMXK00',
     'fc24': 'YMXE00', 'fc06': 'YMXB00', 'fc66': 'YMXX02'}
    
    for gf in gfiles:
        if not 'MFICODE' in gf: continue
        mfi_g2out = gf.split('.grb2')[0]
        mfi_g2out_name = mfi_g2out.split('_')
        hr = mfi_g2out_name[-1].split('.')[0]
        if not hr in mfi_name_structure: 
            print "Couldnt get the IMD MFI name for hour ", hr 
            continue
        rval = mfi_name_structure.get(hr, 'MFICODE')
        # updating RANDOM
        mfi_g2out_name[1] = rval
        mfi_g2out_name = '_'.join(mfi_g2out_name)
        # lets rename the grib2 files as per IMD MFI requirements !
        os.rename(gf, mfi_g2out_name)        
        print "created IMD MFI standard grib2 file from %s to %s" % (gf, mfi_g2out_name)        
    # end of for gf in gfiles:
    os.chdir(cdir)
# end of def renameFiles(outpath):

if __name__ == '__main__':
    
    nkn_server="imd@nkn"
    ftp_server="prod@ftp"
    date = None
    outpath = None
    oftype = None
    utc = None
    helpmsg = 'imd_mfi_rename_g2files_and_put_into_ftp.py --date=20160302 --outpath=path --oftype=analysis --utc=00'
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
    renameFiles(outpath)

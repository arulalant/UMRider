#!/usr/bin/env python

## This is call back script which will be executed after created the grib2 files.
## 
## After created NCUM_EPS grib2 files (which contain all 45 members in it), lets
## do calculate mean over all 45 members and then convert to 2.5X2.5 degree 
## resolution, finally convert to grib1 format for vsdb input purpose.
##
## Arulalan.T
## 25-July-2016.

import os, subprocess, datetime, getopt, sys, iris, numpy, time 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../g2utils')))
from cubeutils import cubeRealizationAverager
from um2grb2 import tweaked_messages, getCubeData


neededVars = [
    ## Pressure Level Variable names & STASH codes
    ('geopotential_height', 'm01s16i202'),      
    ('x_wind', 'm01s15i243'),
    ('y_wind', 'm01s15i244'), 
    ('air_temperature', 'm01s16i203'),
    ('relative_humidity', 'm01s16i256'),
    ## Non Pressure Level Variable names & STASH codes
    ('air_pressure_at_sea_level', 'm01s16i222'),
    ]

#Define _precipVars_
# The following vars should contains only precipitation, rainfall, snow 
# variables, those whose regrid extrapolate should be only in 'linear' mode
# and not in 'mask' mode, and should not have -ve values.
_precipVars_ = [('precipitation_amount', 'm01s05i226'),]

g2ctl = "/gpfs2/home/umtid/Softwares/grib2ctl/g2ctl.pl"
grib2ctl = "/gpfs2/home/umtid/Softwares/grib2ctl/grib2ctl.pl"
gribmap = "/gpfs1/home/Libs/GNU/GRADS/grads-2.0.2.oga.1/Contents/gribmap"
cnvgrib = "/gpfs1/home/Libs/INTEL/CNVGRIB/CNVGRIB-1.4.1/cnvgrib-1.4.1/cnvgrib"
wgrib2 = "/gpfs1/home/Libs/GNU/WGRIB2/v2.0.1/wgrib2"

_targetGridFile_ = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data/sample_global_2p5X2p5_73X144.grib2'))
_targetGrid_ = iris.load(_targetGridFile_)[0]

__grib1FilesNameSuffix__ = None
__wgrib2Arguments__ = ' -set_bin_prec 12 -set_grib_type complex2 -grib_out '
    
def createENSavg_VSDB_Grib1Files(inpath, outpath, today, utc, start_long_fcst_hour, stephr=24):    


    day = int(start_long_fcst_hour) / int(stephr)
    pfileday = str(day).zfill(2)
    pfilename = 'umeps_prg_1cntl_44ens_24hourly_day' + pfileday
    needed_fname = pfilename + '_' + today + '_' + utc + 'Z.grib2'
    files = [f for f in os.listdir(inpath) if f == needed_fname]
    if not files: return 
    
    tDay = datetime.datetime.strptime(today, "%Y%m%d")        
    tvDay = tDay.strftime('%d%m%y')
    
    opath = os.path.join(outpath, today)
    if not os.path.exists(opath): os.makedirs(opath)
        
    g2filepath = os.path.join(opath, 'prg_' + today + pfileday + '.grib2')
    
    ensfpath = os.path.join(inpath, files[0])
    inf = iris.load(ensfpath)
    for varName, varSTASH  in neededVars:
        # define variable name constraint
        varConstraint = iris.Constraint(name=varName)
        print varName, varSTASH
        ensCube = inf.extract(varConstraint)
        if not ensCube: continue        
        print "doing members average for", varName
        # do average over ensembles (1 control + 44 ensemble members)
        ensAvgCube = cubeRealizationAverager(ensCube[0])
        # make memory free
        del ensCube
        
        if (varName, varSTASH) in _precipVars_:
            # DO NOT APPLY iris.analysis.Linear(extrapolation_mode='mask'), 
            # which writes nan every where for the snowfall_flux,  
            # rainfall_flux, precipitation_flux. So donot apply that.         
            exmode = 'linear'
        else:
            # In general all the other variables should not be 
            # extrapolated over masked grid points.
            exmode = 'mask'
        # end of if (...):
        
        scheme = iris.analysis.Linear(extrapolation_mode=exmode)
        try:
            # This lienar interpolate will do extra polate over ocean even 
            # though original data doesnt have values over ocean and wise versa.
            # So lets be aware of this.            
            regdCube = ensAvgCube.regrid(_targetGrid_, scheme)
        except Exception as e:
            print "ALERT !!! Error while regridding!! %s" % str(e)
            print " So skipping this without saving data"
            continue
        # end of try:      
        
        # make memory free
        del ensAvgCube
        
        if (varName, varSTASH) in _precipVars_:
            # Since we are not using 'mask' option for extrapolate while 
            # doing linear regrid, which bring -ve values after regrid in 
            # extrapolated grids. So lets make it as 0 as minimum value.
            regdCube.data[regdCube.data < 0.0] = 0.0
        # end of if (varName, varSTASH) in _precipVars_:
        
        if exmode == 'mask':            
            regdCube.data = numpy.ma.masked_array(regdCube.data, 
                                dtype=numpy.float64, fill_value=9.999e+20) 
        
        print regdCube
        # save into grib2 file 
        iris.fileformats.grib.save_messages(tweaked_messages([regdCube]), 
                                                     g2filepath, append=True) # save grib2 file        
        
        print "Appending %s to grib2 file" % varName
        # make memory free
        del regdCube       
    # end of for varName, varSTASH  in neededVars:    
                    
    wg2filename = 'prg%d%sz%s.grib2' % (day, utc.zfill(2), tvDay)        
    wg2filepath = os.path.join(opath, wg2filename)
        
    # execute post wgrib2 command compression algorithm
    cmd = "%s %s %s %s" % (wgrib2, g2filepath, __wgrib2Arguments__, wg2filepath)
    print cmd
    subprocess.call(cmd, shell=True)            
    time.sleep(10)
    # remove the grib2 file generated by IRIS
    os.remove(g2filepath)
    # rename g2filepath as wg2filepath
    g2filepath = wg2filepath
    print "Created grib2 file using wgrib2 command with compress arguments " 
    
    # Conver grib2 to grib1 
    g1filepath = '.'.join(g2filepath.split('.')[:-1])
    g1filepath = g1filepath if g1filepath else g2filepath[:-1]
    if __grib1FilesNameSuffix__: g1filepath += str(__grib1FilesNameSuffix__)
    
    if os.path.isfile(g1filepath): os.remove(g1filepath)
    
    cmd = [cnvgrib, '-g21', g2filepath, g1filepath]
    subprocess.call(cmd, shell=False)
    cmd = ['chmod', '644', g1filepath]
    subprocess.call(cmd, shell=False)
    print "Converted grib2 to grib1 file : -", g1filepath
    os.remove(g2filepath)                                
# end of def createTarBalls(path, today, ...):

if __name__ == '__main__':

    
    date = None
    inpath = None
    oftype = None
    utc = None
    outpath = '/gpfs4/home/umtid/um2grb2/ArulTest/NCUM_EPS_VSDB'
    
    helpmsg = './ncumeps_create_memavg_vsdb_input.py --date=20160302 --outpath=path --oftype=forecast --utc=00 --start_long_fcst_hour=24 --end_long_fcst_hour=24 --fcst_step_hour=24'
    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:o:t:z:s:e:i", ["date=",
            "outpath=", "oftype=", "utc=", "start_long_fcst_hour=", "end_long_fcst_hour=", "fcst_step_hour="])
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
            inpath = arg 
        elif opt in ("-t", "--oftype"):
            oftype = arg
        elif opt in ("-z", "--utc"):
            utc = arg
        elif opt in ("-s", "--start_long_fcst_hour"):
            start_long_fcst_hour = arg
        elif opt in ("-e", "--end_long_fcst_hour"):
            end_long_fcst_hour = arg
        elif opt in ("-i", "--fcst_step_hour"):
            fcst_step_hour = arg
    # end of for opt, arg in opts:
    
    # create tar balls only if forecast & utc is 00, otherwise skip it!    
    if oftype == 'forecast' and utc == '00' and fcst_step_hour == '24': 
        # pass the arg to function  
        print "calling functrion"
        createENSavg_VSDB_Grib1Files(inpath, outpath, date, utc, start_long_fcst_hour, stephr=fcst_step_hour)    
    # end of if oftype == 'forecast' and utc == '00': 

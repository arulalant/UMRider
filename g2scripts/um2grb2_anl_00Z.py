"""
This is simple script to invoke parallel conversion function from um2grb2
and pass the assimilated / forecasted hour as argument.

UTC : 00Z
fcst hour : 00hr
Output : It creates analysis - 1 file (um_ana_000hr_date_00Z.grib2) 
Written by : Arulalan.T
Date : 07.Dec.2015
"""

import os, sys, datetime
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from loadconfigure import inPath, outPath, tmpPath, date, loadg2utils, \
                debug, targetGridResolution, overwriteFiles, neededVars, \
                requiredLat, requiredLon, anlOutGrib2FilesNameStructure, \
                createGrib2CtlIdxFiles, createGrib1CtlIdxFiles, \
                convertGrib2FilestoGrib1Files, grib1FilesNameSuffix, \
                removeGrib2FilesAfterGrib1FilesCreated, \
                anl_step_hour, pressureLevels, callBackScript, \
                setGrib2TableParameters, anl_aavars_reference_time, \
                anl_aavars_time_bounds, wgrib2Arguments, \
                soilFirstSecondFixedSurfaceUnit, UMtype, targetGridFile, \
                UMInAnlFiles, UMInShortFcstFiles, fillFullyMaskedVars

if loadg2utils == 'system':
    # Load g2utils from system python which has installed through setup.py
    from g2utils.um2grb2 import convertAnlFiles
    print "INFO : imported g2utils.um2grb2 from system python"
elif loadg2utils == 'local':
    # Load g2utils from previous directory for the operational purpose, 
    # where normal user don't have write permission to change the g2utils!
    g2utils_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                            '../g2utils'))
    sys.path.append(g2utils_path)
    from um2grb2 import convertAnlFiles
    print "INFO : imported g2utils.um2grb2 from local previous directory"
    print "loaded from g2utils_path : ", g2utils_path

if isinstance(date, tuple):
    # Got tuple of date string.
    startdate, enddate = date
    print "Got tuple dates"
    print "So um2grb2 ana 00hr conversion - from %s to %s" % (startdate, enddate)
elif isinstance(date, str):
    # only single date
    startdate, enddate = date, date
    print "Got single string date"
    print "So um2grb2 ana 00hr conversion - on %s" % date
# end of if isinstance(date, tuple):

sDay = datetime.datetime.strptime(startdate, "%Y%m%d")
eDay = datetime.datetime.strptime(enddate, "%Y%m%d")
lag = datetime.timedelta(days=1)
while sDay <= eDay:
    # loop through until startdate incremented upto enddate
    print "Going to start progress on", startdate
    # call analysis conversion function w.r.t data assimilated 
    # during short forecast hour - 00UTC.
    convertAnlFiles(inPath, outPath, tmpPath, 
                     UMtype=UMtype, UMInAnlFiles=UMInAnlFiles,
                        UMInShortFcstFiles=UMInShortFcstFiles,
                                targetGridFile=targetGridFile,
                    targetGridResolution=targetGridResolution, 
             date=startdate, utc='00', convertVars=neededVars, 
                  latitude=requiredLat, longitude=requiredLon,
                                pressureLevels=pressureLevels, 
                       fillFullyMaskedVars=fillFullyMaskedVars,
soilFirstSecondFixedSurfaceUnit=soilFirstSecondFixedSurfaceUnit,
            anlFileNameStructure=anlOutGrib2FilesNameStructure, 
                 createGrib2CtlIdxFiles=createGrib2CtlIdxFiles,
                 createGrib1CtlIdxFiles=createGrib1CtlIdxFiles,
   convertGrib2FilestoGrib1Files=convertGrib2FilestoGrib1Files,
                     grib1FilesNameSuffix=grib1FilesNameSuffix,          
   removeGrib2FilesAfterGrib1FilesCreated=removeGrib2FilesAfterGrib1FilesCreated,  
                                  anl_step_hour=anl_step_hour,
          anl_aavars_reference_time=anl_aavars_reference_time,
                anl_aavars_time_bounds=anl_aavars_time_bounds,
                       overwrite=overwriteFiles, lprint=debug,
              setGrib2TableParameters=setGrib2TableParameters,
                              wgrib2Arguments=wgrib2Arguments,
                                callBackScript=callBackScript)
    print "Time lag incremented by 1"
    sDay += lag
    startdate = sDay.strftime('%Y%m%d')
# end of while sDay <= eDay:
print "Successfully completed all the dates till", enddate

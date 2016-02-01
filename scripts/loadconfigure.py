"""
This is script used to load all the parameters from configure text file 
and cross check either all the paths are valid or not.

Written by : Arulalan.T
Date : 07.Dec.2015
"""

import os, sys, time, datetime  
# get this script abspath
scriptPath = os.path.dirname(os.path.abspath(__file__))

# get environment variable for setup file otherwise load default local path 
setupfile = os.environ.get('UMRIDER_SETUP', os.path.join(scriptPath, 'um2grb2_setup.cfg')).strip()
if not os.path.isfile(setupfile):
    raise ValueError("UMRIDER_SETUP file doesnot exists '%s'" % setupfile)

print "loaded UMRIDER_SETUP configure file from ", setupfile
print "Reading configure file to load the paths"
# get the configure lines
clines = [l.strip() for l in open(setupfile).readlines() \
                if not l.startswith(('#', '/', '!', '\n', '%'))]
if not clines:
    raise ValueError("Empty setup list loaded from %s" % setupfile)

# get the dictionary keys, values
cdic = {k.strip(): v.strip() for k,v in [l.split('=') for l in clines]}

# store into local variables
inPath = cdic.get('inPath', None)  
outPath = cdic.get('outPath', None)
tmpPath = cdic.get('tmpPath', None)

startdate = cdic.get('startdate', 'YYYYMMDD')
enddate = cdic.get('enddate', None)
loadg2utils = cdic.get('loadg2utils', 'system')
overwriteFiles = eval(cdic.get('overwriteFiles', 'True'))
debug = eval(cdic.get('debug', 'False'))

targetGridResolution = cdic.get('targetGridResolution', None)
targetGridResolution = None if targetGridResolution in ['None', None] else float(targetGridResolution)

# check the variable's path 
for name, path in [('inPath', inPath), ('outPath', outPath), ('tmpPath', tmpPath)]:
    if path is None:
        raise ValueError("In configure file, '%s' path is not defined !" % name)
    if not os.path.exists(path):
        raise ValueError("In configure file, '%s = %s' path does not exists" % (name, path))
    print name, " = ", path
# end of for name, path in [...]:

# get the current date if not specified
if startdate == 'YYYYMMDD': startdate = time.strftime('%Y%m%d')
if enddate == 'YYYYMMDD': enddate = time.strftime('%Y%m%d')
if startdate in ['None', None]: raise ValueError("Start date can not be None")
if enddate not in ['None', None]:
    sDay = datetime.datetime.strptime(startdate, "%Y%m%d")
    eDay = datetime.datetime.strptime(enddate, "%Y%m%d")
    if sDay > eDay:
        raise ValueError("Start date must be less than End date or wise-versa")
    if sDay == eDay:
        raise ValueError("Both start date and end date are same")     
    
    date = (startdate, enddate)   
else:
    date = startdate
# end of if enddate not in ['None', None]:

# get environment variable for vars file otherwise load default local path 
varfile = os.environ.get('UMRIDER_VARS', os.path.join(scriptPath, 'um2grb2_vars.cfg')).strip()
if not os.path.isfile(varfile):
    raise ValueError("UMRIDER_VARS file doesnot exists '%s'" % varfile)
print "loaded UMRIDER_VARS configure file from ", varfile
print "Reading configure file to load the paths"

# clean it up and store as list of tuples contains both varName and varSTASH
vl = [i.strip() for i in open(varfile).readlines() if i]
neededVars = [j if len(j) == 2 else j[0] for j in 
                    [eval(i) for i in vl if i and not i.startswith('#')]]
if not neededVars:
    raise ValueError("Empty variables list loaded from %s" % varfile)

print "\n" * 4
print "*" * 80
print "date = ", date
print "targetGridResolution = ", targetGridResolution
print "loadg2utils = ", loadg2utils
print "overwriteFiles = ", overwriteFiles
print "debug = ", debug
print "Successfully loaded the above params from UMRIDER_SETUP configure file!", setupfile
print "*" * 80
print "Successfully loaded the below variables from UMRIDER_VARS configure file!", varfile 
print "\n".join([str(i+1)+' : ' + str(tu) for i, tu in enumerate(neededVars)])
print "The above %d variables will be passed to UMRider" % len(neededVars)
print "*" * 80
print "\n" * 4

"""
This is script used to load all the parameters from configure text file 
and cross check either all the paths are valid or not.

Written by : Arulalan.T
Date : 07.Dec.2015
"""

import os, sys, time, datetime  
# get this script abspath
scriptPath = os.path.dirname(os.path.abspath(__file__))

print "Reading configure file to load the paths"
# get the configure lines
clines = [l.strip() for l in open(os.path.join(scriptPath, 'configure')).readlines() \
                if not l.startswith(('#', '/', '!', '\n', '%'))]
# get the dictionary keys, values
cdic = {k.strip(): v.strip() for k,v in [l.split('=') for l in clines]}

# store into local variables
inPath = cdic.get('inPath', None)  
outPath = cdic.get('outPath', None)
tmpPath = cdic.get('tmpPath', None)

startdate = cdic.get('startdate', 'YYYYMMDD')
enddate = cdic.get('enddate', None)
loadg2utils = cdic.get('loadg2utils', 'system')
debug = cdic.get('debug', False)
debug = True if debug == 'True' else False

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
    
print "date = ", date
print "loadg2utils = ", loadg2utils
print "debug = ", debug
print "Successfully loaded the above params from configure file!"

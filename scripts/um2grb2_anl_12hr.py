"""
This is simple script to invoke parallel conversion function from um2grb2
and pass the assimilated / forecasted hour as argument.

hour : 12
Output : It creates analysis - 1 file (um_ana_012hr_date.grib2).

Written by : Arulalan.T
Date : 07.Dec.2015
"""

import os, sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from loadconfigure import inPath, outPath, tmpPath, date, loadg2utils, debug

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

# call analysis conversion function w.r.t data assimilated during short forecast hour.
convertAnlFiles(inPath, outPath, tmpPath, date, hr='12', lprint=debug)


We have to copy the following files of iris package which are updated by me for the purpose of NCMRWF grib2 write/read conversion.

__init__.py
_grib_cf_map.py
_load_convert.py
_save_rules.py
load_rules.py

Copy the above scripts into the following path
~/miniconda2/envs/iris-1.10.dev1/lib/python2.7/site-packages/iris/fileformats/grib/
 

Note : Dont forget to take copy of original scripts in the installed path. Also check the diff between these files with sys installed iris path scripts.

ARULALAN.T <arulalan@ncmrwf.gov.in>
12-Sep-2016

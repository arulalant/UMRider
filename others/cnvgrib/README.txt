Get the required / supporting libraries from here
http://www.nco.ncep.noaa.gov/pmb/codes/GRIB2/ and read the installation
instruction from Build_cnvgrib.pdf
<http://www.nco.ncep.noaa.gov/pmb/codes/GRIB2/Build_cnvgrib.pdf> in the
same link.

STEP 1:

$ cd g2lib-1.4.0

# vim params.f   ---> end of file, I included NCUM grib2 to grib1 param table.
At the moment (11-04-2016) supported for 5 variables. Infuture, we need to
include all the NCUM grib2 param codes.

$ cp ~/params.f .  ----> Just copy param.f which we already edited and kept in this folder.

$ make clean
$ make


STEP 2:

$ cd w3lib-2.0.2
$ make clean
$ make

STEP 3:

$ cd cnvgrib-1.4.1
$ make clean
$ make

Regards,
Arulalan.T

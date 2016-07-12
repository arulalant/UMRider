#!/bin/bash
set -x

date=$1
cyc=$2

eleventhDay=$3


OUTDIR=/gpfs3/home/umfcst/ShortJobs/NCUMReg_WIND_ENERGY/0.04/CSVFiles/$date
mkdir -p $OUTDIR

cp /gpfs3/home/umfcst/ShortJobs/NCUMReg_WIND_ENERGY/0.04/TarFiles/WindEnergy.India.0.04.$date$cyc.tar.gz  $OUTDIR/
cp -r timeseries_reg.gs latlon.csv ${OUTDIR}/

cd $OUTDIR

rm -f *.nc
rm -f *tar.gz
rm -f wind*csv

tar -zxf WindEnergy.India.0.04.$date$cyc.tar.gz
for tt in {01..72..01} ; do
wgrib2 -netcdf ncum_reg_fcst_${tt}h${date}.nc ncum_reg_fcst_${tt}h${date}.grb2
done
cdo cat *.nc final.nc
rm -f ncum*
grads -blc timeseries_reg.gs

OUTFILE=wind50reg${date}${cyc}.zip

zip $OUTFILE wind50*csv
rm -f wind50*csv
rm timeseries_reg.gs

ssh ncmlogin3 "scp -p $OUTDIR/$OUTFILE  prod@ftp:/data/energon/4km/"

OUTFILE11Day=wind50reg${eleventhDay}${cyc}.zip

ssh ncmlogin3 "ssh prod@ftp rm -rf /data/energon/4km/$OUTFILE11Day"

exit



#!/bin/bash
#
#BSUB -a poe                  # set parallel operating environment
#BSUB -J umeps2g2anl          # job name
#BSUB -W 06:00                # wall-clock time (hrs:mins)
#BSUB -n 16                  # number of tasks in job
#BSUB -x
#BSUB -R span[ptile=16]
#BSUB -R rusage[mem=61440]
#BSUB -q small             	  # queue
#BSUB -e /gpfs3/home/umeps/UMRiderLogs/post/bsub/um2grb2.anl.00hr.err.%J.%I.hybrid     # error file name in which %J is replaced by the job ID
#BSUB -o /gpfs3/home/umeps/UMRiderLogs/post/bsub/um2grb2.anl.00hr.out.%J.%I.hybrid     # output file name in which %J is replaced by the job ID


# find out the directory of this bash script after submitted to bsub
DIR="$( cd "$( dirname "${BASH_SOURCE[1]}" )" && pwd )"

# get the absolute path of the local table 
localTable_relative_dir="$DIR/../../tables/local/ncmr/v1/"
localTable_absolute_dir="$( cd "$localTable_relative_dir" && pwd )"
localTable=$localTable_absolute_dir/ncmr_grib2_local_table

# get the absolute path of the script for forecast 00utc
g2scripts_relative_dir="$DIR/../../g2scripts/"
g2scripts_absolute_dir="$( cd "$g2scripts_relative_dir" && pwd )"
g2script=$g2scripts_absolute_dir/umeps2grb2_fcst_00Z.py

# export the configure paths to needed variables
export UMRIDER_SETUP=$DIR/ncumeps_global_post_um2grb2_24hourly_setup.cfg  ## yes, this is correct program.
export UMRIDER_VARS=$DIR/ncumeps_global_post_um2grb2_anl_vars.cfg
export GRIB2TABLE=$localTable

echo "export UMRIDER_SETUP="$UMRIDER_SETUP
echo "export UMRIDER_VARS="$UMRIDER_VARS
echo "export GRIB2TABLE="$GRIB2TABLE

# sourcing umtid_bashrc to load module python-uvcdat-iris!
#source "$DIR/../umtid_bashrc"
# execute the script
#srun -n 1 -c 64 python $g2script

export SHELL=/bin/bash

hour=0  # 0 will produce 0th hour prognostic data (equivalent to analysis)
echo "hour="${hour}

python $g2script --start_long_fcst_hour=${hour} --end_long_fcst_hour=${hour}

#!/bin/bash
#
#BSUB -L /bin/bash
#BSUB -a poe                  # set parallel operating environment
#BSUB -J um2grb2              # job name
#BSUB -W 06:00                # wall-clock time (hrs:mins)
#BSUB -n 6                    # number of tasks in job
#BSUB -q small                # queue
#BSUB -e um2grb2.anl.00hr.err.%J.hybrid     # error file name in which %J is replaced by the job ID
#BSUB -o um2grb2.anl.00hr.out.%J.hybrid     # output file name in which %J is replaced by the job ID

/gpfs2/home/umtid/Pythons/Python-2.7.9/bin/python /gpfs2/home/umtid/UMRider/scripts/um2grb2_anl_00Z.py


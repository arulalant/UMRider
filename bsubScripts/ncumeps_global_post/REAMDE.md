This will produce forecast ensembles 00utc files.

$ bsub -J umeps2g2[$i] < ncumeps_global_post_fcst_00Z_6hourly.bash

This will produce long forecast 00utc files with 6-hourly intervals by packing all 45 members together.

$ bsub -J umeps2g2[$i] < ncumeps_global_post_fcst_00Z_24hourly.bash

This will produce long forecast 00utc files with 24-hourly intervals by packing all 45 members together.

$ bsub -J umeps2vsdbanl < ncumeps_global_vsdb_anl.bash 
$ bsub -J umeps2vsdbfcst[$i] < ncumeps_global_vsdb_fcst.bash 

This will produce long forecast 00utc files with 24-hourly intervals by averaging all 45 members and convert to grib1 standard for vsdb input.

Arulalan.T
